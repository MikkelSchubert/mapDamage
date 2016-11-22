/* The MIT License

   Copyright (c) 20082-2012 by Heng Li <lh3@me.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <Python.h>

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <stdarg.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


static PyObject* SeqtkError;


typedef struct {
	int n, m;
	uint64_t *a;
} reglist_t;

unsigned char seq_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15 /*'-'*/,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

unsigned char seq_nt16to4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
int bitcnt_table[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };


int strbuff_init(kstring_t* str)
{
	str->l = 0;
	str->m = 1024;
	str->s = (char*)malloc(str->m);

	if (str->s) {
		str->s[str->l] = 0;
	}

	return (int)str->s;
}


PyObject* strbuff_finalize(kstring_t* str)
{
	PyObject* result = Py_BuildValue("s", str->s);
	free(str->s);

	return result;
}


int strbuff_printf(kstring_t* str, const char* fmt, ...)
{
	va_list args1, args2;
	va_start(args1, fmt);
	va_copy(args2, args1);

	const int len = vsnprintf(NULL, 0, fmt, args1);
	if (len + 1 > str->m - str->l) {
		do {
			str->m *= 2;
		} while (len + 1 > str->m - str->l);

		str->s = realloc(str->s, str->m);
		if (!str->s) {
			va_end(args1);
			va_end(args2);
			PyErr_NoMemory();
			return 1;
		}
	}

	vsprintf(str->s + str->l, fmt, args2);

	str->l += len;
	str->s[str->l] = '\0';

	va_end(args1);
	va_end(args2);
	return 0;
}


int _comp(kstring_t* str, gzFile fp)
{
	kseq_t *seq;
	long l = 0;
	int err = 0;
	reglist_t dummy;

	seq = kseq_init(fp);
	dummy.n= dummy.m = 1; dummy.a = calloc(1, 8);
	while ((l = kseq_read(seq)) >= 0) {
		int i, k;
		reglist_t *p = 0;
		p = &dummy;
		dummy.a[0] = l;

		for (k = 0; p && k < p->n; ++k) {
			int beg = p->a[k]>>32, end = p->a[k]&0xffffffff;
			int la, lb, lc, na, nb, nc, cnt[11];
			if (beg > 0) la = seq->seq.s[beg-1], lb = seq_nt16_table[la], lc = bitcnt_table[lb];
			else la = 'a', lb = -1, lc = 0;
			na = seq->seq.s[beg]; nb = seq_nt16_table[na]; nc = bitcnt_table[nb];
			memset(cnt, 0, 11 * sizeof(int));
			for (i = beg; i < end; ++i) {
				int is_CpG = 0, a, b, c;
				a = na; b = nb; c = nc;
				na = seq->seq.s[i+1]; nb = seq_nt16_table[na]; nc = bitcnt_table[nb];
				if (b == 2 || b == 10) { // C or Y
					if (nb == 4 || nb == 5) is_CpG = 1;
				} else if (b == 4 || b == 5) { // G or R
					if (lb == 2 || lb == 10) is_CpG = 1;
				}

				if (c > 1) ++cnt[c+2];
				if (c == 1) ++cnt[seq_nt16to4_table[b]];
				if (b == 10 || b == 5) ++cnt[9];
				else if (c == 2) {
					++cnt[8];
				}
				if (is_CpG) {
					++cnt[7];
					if (b == 10 || b == 5) ++cnt[10];
				}

				la = a; lb = b; lc = c;
			}

			if (strbuff_printf(str, "%s\t%ld", seq->name.s, l)) {
				err = 1;
				break;
			}

			for (i = 0; i < 11; ++i) {
				if (strbuff_printf(str, "\t%d", cnt[i])) {
					err = 1;
					break;
				}
			}

			if (err || strbuff_printf(str, "\n")) {
				err = 1;
				break;
			}
		}
	}

	// Errors only relevant for FASTQ files
	if (l == -3) {
        PyErr_SetString(SeqtkError, "qual string has different length than sequence");
	} else if (l == -2) {
        PyErr_SetString(SeqtkError, "missing quality string in file");
	}

	free(dummy.a);
	kseq_destroy(seq);

	return (l < -1) || err;
}


static PyObject*
seqtk_comp(PyObject* self, PyObject* args)
{
	errno = 0;

	const char* filename = NULL;
	if (!PyArg_ParseTuple(args, "s", &filename)) {
		return NULL;
	}

	kstring_t buffer;
	if (!strbuff_init(&buffer)) {
		return PyErr_NoMemory();
	}

	printf("Opening '%s'\n", filename);
	gzFile fp = gzopen(filename, "r");
	if (!fp) {
		return PyErr_SetFromErrno(PyExc_OSError);
	}

	if (_comp(&buffer, fp)) {
		gzclose(fp);
		return NULL;
	}

	gzclose(fp);

	return strbuff_finalize(&buffer);
}


static PyMethodDef SeqtkMethods[] = {
	{"comp",  seqtk_comp, METH_VARARGS,
	 "Call the 'seqk comp' function."},
	{NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
initseqtk(void)
{
	PyObject* module = Py_InitModule("seqtk", SeqtkMethods);
	if (!module) {
		return;
	}

    SeqtkError = PyErr_NewException("mapdamage.seqtk.error", NULL, NULL);
    Py_INCREF(SeqtkError);
    PyModule_AddObject(module, "error", SeqtkError);
}
