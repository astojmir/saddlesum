/*
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* Code author:  Aleksandar Stojmirovic
*
* Reference: A. Stojmirovic and Y-K Yu. Robust and accurate data enrichment
*            statistics via distribution function of sum of weights.
*            Bioinformatics, 26(21):2752-2759, 2010.
*
*/

#ifndef _UTILS_H
#define _UTILS_H
#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdio.h>
#include <inttypes.h>

/* Inline macros */
#define max(a, b)  ((a) > (b) ? (a) : (b))
#define min(a, b)  ((a) < (b) ? (a) : (b))

/* Memory allocation with error checking */
void *calloc_(size_t nmemb, size_t size);
void *malloc_(size_t size);
void *realloc_(void *ptr, size_t size);
char *strdup_(const char *s);

/* Hash functions to be used with hashtable */
unsigned int hash_from_string(void *str);
int str_equal(void *str1, void *str2);
unsigned int hash_from_int(void *i_);
int int_equal(void *i_, void *j_);

/* Printing into buffers that can grow */
typedef struct _PrintBuf_s {
        char *buf;
        size_t size;
        size_t len;
        uint8_t owned;
} PrintBuf;

PrintBuf *PrintBuf_init(PrintBuf *pbuf);
void PrintBuf_delete(PrintBuf *pbuf);
void PrintBuf_printf(PrintBuf *pbuf, int pos, const char *fmt, ...);

/* Binary file helper routines */
void fread_(void *ptr, size_t size, size_t nmemb, FILE *stream);
void fread_uint32(uint32_t *ptr, size_t nmemb, FILE *stream);
char *fread_buf(FILE *stream, uint32_t *n);
void fskip_buf(FILE *stream);
void strbuf2array(char **dest, char *buf, int n);

#ifdef __cplusplus
}
#endif
#endif /* !_UTILS_H */
