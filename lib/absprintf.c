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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "miscutils.h"

#define ABSPRINTF_BUF 256


PrintBuf *PrintBuf_init(PrintBuf *pbuf)
{
        if (pbuf == NULL) {
                pbuf = malloc_(sizeof(PrintBuf));
                pbuf->owned = 1;
        }
        else {
                pbuf->owned = 0;
        }
        pbuf->buf = malloc_(ABSPRINTF_BUF);
        pbuf->size = ABSPRINTF_BUF;
        pbuf->len = 0;
        return pbuf;
}


void PrintBuf_delete(PrintBuf *pbuf)
{
        free(pbuf->buf);
        if (pbuf->owned) {
                free(pbuf);
        }
}


void PrintBuf_printf(PrintBuf *pbuf, int pos, const char *fmt, ...)
{
	va_list ap;
	size_t diff;
	int n;
        pbuf->len = pos >= 0 ? pos : pbuf->len;
	while (1) {
		/* Try to print in the allocated space. */
		diff = pbuf->size - pbuf->len;
		va_start(ap, fmt);
		n = vsnprintf (pbuf->buf + pbuf->len, diff, fmt, ap);
		va_end(ap);
		/* If that worked, get out. */
		if (n > -1 && n < diff) {
			pbuf->len += n;
			break;
		}
		/* Else try again with more space. */
		if (n > -1) {   /* glibc 2.1 */
			pbuf->size += max(ABSPRINTF_BUF, n+1-diff);
                }
		else {           /* glibc 2.0 */
			pbuf->size += ABSPRINTF_BUF;
                }
		pbuf->buf = realloc_(pbuf->buf, pbuf->size);
	}
}

