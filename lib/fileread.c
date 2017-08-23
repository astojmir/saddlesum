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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "miscutils.h"


void fread_(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
        size_t count = fread(ptr, size, nmemb, stream);
        if (count != nmemb) {
                fprintf(stderr, "Problem reading file.");
                exit(EXIT_FAILURE);
        }
}

void fread_uint32(uint32_t *ptr, size_t nmemb, FILE *stream)
{
        uint32_t *x = ptr;
        uint32_t *y = x + nmemb;
        for (; x < y; x++) {
                *x = fgetc(stream) & 0xFF;            /* read lowest byte */
                *x |= (fgetc(stream) & 0xFF) << 8;    /* read 2nd byte */
                *x |= (fgetc(stream) & 0xFF) << 16;   /* read 3rd byte */
                *x |= (fgetc(stream) & 0xFF) << 24;   /* read highest byte */
        }
}

char *fread_buf(FILE *stream, uint32_t *n)
{
        char *buf;
        fread_uint32(n, 1, stream);
        buf = malloc_(*n);
        fread_(buf, 1, *n, stream);
        return buf;
}

void fskip_buf(FILE *stream)
{
        uint32_t n;
        fread_uint32(&n, 1, stream);
        fseek(stream, n, SEEK_CUR);
}


void strbuf2array(char **dest, char *buf, int n)
{
        int i;
        for (i=0; i < n; i++, dest++) {
                *dest = buf;
                buf += strlen(buf) + 1;
        }
}
