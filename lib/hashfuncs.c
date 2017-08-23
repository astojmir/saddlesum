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

#include <string.h>
#include "miscutils.h"

/* sdbm function from http://www.cse.yorku.ca/~oz/hash.html */
unsigned int hash_from_string(void *str)
{
	unsigned char *s = (unsigned char *) str;
        unsigned int hash = 0;
        int c;

        while ((c = *s++)) {
		hash = c + (hash << 6) + (hash << 16) - hash;
	}
        return hash;
}

int str_equal(void *str1, void *str2)
{
	if (strcmp((char *) str1, (char *) str2) == 0) {
		return 1;
	}
	return 0;
}

unsigned int hash_from_int(void *i_)
{
        uint32_t *i = (uint32_t *) i_;
        return *i;
 }

int int_equal(void *i_, void *j_)
{
        uint32_t *i = (uint32_t *) i_;
        uint32_t *j = (uint32_t *) j_;
        if (*i == *j) {
                return 1;
        }
	return 0;
}
