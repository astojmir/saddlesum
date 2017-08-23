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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "miscutils.h"
#include "hypergeom.h"


HypergeomStats *HypergeomStats_init(unsigned int n, unsigned int c)
{
        HypergeomStats *data;
        register double x = 0.0;
        register unsigned int i;
        register double *F;

        data = calloc_(1, sizeof(HypergeomStats));
        data->logFact = malloc_((n+1) * sizeof(double));
        data->n = n;
        data->c = c;

        F = data->logFact;
        *F++ = 0.0;
        for (i=1; i < n+1; i++, F++) {
                x += log((double) i);
                *F = x;
        }
        return data;
}


void  HypergeomStats_del(HypergeomStats *data)
{
        free(data->logFact);
        free(data);
}


double HypergeomStats_pvalue(HypergeomStats *data, unsigned int S, unsigned int m)
{
        register double *F = data->logFact;
        unsigned int n = data->n;
        unsigned int c = data->c;
        unsigned int K = min(c, m)+1;
        double neg_log_nCc = F[c] + F[n-c] - F[n];
        register unsigned int i;
        double x;
        double pval = 0.0;
        double old_pval = -1.0;
        

        for (i=S; i < K; i++, old_pval=pval) {
                x = F[m] + F[n-m] - F[i] - F[m-i] - F[c-i] - F[n-m-c+i] + neg_log_nCc;
                pval += exp(x);
		if (old_pval == pval) {
			break;
		}
        }
        return pval;
}
