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

#include <R.h>
#include <Rdefines.h>
#include "saddlesum.h"

static void saddleSumFinalizer(SEXP ptr)
{
    if(!R_ExternalPtrAddr(ptr)) return;
    SADDLE_SUM_del(R_ExternalPtrAddr(ptr));
    R_ClearExternalPtr(ptr);
}

SEXP saddleSumCreate(SEXP bw)
{
    SEXP ans, ptr;
    double *bkgrnd_weights = NUMERIC_POINTER(bw);
    int num_weights = GET_LENGTH(bw);
    SDDLSUM *data;

    data = SADDLE_SUM_init(bkgrnd_weights, num_weights);

    if (data != NULL) {
	PROTECT(ans = NEW_LIST(0));
	PROTECT(ptr = R_MakeExternalPtr(data, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(ptr, saddleSumFinalizer, TRUE);
	SET_ATTR(ans, install("saddleSum_ptr"), ptr);
	UNPROTECT(2);
    }
    else {
	ans = R_NilValue;
    }
    return ans;
}

SEXP saddleSumPvalue(SEXP pdt, SEXP psc, SEXP pnh, SEXP pcp, SEXP pmi, SEXP ptl)
{

    double *score = NUMERIC_POINTER(psc);
    int *num_hits = INTEGER_POINTER(pnh);
    double *cutoff_pvalue = NUMERIC_POINTER(pcp);
    int *maxiter = INTEGER_POINTER(pmi); 
    double *tol = NUMERIC_POINTER(ptl);
    SDDLSUM *data;
    double pval;
    SEXP ans, ptr;

    ptr = GET_ATTR(pdt, install("saddleSum_ptr"));
    data = R_ExternalPtrAddr(ptr);
    if (data) {
	pval = SADDLE_SUM_pvalue(data, *score, *num_hits, *cutoff_pvalue, *maxiter, *tol);
	PROTECT(ans = NEW_NUMERIC(1));
	NUMERIC_POINTER(ans)[0] = pval;
	UNPROTECT(1);
    }
    else {
	ans = R_NilValue;
    }
    return ans;
}
