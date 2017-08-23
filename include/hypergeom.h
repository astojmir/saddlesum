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

#ifndef _HYPERGEOM_H
#define _HYPERGEOM_H
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
        unsigned int n;
        unsigned int c;
        double *logFact;
} HypergeomStats;

HypergeomStats *HypergeomStats_init(unsigned int n, unsigned int c);
void  HypergeomStats_del(HypergeomStats *data);
double HypergeomStats_pvalue(HypergeomStats *data, unsigned int S, unsigned int m);


#ifdef __cplusplus
}
#endif
#endif /* !_HYPERGEOM_H */
