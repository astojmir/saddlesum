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
#include <string.h>
#include <math.h>
#include "saddlesum.h"

const int INITIAL_LAMBDAS = 2;

extern double SQRT2;
extern double SQ2OPI;
extern double ndtr(double Z);


static int dbl_compare (const void * M1, const void * M2)
{
  const double *i1 = M1;
  const double *i2 = M2;
  if (*i1 > *i2)
    return 1;
  else if (*i1 < *i2)
    return -1;
  else
    return 0;
}

static int SADDLE_SUM_bisect(SDDLSUM *data, double x)
{
    int a, b, c;
    LMDBITEM **lmbd = data->sorted_lambdas;

    if (data->num_lambdas == 0) {
	return 0;
    }

    if (x > lmbd[data->num_lambdas-1]->mean)
	return data->num_lambdas;
    if (x < lmbd[0]->mean)
	return 0;

    a = 0;
    b = data->num_lambdas - 1;
    while (a < b-1) {
	c = a + ((b - a) / 2);
	if (lmbd[c]->mean < x)
	    a = c;
	else
	    b = c;
    }
    return b;
}

static LMDBITEM *SADDLE_SUM_add_item(SDDLSUM *data, double lmbd)
{
    int i;
    int N = data->num_weights;
    LMDBITEM *item;
    double tmp, w;
    double Nrho = 0.0;
    double Nrho1 = 0.0;
    double Nrho2 = 0.0;
    double D1K, D2K;
    const double wmax = data->max_weight;

    if (data->cur_item == NULL) {
	item = malloc(sizeof(LMDBITEM));
	if (item  == NULL) {
	    return NULL;
	}
	data->cur_item = item;
    }
    else {
	item = data->cur_item;
    }


    /* Here factor out wmax for improved numerical stability. */
    /* For the same reason, sum smallest to largest (bkgrnd_weights are sorted */
    /* in  SADDLE_SUM_init) */ 
    for (i=0; i < N; i++) {
	w = data->bkgrnd_weights[i];
	tmp = exp(lmbd*(w-wmax));
	Nrho += tmp;
	tmp *= w;
	Nrho1 += tmp;
	tmp *= w;
	Nrho2 += tmp;
    }
    D1K = Nrho1 / Nrho;
    D2K = Nrho2/Nrho - D1K*D1K;

    item->mean = D1K;
    item->lambda = lmbd;
    item->D2K = D2K;
    item->expH = Nrho * exp(lmbd*(wmax-D1K)) / N;

    item->C = 2*lmbd*sqrt(D2K);
    /* By omitting sgn(lmbd) here we explicitely assume lmbd > 0 */
    item->D = -SQRT2 * sqrt(lmbd*(D1K-wmax) - log(Nrho) + log(N));

    return item;
}

static double LMBD_ITEM_pvalue(LMDBITEM *item, int m)
{
    double sqrtm, phi;
    sqrtm = sqrt(m);

    /* Exit if we are too close to the mean */
    if (item->D * sqrtm > -1.0) 
	return 1.0;

    /* Lugannani-Rice formula */
    phi = SQ2OPI * pow(item->expH, m);
    return ndtr(item->D*sqrtm) + (phi / item->C / sqrtm) + (phi / item->D / sqrtm / 2);
}


static LMDBITEM *SADDLE_SUM_insert_item(SDDLSUM *data, LMDBITEM *item)
{
    int i, n;
    LMDBITEM **ppl;

    if (data->num_lambdas+1 >= data->size_lambdas) {
	ppl = realloc(data->sorted_lambdas, 2*data->size_lambdas*sizeof(LMDBITEM *));
	if (ppl  == NULL) {
	    return NULL;
	}
	data->size_lambdas *= 2;
	data->sorted_lambdas = ppl;
    }

    i = SADDLE_SUM_bisect(data, item->mean);
    if (i < data->num_lambdas) {
	n = data->num_lambdas - i;
	memmove(data->sorted_lambdas+i+1, data->sorted_lambdas+i, n*sizeof(LMDBITEM *));
    }
    data->num_lambdas++;
    data->sorted_lambdas[i] = item;
    data->cur_item = NULL;
    return item;
}


SDDLSUM *SADDLE_SUM_init(double *bkgrnd_weights, int num_weights)
{
    SDDLSUM *data;
    int i;
    double sum = 0;


    data = calloc(sizeof(SDDLSUM),1);
    if (data  == NULL)
	return NULL;

    data->bkgrnd_weights = malloc(num_weights * sizeof(double));
    if (data->bkgrnd_weights == NULL) {
	SADDLE_SUM_del(data);
	return NULL;
    }

    data->sorted_lambdas = calloc(INITIAL_LAMBDAS, sizeof(LMDBITEM *));
    if (data->sorted_lambdas == NULL) {
	SADDLE_SUM_del(data);
	return NULL;
    }

    data->cur_item = NULL;
    data->size_lambdas = INITIAL_LAMBDAS;
    data->num_weights = num_weights;

    /* Copy and sort background weights */
    for (i=0; i < num_weights; i++) {
	data->bkgrnd_weights[i] = bkgrnd_weights[i];
	sum += bkgrnd_weights[i];
    }
    data->bkgrnd_mean = sum / num_weights;
    qsort(data->bkgrnd_weights, num_weights, sizeof(double), dbl_compare);
    data->max_weight = data->bkgrnd_weights[num_weights-1];

    return data;
}

void SADDLE_SUM_del(SDDLSUM *data)
{
    int i;
    LMDBITEM **ppl;

    free(data->cur_item);
    free(data->bkgrnd_weights);
    for (ppl=data->sorted_lambdas, i=0; i < data->num_lambdas; i++,ppl++) {
	free(*ppl);
    }
    free(data->sorted_lambdas);
    free(data);
}

double SADDLE_SUM_pvalue(SDDLSUM *data, double score, int num_hits,
			 double cutoff_pvalue, int maxiter, double tol) 
{
    int i, iter;
    LMDBITEM *item;
    double x = score / num_hits;
    double y;
    double ya = 0.0;
    double yb = -1.0;
    double yc = 0.0;
    double pval = 1.0;
    double min_pval = pow(1.0/data->num_weights, num_hits);
    double diff_means;

    if (x <= data->bkgrnd_mean)
	return 1.0;

    /* Make initial guess of lambda and bracket it */
    i = SADDLE_SUM_bisect(data, x);
    if (i > 0) {
	ya = data->sorted_lambdas[i-1]->lambda;
    }
    if (i < data->num_lambdas) {
	item = data->sorted_lambdas[i];
	yb = item->lambda;
	yc = 0.5*(ya+yb);
	pval = LMBD_ITEM_pvalue(item, num_hits); 
	if ((pval > cutoff_pvalue) ||  (yb-ya) < tol){
	    /* Skip Newton's method */
	    maxiter = 0;
	}
    }
    else {
	/* We need to establish the right bracket */
	yb = 0.5;
	do {
	    yb *= 2.0;
	    item = SADDLE_SUM_add_item(data, yb);
	    item = SADDLE_SUM_insert_item(data, item);
	} while (fabs(item->mean - data->max_weight) >= tol);
	yc = 0.5*(ya+yb);
    }

    /* Now iterate Newton's method */
    for (iter=0; iter < maxiter; iter++) {
	
	item = SADDLE_SUM_add_item(data, yc);
	if (item == NULL) {
	    /* If we cannot allocate memory, we can at least report the previous
	       approximation.*/ 
	    break;
	}
	pval = LMBD_ITEM_pvalue(item, num_hits); 

	diff_means = item->mean - x;
	if (diff_means < 0.0) {
	    ya = yc;
	}
	else {
	    yb = yc;
	    if (pval > cutoff_pvalue) {
		/* Here x < item->mean => pval is an underestimate. */
		/* By monotonicity of the saddlepoint function, this implies that the */
		/* true value of lambda will give even larger p-value. Hence we can */
		/* terminate the iteration right here without computing the exact */
		/* solution.  */
		break;
	    }
	}
	
	/* Newton's approximation */
	y = yc - diff_means / item->D2K;

	/* Don't allow the approximations to leave the bracket. */
	if ((y < ya) || (y > yb)) {
	    y = 0.5*(ya+yb);
	}

	/* Termination condition */
	if ((fabs(y-yc) < tol) || (fabs(diff_means) < tol)) {
	    break;
	}

	item = SADDLE_SUM_insert_item(data, item);
	if (item == NULL) {
	    break;
	}
	yc = y;
    }
    
    /* Ensure the pvalue is not lower or higher than reasonable */
    pval = pval > min_pval ? pval : min_pval;
    pval = pval < 1.0 ? pval : 1.0;
    
    return pval;
}
