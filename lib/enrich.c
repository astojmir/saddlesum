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
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "miscutils.h"
#include "fsfile.h"
#include "hashtable.h"
#include "enrich.h"
#include "saddlesum.h"
#include "hypergeom.h"

#ifndef WEIGHTS_MAX_FIELD_WIDTH
#define WEIGHTS_MAX_FIELD_WIDTH 128
#endif

#define INITIAL_TERM_HITS 16

#ifndef SADDLESUM_MAX_ITERS
#define SADDLESUM_MAX_ITERS 50
#endif
#ifndef SADDLESUM_TOLERANCE
#define SADDLESUM_TOLERANCE 1.0e-11
#endif



EnrichContext *EnrichContext_init(const char *db_name,
				  uint32_t min_term_size,
                                  double Evalue_cutoff,
                                  double effective_db_size,
                                  EnrichStats statistics_type,
				  TransformType transform_type,
                                  uint8_t discretized_weights,
				  CutoffType cutoff_type,
				  uint32_t rank_cutoff,
				  double weight_cutoff,
                                  uint8_t use_all_weights)
{

        EnrichContext *cntxt = calloc_(1, sizeof(EnrichContext));
        cntxt->db_name = db_name;
        cntxt->min_term_size = min_term_size;
        cntxt->Evalue_cutoff = Evalue_cutoff;
        cntxt->effective_db_size = effective_db_size;
        cntxt->statistics_type = statistics_type;
	cntxt->transform_type = transform_type;
        cntxt->discretized_weights = discretized_weights;
	cntxt->cutoff_type = cutoff_type;
        cntxt->rank_cutoff = rank_cutoff;
        cntxt->weight_cutoff = weight_cutoff;
        cntxt->use_all_weights = use_all_weights;
	cntxt->term_hits = calloc_(INITIAL_TERM_HITS, sizeof(TermHit));
	cntxt->max_term_hits = INITIAL_TERM_HITS;
        return cntxt;
}

void EnrichContext_delete(EnrichContext *cntxt)
{
        int i;
        while (cntxt->first_warning != NULL) {
                cntxt->last_warning = cntxt->first_warning;
                cntxt->first_warning = cntxt->first_warning->next;
                free(cntxt->last_warning->msg);
                free(cntxt->last_warning);
        }
        cntxt->first_warning = NULL;
        cntxt->last_warning = NULL;

        if (cntxt->weights != NULL) {
                for (i=0; i < cntxt->num_entities; i++) {
                        if (cntxt->input_symbols[i] != NULL) {
                                free(cntxt->input_symbols[i]);
                        }
                }
                free(cntxt->input_symbols);
                free(cntxt->weights);
                free(cntxt->used_indices);
                cntxt->weights = NULL;
                cntxt->used_indices = NULL;
        }
        if (cntxt->term_hits != NULL) {
                free(cntxt->term_hits);
                cntxt->term_hits = NULL;
        }
        free(cntxt);
}


static
void EnrichResults_add_warning(EnrichContext *cntxt, EntityWarning *warning)
{
        if (warning != NULL) {
                if (cntxt->first_warning == NULL) {
                        cntxt->first_warning = warning;
                        cntxt->last_warning = warning;
                }
                else {
                        cntxt->last_warning->next = warning;
                        cntxt->last_warning = warning;
                }
                switch (warning->code) {
                case UNKNOWN_ID:
                        cntxt->num_unknown_ids++;
                        break;
                case UNRESOLVABLE_CONFLICT:
                        cntxt->num_conflicting_ids++;
                        break;
                case DUPLICATE_ID:
                        cntxt->num_duplicate_ids++;
                        break;
                case RESOLVABLE_CONFLICT:
                        cntxt->num_resolvable_ids++;
                        break;
                }
        }
}


static int dbl_reverse_compare (const void *M1, const void *M2)
{
        const double *v1 = (const double *) M1;
        const double *v2 = (const double *) M2;
        return (*v1 < *v2) - (*v1 > *v2);
}



static int TermHit_compare (const void * M1, const void * M2)
{
        const TermHit *th1 = (const TermHit *) M1;
        const TermHit *th2 = (const TermHit *) M2;
        register int retval;

        /* First compare by namespace - compare pointers since they ought to indicate
           the order of entry */
        retval = (th1->term->namespace > th2->term->namespace) - (th1->term->namespace < th2->term->namespace);

        if (!retval) {
                retval = (th1->Evalue > th2->Evalue) - (th1->Evalue < th2->Evalue);
                retval = retval ? retval : (th1 > th2) - (th1 < th2);
        }
        return retval;
}

static void EnrichResults_insert_term_hit(EnrichContext *cntxt, uint32_t term_index,
					  double score, uint32_t num_entities, double Pvalue)
{
	TermHit *term_hit;
	if (cntxt->num_term_hits >= cntxt->max_term_hits) {
		cntxt->max_term_hits *= 2;
		cntxt->term_hits = realloc_(cntxt->term_hits, cntxt->max_term_hits*sizeof(TermHit));
	}
	term_hit = cntxt->term_hits + cntxt->num_term_hits++;
	term_hit->term = NULL;
	term_hit->term_index = term_index;
	term_hit->score = score;
	term_hit->num_entities = num_entities;
	term_hit->Evalue = Pvalue;
	term_hit->Pvalue = Pvalue;
}

void EnrichResults_load_weights(EnrichContext *cntxt, const char *weights_filename,
				EntityDb *entity_db, TermMappingDb *mapping_db)
{

	unsigned int line_num = 1;
        uint32_t entity_index;
        double wght;
        char *endptr;
        char buf1[WEIGHTS_MAX_FIELD_WIDTH];
        char buf2[WEIGHTS_MAX_FIELD_WIDTH];

        PrintBuf *pbuf = PrintBuf_init(NULL);;
	FSFile_proc_code_type retcode;
        EntityWarning *warning;

	FILE *fp = stdin;

        uint8_t *mapped_indices;
	uint32_t term_index;
	uint32_t *hits;
	uint32_t *end_hits;
	uint32_t num_hits;


        if (strcmp(weights_filename, "-")) {
                fp = fopen(weights_filename, "r");
                if (fp == NULL) {
                        fprintf(stderr, "Could not open file %s.", weights_filename);
                        exit(EXIT_FAILURE);
                }
        }

        cntxt->num_entities = entity_db->num_entities;
        cntxt->weights = calloc_(entity_db->num_entities, sizeof(double));
        cntxt->used_indices = calloc_(entity_db->num_entities, sizeof(uint8_t));
        cntxt->input_symbols = calloc_(entity_db->num_entities, sizeof(char *));

        /* For background, we either use all reckognised weights or only those */
        /* weights that have hits mapping onto them. */
        mapped_indices = malloc_(entity_db->num_entities * sizeof(uint8_t));
        mapping_db->reset(mapping_db);
        if (cntxt->use_all_weights) {
                memset(mapped_indices, 1, entity_db->num_entities);
        }
        else {
                memset(mapped_indices, 0, entity_db->num_entities);
                while (mapping_db->get_next_mapping(mapping_db, &term_index, &hits, &num_hits)) {
                        for (end_hits=hits+num_hits; hits < end_hits; hits++) {
                                mapped_indices[*hits] = 1;
                        }
                }
        }

	while (END_OF_FILE != (retcode = SSFile_next_field(fp, buf1, WEIGHTS_MAX_FIELD_WIDTH))) {

 		if (retcode != END_OF_FIELD) {
			fprintf(stderr, "Invalid weight file format (line %d, field #1).\n",
				line_num);
			exit(EXIT_FAILURE);
		}

                warning = entity_db->map_symbol(entity_db, buf1, &entity_index);
                EnrichResults_add_warning(cntxt, warning);

		retcode = SSFile_next_field(fp, buf2, WEIGHTS_MAX_FIELD_WIDTH);
		if (retcode != END_OF_LINE) {
			fprintf(stderr, "Invalid weight file format (line %d, field #2).\n",
				line_num);
			exit(EXIT_FAILURE);
		}

                wght = strtod(buf2, &endptr);
                if (*endptr != '\0') {
			fprintf(stderr, "Invalid weight %s (line %d, field #2).\n",
				buf2, line_num);
			exit(EXIT_FAILURE);
                }

                if (warning == NULL || warning->code == RESOLVABLE_CONFLICT) {

                        if (cntxt->used_indices[entity_index]) {
                                warning = calloc_(1, sizeof(EntityWarning));
                                PrintBuf_printf(pbuf, 0,
                                                "Duplicate weight for %s (line %d)"
                                                " - additional instance IGNORED.",
                                                buf1, line_num);
                                warning->msg = strdup_(pbuf->buf);
                                warning->code = DUPLICATE_ID;
                                EnrichResults_add_warning(cntxt, warning);
                        }
                        else if (mapped_indices[entity_index]) {
                                cntxt->weights[entity_index] = wght;
                                cntxt->used_indices[entity_index] = 1;
                                cntxt->num_valid_ids++;
                                cntxt->input_symbols[entity_index] = strdup_(buf1);
                        }
                }

		line_num++;
                cntxt->num_raw_weights++;
	}
        PrintBuf_delete(pbuf);
        free(mapped_indices);
        fclose(fp);
        cntxt->num_unused_entities = cntxt->num_entities - cntxt->num_valid_ids;
}


void EnrichResults_process_weights(EnrichContext *cntxt)
{

        double *used_weights;
        uint32_t i;
        uint32_t j;

	/* Apply transformations */
	switch (cntxt->transform_type) {
	case NO_TRANSFORM:
		break;
	case FLIP:
		for (i=0; i < cntxt->num_entities; i++) {
			cntxt->weights[i] = -cntxt->weights[i];
		}
		break;
	case ABS:
		for (i=0; i < cntxt->num_entities; i++) {
			cntxt->weights[i] = fabs(cntxt->weights[i]);
		}
		break;
	}

        /* Apply cutoffs */
        /* All top weights up to and including weights cutoff are left */
        if (cntxt->cutoff_type == RANK) {
                used_weights = malloc_(cntxt->num_valid_ids * sizeof(double));
                for (i=0,j=0; i < cntxt->num_entities; i++) {
                        if (cntxt->used_indices[i]) {
                                used_weights[j++] = cntxt->weights[i];
                        }
                }
                qsort(used_weights, cntxt->num_valid_ids, sizeof(double), dbl_reverse_compare);
		cntxt->weight_cutoff = used_weights[cntxt->rank_cutoff - 1];
		free(used_weights);
        }
        if (cntxt->cutoff_type != NONE) {
                for (i=0,j=0; i < cntxt->num_entities; i++) {
                        if (cntxt->used_indices[i]) {
				if (cntxt->weights[i] < cntxt->weight_cutoff) {
					cntxt->weights[i] = 0.0;
				}
				else {
					j++;
				}
			}
                }
		cntxt->rank_cutoff = j;
	}

	/* Discretize weights */
        if (cntxt->discretized_weights) {
                for (i=0; i < cntxt->num_entities; i++) {
                        if (cntxt->used_indices[i] && cntxt->weights[i] > 0.0) {
                                cntxt->weights[i] = 1.0;
                        }
			else {
                                cntxt->weights[i] = 0.0;
			}
                }
        }

	/* Must count here for non-zero weights since any of previous
	   transformations could have affected them */
	for (i=0; i < cntxt->num_entities; i++) {
		if (cntxt->used_indices[i] && cntxt->weights[i] != 0.0) {
			cntxt->num_nonzero_valid_ids++;
		}
	}
}


static void EnrichResults_wsum_pvalues(EnrichContext *cntxt, CVTermDb *term_db,
				       TermMappingDb *mapping_db)
{
	uint32_t term_index;
	uint32_t *hits;
	uint32_t *end_hits;
	uint32_t num_hits;
	uint32_t num_used_hits;

        SDDLSUM *sddlsum;
        double score;
        double Pvalue;
        double *used_weights;
        unsigned int i;
        unsigned int j;

	/* Copy used weights to obtain background distribution */
	used_weights = malloc_(cntxt->num_valid_ids * sizeof(double));
	for (i=0,j=0; i < cntxt->num_entities; i++) {
		if (cntxt->used_indices[i]) {
			used_weights[j++] = cntxt->weights[i];
		}
	}
	sddlsum = SADDLE_SUM_init(used_weights, cntxt->num_valid_ids);
	if (sddlsum == NULL) {
		fprintf(stderr, "Could not allocate saddlesum context.\n");
		exit(EXIT_FAILURE);
	}
	free(used_weights);

	/* SADDLESUM - main loop */
	mapping_db->reset(mapping_db);
	while (mapping_db->get_next_mapping(mapping_db, &term_index, &hits, &num_hits)) {
		score = 0.0;
		for (num_used_hits=0, end_hits=hits+num_hits; hits < end_hits; hits++) {
			score += cntxt->weights[*hits];
			num_used_hits += cntxt->used_indices[*hits];
		}
		if (num_used_hits < cntxt->min_term_size) {
			continue;
		}
		Pvalue = SADDLE_SUM_pvalue(sddlsum, score, num_used_hits,
					   cntxt->Pvalue_cutoff,
					   SADDLESUM_MAX_ITERS,
					   SADDLESUM_TOLERANCE);

		if (Pvalue <= cntxt->Pvalue_cutoff) {
			EnrichResults_insert_term_hit(cntxt, term_index,
						      score, num_used_hits, Pvalue);
		}
	}
	SADDLE_SUM_del(sddlsum);
}


static
void EnrichResults_wsum_single_pvalue(EnrichContext *cntxt, CVTermDb *term_db,
                                      TermMappingDb *mapping_db, uint32_t term_index)
{
	uint32_t *hits;
	uint32_t *end_hits;
	uint32_t num_hits;
	uint32_t num_used_hits;

        SDDLSUM *sddlsum;
        double score;
        double Pvalue = -1.0;
        double *used_weights;
        unsigned int i;
        unsigned int j;

	/* Copy used weights to obtain background distribution */
	used_weights = malloc_(cntxt->num_valid_ids * sizeof(double));
	for (i=0,j=0; i < cntxt->num_entities; i++) {
		if (cntxt->used_indices[i]) {
			used_weights[j++] = cntxt->weights[i];
		}
	}
	sddlsum = SADDLE_SUM_init(used_weights, cntxt->num_valid_ids);
	if (sddlsum == NULL) {
		fprintf(stderr, "Could not allocate saddlesum context.\n");
		exit(EXIT_FAILURE);
	}
	free(used_weights);

	mapping_db->current_term = term_index;
	(void) mapping_db->get_next_mapping(mapping_db, &term_index, &hits, &num_hits);
        score = 0.0;
        for (num_used_hits=0, end_hits=hits+num_hits; hits < end_hits; hits++) {
                score += cntxt->weights[*hits];
                num_used_hits += cntxt->used_indices[*hits];
        }
        if (num_used_hits >= cntxt->min_term_size) {
                Pvalue = SADDLE_SUM_pvalue(sddlsum, score, num_used_hits,
                                           cntxt->Pvalue_cutoff,
                                           SADDLESUM_MAX_ITERS,
                                           SADDLESUM_TOLERANCE);
        }
        EnrichResults_insert_term_hit(cntxt, term_index, score, num_used_hits, Pvalue);
	SADDLE_SUM_del(sddlsum);
}


static
void EnrichResults_hgem_pvalues(EnrichContext *cntxt, CVTermDb *term_db,
				       TermMappingDb *mapping_db)
{
	uint32_t term_index;
	uint32_t *hits;
	uint32_t *end_hits;
	uint32_t num_hits;
	uint32_t num_used_hits;

        double Pvalue;
        unsigned int i;
        HypergeomStats *hgeom;

	hgeom = HypergeomStats_init(cntxt->num_valid_ids, cntxt->num_nonzero_valid_ids);
	/* FISHER_EXACT - main loop */
	mapping_db->reset(mapping_db);
	while (mapping_db->get_next_mapping(mapping_db, &term_index, &hits, &num_hits)) {
		i = 0;
		for (num_used_hits=0, end_hits=hits+num_hits; hits < end_hits; hits++) {
			i += (cntxt->weights[*hits] > 0.0);
			num_used_hits += cntxt->used_indices[*hits];
		}
		if (num_used_hits < cntxt->min_term_size) {
			continue;
		}
		Pvalue = HypergeomStats_pvalue(hgeom, i, num_used_hits);

		if (Pvalue <= cntxt->Pvalue_cutoff) {
			EnrichResults_insert_term_hit(cntxt, term_index,
						      (double) i, num_used_hits, Pvalue);
		}
	}
	HypergeomStats_del(hgeom);
}


static
void EnrichResults_hgem_single_pvalue(EnrichContext *cntxt, CVTermDb *term_db,
                                      TermMappingDb *mapping_db, uint32_t term_index)
{
	uint32_t *hits;
	uint32_t *end_hits;
	uint32_t num_hits;
	uint32_t num_used_hits;
        double Pvalue = -1.0;
        unsigned int i;
        HypergeomStats *hgeom;

	hgeom = HypergeomStats_init(cntxt->num_valid_ids, cntxt->num_nonzero_valid_ids);
	mapping_db->current_term = term_index;
        (void) mapping_db->get_next_mapping(mapping_db, &term_index, &hits, &num_hits);
        i = 0;
        for (num_used_hits=0, end_hits=hits+num_hits; hits < end_hits; hits++) {
                i += (cntxt->weights[*hits] > 0.0);
                num_used_hits += cntxt->used_indices[*hits];
        }
        if (num_used_hits >= cntxt->min_term_size) {
                Pvalue = HypergeomStats_pvalue(hgeom, i, num_used_hits);
        }
        EnrichResults_insert_term_hit(cntxt, term_index, (double) i, num_used_hits, Pvalue);
	HypergeomStats_del(hgeom);
}


static
void EnrichResults_count_used_terms(EnrichContext *cntxt, CVTermDb *term_db,
                                    TermMappingDb *mapping_db)
{
	/* Get effective sample size and hence a Pvalue_cutoff */
	/* Here we need to do one more scan of term mappings to count the terms
	   that are relevant */

	uint32_t term_index;
	uint32_t *hits;
	uint32_t *end_hits;
	uint32_t num_hits;
	uint32_t num_used_hits;
	cntxt->num_terms = term_db->num_terms;
	mapping_db->reset(mapping_db);
	while (mapping_db->get_next_mapping(mapping_db, &term_index, &hits, &num_hits)) {
		for (num_used_hits=0, end_hits=hits+num_hits; hits < end_hits; hits++) {
			if (cntxt->used_indices[*hits]) {
				num_used_hits++;
			}
		}
		if (num_used_hits >= cntxt->min_term_size) {
			cntxt->num_used_terms++;
		}
	}
	if (cntxt->effective_db_size <= 0.0) {
		cntxt->effective_db_size = (double) cntxt->num_used_terms;
	}
        cntxt->Pvalue_cutoff = cntxt->Evalue_cutoff / cntxt->effective_db_size;
}


void EnrichResults_calc_pvalues(EnrichContext *cntxt, CVTermDb *term_db,
				TermMappingDb *mapping_db)
{
        TermHit *term_hits;
        TermHit *end_term_hits;

        EnrichResults_count_used_terms(cntxt, term_db, mapping_db);

	/* Main run */
        switch (cntxt->statistics_type) {
        case SADDLESUM:
		EnrichResults_wsum_pvalues(cntxt, term_db, mapping_db);
                break;
        case FISHER_EXACT:
		EnrichResults_hgem_pvalues(cntxt, term_db, mapping_db);
                break;
        }

	/* Update term_hit data */
        term_hits = cntxt->term_hits;
        end_term_hits = cntxt->term_hits + cntxt->num_term_hits;
        for (; term_hits < end_term_hits; term_hits++) {
                term_hits->term = term_db->get_term_from_index(term_db, term_hits->term_index);
                term_hits->Evalue = term_hits->Pvalue * cntxt->effective_db_size;
        }

	/* Sort term_hits */
        qsort(cntxt->term_hits, cntxt->num_term_hits, sizeof(TermHit), TermHit_compare);
}


void EnrichResults_calc_single_pvalue(EnrichContext *cntxt, CVTermDb *term_db,
                                      TermMappingDb *mapping_db, uint32_t term_index)
{
        EnrichResults_count_used_terms(cntxt, term_db, mapping_db);
        switch (cntxt->statistics_type) {
        case SADDLESUM:
		EnrichResults_wsum_single_pvalue(cntxt, term_db, mapping_db, term_index);
                break;
        case FISHER_EXACT:
		EnrichResults_hgem_single_pvalue(cntxt, term_db, mapping_db, term_index);
                break;
        }
	/* Update term_hit data */
        cntxt->term_hits->term = term_db->get_term_from_index(term_db, cntxt->term_hits->term_index);
        cntxt->term_hits->Evalue = cntxt->term_hits->Pvalue * cntxt->effective_db_size;
}
