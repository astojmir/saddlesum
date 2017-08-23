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

#ifndef _ENRICHCONTEXT_H
#define _ENRICHCONTEXT_H
#ifdef __cplusplus
extern "C" {
#endif

#include "entity.h"
#include "cvterm.h"
#include "termdb2entities.h"

#define INIT_MAX_ENTITIES 256
#define INIT_MAX_TERMS 64

typedef enum {TEXT, TAB} OutputType;

typedef enum {
        SADDLESUM,
        FISHER_EXACT
} EnrichStats;

typedef enum {
        NO_TRANSFORM,
	FLIP,
        ABS
} TransformType;

typedef enum {
        NONE,
	RANK,
        MIN_VALUE
} CutoffType;

typedef struct _TermHit_s {
        CVTerm *term;
        uint32_t term_index;
        double score;
        uint32_t num_entities;
        double Evalue;
        double Pvalue;
} TermHit;


typedef struct _EnrichContext_s {
        const char *db_name;
        uint32_t num_terms;
        uint32_t num_used_terms;
	uint32_t min_term_size;
        uint32_t num_entities;
        uint32_t num_raw_weights;
        uint32_t num_valid_ids;
        uint32_t num_nonzero_valid_ids;
        uint32_t num_unknown_ids;
        uint32_t num_duplicate_ids;
        uint32_t num_conflicting_ids;
        uint32_t num_resolvable_ids;
        uint32_t num_unused_entities;
        double Evalue_cutoff;
        double Pvalue_cutoff;
        double effective_db_size;
        EnrichStats statistics_type;
        uint8_t num_namespaces;
	TransformType transform_type;
        uint8_t discretized_weights;
	CutoffType cutoff_type;
        uint32_t rank_cutoff;
	double weight_cutoff;
        uint8_t use_all_weights;
        EntityWarning *first_warning;
        EntityWarning *last_warning;
        double *weights;
        uint8_t *used_indices;
        uint32_t num_term_hits;
        uint32_t max_term_hits;
        TermHit *term_hits;
        char **input_symbols;
} EnrichContext;


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
                                  uint8_t use_all_weights);

void EnrichContext_delete(EnrichContext *cntxt);

void EnrichResults_load_weights(EnrichContext *cntxt, const char *weights_filename,
				EntityDb *entity_db, TermMappingDb *mapping_db);

void EnrichResults_process_weights(EnrichContext *cntxt);

void EnrichResults_calc_pvalues(EnrichContext *cntxt, CVTermDb *term_db,
				TermMappingDb *mapping_db);

void EnrichResults_calc_single_pvalue(EnrichContext *cntxt, CVTermDb *term_db,
                                      TermMappingDb *mapping_db, uint32_t term_index);

void EnrichResults_print_all_text(EnrichContext *cntxt, FILE *fp,
                                  uint8_t print_warnings,
                                  uint8_t print_unknown_ids);

void EnrichResults_print_all_tabsep(EnrichContext *cntxt, FILE *fp, CVTermDb *term_db,
                                    TermMappingDb *mapping_db);



int GMT_enrichment_context(const char *gmt_filename, const char *namespace, EntityDb **entity_db,
                           CVTermDb **term_db, TermMappingDb **mapping_db);

int ETD_enrichment_context(const char *etd_filename, EntityDb **entity_db_,
			   CVTermDb **term_db_, TermMappingDb **mapping_db_,
                           const char **excluded_namespaces, int num_excluded);

void ETDTermDb_print_info(const char *etd_filename, FILE *fp, OutputType output_type);

void ETDTermDb_print_namespaces(const char *etd_filename, FILE *fp_out,
                                OutputType output_type);

void EnrichResults_print_term_text(EnrichContext *cntxt, FILE *fp, EntityDb *entity_db,
                                   TermMappingDb *mapping_db);

void EnrichResults_print_term_tabsep(EnrichContext *cntxt, FILE *fp, EntityDb *entity_db,
                                     TermMappingDb *mapping_db);


#ifdef __cplusplus
}
#endif
#endif /* !_ENRICHCONTEXT_H */
