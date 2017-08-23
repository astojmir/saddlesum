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
#include "stack.h"
#include "hashtable.h"
#include "hashtable_itr.h"
#include "enrich.h"
#include "termdb2entities.h"

#define DASHED_LINE "--------------------------------------------------" \
                    "---------------------------------------\n"



static
void EnrichResults_print_header(EnrichContext *cntxt, FILE *fp, const char *fmt,
                                const char *heading_fmt)
{
        PrintBuf *pbuf = PrintBuf_init(NULL);

        fprintf(fp, heading_fmt, "QUERY AND DATABASE SUMMARY");

        if (cntxt->db_name != NULL) {
                fprintf(fp, fmt, "Database name", cntxt->db_name);
        }

        PrintBuf_printf(pbuf, 0, "%d", cntxt->num_terms);
        fprintf(fp, fmt, "Total database terms", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%d", cntxt->num_entities);
        fprintf(fp, fmt, "Total database entities", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%d", cntxt->num_raw_weights);
        fprintf(fp, fmt, "Submitted weights", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%d", cntxt->num_valid_ids);
        fprintf(fp, fmt, "Valid submitted entity ids", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%d", cntxt->min_term_size);
        fprintf(fp, fmt, "Minimum term size (weighted entities per term)", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%d", cntxt->num_used_terms);
        fprintf(fp, fmt, "Used database terms", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%d", cntxt->num_nonzero_valid_ids);
        fprintf(fp, fmt, "Non-zero weight entities", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%d", cntxt->num_unknown_ids);
        fprintf(fp, fmt, "Unknown submitted entity ids", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%d", cntxt->num_duplicate_ids);
        fprintf(fp, fmt, "Duplicate submitted entity ids", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%d", cntxt->num_conflicting_ids);
        fprintf(fp, fmt, "Unresolvable (ignored) conflicting entity ids", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%d", cntxt->num_resolvable_ids);
        fprintf(fp, fmt, "Resolvable (accepted) conflicting entity ids", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%d", cntxt->num_unused_entities);
        fprintf(fp, fmt, "Entities without submitted weight", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%8.2e", cntxt->Evalue_cutoff);
        fprintf(fp, fmt, "E-value cutoff", pbuf->buf);

        PrintBuf_printf(pbuf, 0, "%8.2e", cntxt->effective_db_size);
        fprintf(fp, fmt, "Effective database size", pbuf->buf);

        switch (cntxt->statistics_type) {
        case SADDLESUM:
                PrintBuf_printf(pbuf, 0, "%s", "Lugananni-Rice (sum of weights)");
                break;
        case FISHER_EXACT:
                PrintBuf_printf(pbuf, 0, "%s", "One-sided Fisher's Exact test");
                break;
        }
        fprintf(fp, fmt, "Statistics", pbuf->buf);

	switch (cntxt->cutoff_type) {
	case NONE:
                PrintBuf_printf(pbuf, 0, "%s", "All");
		fprintf(fp, fmt, "Top-ranked weights selected", pbuf->buf);
                PrintBuf_printf(pbuf, 0, "%s", "N/A");
		fprintf(fp, fmt, "Minimum weight selected", pbuf->buf);
		break;
	default:
                PrintBuf_printf(pbuf, 0, "%d", cntxt->rank_cutoff);
		fprintf(fp, fmt, "Top-ranked weights selected", pbuf->buf);
                PrintBuf_printf(pbuf, 0, "%.4f", cntxt->weight_cutoff);
		fprintf(fp, fmt, "Minimum weight selected", pbuf->buf);
		break;
	}

        if (cntxt->discretized_weights) {
                PrintBuf_printf(pbuf, 0, "%s", "Yes");
        }
        else {
                PrintBuf_printf(pbuf, 0, "%s", "No");
        }
        fprintf(fp, fmt, "Discretized weights", pbuf->buf);
        PrintBuf_delete(pbuf);
}


static
void EnrichResults_print_warnings(EnrichContext *cntxt, FILE *fp, const char *fmt,
                                  const char *heading_fmt)
{
        EntityWarning *warning = cntxt->first_warning;

        fprintf(fp, heading_fmt, "NOMENCLATURE WARNINGS");

        while (warning != NULL) {
                if (warning->code != UNKNOWN_ID) {
                        fprintf(fp, fmt, warning->msg);
                }
                warning = warning->next;
        }
}


static
void EnrichResults_print_unknown_ids(EnrichContext *cntxt, FILE *fp, const char *sep,
                                     const char *heading_fmt, uint16_t break_after)
{
        EntityWarning *warning = cntxt->first_warning;
        PrintBuf *pbuf = PrintBuf_init(NULL);
        uint16_t line_chars = 0;
        uint16_t seplen = strlen(sep);
        uint16_t msglen;

        fprintf(fp, heading_fmt, "UNKNOWN IDS");
        while (warning != NULL) {
                if (warning->code == UNKNOWN_ID) {
                        msglen = strlen(warning->msg);
                        if (line_chars + msglen + seplen > break_after) {
                                PrintBuf_printf(pbuf, -1, "\n");
                                line_chars = 0;
                        }
                        PrintBuf_printf(pbuf, -1, "%s%s", warning->msg, sep);
                        line_chars += msglen + seplen;
                }
                warning = warning->next;
        }
        if (pbuf->len) { /* only print if there is something to print */
                /* Delete last separator and add newline*/
                if (pbuf->len > seplen) {
                        pbuf->len -= seplen;
                }
                PrintBuf_printf(pbuf, -1, "\n");
                fprintf(fp, "%s", pbuf->buf);
        }
        PrintBuf_delete(pbuf);
}


static
void EnrichResults_print_term_hits(EnrichContext *cntxt, FILE *fp, const char *fmt,
                                   const char *heading_fmt, uint8_t print_table_headings)
{
        PrintBuf *pbuf1 = PrintBuf_init(NULL);
        PrintBuf *pbuf2 = PrintBuf_init(NULL);
        PrintBuf *pbuf3 = PrintBuf_init(NULL);
        PrintBuf *pbuf4 = PrintBuf_init(NULL);
        const char *old_namespace = NULL;
        TermHit *term_hits;
        TermHit *end_term_hits;

        term_hits = cntxt->term_hits;
        end_term_hits = cntxt->term_hits + cntxt->num_term_hits;
        for (; term_hits < end_term_hits; term_hits++) {

                if (old_namespace != term_hits->term->namespace) {
                        if (old_namespace != NULL) {
                                fprintf(fp, DASHED_LINE);
                        }
                        fprintf(fp, heading_fmt, term_hits->term->namespace);
                        if (print_table_headings) {
                                fprintf(fp, fmt, "Term ID", "Name", "Associations",
                                        "Score", "E-value", "URL");
                                fprintf(fp, DASHED_LINE);
                        }
                }
                old_namespace = term_hits->term->namespace;

                PrintBuf_printf(pbuf1, 0, "%d",  term_hits->num_entities);
                PrintBuf_printf(pbuf2, 0, "%.4f",  term_hits->score);
                PrintBuf_printf(pbuf3, 0, "%.2e",  term_hits->Evalue);

                if (term_hits->term->print_url != NULL) {
                        term_hits->term->print_url(term_hits->term, pbuf4);
                }
                else {
                        PrintBuf_printf(pbuf4, 0, "");
                }
                fprintf(fp, fmt, term_hits->term->term_id, term_hits->term->description,
                        pbuf1->buf, pbuf2->buf, pbuf3->buf, pbuf4->buf);
        }
        if (cntxt->num_term_hits) {
                fprintf(fp, DASHED_LINE);
        }
        PrintBuf_delete(pbuf1);
        PrintBuf_delete(pbuf2);
        PrintBuf_delete(pbuf3);
        PrintBuf_delete(pbuf4);
}


static
unsigned int hash_from_pCVTerm(void *ptr)
{
        CVTerm **t = (CVTerm **) ptr;
        return hash_from_string((*t)->term_id);
}


static
int pCVTerm_equal(void *p1, void *p2)
{
        CVTerm **t1 = (CVTerm **) p1;
        CVTerm **t2 = (CVTerm **) p2;
        return str_equal((*t1)->term_id, (*t2)->term_id);
}


static inline
void insert_term_(struct hashtable *visited, CVTerm *term, intptr_t level)
{
        CVTerm **pterm = malloc_(sizeof(CVTerm *));
        *pterm = term;
        if (! hashtable_insert(visited, pterm, (void *) level)) {
                fprintf(stderr, "Could not insert a term into 'visited'.\n");
                exit(EXIT_FAILURE);
        }
}

static
void print_input_symbols(EnrichContext *cntxt, PrintBuf *pbuf, CVTermDb *term_db,
                         TermMappingDb *mapping_db, CVTerm *term)
{
        uint32_t term_index;
        uint32_t *hits;
        uint32_t num_hits;
        int i;

        term_index = term_db->get_index_from_term_id(term_db, term->term_id);
        mapping_db->current_term = term_index;
        (void) mapping_db->get_next_mapping(mapping_db, &term_index, &hits, &num_hits);

        PrintBuf_printf(pbuf, 0, "");
        if (num_hits == 0) {
                return;
        }

        for (i=0; i < num_hits; i++) {
                if ( cntxt->used_indices[hits[i]] ) {
                        PrintBuf_printf(pbuf, -1, "%s,", cntxt->input_symbols[hits[i]]);
                }
        }
        pbuf->len--;
        PrintBuf_printf(pbuf, -1, "");
}



static
void EnrichResults_print_term_relationship_graph(EnrichContext *cntxt, FILE *fp,
                                                 const char *edge_fmt,
                                                 const char *node_fmt,
                                                 const char *heading_fmt,
                                                 CVTermDb *term_db,
                                                 TermMappingDb *mapping_db)
{
        struct hashtable *visited;
        struct hashtable_itr *itr;
        Stack *stack = Stack_init();
        PrintBuf *pbuf = PrintBuf_init(NULL);
        PrintBuf *pbuf2 = PrintBuf_init(NULL);
        CVTerm **pterm;
        intptr_t level;

        int i;
        uint32_t num_parents;
        CVTerm **parents;
        char **edgetypes;
        CVTerm *term;

        TermHit *term_hits;
        TermHit *end_term_hits;

        visited = create_hashtable(cntxt->num_terms, hash_from_pCVTerm, pCVTerm_equal);
        fprintf(fp, heading_fmt, "TERM RELATIONSHIPS");
        /* Transitive closure using stack */

        /* Insert significant terms into visited */
        term_hits = cntxt->term_hits;
        end_term_hits = cntxt->term_hits + cntxt->num_term_hits;
        for (; term_hits < end_term_hits; term_hits++) {
                if ((term_hits->term->get_parents == NULL) ||   \
                    ((num_parents = term_hits->term->get_parents(term_hits->term,
                                                                 &parents,
                                                                 &edgetypes)) == 0)) {
                        /* Do not insert terms not related to any other term */
                        continue;
                }
                insert_term_(visited, term_hits->term, floor(-log10(term_hits->Evalue)));
                Stack_push(stack, term_hits->term);
        }

        /* Depth first traversal of relationship graph */
        while (!Stack_is_empty(stack)) {
                term = (CVTerm *) Stack_pop(stack);
                if (term->get_parents == NULL) {
                        continue;
                }
                num_parents = term->get_parents(term, &parents, &edgetypes);
                for (i=0; i < num_parents; i++, parents++, edgetypes++) {
                        fprintf(fp, edge_fmt,term->term_id, *edgetypes,
                                (*parents)->term_id);
                        if (hashtable_search(visited, parents) == NULL) {
                                Stack_push(stack, *parents);
                                insert_term_(visited, *parents, -1);
                        }
                }
        }

        fprintf(fp, heading_fmt, "NODE PROPERTIES");
        if (hashtable_count(visited) > 0) {
                itr = hashtable_iterator(visited);
                do {
                        pterm = (CVTerm **) hashtable_iterator_key(itr);
                        level = (intptr_t) hashtable_iterator_value(itr);

                        if (level != -1) {
                                print_input_symbols(cntxt, pbuf2, term_db, mapping_db, *pterm);
                        }
                        else {
                                PrintBuf_printf(pbuf2, 0, "");
                        }

                        level = (level == -1) ? 0 : level;
                        if ( (*pterm)->print_url != NULL) {
                                (*pterm)->print_url(*pterm, pbuf);
                        }
                        else {
                                PrintBuf_printf(pbuf, 0, "");
                        }
                        fprintf(fp, node_fmt, (*pterm)->term_id,
                                (*pterm)->description, (int) level, pbuf->buf, pbuf2->buf);
                } while (hashtable_iterator_advance(itr));
        }
        hashtable_destroy(visited, 0);
        Stack_delete(stack);
        PrintBuf_delete(pbuf);
        PrintBuf_delete(pbuf2);
}


void EnrichResults_print_all_text(EnrichContext *cntxt, FILE *fp,
                                  uint8_t print_warnings,
                                  uint8_t print_unknown_ids)
{
        fprintf(fp, DASHED_LINE);
        fprintf(fp, "               SADDLESUM RESULTS\n");
        fprintf(fp, DASHED_LINE);

        EnrichResults_print_header(cntxt, fp, "%-48.48s %s\n", "\n**** %s ****\n");

        if (print_warnings) {
                EnrichResults_print_warnings(cntxt, fp, "%s\n",  "\n**** %s ****\n");
        }
        if (print_unknown_ids) {
                EnrichResults_print_unknown_ids(cntxt, fp, ", ", "\n**** %s ****\n", 80);
        }

        EnrichResults_print_term_hits(cntxt, fp, "%-15.15s %-40.40s %6.6s %12.12s %10.10s%.0s\n",
                                      "\n**** %s ****\n", 1);
        fprintf(fp, DASHED_LINE);
}


void EnrichResults_print_all_tabsep(EnrichContext *cntxt, FILE *fp, CVTermDb *term_db,
                                    TermMappingDb *mapping_db)
{
        EnrichResults_print_header(cntxt, fp, "%s\t%s\n", "#\n# %s\n#\n");
        EnrichResults_print_warnings(cntxt, fp, "%s\n", "#\n# %s\n#\n");
        EnrichResults_print_unknown_ids(cntxt, fp, ",", "#\n# %s\n#\n", 30000);
        EnrichResults_print_term_relationship_graph(cntxt, fp,
                                                    "%s\t%s\t%s\n",
                                                    "%s\t%s\t%d\t%s\t%s\n",
                                                    "#\n# %s\n#\n",
                                                    term_db,
                                                    mapping_db);
        EnrichResults_print_term_hits(cntxt, fp, "%s\t%s\t%s\t%s\t%s\t%s\n",
                                      "#\n# %s\n#\n", 0);
}


typedef struct _EntityWeight_s {
        Entity *entity;
        double weight;
        uint8_t used_flag;
} EntityWeight;


static int EntityWeight_compare (const void * M1, const void * M2)
{
        const EntityWeight *ew1 = (const EntityWeight *) M1;
        const EntityWeight *ew2 = (const EntityWeight *) M2;
        register int retval;

        /* Compare by weight, then by ID */
        retval = ew2->used_flag - ew1->used_flag;
        if (!retval) {
                retval = (ew1->weight < ew2->weight) - (ew1->weight > ew2->weight);
        }
        if (!retval) {
                retval = strcmp(ew1->entity->symbol, ew2->entity->symbol);
        }
        return retval;
}


static
void EnrichResults_print_term(EnrichContext *cntxt, FILE *fp, EntityDb *entity_db,
                              TermMappingDb *mapping_db, const char *fmt,
                              const char *summary_fmt, const char *heading_fmt,
                              uint8_t print_table_headings)
{
        CVTerm *term;
        PrintBuf *pbuf1;
        PrintBuf *pbuf2;
        PrintBuf *pbuf3;
        uint32_t term_index;
        uint32_t *hits;
        uint32_t num_hits;
        uint32_t num_used_hits = 0;
        uint32_t num_positive_hits = 0;
        uint32_t num_negative_hits = 0;
        uint32_t num_zero_hits = 0;
        size_t i;
        EntityWeight *ranked;
        Entity *entity;
        const char *desc;

        if (!cntxt->num_term_hits) {
                return;
        }
        term = cntxt->term_hits->term;
        pbuf1 = PrintBuf_init(NULL);
        pbuf2 = PrintBuf_init(NULL);
        pbuf3 = PrintBuf_init(NULL);

        mapping_db->current_term = cntxt->term_hits->term_index;
        (void) mapping_db->get_next_mapping(mapping_db, &term_index,
                                              &hits, &num_hits);
        ranked = malloc_(num_hits * sizeof(EntityWeight));
        for (i=0; i < num_hits; i++) {
                ranked[i].entity = entity_db->get_entity_from_index(entity_db, hits[i]);
                ranked[i].weight = cntxt->weights[hits[i]];
                ranked[i].used_flag = cntxt->used_indices[hits[i]];
                if (ranked[i].used_flag) {
                        num_zero_hits += (cntxt->weights[hits[i]] == 0.0);
                        num_used_hits += 1;
                        num_positive_hits += (cntxt->weights[hits[i]] > 0.0);
                        num_negative_hits += (cntxt->weights[hits[i]] < 0.0);
                }
        }
        qsort(ranked, num_hits, sizeof(EntityWeight), EntityWeight_compare);

        fprintf(fp, heading_fmt, "TERM SUMMARY");
        if (cntxt->db_name != NULL) {
                fprintf(fp, summary_fmt, "Database name", cntxt->db_name);
        }
        fprintf(fp, summary_fmt, "Term ID", term->term_id);

        fprintf(fp, summary_fmt, "Term definition", term->description);

        fprintf(fp, summary_fmt, "Namespace", term->namespace);

        PrintBuf_printf(pbuf1, 0, "%d", num_hits);
        fprintf(fp, summary_fmt, "Total associations", pbuf1->buf);

        PrintBuf_printf(pbuf1, 0, "%d", num_used_hits);
        fprintf(fp, summary_fmt, "Weighted associations", pbuf1->buf);

        PrintBuf_printf(pbuf1, 0, "%d", num_positive_hits);
        fprintf(fp, summary_fmt, "Positive-weighted associations", pbuf1->buf);

        PrintBuf_printf(pbuf1, 0, "%d", num_negative_hits);
        fprintf(fp, summary_fmt, "Negative-weighted associations", pbuf1->buf);

        PrintBuf_printf(pbuf1, 0, "%d", num_zero_hits);
        fprintf(fp, summary_fmt, "Zero-weighted associations", pbuf1->buf);

        PrintBuf_printf(pbuf1, 0, "%.4f", cntxt->term_hits->score);
        fprintf(fp, summary_fmt, "Score", pbuf1->buf);

        if (cntxt->term_hits->Pvalue >= 0.0) {
                PrintBuf_printf(pbuf1, 0, "%.4e", cntxt->term_hits->Pvalue);
        }
        else {
                PrintBuf_printf(pbuf1, 0, "N/A");
        }
        fprintf(fp, summary_fmt, "P-value", pbuf1->buf);

        if (cntxt->term_hits->Evalue >= 0.0) {
                PrintBuf_printf(pbuf1, 0, "%.4e", cntxt->term_hits->Evalue);
        }
        else {
                PrintBuf_printf(pbuf1, 0, "N/A");
        }
        fprintf(fp, summary_fmt, "E-value", pbuf1->buf);

        fprintf(fp, heading_fmt, "ENTITIES");
        if (print_table_headings) {
                fprintf(fp, fmt, "Rank", "Identifier", "Description", "Weight", "Links");
                fprintf(fp, DASHED_LINE);
        }
        for (i=0; i < num_hits; i++) {
                entity = ranked[i].entity;
                desc = entity->description == NULL ? "" : entity->description;
                PrintBuf_printf(pbuf3, 0, "%d.", (int) i+1);
                if (ranked[i].used_flag) {
                        PrintBuf_printf(pbuf1, 0, "%.4f", ranked[i].weight);
                }
                else {
                        PrintBuf_printf(pbuf1, 0, "None");
                }
                if (entity->print_url != NULL) {
                        entity->print_url(entity, pbuf2);
                }
                else {
                        PrintBuf_printf(pbuf2, 0, "");
                }
                fprintf(fp, fmt, pbuf3->buf, entity->symbol, desc,
                        pbuf1->buf, pbuf2->buf);
        }
        free(ranked);
        PrintBuf_delete(pbuf1);
        PrintBuf_delete(pbuf2);
        PrintBuf_delete(pbuf3);
}


void EnrichResults_print_term_text(EnrichContext *cntxt, FILE *fp, EntityDb *entity_db,
                                   TermMappingDb *mapping_db)
{
        EnrichResults_print_term(cntxt, fp,entity_db, mapping_db,
                                 "%6.6s %-15.15s %-40.40s %12.12s %.0s\n",
                                 "%-48.48s %s\n", "\n**** %s ****\n", 0);
}


void EnrichResults_print_term_tabsep(EnrichContext *cntxt, FILE *fp, EntityDb *entity_db,
                                     TermMappingDb *mapping_db)
{
        EnrichResults_print_term(cntxt, fp,entity_db, mapping_db,
                                 "%s\t%s\t%s\t%s\t%s\n", "%s\t%s\n",
                                 "#\n# %s\n#\n", 0);
}
