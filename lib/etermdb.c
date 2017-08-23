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
#include "miscutils.h"
#include "hashtable.h"
#include "cvterm.h"
#include "enrich.h"

#define EXTERMDB_MAGIC 1644632861U
#define KEGGTERMDB_MAGIC 2264738403U
#define GOTERMDB_MAGIC 2187050528U
#define EXTERMDB_BUF_INCR 2048


/* Term database */
static
void ETDTermDb_delete(CVTermDb *term_db_)
{
	int i;
	ETDTermDb *term_db = (ETDTermDb *) term_db_;
	for (i=0; i < term_db->num_terms; i++) {
		term_db->terms[i]->delete(term_db->terms[i]);
	}
	free(term_db->terms);
	term_db->terms = NULL;
	term_db->num_terms = 0;
	term_db->max_terms = 0;
        free(term_db->db_name);
        free(term_db->namespaces);
        free(term_db->namespace_buf);
	for (i=0; i < term_db->num_namespaces; i++) {
                free(term_db->termid_bufs[i]);
                free(term_db->desc_bufs[i]);
                free(term_db->edgetype_bufs[i]);
                free(term_db->metadata_bufs[i]);
	}
        free(term_db->termid_bufs);
        free(term_db->desc_bufs);
        free(term_db->edgetype_bufs);
        free(term_db->metadata_bufs);
	hashtable_destroy(term_db->termid2index, 0);
	free(term_db);
}


static
CVTerm *ETDTermDb_get_term_from_index(CVTermDb *term_db_, uint32_t i)
{
	ETDTermDb *term_db = (ETDTermDb *) term_db_;
	if (i < term_db->num_terms) {
		return term_db->terms[i];
	}
	return NULL;
}


static
int ETDTermDb_get_index_from_term_id(CVTermDb *term_db_, char *term_id)
{
	ETDTermDb *term_db = (ETDTermDb *) term_db_;
	void *indx = hashtable_search(term_db->termid2index, term_id);
	if (indx == NULL) {
                return -1;
        }
	return  (int) ((intptr_t) indx) - 1;
}


static
CVTerm *ETDTermDb_get_term_from_term_id(CVTermDb *term_db_, char *term_id)
{
	ETDTermDb *term_db = (ETDTermDb *) term_db_;
        int i = ETDTermDb_get_index_from_term_id(term_db_, (char *) term_id);
        CVTerm *term;
        if (i > 0) {
                term = *(term_db->terms + i);
                return term;
        }
        return NULL;
}


static
int ETDTermDb_insert_term(CVTermDb *term_db_, const char *term_id,
                          const char *namespace, const char *description)
{
	ETDTermDb *term_db = (ETDTermDb *) term_db_;
	CVTerm *term = ETDTermDb_get_term_from_term_id(term_db_, (char *) term_id);

	if (term != NULL) {
		return 0; /* term_id already exists */
	}

	if (term_db->num_terms  >= term_db->max_terms) {
		term_db->max_terms += EXTERMDB_BUF_INCR;
		term_db->terms = realloc_(term_db->terms,
					  term_db->max_terms*sizeof(CVTerm *));
	}
        term = CVTerm_init(NULL);
	term->term_id = strdup_(term_id);
	term->description = strdup_(description);
	term->namespace = namespace;
        term_db->terms[term_db->num_terms++] = term;

	if (! hashtable_insert(term_db->termid2index, strdup_(term_id),
                               (void *) ((intptr_t) term_db->num_terms + 1))) {
		fprintf(stderr, "Could not insert term %s.\n", term_id);
		exit(EXIT_FAILURE);
	}
	return 1;
}


static
FILE *ETDTermDb_open(const char *etd_filename)
{
        char header[8];
        uint32_t magic;
        FILE *fp = fopen(etd_filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "Could not open file %s.\n", etd_filename);
		exit(EXIT_FAILURE);
	}
        fread_(header, 1, 8, fp);
        fread_uint32(&magic, 1, fp);
        if (memcmp(header, "EXTERMDB", 8) || (magic != EXTERMDB_MAGIC)) {
                fprintf(stderr, "%s is not a valid ETD format database.\n",
                        etd_filename);
                exit(EXIT_FAILURE);
        }
        return fp;
}


static
ETDTermDb *ETDTermDb_init(FILE *fp)
{
	ETDTermDb *term_db;
        uint32_t bufsize;

        term_db = calloc_(1, sizeof(ETDTermDb));
	term_db->delete = ETDTermDb_delete;
	term_db->get_term_from_index = ETDTermDb_get_term_from_index;
	term_db->get_index_from_term_id = ETDTermDb_get_index_from_term_id;
	term_db->get_term_from_term_id = ETDTermDb_get_term_from_term_id;
        term_db->insert_term = ETDTermDb_insert_term;
	term_db->termid2index = create_hashtable(INIT_MAX_TERMS, hash_from_string, str_equal);

        term_db->db_name = fread_buf(fp, &bufsize);
        fread_uint32(&term_db->num_namespaces, 1, fp);
        term_db->namespaces = malloc_(term_db->num_namespaces * sizeof(char *));
        term_db->namespace_buf = fread_buf(fp, &bufsize);
        strbuf2array(term_db->namespaces, term_db->namespace_buf,
                     term_db->num_namespaces);

        term_db->num_ns_terms = calloc_(1, term_db->num_namespaces * sizeof(uint32_t));
        term_db->termid_bufs = calloc_(1, term_db->num_namespaces * sizeof(char *));
        term_db->desc_bufs = calloc_(1, term_db->num_namespaces * sizeof(char *));
        term_db->edgetype_bufs = calloc_(1, term_db->num_namespaces * sizeof(char *));
        term_db->metadata_bufs = calloc_(1, term_db->num_namespaces * sizeof(char *));

	return term_db;
}


static
uint32_t ETDTermDb_read_ns_header(FILE *fp)
{
        char header[8];
        uint32_t magic;
        fread_(header, 1, 8, fp);
        fread_uint32(&magic, 1, fp);
        if (memcmp(header, "TERMDBNS", 8)
            || ((magic != KEGGTERMDB_MAGIC)
                && (magic != GOTERMDB_MAGIC)) ) {
                fprintf(stderr, "Invalid ETD format database.\n");
                exit(EXIT_FAILURE);
        }
        return magic;
}


static
int ETDTermDb_read_ns_data(FILE *fp, ETDTermDb *term_db, TermMappingDb *mapping_db,
                           CVTerm *(*term_init_func) (CVTerm *), int cns)
{
        int term_offset = term_db->num_terms;
        uint32_t num_edgetypes;
        char **tmp_edgetypes;
        uint32_t bufsize;
        uint32_t M;
        uint32_t *tmp_counts;
        uint32_t tmp;
        int i;
        int j;
        char *cur;
        ETDTerm *term;

        fread_uint32(&num_edgetypes, 1, fp);
        term_db->edgetype_bufs[cns] = fread_buf(fp, &bufsize);
        tmp_edgetypes = malloc_(num_edgetypes * sizeof(char *));
        strbuf2array(tmp_edgetypes, term_db->edgetype_bufs[cns], num_edgetypes);

        fread_uint32(&M, 1, fp);
        tmp_counts = malloc_(M * sizeof(uint32_t));
        term_db->num_ns_terms[cns] = M;

	if (M + term_db->num_terms  >= term_db->max_terms) {
		term_db->max_terms += M;
		term_db->terms = realloc_(term_db->terms,
					  term_db->max_terms*sizeof(CVTerm *));
	}
        term_db->num_terms += M;
        for (i=term_offset; i < term_db->num_terms; i++) {
                term_db->terms[i] = term_init_func(NULL);
                term = (ETDTerm *) term_db->terms[i];
                term->namespace = term_db->namespaces[cns];
                fread_uint32(&term->flag, 1, fp);
        }

        /* If this is too slow, we can later replace it with inlined code for bulk insertion */
        fread_uint32(tmp_counts, M, fp);
        for (i=0; i < M; i++) {
                mapping_db->insert_new_mapping(mapping_db);
                for (j=0; j < tmp_counts[i]; j++) {
                        fread_uint32(&tmp, 1, fp);
                        mapping_db->insert_hit(mapping_db, tmp);
                }
        }
        free(tmp_counts);

        term_db->termid_bufs[cns] = fread_buf(fp, &bufsize);
        cur = term_db->termid_bufs[cns];
        for (i=term_offset; i < term_db->num_terms; i++) {
                term_db->terms[i]->term_id = cur;
                cur += (strlen(cur) + 1);
                if (! hashtable_insert(term_db->termid2index,
                                       strdup_(term_db->terms[i]->term_id),
                                       (void *) ((intptr_t) i + 1))) {
                        fprintf(stderr, "Could not insert term %s.\n",
                                term_db->terms[i]->term_id);
                        exit(EXIT_FAILURE);
                }
        }

        term_db->desc_bufs[cns] = fread_buf(fp, &bufsize);
        cur = term_db->desc_bufs[cns];
        for (i=term_offset; i < term_db->num_terms; i++) {
                term_db->terms[i]->description = cur;
                cur += (strlen(cur) + 1);
        }

        for (i=term_offset; i < term_db->num_terms; i++) {
                term = (ETDTerm *) term_db->terms[i];
                fread_uint32(&tmp, 1, fp);
                term->num_parents = tmp;
                if (tmp > 0) {
                        term->parents = malloc_(tmp*sizeof(CVTerm *));
                        term->edgetypes = malloc_(tmp*sizeof(const char *));
                }
        }

        for (i=term_offset; i < term_db->num_terms; i++) {
                term = (ETDTerm *) term_db->terms[i];
                for (j=0; j < term->num_parents; j++) {
                        fread_uint32(&tmp, 1, fp);
                        term->parents[j] = term_db->terms[term_offset + tmp];
                }
        }

        for (i=term_offset; i < term_db->num_terms; i++) {
                term = (ETDTerm *) term_db->terms[i];
                for (j=0; j < term->num_parents; j++) {
                        fread_uint32(&tmp, 1, fp);
                        term->edgetypes[j] = tmp_edgetypes[tmp];
                }
        }
        free(tmp_edgetypes);

        term_db->metadata_bufs[cns] = fread_buf(fp, &bufsize);
        return term_offset;
}


static
void ETDTermDb_skip_ns_data(FILE *fp)
{
        uint32_t tmp;
        uint32_t M;
        uint32_t n;
        int i;
        fseek(fp, sizeof(uint32_t), SEEK_CUR);
        fskip_buf(fp);
        fread_uint32(&M, 1, fp);
        fseek(fp, M * sizeof(uint32_t), SEEK_CUR);
        for (n=0, i=0; i < M; i++) {
                fread_uint32(&tmp, 1, fp);
                n += tmp;
        }
        fseek(fp, n * sizeof(uint32_t), SEEK_CUR);
        fskip_buf(fp);
        fskip_buf(fp);
        for (n=0, i=0; i < M; i++) {
                fread_uint32(&tmp, 1, fp);
                n += tmp;
        }
        fseek(fp, 2 * n * sizeof(uint32_t), SEEK_CUR);
        fskip_buf(fp);
}

static
int is_excluded(const char *namespace, const char **excluded_namespaces, int n)
{
        int i;
        for (i=0; i < n; i++) {
                if (!strcmp(namespace, excluded_namespaces[i])) {
                        return 1;
                }
        }
        return 0;
}

static
void ETDTermDb_update_KEGG_terms(ETDTermDb *term_db, int cns, int term_offset)
{
        int i;
        KEGGTerm *term;
        for (i=term_offset; i < term_db->num_terms; i++) {
                term = (KEGGTerm *) term_db->terms[i];
                term->org_prefix = term_db->metadata_bufs[cns];
        }
}


int ETD_enrichment_context(const char *etd_filename, EntityDb **entity_db_,
			   CVTermDb **term_db_, TermMappingDb **mapping_db_,
                           const char **excluded_namespaces, int num_excluded)
{
        FILE *fp;
        NCBIGenesDb *entity_db;
        ETDTermDb *term_db;
        TermMappingDb *mapping_db = *mapping_db_;
        uint32_t magic;
        int i;
        int x;  /* term_offset */

        fp = ETDTermDb_open(etd_filename);

        /* read ETD file header */
        term_db = ETDTermDb_init(fp);

        /* read NCBI Genes database */
        entity_db = NCBIGenesDb_init(fp);

        /* read namespaces */
	if (mapping_db == NULL) {
		mapping_db = TermMappingDb_init();
	}
        for (i=0; i < term_db->num_namespaces; i++) {
                magic = ETDTermDb_read_ns_header(fp);
                if (is_excluded(term_db->namespaces[i], excluded_namespaces,
                                num_excluded)) {
                        ETDTermDb_skip_ns_data(fp);
                }
                else if (magic == KEGGTERMDB_MAGIC) {
                        x = ETDTermDb_read_ns_data(fp, term_db, mapping_db,
                                                   KEGGTerm_init, i);
                        /* need to add org_prefix to terms */
                        ETDTermDb_update_KEGG_terms(term_db, i, x);
                }
                else {  /* GOTERMDB_MAGIC */
                        (void) ETDTermDb_read_ns_data(fp, term_db, mapping_db,
                                                      GOTerm_init, i);
                        /* No need for further processing */
                }
        }
	fclose(fp);

	*entity_db_ = (EntityDb *) entity_db;
	*term_db_ = (CVTermDb *) term_db;
	*mapping_db_ = mapping_db;

	return 1;
}


void ETDTermDb_print_info(const char *etd_filename, FILE *fp,
                          OutputType output_type)
{
        PrintBuf *pbuf1 = PrintBuf_init(NULL);
        PrintBuf *pbuf2 = PrintBuf_init(NULL);
        NCBIGenesDb *entity_db = NULL;
        ETDTermDb *term_db = NULL;
        TermMappingDb *mapping_db = NULL;
        int i;
        const char *fmt;
        const char *heading_fmt;

        switch (output_type) {
        case TEXT:
                fmt = "%-48.48s %s\n";
                heading_fmt = "\n**** %s ****\n";
                break;
        case TAB:
        default:
                fmt = "%s\t%s\n";
                heading_fmt = "#\n# %s\n#\n";
                break;
        }

        (void) ETD_enrichment_context(etd_filename,
                                      (EntityDb **) &entity_db,
                                      (CVTermDb **) &term_db,
                                      (TermMappingDb **) &mapping_db,
                                      NULL, 0);


        fprintf(fp, heading_fmt, "EXTENDED TERM DATABASE GENERAL INFO");

        fprintf(fp, fmt, "Database name", term_db->db_name);
        fprintf(fp, fmt, "Database file", etd_filename);
        fprintf(fp, fmt, "NCBI Gene info file", entity_db->gene_info_file);

        PrintBuf_printf(pbuf1, 0, "%d", entity_db->tax_id);
        fprintf(fp, fmt, "NCBI Taxonomy ID", pbuf1->buf);

        PrintBuf_printf(pbuf1, 0, "%d", entity_db->num_entities);
        fprintf(fp, fmt, "Total Gene entries", pbuf1->buf);

        PrintBuf_printf(pbuf1, 0, "%d", term_db->num_terms);
        fprintf(fp, fmt, "Total terms", pbuf1->buf);

        PrintBuf_printf(pbuf1, 0, "%d", term_db->num_namespaces);
        fprintf(fp, fmt, "Term namespaces", pbuf1->buf);

        fprintf(fp, heading_fmt, "NAMESPACES (TERM COUNTS)");
        for (i=0; i < term_db->num_namespaces; i++) {
                PrintBuf_printf(pbuf1, 0, "%d", term_db->num_ns_terms[i]);
                PrintBuf_printf(pbuf2, 0, "%s",
                                term_db->namespaces[i]);
                fprintf(fp, fmt, pbuf2->buf, pbuf1->buf);
        }
        PrintBuf_delete(pbuf1);
        PrintBuf_delete(pbuf2);
}


void ETDTermDb_print_namespaces(const char *etd_filename, FILE *fp_out,
                                OutputType output_type)
{
        FILE *fp_in;
        ETDTermDb *term_db;
        int i;

        fp_in = ETDTermDb_open(etd_filename);
        term_db = ETDTermDb_init(fp_in);
	fclose(fp_in);

        if (output_type == TEXT) {
                fprintf(fp_out, "Term database contains the following namespaces:\n");
        }
        for (i=0; i < term_db->num_namespaces; i++) {
                fprintf(fp_out, "%s\n", term_db->namespaces[i]);
        }
}
