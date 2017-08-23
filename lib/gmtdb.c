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
#include "fsfile.h"
#include "hashtable.h"
#include "enrich.h"


#ifndef GMT_MAX_FIELD_WIDTH
#define GMT_MAX_FIELD_WIDTH 512
#endif


/* Entity database */
typedef struct _GMTEntityDb_s {
	EntityDb_HEAD
        Entity *entities;
        uint32_t max_entities;
        struct hashtable *symbol2index;
} GMTEntityDb;

static
void GMTEntityDb_delete(EntityDb *entity_db_)
{
	int i;
	GMTEntityDb *entity_db = (GMTEntityDb *) entity_db_;
	for (i=0; i < entity_db->num_entities; i++) {
		entity_db->entities[i].delete(entity_db->entities + i);
	}
	free(entity_db->entities);
	entity_db->entities = NULL;
	entity_db->num_entities = 0;
	entity_db->max_entities = 0;
	hashtable_destroy(entity_db->symbol2index, 0);
	free(entity_db);
}

static
uint32_t GMTEntityDb_insert_item(EntityDb *entity_db_, const char *symbol,
                                 const char *description)
{
	GMTEntityDb *entity_db = (GMTEntityDb *) entity_db_;
        Entity *item;
	uint32_t i;
	void *indx = hashtable_search(entity_db->symbol2index, (char *) symbol);
	if (indx == NULL) {
		if (entity_db->num_entities >= entity_db->max_entities) {
			entity_db->max_entities *= 2;
			entity_db->entities = realloc_(entity_db->entities,
						       entity_db->max_entities*sizeof(Entity));
		}
		i = entity_db->num_entities;
		item = entity_db->entities + entity_db->num_entities;
		entity_db->num_entities++;
		item->delete = Entity_delete_data;
                item->print_url = NULL;
		item->symbol = strdup_(symbol);
                item->description = NULL;

		if (! hashtable_insert(entity_db->symbol2index, strdup_(symbol),
                                       (void *) ((intptr_t) i)+1) ) {
			fprintf(stderr, "Could not insert entity item %s.\n", symbol);
			exit(EXIT_FAILURE);
		}
	}
	else {
		i = (uint32_t) ((size_t) indx) - 1;
	}
	return (uint32_t) i;
}

static
EntityWarning *GMTEntityDb_map_symbol(EntityDb *entity_db_, char *symbol, uint32_t *i)
{
	PrintBuf pbuf;
	EntityWarning *warning;
	GMTEntityDb *entity_db = (GMTEntityDb *) entity_db_;
	void *indx = hashtable_search(entity_db->symbol2index, symbol);

	if (indx == NULL) {
                /* Initialize PrintBuf on the stack. EntityWarning instance
                   will take ownership of the actual buffer so PrintBuf_delete
                   should not be called. */
                (void) PrintBuf_init(&pbuf);

		warning = calloc_(1, sizeof(EntityWarning));
		PrintBuf_printf(&pbuf, 0, "%s", symbol);
		warning->msg = pbuf.buf;
                warning->code = UNKNOWN_ID;
		return warning;
	}
	*i = (uint32_t) ((size_t) indx) - 1;
	return NULL;
}

static
Entity * GMTEntityDb_get_entity_from_index(EntityDb *entity_db_, uint32_t i)
{
	GMTEntityDb *entity_db = (GMTEntityDb *) entity_db_;

	if (i < entity_db->num_entities) {
		return entity_db->entities + i;
	}
	return NULL;
}

static
GMTEntityDb *GMTEntityDb_init(void)
{
	GMTEntityDb *entity_db = calloc_(1, sizeof(GMTEntityDb));

	entity_db->delete = GMTEntityDb_delete;
	entity_db->map_symbol = GMTEntityDb_map_symbol;
	entity_db->get_entity_from_index = GMTEntityDb_get_entity_from_index;
        entity_db->insert_item = GMTEntityDb_insert_item;
	entity_db->entities = calloc_(INIT_MAX_ENTITIES, sizeof(Entity));
	entity_db->num_entities = 0;
	entity_db->max_entities = INIT_MAX_ENTITIES;
	entity_db->symbol2index = create_hashtable(INIT_MAX_ENTITIES, hash_from_string, str_equal);
	return entity_db;
}

/* Term database */
typedef struct _GMTTermDb_s {
	CVTermDb_HEAD
        CVTerm *terms;
        uint32_t max_terms;
        struct hashtable *termid2index;
} GMTTermDb;

static
void GMTTermDb_delete(CVTermDb *term_db_)
{
	int i;
	GMTTermDb *term_db = (GMTTermDb *) term_db_;
	for (i=0; i < term_db->num_terms; i++) {
		term_db->terms[i].delete(term_db->terms + i);
	}
	free(term_db->terms);
	term_db->terms = NULL;
	term_db->num_terms = 0;
	term_db->max_terms = 0;
	hashtable_destroy(term_db->termid2index, 0);
	free(term_db);
}

static
CVTerm *GMTTermDb_get_term_from_index(CVTermDb *term_db_, uint32_t i)
{
	GMTTermDb *term_db = (GMTTermDb *) term_db_;

	if (i < term_db->num_terms) {
		return term_db->terms + i;
	}
	return NULL;
}

static
int GMTTermDb_get_index_from_term_id(CVTermDb *term_db_, char *term_id)
{
	GMTTermDb *term_db = (GMTTermDb *) term_db_;
	void *indx = hashtable_search(term_db->termid2index, term_id);
	if (indx == NULL) {
                return -1;
        }
	return  (int) ((intptr_t) indx) - 1;
}

static
CVTerm *GMTTermDb_get_term_from_term_id(CVTermDb *term_db_, char *term_id)
{
	GMTTermDb *term_db = (GMTTermDb *) term_db_;
        int i = GMTTermDb_get_index_from_term_id(term_db_, (char *) term_id);
        CVTerm *term;
        if (i > 0) {
                term = term_db->terms + i;
                return term;
        }
        return NULL;
}


static
int GMTTermDb_insert_term(CVTermDb *term_db_, const char *term_id,
                          const char *namespace, const char *description)
{
	GMTTermDb *term_db = (GMTTermDb *) term_db_;
	CVTerm *term = GMTTermDb_get_term_from_term_id(term_db_, (char *) term_id);

	if (term != NULL) {
		return 0; /* term_id already exists */
	}

	if (term_db->num_terms  >= term_db->max_terms) {
		term_db->max_terms *= 2;
		term_db->terms = realloc_(term_db->terms,
					  term_db->max_terms*sizeof(CVTerm));
	}
	term = term_db->terms + term_db->num_terms;
	term_db->num_terms++;
	term = CVTerm_init(term);
	term->term_id = strdup_(term_id);
	term->description = strdup_(description);
	term->namespace = namespace;

	if (! hashtable_insert(term_db->termid2index, strdup_(term_id),
                               (void *) ((intptr_t) term_db->num_terms))) {
		fprintf(stderr, "Could not insert term %s.\n", term_id);
		exit(EXIT_FAILURE);
	}
	return 1;
}

static
GMTTermDb *GMTTermDb_init(void)
{
	GMTTermDb *term_db = calloc_(1, sizeof(GMTTermDb));
	term_db->delete = GMTTermDb_delete;
	term_db->get_term_from_index = GMTTermDb_get_term_from_index;
	term_db->get_term_from_term_id = GMTTermDb_get_term_from_term_id;
        term_db->get_index_from_term_id = GMTTermDb_get_index_from_term_id;
        term_db->insert_term = GMTTermDb_insert_term;
	term_db->terms = calloc_(INIT_MAX_TERMS, sizeof(CVTerm));
	term_db->num_terms = 0;
	term_db->max_terms = INIT_MAX_TERMS;
	term_db->termid2index = create_hashtable(INIT_MAX_TERMS, hash_from_string, str_equal);
	return term_db;
}




int GMT_enrichment_context(const char *gmt_filename, const char *namespace, EntityDb **entity_db_,
			   CVTermDb **term_db_, TermMappingDb **mapping_db_)
{
	unsigned int line_num = 1;
	unsigned int current_field;
	char buf1[GMT_MAX_FIELD_WIDTH];
	char buf2[GMT_MAX_FIELD_WIDTH];
	FSFile_proc_code_type retcode;
	uint32_t entity_index;
	int active_term;

	EntityDb *entity_db = *entity_db_;
	CVTermDb *term_db = *term_db_;
        TermMappingDb *mapping_db = *mapping_db_;
	FILE *fp;

	if (entity_db == NULL) {
		entity_db = (EntityDb *) GMTEntityDb_init();
	}
	if (term_db == NULL) {
		term_db = (CVTermDb *) GMTTermDb_init();
	}
	if (mapping_db == NULL) {
		mapping_db = TermMappingDb_init();
	}

	fp = fopen(gmt_filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "Could not open file %s.\n", gmt_filename);
		exit(EXIT_FAILURE);
	}


	while (END_OF_FILE != (retcode = FSFile_next_field(fp, '\t', buf1, GMT_MAX_FIELD_WIDTH))) {
		/* First read term_id and description */

 		if (retcode != END_OF_FIELD) {
			fprintf(stderr, "Invalid gmt file format (file %s, line %d, field #1).\n",
				gmt_filename, line_num);
			exit(EXIT_FAILURE);
		}
		retcode = FSFile_next_field(fp, '\t', buf2, GMT_MAX_FIELD_WIDTH);
		if (retcode != END_OF_FIELD) {
			fprintf(stderr, "Invalid gmt file format (file %s, line %d, field #2).\n",
				gmt_filename, line_num);
			exit(EXIT_FAILURE);
		}

		active_term = term_db->insert_term(term_db, buf1, namespace, buf2);

		if (active_term) {
			mapping_db->insert_new_mapping(mapping_db);
		}

		/* Now process all entities (genes) for that term */
		current_field = 3;
		do {
			retcode = FSFile_next_field(fp, '\t', buf1, GMT_MAX_FIELD_WIDTH);
			if (retcode == TRUNCATED_FIELD) {
				fprintf(stderr, "Maximum field width (%d chars) exceeded (file %s, line %d, field #%d).\n",
					GMT_MAX_FIELD_WIDTH-1, gmt_filename, line_num, current_field);
				exit(EXIT_FAILURE);
			}
                        if (strlen(buf1) < 1) {
				fprintf(stderr, "Empty entity field (file %s, line %d, field #%d).\n",
					gmt_filename, line_num, current_field);
				exit(EXIT_FAILURE);
                        }

			if (active_term) {
				entity_index = entity_db->insert_item(entity_db, buf1, NULL);
                                mapping_db->insert_hit(mapping_db, entity_index);
			}
			current_field++;
		} while (retcode == END_OF_FIELD);

		line_num++;
	}

	fclose(fp);

	*entity_db_ = entity_db;
	*term_db_ = term_db;
	*mapping_db_ = mapping_db;

	return 1;
}
