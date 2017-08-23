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
#include "enrich.h"

#define NCBIGENES_ITEM_INCR 4096
#define NCBIGENES_SYMB_BUF_INCR 4096
#define NCBIGENES_DESC_BUF_INCR 16384

#define NCBIGENE_MAGIC 1200900292U  /* from ncbi_gene.py */
#define NCBIGENE_URL_FMT "http://www.ncbi.nlm.nih.gov/sites/entrez?" \
                         "db=gene&cmd=Retrieve&dopt=Graphics"        \
                         "&list_uids=%ld"

/* GeneEntity */
typedef struct _GeneEntity_s {
	Entity_HEAD
        uint32_t gene_id;
} GeneEntity;


static
void GeneEntity_print_url(Entity *entity_, PrintBuf *pbuf)
{
        GeneEntity *entity = (GeneEntity *) entity_;
        PrintBuf_printf(pbuf, 0, NCBIGENE_URL_FMT, (long) entity->gene_id);
}


/* Symbol conflicts */
typedef struct _SymbolConflict_s {
        EntityWarning_code_type wtype;
        uint32_t num_aliases;
        char *aliases; /* pointer to first symbol */
} SymbolConflict;


static
void print_id_list(char *src, int num_items, PrintBuf *pbuf)
{
        int i;
        if (num_items < 1) {
                return; /* just in case */
        }
        PrintBuf_printf(pbuf, -1, "%s", src);
        src += (strlen(src) + 1);
        for (i=1; i < (num_items-1); i++) {
                PrintBuf_printf(pbuf, -1, ", %s", src);
                src += (strlen(src) + 1);
        }
        if (num_items > 1) {
                PrintBuf_printf(pbuf, -1, " and %s", src);
        }
}


static
EntityWarning *SymbolConflict_warning(SymbolConflict *conflict)
{
        PrintBuf pbuf;
	//	EntityWarning *warning = malloc_(sizeof(EntityWarning));
	EntityWarning *warning = calloc_(1, sizeof(EntityWarning));
        warning->code = conflict->wtype;
        char *cur = conflict->aliases;

        /* Initialize PrintBuf on the stack. EntityWarning instance
           will take ownership of the actual buffer so PrintBuf_delete
           should not be called. */
        (void) PrintBuf_init(&pbuf);

        if (conflict->wtype == RESOLVABLE_CONFLICT) {
                PrintBuf_printf(&pbuf, -1,  "Identifier %s is also a synonym for ",
                                cur);
                cur += (strlen(cur) + 1);
                print_id_list(cur, conflict->num_aliases - 1, &pbuf);
                PrintBuf_printf(&pbuf, -1,  ".");
        }
        else {
                PrintBuf_printf(&pbuf, -1,
                                "%s is not a primary identifier for any"
                                " gene while being an alias for ",
                                cur);
                cur += (strlen(cur) + 1);
                print_id_list(cur, conflict->num_aliases - 1, &pbuf);
                PrintBuf_printf(&pbuf, -1,  ".");
        }
        warning->msg = pbuf.buf;
        return warning;
}


/* Entity database */
static
void NCBIGenesDb_delete(EntityDb *entity_db_)
{
	NCBIGenesDb *entity_db = (NCBIGenesDb *) entity_db_;
	hashtable_destroy(entity_db->geneid2index, 0);
	hashtable_destroy(entity_db->alias2index, 0);
	hashtable_destroy(entity_db->conflicts, 1);
	free(entity_db->gene_ids);
	free(entity_db->symbols);
        free(entity_db->descriptions);
	free(entity_db->metadata_buf);
	free(entity_db->symbols_buf);
	free(entity_db->desc_buf);
        free(entity_db->conflicts1_buf);
        free(entity_db->conflicts2_buf);
	free(entity_db);
}


static
int NCBIGenesDb_search_symbol(NCBIGenesDb *entity_db, const char *symbol, uint32_t *i)
{
        char *endptr;
        uint32_t gene_id;
        gene_id = (uint32_t) strtoul(symbol, &endptr, 0);
        void *indx = NULL;
        if (*endptr == '\0') {  /* string is valid */
                indx = hashtable_search(entity_db->geneid2index, &gene_id);
        }
	if (indx == NULL) {
                /* The cast is OK since the key is not modified by hashtable_search */
                indx = hashtable_search(entity_db->alias2index, (char *) symbol);
        }
	if (indx == NULL) {
                return 0;
        }
        *i = (uint32_t) ((intptr_t) indx) - 1;
        return 1;
}


static
uint32_t NCBIGenesDb_insert_item(EntityDb *entity_db_, const char *symbol,
                                 const char *description)
{
        int n;
	NCBIGenesDb *entity_db = (NCBIGenesDb *) entity_db_;
	uint32_t i;
        int found = NCBIGenesDb_search_symbol(entity_db, symbol, &i);

	if (!found) {
		if (entity_db->num_entities >= entity_db->max_entities) {
			entity_db->max_entities += NCBIGENES_ITEM_INCR;
			entity_db->gene_ids = realloc_(entity_db->gene_ids,
						       entity_db->max_entities*sizeof(uint32_t));
			entity_db->symbols = realloc_(entity_db->symbols,
                                                      entity_db->max_entities*sizeof(char *));
			entity_db->descriptions = realloc_(entity_db->descriptions,
                                                           entity_db->max_entities*sizeof(char *));
		}
		i = entity_db->num_entities++;
                entity_db->gene_ids[i] = 0;

                n = strlen(symbol) + 1;
                if (entity_db->symbols_buf_len + n >= entity_db->symbols_buf_size) {
                        entity_db->symbols_buf_size += (n + NCBIGENES_SYMB_BUF_INCR);
                        entity_db->symbols_buf = realloc_(entity_db->symbols_buf,
                                                          entity_db->symbols_buf_size);
                }
                memcpy(entity_db->symbols_buf + entity_db->symbols_buf_len, symbol, n);
                entity_db->symbols[i] = entity_db->symbols_buf + entity_db->symbols_buf_len;
                entity_db->symbols_buf_len += n;

                if (description == NULL) {
                        description = "";
                }

                n = strlen(description) + 1;
                if (entity_db->desc_buf_len + n >= entity_db->desc_buf_size) {
                        entity_db->desc_buf_size += (n + NCBIGENES_DESC_BUF_INCR);
                        entity_db->desc_buf = realloc_(entity_db->desc_buf,
                                                          entity_db->desc_buf_size);
                }
                memcpy(entity_db->desc_buf + entity_db->desc_buf_len, description, n);
                entity_db->descriptions[i] = entity_db->desc_buf + entity_db->desc_buf_len;
                entity_db->desc_buf_len += n;

		if (! hashtable_insert(entity_db->alias2index, strdup_(symbol),
                                       (void *) ((intptr_t) (i+1))) ) {
			fprintf(stderr, "Could not insert entity item %s.\n", symbol);
			exit(EXIT_FAILURE);
		}
	}
	return i;
}


static
EntityWarning *NCBIGenesDb_map_symbol(EntityDb *entity_db_, char *symbol, uint32_t *i)
{
	PrintBuf pbuf;
        SymbolConflict *conflict;
	EntityWarning *warning = NULL;
	NCBIGenesDb *entity_db = (NCBIGenesDb *) entity_db_;
        int found = NCBIGenesDb_search_symbol(entity_db, symbol, i);

	if (!found) {
                conflict = (SymbolConflict *) hashtable_search(entity_db->conflicts, symbol);
                if (conflict != NULL) {
                        /* unresolvable conflict */
                        warning = SymbolConflict_warning(conflict);
                }
                else {
                        /* Initialize PrintBuf on the stack. EntityWarning instance
                           will take ownership of the actual buffer so PrintBuf_delete
                           should not be called. */
                        (void) PrintBuf_init(&pbuf);

                        warning = calloc_(1, sizeof(EntityWarning));
                        PrintBuf_printf(&pbuf, 0, "%s", symbol);
                        warning->msg = pbuf.buf;
                        warning->code = UNKNOWN_ID;
                }
	}
        else {
                conflict = (SymbolConflict *) hashtable_search(entity_db->conflicts, symbol);
                if (conflict != NULL) {
                        /* resolvable conflict */
                        warning = SymbolConflict_warning(conflict);
                }
        }
	return warning;
}


static
Entity *NCBIGenesDb_get_entity_from_index(EntityDb *entity_db_, uint32_t i)
{
	NCBIGenesDb *entity_db = (NCBIGenesDb *) entity_db_;
        GeneEntity *entity;
	if (i < entity_db->num_entities) {
                entity = calloc_(1, sizeof(GeneEntity));
                entity->delete = Entity_delete_struct_only;
                entity->symbol = entity_db->symbols[i];
                entity->description = entity_db->descriptions[i];
                entity->gene_id = entity_db->gene_ids[i];
                if (entity->gene_id) {
                        entity->print_url = GeneEntity_print_url;
                }
                else {
                        entity->print_url = NULL;
                }
                return (Entity *) entity;
	}
	return NULL;
}

static
void NCBIGenesDb_insert_conflicts(FILE *fp, NCBIGenesDb *entity_db,
                                  EntityWarning_code_type wtype)
{
        uint32_t n;
        uint32_t tmp;
        uint32_t *counts;
        char *cur;
        SymbolConflict *cnf;
        int i;
        int j;

        fread_uint32(&n, 1, fp);
        counts = malloc_(n * sizeof(uint32_t));
        fread_uint32(counts, n, fp);

        cur = fread_buf(fp, &tmp);
        if (wtype == RESOLVABLE_CONFLICT) {
                entity_db->conflicts1_buf = cur;
        }
        else {
                entity_db->conflicts2_buf = cur;
        }

        for (i=0; i < n; i++) {
                cnf = malloc_(sizeof(SymbolConflict));
                cnf->wtype = wtype;
                cnf->num_aliases = counts[i];
                cnf->aliases = cur;
                if (! hashtable_insert(entity_db->conflicts, strdup_(cur), (void *) cnf) ) {
                        fprintf(stderr, "Could not insert conflict for %s.\n", cur);
                        exit(EXIT_FAILURE);
                }
                for (j=0; j < counts[i]; j++) {
                        cur += (strlen(cur) + 1);
                }
        }
        free(counts);
}

NCBIGenesDb *NCBIGenesDb_init(FILE *fp)
{
        char header[8];
        uint32_t tmp;
        uint32_t n;
        uint32_t *counts;
        uint32_t *gene_id_;
        char *buf;
        char *cur;
        char *dest;
        intptr_t i;
        int j;
	NCBIGenesDb *entity_db = calloc_(1, sizeof(NCBIGenesDb));

        /* header */
        fread_(header, 1, 8, fp);
        fread_uint32(&tmp, 1, fp);
        if (memcmp(header, "NCBIGENE", 8) || (tmp != NCBIGENE_MAGIC)) {
                fprintf(stderr, "Invalid binary file format for gene_info index.\n");
                exit(EXIT_FAILURE);
        }

        /* metadata */
        entity_db->metadata_buf = fread_buf(fp, &tmp);
        entity_db->gene_info_file = entity_db->metadata_buf;
        entity_db->url_fmt = entity_db->metadata_buf + strlen(entity_db->metadata_buf) + 1;

        /* skip checksum */
        fread_uint32(&tmp, 1, fp);

        /* tax_id */
        fread_uint32(&entity_db->tax_id, 1, fp);

        /* gene count - allocate all related arrays */
        fread_uint32(&entity_db->num_entities, 1, fp);
	entity_db->max_entities = entity_db->num_entities;
        entity_db->gene_ids = malloc_(entity_db->num_entities * sizeof(uint32_t));
        entity_db->symbols = malloc_(entity_db->num_entities * sizeof(char *));
        entity_db->descriptions = malloc_(entity_db->num_entities * sizeof(char *));
	entity_db->alias2index = create_hashtable(entity_db->num_entities,
                                                  hash_from_string, str_equal);
	entity_db->geneid2index = create_hashtable(entity_db->num_entities,
                                                   hash_from_int, int_equal);
	entity_db->conflicts = create_hashtable(entity_db->num_entities / 2,
                                                hash_from_string, str_equal);

        /* gene_ids */
        fread_uint32(entity_db->gene_ids, entity_db->num_entities, fp);
        for (i=0; i < entity_db->num_entities; i++) {
                gene_id_ = malloc_(sizeof(uint32_t));
                *gene_id_ = entity_db->gene_ids[i];
		if (! hashtable_insert(entity_db->geneid2index,
                                       gene_id_, (void *) (i+1)) ) {
			fprintf(stderr, "Could not insert gene id %ld.\n",
                                (long int) entity_db->gene_ids[i]);
			exit(EXIT_FAILURE);
		}
        }

        /* skip offsets */
        fseek(fp, (entity_db->num_entities+1)*sizeof(uint32_t), SEEK_CUR);

        /* symbols */
        counts = malloc_(entity_db->num_entities * sizeof(uint32_t));
        fread_uint32(counts, entity_db->num_entities, fp);
        buf = fread_buf(fp, &tmp); /* temporary buffer */

        /* insert all symbols and aliases into entity_db->alias2index */
        /*    (also count cannonical symbols) */
        cur = buf;
        for (i=0; i < entity_db->num_entities; i++) {
                entity_db->symbols_buf_len += (strlen(cur) + 1);
                for (j=0; j < counts[i]; j++) {
                        if (! hashtable_insert(entity_db->alias2index,
                                               strdup_(cur),
                                               (void *) (i+1)) ) {
                                fprintf(stderr, "Could not insert entity item %s.\n", cur);
                                exit(EXIT_FAILURE);
                        }
                        cur += (strlen(cur) + 1);
                }
        }

        /* Insert cannonical symbols into entity_db->symbols_buf and make
           entity_db->symbols items point to them */
        entity_db->symbols_buf_size = entity_db->symbols_buf_len;
        entity_db->symbols_buf = malloc_(entity_db->symbols_buf_size);
        cur = buf;
        dest = entity_db->symbols_buf;
        for (i=0; i < entity_db->num_entities; i++) {
                n = strlen(cur) + 1;
                memcpy(dest, cur, n);
                entity_db->symbols[i] = dest;
                dest += n;
                cur += n;
                for (j=1; j < counts[i]; j++) {
                        n = strlen(cur) + 1;
                        cur += n;
                }
        }
        free(buf);
        free(counts);

        /* descriptions */
        entity_db->desc_buf = fread_buf(fp, &entity_db->desc_buf_size);
        entity_db->desc_buf_len = entity_db->desc_buf_size;
        strbuf2array(entity_db->descriptions, entity_db->desc_buf, entity_db->num_entities);

        /* conflicts */
        NCBIGenesDb_insert_conflicts(fp, entity_db, RESOLVABLE_CONFLICT);
        NCBIGenesDb_insert_conflicts(fp, entity_db, UNRESOLVABLE_CONFLICT);

        /* function pointers */
	entity_db->delete = NCBIGenesDb_delete;
	entity_db->map_symbol = NCBIGenesDb_map_symbol;
	entity_db->get_entity_from_index = NCBIGenesDb_get_entity_from_index;
        entity_db->insert_item = NCBIGenesDb_insert_item;
	return entity_db;
}







