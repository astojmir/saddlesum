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

#ifndef _ENTITY_H
#define _ENTITY_H
#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include <stdio.h>
#include "miscutils.h"

/* Basic entity object interface */
struct _Entity_s;

#define Entity_HEAD			              \
        void (*delete) (struct _Entity_s *);          \
        void (*print_url) (struct _Entity_s *,        \
                           PrintBuf *);               \
	char *symbol;                                 \
        char *description;


typedef struct _Entity_s {
	Entity_HEAD
} Entity;

void Entity_delete_data(Entity *entity);
void Entity_delete_all(Entity *entity);
void Entity_delete_struct_only(Entity *entity);


typedef enum {UNKNOWN_ID,
	      RESOLVABLE_CONFLICT,
	      UNRESOLVABLE_CONFLICT,
              DUPLICATE_ID} EntityWarning_code_type;


typedef struct _EntityWarning_s {
	struct _EntityWarning_s *next;
	char *msg;
        int code;
} EntityWarning;


/* Basic entity database object interface */
struct _EntityDb_s;

#define EntityDb_HEAD                                                              \
        uint32_t num_entities;                                                     \
        void (*delete) (struct _EntityDb_s *);                                     \
        EntityWarning *(*map_symbol) (struct _EntityDb_s *, char *, uint32_t *);   \
        Entity * (*get_entity_from_index) (struct _EntityDb_s *, uint32_t);        \
        uint32_t (*insert_item) (struct _EntityDb_s *, const char *, const char *);


typedef struct _EntityDb_s {
	EntityDb_HEAD
} EntityDb;


typedef struct _NCBIGenesDb_s {
	EntityDb_HEAD
        uint32_t max_entities;
        const char *gene_info_file;
        const char *url_fmt;  /* ignored in favor of NCBIGENE_URL_FMT */
        uint32_t tax_id;
        uint32_t *gene_ids;
        char **symbols;
        char **descriptions;
        /* buffers: symbols_buf and desc_buf can grow, others can't */
        char *metadata_buf;
        uint32_t symbols_buf_size;
        uint32_t symbols_buf_len;
        char *symbols_buf;
        uint32_t desc_buf_size;
        uint32_t desc_buf_len;
        char *desc_buf;
        char *conflicts1_buf;
        char *conflicts2_buf;
        struct hashtable *geneid2index;
        struct hashtable *alias2index;
        struct hashtable *conflicts;
} NCBIGenesDb;

NCBIGenesDb *NCBIGenesDb_init(FILE *fp);


#ifdef __cplusplus
}
#endif
#endif /* !_ENTITY_H */
