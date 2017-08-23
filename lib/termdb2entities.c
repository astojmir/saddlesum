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
#include "miscutils.h"
#include "enrich.h"


static
void TermMappingDb_delete(TermMappingDb *mapping_db)
{
	free(mapping_db->hits);
	mapping_db->hits = NULL;
	mapping_db->max_hits = 0;
	mapping_db->num_hits = 0;
	free(mapping_db->offsets);
	mapping_db->offsets = NULL;
	mapping_db->max_mappings = 0;
	mapping_db->num_mappings = 0;
	free(mapping_db);
}


static
void TermMappingDb_reset(TermMappingDb *mapping_db)
{
	mapping_db->current_term = 0;
}


static
int TermMappingDb_get_next_mapping(TermMappingDb *mapping_db, uint32_t *term_index,
                                   uint32_t **hits, uint32_t *num_hits)
{
	uint32_t i = mapping_db->current_term;
	if (i >= mapping_db->num_mappings) {
		*hits = NULL;
		return 0;
	}
	*term_index = i;
	*hits = mapping_db->hits + mapping_db->offsets[i];
	*num_hits = mapping_db->offsets[i+1] - mapping_db->offsets[i];
	mapping_db->current_term++;
	return 1;
}


static
void TermMappingDb_insert_new_mapping(TermMappingDb *mapping_db)
{
        mapping_db->num_mappings++;
        if (mapping_db->num_mappings >= mapping_db->max_mappings) {
                mapping_db->max_mappings *= 2;
                mapping_db->offsets = realloc_(mapping_db->offsets,
                                               mapping_db->max_mappings*sizeof(uint32_t));
        }
        mapping_db->offsets[mapping_db->num_mappings] = mapping_db->num_hits;
}


static
void TermMappingDb_insert_hit(TermMappingDb *mapping_db, uint32_t entity_index)
{
        if (mapping_db->num_hits >= mapping_db->max_hits) {
                mapping_db->max_hits *= 2;
                mapping_db->hits = realloc_(mapping_db->hits,
                                            mapping_db->max_hits*sizeof(uint32_t));
        }
        mapping_db->hits[mapping_db->num_hits++] = entity_index;
        mapping_db->offsets[mapping_db->num_mappings] = mapping_db->num_hits;
}


TermMappingDb *TermMappingDb_init(void)
{
	TermMappingDb *mapping_db = calloc_(1, sizeof(TermMappingDb));
	mapping_db->delete = TermMappingDb_delete;
	mapping_db->reset = TermMappingDb_reset;
	mapping_db->get_next_mapping = TermMappingDb_get_next_mapping;
        mapping_db->insert_new_mapping = TermMappingDb_insert_new_mapping;
        mapping_db->insert_hit = TermMappingDb_insert_hit;
	mapping_db->hits = calloc_(INIT_MAX_ENTITIES, sizeof(uint32_t));
	mapping_db->max_hits = INIT_MAX_ENTITIES;
	mapping_db->num_hits = 0;
	mapping_db->offsets = calloc_(INIT_MAX_TERMS + 1, sizeof(uint32_t));
	mapping_db->max_mappings = INIT_MAX_TERMS + 1;
	mapping_db->num_mappings = 0;
	mapping_db->current_term = 0;
	return mapping_db;
}

