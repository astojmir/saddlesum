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

#ifndef _TERMDB2ENTITIES_H
#define _TERMDB2ENTITIES_H
#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>

/* Basic interface for a collection of mappings from terms to sets of entities */
struct _TermMappingDb_s;

#define TermMappingDb_HEAD                                              \
        uint32_t num_mappings;                                          \
        void (*delete) (struct _TermMappingDb_s *);                     \
        void (*reset) (struct _TermMappingDb_s *);                      \
        int (*get_next_mapping) (struct _TermMappingDb_s *, uint32_t *, \
                                 uint32_t **, uint32_t *);              \
        void (*insert_new_mapping) (struct _TermMappingDb_s *);         \
        void (*insert_hit) (struct _TermMappingDb_s *, uint32_t);       \
        uint32_t *hits;                                                 \
        uint32_t num_hits;                                              \
        uint32_t max_hits;                                              \
        uint32_t *offsets;                                              \
        uint32_t max_mappings;                                          \
        uint32_t current_term;

typedef struct _TermMappingDb_s {
	TermMappingDb_HEAD
} TermMappingDb;

TermMappingDb *TermMappingDb_init(void);


#ifdef __cplusplus
}
#endif
#endif /* !_TERMDB2ENTITIES_H */
