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

#ifndef _CVTERM_H
#define _CVTERM_H
#ifdef __cplusplus
extern "C" {
#endif


#include <inttypes.h>
#include "miscutils.h"

/* Basic cvterm object interface */
struct _CVTerm_s;

#define CVTerm_HEAD		                                          \
        void (*delete) (struct _CVTerm_s *);                              \
        void (*print_url) (struct _CVTerm_s *, PrintBuf *);               \
        uint32_t (*get_parents) (struct _CVTerm_s *,                      \
                                 struct _CVTerm_s ***, char ***);         \
	char *term_id;                                                    \
        const char *namespace;                                            \
        char *description;


typedef struct _CVTerm_s {
	CVTerm_HEAD
} CVTerm;

/* ETDTerm extends CVTerm */
#define ETDTerm_HEAD		                                          \
        uint32_t flag;                                                    \
        uint32_t num_parents;                                             \
        CVTerm **parents;                                                 \
        char **edgetypes;

typedef struct _ETDTerm_s {
	CVTerm_HEAD
        ETDTerm_HEAD
} ETDTerm;

/* KEGGTerm extends ETDTerm */
typedef struct _KEGGTerm_s {
	CVTerm_HEAD
        ETDTerm_HEAD
        const char *org_prefix;
} KEGGTerm;

CVTerm *CVTerm_init(CVTerm *);
CVTerm *GOTerm_init(CVTerm *);
CVTerm *KEGGTerm_init(CVTerm *);



/* Basic term database object interface */
struct _CVTermDb_s;

#define CVTermDb_HEAD                                                     \
        uint32_t num_terms;                                               \
        void (*delete) (struct _CVTermDb_s *);                            \
        CVTerm * (*get_term_from_index) (struct _CVTermDb_s *, uint32_t); \
        CVTerm * (*get_term_from_term_id) (struct _CVTermDb_s *, char *); \
        int (*get_index_from_term_id) (struct _CVTermDb_s *, char *);  \
        int (*insert_term) (struct _CVTermDb_s *, const char *,           \
                            const char *, const char *);

typedef struct _CVTermDb_s {
	CVTermDb_HEAD
} CVTermDb;


typedef struct _ETDTermDb_s {
	CVTermDb_HEAD
        uint32_t max_terms;
        CVTerm **terms;
        char *db_name;
        uint32_t num_namespaces;
        uint32_t *num_ns_terms;
        char **namespaces;
        char *namespace_buf;
        char **termid_bufs;
        char **desc_bufs;
        char **edgetype_bufs;
        char **metadata_bufs;
        struct hashtable *termid2index;
} ETDTermDb;

#ifdef __cplusplus
}
#endif
#endif /* !_CVTERM_H */
