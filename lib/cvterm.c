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
#include "cvterm.h"
#include "miscutils.h"

#define GOTERM_URL_FMT "http://amigo.geneontology.org/cgi-bin/amigo/"     \
                       "term-details.cgi?term=%s"

#define KEGG_LEAF_URL_FMT "http://www.genome.jp/kegg-bin/show_pathway?"   \
                          "org_name=%s&amp;mapno=%s&amp;mapscale=1.0"     \
                          "&amp;show_description=show"

#define KEGG_HIGHER_URL_FMT "http://www.genome.jp/kegg-bin/get_htext?"    \
                            "htext=br08901.keg&amp;filedir=%%2ffiles"     \
                            "&amp;extend=&amp;open=%s"

#define KEGG_ROOT_URL_FMT "http://www.genome.jp/kegg/pathway.html"


/* Plain CVTerm owns term_id and description members */
static
void CVTerm_delete_data(CVTerm *term)
{
        free(term->term_id);
        free(term->description);
	/* Do not free namespace - it will be owned by the term database */
}


static
void CVTerm_delete_all(CVTerm *term)
{
        CVTerm_delete_data(term);
        free(term);
}


CVTerm *CVTerm_init(CVTerm *term)
{
        if (term == NULL) {
                term = calloc_(1, sizeof(CVTerm));
                term->delete = CVTerm_delete_all;
        }
        else {
                term->delete = CVTerm_delete_data;
        }
        term->get_parents = NULL;
        term->print_url = NULL;
        return term;
}


/* ETDTerm does not own most of its data */
static
void ETDTerm_delete_data(CVTerm *term_)
{
        ETDTerm *term = (ETDTerm *) term_;
        if (term->num_parents > 0) {
                free(term->parents);
                free(term->edgetypes);
        }
        free(term);
}


static
void ETDTerm_delete_all(CVTerm *term)
{
        ETDTerm_delete_data(term);
        free(term);
}


static
uint32_t ETDTerm_get_parents(CVTerm *term_, CVTerm ***parents,
                             char ***edgetypes)
{
        ETDTerm *term = (ETDTerm *) term_;
        if (term->num_parents > 0) {
                *parents = term->parents;
                *edgetypes = term->edgetypes;
        }
        return term->num_parents;
}


static
void GOTerm_print_url(CVTerm *term_, PrintBuf *pbuf)
{
        ETDTerm *term = (ETDTerm *) term_;
        PrintBuf_printf(pbuf, 0, GOTERM_URL_FMT, term->term_id);
}


CVTerm *GOTerm_init(CVTerm *term_)
{
        ETDTerm *term = (ETDTerm *) term_;
        if (term == NULL) {
                term = calloc_(1, sizeof(KEGGTerm));
                term->delete = ETDTerm_delete_all;
        }
        else {
                term->delete = ETDTerm_delete_data;
        }
        term->get_parents = ETDTerm_get_parents;
        term->print_url = GOTerm_print_url;
        return (CVTerm *) term;
}


static
void KEGGTerm_print_url(CVTerm *term_, PrintBuf *pbuf)
{
        KEGGTerm *term = (KEGGTerm *) term_;
        const char *term_id;
        switch (term->flag) {
        case 0:
                term_id = term->term_id + 8; /* skip 'KEGG:xxx' */
                PrintBuf_printf(pbuf, 0, KEGG_LEAF_URL_FMT,
                                term->org_prefix, term_id);
                break;
        case 1:
                term_id = term->term_id + 5; /* skip 'KEGG:' */
                PrintBuf_printf(pbuf, 0, KEGG_HIGHER_URL_FMT, term_id);
                break;
        case 2:
        default:
                PrintBuf_printf(pbuf, 0, KEGG_ROOT_URL_FMT);
                break;
        }
}


CVTerm *KEGGTerm_init(CVTerm *term_)
{
        KEGGTerm *term = (KEGGTerm *) term_;
        if (term == NULL) {
                term = calloc_(1, sizeof(KEGGTerm));
                term->delete = ETDTerm_delete_all;
        }
        else {
                term->delete = ETDTerm_delete_data;
        }
        term->get_parents = ETDTerm_get_parents;
        term->print_url = KEGGTerm_print_url;
        return (CVTerm *) term;
}



