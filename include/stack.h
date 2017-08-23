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

#ifndef _ENRICHSTACK_H
#define _ENRICHSTACK_H
#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

typedef struct _Stack_s {
        unsigned int data_len;
        unsigned int data_size;
        void **data;
} Stack;

Stack *Stack_init(void);
void Stack_delete(Stack *stack);
int Stack_is_empty(Stack *stack);
void Stack_push(Stack *stack, void * item);
void * Stack_pop(Stack *stack);

#ifdef __cplusplus
}
#endif
#endif /* !_ENRICHSTACK_H */
