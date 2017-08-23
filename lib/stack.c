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

#include <stdint.h>
#include <stdlib.h>
#include "miscutils.h"
#include "stack.h"

#define INITIAL_TERM_STACK_SIZE 1024


Stack *Stack_init(void)
{
        Stack *stack = malloc_(sizeof(Stack));
        stack->data_len = 0;
        stack->data_size = INITIAL_TERM_STACK_SIZE;
        stack->data = malloc_(INITIAL_TERM_STACK_SIZE * sizeof(void *));
        return stack;
}


void Stack_delete(Stack *stack)
{
        free(stack->data);
        free(stack);
}


int Stack_is_empty(Stack *stack)
{
        return (stack->data_len == 0);
}


void Stack_push(Stack *stack, void * item)
{
        if (stack->data_len >= stack->data_size) {
                stack->data_size *= 2;
                stack->data = realloc_(stack->data,
                                       stack->data_size * sizeof(void *));
        }
        stack->data[stack->data_len++] = item;
}


void * Stack_pop(Stack *stack)
{
        if (stack->data_len == 0) {
                return NULL;
        }
        return stack->data[--stack->data_len];
}
