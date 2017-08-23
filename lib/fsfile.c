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

#include <stdio.h>
#include <ctype.h>
#include "fsfile.h"

FSFile_proc_code_type FSFile_next_field(FILE *fp, char fs, char *buf,
					uint32_t max_chars)
{
	int c;
	int i;
	char *p;

	for (i=0, p=buf; i < max_chars-1; p++,i++) {
		c = getc(fp);
		if ((c == EOF) || (c == fs) || (c == '\n')) {
			break;
		}
		*p = c;
	}
	*p = '\0';

	if (i >= (max_chars-1)) {
		return TRUNCATED_FIELD;
	}
	if (c == '\n') {
		return END_OF_LINE;
	}
	if (c == EOF) {
		return END_OF_FILE;
	}
	return END_OF_FIELD;
}

FSFile_proc_code_type SSFile_next_field(FILE *fp, char *buf,
					uint32_t max_chars)
{
	int c;
	int i = 0;
	char *p;

	do {
		c = getc(fp);
	} while (isspace(c));
	ungetc(c, fp);

	for (i=0, p=buf; i < max_chars-1; p++,i++) {
		c = getc(fp);
		if ((c == EOF) || (c == '\n')) {
			break;
		}
		if (isspace(c)) {
			do {
				c = getc(fp);
				if (c == '\n') {
					break;
				}
				if (!isspace(c)) {
					ungetc(c, fp);
					break;
				}
			} while (c != EOF);
			break;
		}
		*p = c;
	}
	*p = '\0';

	if (i >= (max_chars-1)) {
		return TRUNCATED_FIELD;
	}
	if (c == '\n') {
		return END_OF_LINE;
	}
	if (c == EOF) {
		return END_OF_FILE;
	}
	return END_OF_FIELD;
}
