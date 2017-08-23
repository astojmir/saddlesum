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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include "enrich.h"
#include "miscutils.h"
#include "help.h"

#define MIN_ARGS 1
#define OUTPUT_FMT_OPTS 2
#define FULL_VERSION VERSION " (qmbpmn-tools-" VERSION ")"

#define option_err_msg(msg) fprintf(stderr, "%s: %s\nFor help type %s -h\n", \
                                    argv[0], msg, argv[0]); \
                            exit(EXIT_FAILURE)


char help_msg[] = HELP_SHOW_ETD;



int main(int argc, char **argv)
{
        /* Option parsing */
        int c;
        const char *output_fmt_types[OUTPUT_FMT_OPTS] = {"txt", "tab"};
        int i;

        /* Output arguments */
        const char *output_filename = NULL;
        OutputType output_type = TEXT;
        FILE *fp = stdout;
        uint8_t show_namespaces_only = 0;

	const char *etd_filename;

        opterr = 0;
        while ( (c = getopt(argc, argv, "VhNO:F:")) != -1) {
                switch (c) {
                case 'V':
                        printf("%s: version %s\n", argv[0], FULL_VERSION);
                        exit(EXIT_SUCCESS);
                        break;
                case 'h':
                        printf("%s", help_msg);
                        exit(EXIT_SUCCESS);
                        break;
                case 'N':
                        show_namespaces_only = 1;
                        break;
                case 'O':
                        output_filename = optarg;
                        fp = fopen(output_filename, "w");
                        if (fp == NULL) {
                                fprintf(stderr, "%s: Could not open output file %s.\n",
                                        argv[0], output_filename);
                                exit(EXIT_FAILURE);
                        }
                        break;
                case 'F':
                        for (i=0; i < OUTPUT_FMT_OPTS; i++) {
                                if (!strcmp(optarg, output_fmt_types[i])) {
                                        output_type = i;
                                        break;
                                }
                        }
                        if (i >= OUTPUT_FMT_OPTS) {
                                option_err_msg("Invalid argument for option -F.");
                        }
                        break;
                case '?':
                        fprintf(stderr, "%s: Invalid option -- %c.\nFor help type %s -h.\n",
                                argv[0], optopt, argv[0]);
                        exit(EXIT_FAILURE);
                        break;
                default:
                        option_err_msg("Invalid arguments.");
                }
        }
        if (argc < optind + MIN_ARGS) {
                option_err_msg("Insufficient arguments.");
        }

        etd_filename = argv[optind];

        if (show_namespaces_only) {
                ETDTermDb_print_namespaces(etd_filename, fp, output_type);
        }
        else {
                ETDTermDb_print_info(etd_filename, fp, output_type);
        }

        if (output_filename != NULL) {
                fclose(fp);
        }
	return EXIT_SUCCESS;
}
