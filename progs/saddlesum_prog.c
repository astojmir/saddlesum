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

#define MAX_EXCLUDED 4
#define MIN_ARGS 2
#define STAT_OPTS 2
#define TRANSFORM_OPTS 2
#define OUTPUT_FMT_OPTS 2
#define FULL_VERSION VERSION " (qmbpmn-tools-" VERSION ")"

#define option_err_msg(msg) fprintf(stderr, "%s: %s\nFor help type %s -h\n", \
                                    argv[0], msg, argv[0]); \
                            exit(EXIT_FAILURE)


char help_msg[] = HELP_SADDLESUM;

static void get_namespace_and_file(char *arg, char **namespace, char **filename,
				   char sepchar)
{
	*namespace = NULL;
	*filename = arg;

	for (; *arg; arg++) {
		if (*arg == sepchar) {
			*arg = '\0';
			*namespace = *filename;
			*filename = arg+1;
			break;
		}
	}
}

int main(int argc, char **argv)
{
        /* Option parsing */
        int c;
        long int tmp_long;
        double tmp_dbl;
        char *tailptr;
        const char *stats_types[STAT_OPTS] = {"wsum", "hgem"};
        const char *transform_types[TRANSFORM_OPTS] = {"flip", "abs"};
        const char *output_fmt_types[OUTPUT_FMT_OPTS] = {"txt", "tab"};
        int i;

        /* Enrichment context arguments (set to defaults) */
        const char *db_name = NULL;
        uint32_t min_term_size = 3;
        double Evalue_cutoff = 1e-02;
        double effective_db_size = -1.0;
        EnrichStats statistics_type = SADDLESUM;
        uint8_t discretized_weights = 0;
	TransformType transform_type = NO_TRANSFORM;
        CutoffType cutoff_type = NONE;
        uint32_t rank_cutoff = 0;
        double weight_cutoff = 0.0;
        uint32_t use_all_weights = 0;
        EnrichContext *cntxt;


        /* Output arguments */
        const char *output_filename = NULL;
        OutputType output_type = TEXT;
        uint8_t print_warnings = 0;
        uint8_t print_unknown_ids = 0;
        FILE *fp = stdout;

        /* Database and weights variables */
	char *tdb_filename;
	const char *weights_filename;
	char *namespace;
	EntityDb *entity_db = NULL;
	CVTermDb *term_db = NULL;
	TermMappingDb *mapping_db = NULL;

        /* Excluded namespaces (ETD term_db only) */
        int num_excluded = 0;
        int max_excluded = MAX_EXCLUDED;
        const char **excluded_namespaces = malloc_(MAX_EXCLUDED * sizeof(const char *));

        /* Option to print specific term instead of doing search */
        const char *term_id = NULL;
        int term_index;

        opterr = 0;
        while ( (c = getopt(argc, argv, "Vhm:e:n:s:t:dr:w:x:aT:O:F:WU")) != -1) {
                switch (c) {
                case 'V':
                        printf("%s: standalone SaddleSum, version %s\n", argv[0], FULL_VERSION);
                        exit(EXIT_SUCCESS);
                        break;
                case 'h':
                        printf("%s", help_msg);
                        exit(EXIT_SUCCESS);
                        break;
                case 'm':
                        tmp_long = strtol(optarg, &tailptr, 10);
                        if (tailptr == optarg || tmp_long < 1 || tmp_long > UINT32_MAX) {
                                option_err_msg("Invalid argument for option -m.");
                        }
                        min_term_size = tmp_long;
                        break;
                case 'e':
                        tmp_dbl = strtod(optarg, &tailptr);
                        if (tailptr == optarg || tmp_dbl < 0.0) {
                                option_err_msg("Invalid argument for option -e.");
                        }
                        Evalue_cutoff = tmp_dbl;
                        break;
                case 'n':
                        tmp_dbl = strtod(optarg, &tailptr);
                        if (tailptr == optarg || tmp_dbl < 0.0) {
                                option_err_msg("Invalid argument for option -n.");
                        }
                        effective_db_size = tmp_dbl;
                        break;
                case 's':
                        for (i=0; i < STAT_OPTS; i++) {
                                if (!strcmp(optarg, stats_types[i])) {
                                        statistics_type = i;
                                        break;
                                }
                        }
                        if (i >= STAT_OPTS) {
                                option_err_msg("Invalid argument for option -s.");
                        }
                        break;
                case 't':
                        for (i=0; i < TRANSFORM_OPTS; i++) {
                                if (!strcmp(optarg, transform_types[i])) {
                                        transform_type = i+1;
                                        break;
                                }
                        }
                        if (i >= TRANSFORM_OPTS) {
                                option_err_msg("Invalid argument for option -t.");
                        }
                        break;
                case 'd':
                        discretized_weights = 1;
                        break;
                case 'r':
                        tmp_long = strtol(optarg, &tailptr, 10);
                        if (tailptr == optarg || tmp_long < 1 || tmp_long > UINT32_MAX) {
                                option_err_msg("Invalid argument for option -r.");
                        }
                        if (cutoff_type != NONE) {
                                option_err_msg("Cannot specify both -r and -w options.");
                        }
                        cutoff_type = RANK;
                        rank_cutoff = tmp_long;
                        break;
                case 'w':
                        tmp_dbl = strtod(optarg, &tailptr);
                        if (tailptr == optarg || tmp_dbl < 0.0) {
                                option_err_msg("Invalid argument for option -w.");
                        }
                        if (cutoff_type != NONE) {
                                option_err_msg("Cannot specify both -r and -w options.");
                        }
                        cutoff_type = MIN_VALUE;
                        weight_cutoff = tmp_dbl;
                        break;
                case 'x':
                        if (num_excluded >= max_excluded) {
                                max_excluded *= 2;
                                excluded_namespaces = realloc_(excluded_namespaces,
                                                               max_excluded * sizeof(const char *));
                        }
                        excluded_namespaces[num_excluded++] = optarg;
                        break;
                case 'a':
                        use_all_weights = 1;
                        break;
                case 'T':
                        term_id = optarg;
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
                case 'W':
                        print_warnings = 1;
                        break;
                case 'U':
                        print_unknown_ids = 1;
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
        if (statistics_type == FISHER_EXACT && cutoff_type == NONE) {
                option_err_msg("Must choose a cutoff using -r or -w when requesting"
                               " Fisher's Exact test.\n");
        }

	/* Get weights */
	weights_filename = argv[optind++];

	/* Get first database (namespace is optional) */
	get_namespace_and_file(argv[optind++], &namespace, &tdb_filename, ':');
	if (namespace == NULL) {
                /* Assume first termdb file is ETD */
                ETD_enrichment_context(tdb_filename, &entity_db, &term_db, &mapping_db,
                                       excluded_namespaces, num_excluded);
                db_name = ((ETDTermDb *)term_db)->db_name;
	}
        else {
                /* Assume first termdb file is GMT */
                GMT_enrichment_context(tdb_filename, namespace, &entity_db, &term_db, &mapping_db);
        }

	/* From second database, all termdbs are GMT and namespace is compulsory */
	for (;optind < argc; optind++) {
		get_namespace_and_file(argv[optind], &namespace, &tdb_filename, ':');
		if (namespace == NULL) {
			option_err_msg("Specifying namespace is mandatory for second "
				       "and subsequent term databases.");
		}
		GMT_enrichment_context(tdb_filename, namespace, &entity_db, &term_db, &mapping_db);
	}


        cntxt = EnrichContext_init(db_name, min_term_size, Evalue_cutoff,
                                   effective_db_size, statistics_type,
				   transform_type, discretized_weights,
				   cutoff_type, rank_cutoff, weight_cutoff,
                                   use_all_weights);

        EnrichResults_load_weights(cntxt, weights_filename, entity_db, mapping_db);

        EnrichResults_process_weights(cntxt);

        if (term_id != NULL) {
                term_index = term_db->get_index_from_term_id(term_db, (char *) term_id);
                if (term_index < 0) {
                        fprintf(stderr, "Could not retrieve term with ID %s.\n",
                                term_id);
                        exit(EXIT_FAILURE);
                }
                EnrichResults_calc_single_pvalue(cntxt, term_db, mapping_db, term_index);
                switch (output_type) {
                case TEXT:
                        EnrichResults_print_term_text(cntxt, fp, entity_db, mapping_db);
                        break;
                case TAB:
                        EnrichResults_print_term_tabsep(cntxt, fp, entity_db, mapping_db);
                        break;
                }
        }
        else {
                EnrichResults_calc_pvalues(cntxt, term_db, mapping_db);
                switch (output_type) {
                case TEXT:
                        EnrichResults_print_all_text(cntxt, fp, print_warnings,
                                                     print_unknown_ids);
                        break;
                case TAB:
                        EnrichResults_print_all_tabsep(cntxt, fp, term_db, mapping_db);
                        break;
                }
        }
        if (output_filename != NULL) {
                fclose(fp);
        }
	EnrichContext_delete(cntxt);
        free(excluded_namespaces);
	return EXIT_SUCCESS;
}
