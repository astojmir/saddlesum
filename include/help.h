/* This file was automatically generated from RST docs */ 

#define HELP_SADDLESUM "" \
"SYNOPSIS\n" \
"\n" \
"   The standalone saddlesum interface takes the following general form:\n" \
"\n" \
" saddlesum [options] <weights_file> [<namespace>:]<term_db> [<namespace>:<term_db> ...]\n" \
"\n" \
"OPTIONS\n" \
"\n" \
"  Arguments\n" \
"\n" \
"   <weights_file>\n" \
"\n" \
"           Tab-delimited file with entity ids and weights. If <weights_file>\n" \
"           is specified as -, use standard input.\n" \
"\n" \
"   <namespace>:<term_db>\n" \
"\n" \
"           The text up to the first column denotes the namespace label while\n" \
"           the rest of the argument denotes a path to a term database. For\n" \
"           the first database, the namespace component can be omitted, in\n" \
"           which case the file is assumed to be in ETD format. All subsequent\n" \
"           databases should be in GMT format and have namespeace specified.\n" \
"           This allows combining multiple databases to obtain joint results.\n" \
"           Significant terms for each namespace are treated separately.\n" \
"\n" \
"  Generic options\n" \
"\n" \
"   -h\n" \
"\n" \
"           Print a description of all command line options.\n" \
"\n" \
"   -V\n" \
"\n" \
"           Print the version number and exit. The version number always\n" \
"           matches that of the qmbpmn-tools.\n" \
"\n" \
"  Statistical options\n" \
"\n" \
"   -m <min_term_size>\n" \
"\n" \
"           Set the minimum number of entities for a term to be considered\n" \
"           (default = 2). Only entites with supplied weights count towards\n" \
"           the term size.\n" \
"\n" \
"   -e <Evalue_cutoff>\n" \
"\n" \
"           Set the largest E-value for a term to be considered significant\n" \
"           (default = 1e-02).\n" \
"\n" \
"           Note\n" \
"\n" \
"           Since saddlesum uses an algorithm that quickly rejects those terms\n" \
"           that cannot have an E-value smaller than <Evalue_cutoff>, the\n" \
"           choice of <Evalue_cutoff> can significantly affect the running\n" \
"           time of the program.\n" \
"\n" \
"   -n <effective_db_size>\n" \
"\n" \
"           Set the effective term database size for applying Bonferroni\n" \
"           correction (i.e. calculating E-values) to P-values output by the\n" \
"           algorithm (default = the total number of terms considered).\n" \
"           Setting <effective_db_size> to 1.0 will result in the original\n" \
"           P-values being returned.\n" \
"\n" \
"   -s <statistics_type>\n" \
"\n" \
"           Set statistical method used to evaluate P-values. The argument\n" \
"           must be one of the following:\n" \
"\n" \
"                wsum (default)\n" \
"                        use SaddleSum statistics based on sum-of-weights\n" \
"                        score for each term\n" \
"\n" \
"                hgem\n" \
"                        use one-sided Fisher’s Exact test.\n" \
"\n" \
"           Note\n" \
"\n" \
"           Option -s hgem implies option -d (whether specified or not). Also,\n" \
"           one of -r or -w options must be used to specify the cutoff for\n" \
"           selecting weights.\n" \
"\n" \
"   -a\n" \
"\n" \
"           When this option is specified, all weights that can be mapped to\n" \
"           valid entities are used as statistical background. Otherwise, only\n" \
"           those weights that both map to a valid entity and a vocabulary\n" \
"           term in the term database are used.\n" \
"\n" \
"   -x <namespace>\n" \
"\n" \
"           Exclude <namespace> from ETD database. Each ETD database may\n" \
"           contain multiple namespaces. This option allows specific\n" \
"           namespaces to be totally excluded from consideration. It affects\n" \
"           the effective database and hence the term E-values. More than one\n" \
"           namespace can be excluded by using this option multiple times.\n" \
"\n" \
"           Note\n" \
"\n" \
"           Use saddlesum-show-etd program to discover the names of all\n" \
"           namespaces present in an ETD file.\n" \
"\n" \
"   -T <term_id>\n" \
"\n" \
"           Compute statistics only for the term with ID <term_id> and display\n" \
"           the list of entities associated with that term, together with\n" \
"           their weights. All other statistical and weight processing options\n" \
"           can be specified in conjunction with -T but -e has no effect.\n" \
"\n" \
"           Note\n" \
"\n" \
"           If specifying -m results in the term to be excluded from\n" \
"           computation of P-value, no statistics will be computed and\n" \
"           displayed.\n" \
"\n" \
"  Weight processing options\n" \
"\n" \
"   -t <weight_transformation>\n" \
"\n" \
"           Apply a transformation to each of the provided weights prior to\n" \
"           other applying other processing options (see below) and\n" \
"           calculating enrichment statistics. The argument must be one of the\n" \
"           following:\n" \
"\n" \
"                flip\n" \
"                        flip the sign of each weight\n" \
"\n" \
"                abs\n" \
"                        take the absolute value of each weight.\n" \
"\n" \
"           When this option is omitted, no transformation will be performed.\n" \
"\n" \
"   -r <rank_cutoff>\n" \
"\n" \
"           Set all weights ranked lower than <rank_cutoff> to 0. If there are\n" \
"           several weights tied at <rank_cutoff>, keep all of them.\n" \
"\n" \
"   -w <weight_cutoff>\n" \
"\n" \
"           Set all weights smaller than <weight_cutoff> to 0.\n" \
"\n" \
"   Note\n" \
"\n" \
"   Only one of the -r and -w options can be set at the same time.\n" \
"\n" \
"   -d\n" \
"\n" \
"           Discretize weights. Set all weights greater than 0 to 1 and all\n" \
"           those smaller than 0 to 0.\n" \
"\n" \
"   Note\n" \
"\n" \
"   The weight processing options are applied in this order: -t, then -r or -w\n" \
"   and finally -d.\n" \
"\n" \
"  Output options\n" \
"\n" \
"   -O <output_file>\n" \
"\n" \
"           Output results to <output_file> instead of to the standard output.\n" \
"\n" \
"   -F <output_format>\n" \
"\n" \
"           Select output format. The argument must be one of the following:\n" \
"\n" \
"                txt (default)\n" \
"                        print results as formatted (pretty) text.\n" \
"\n" \
"                tab\n" \
"                        print results as a tab-delimited file. Different\n" \
"                        sections are separated by heading lines starting with\n" \
"                        # character.\n" \
"\n" \
"   -W\n" \
"\n" \
"           Print warnings about entity identifiers from the <weights_file>.\n" \
"\n" \
"   -U\n" \
"\n" \
"           Print ids from <weights_file> that are not present in the term\n" \
"           databases.\n" \
"\n" \
"   Note\n" \
"\n" \
"   Options -W and -U apply only to text output. Tab-delimited ouput always\n" \
"   contains sections containing warnings and unknown ids.\n" \
"\n" 


#define HELP_SHOW_ETD "" \
"SYNOPSIS\n" \
"\n" \
" saddlesum-show-etd [options] <term_db>\n" \
"\n" \
"OPTIONS\n" \
"\n" \
"  Arguments\n" \
"\n" \
"   <term_db>\n" \
"\n" \
"           A database in ETD format.\n" \
"\n" \
"  Generic options\n" \
"\n" \
"   -h\n" \
"\n" \
"           Print a description of all command line options.\n" \
"\n" \
"   -V\n" \
"\n" \
"           Print the version number and exit. The version number always\n" \
"           matches that of the qmbpmn-tools.\n" \
"\n" \
"  Output options\n" \
"\n" \
"   -N\n" \
"\n" \
"           Only list the namespaces present in the database, one per line.\n" \
"           Otherwise, a more detailed info table is output.\n" \
"\n" \
"           Note\n" \
"\n" \
"           Specify -Ftab to obtain only the list of namespaces, without any\n" \
"           other text.\n" \
"\n" \
"   -O <output_file>\n" \
"\n" \
"           Output results to <output_file> instead of to the standard output.\n" \
"\n" \
"   -F <output_format>\n" \
"\n" \
"           Select output format. The argument must be one of the following:\n" \
"\n" \
"                txt (default)\n" \
"                        print results as formatted (pretty) text.\n" \
"\n" \
"                tab\n" \
"                        print results as a tab-delimited file. Different\n" \
"                        sections are separated by heading lines starting with\n" \
"                        # character.\n" \
"\n" 


