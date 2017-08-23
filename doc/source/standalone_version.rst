.. _standalone-label:

Standalone Program
==================

The standalone version of SaddleSum is written in the C programming
language and requires no additional libraries. It offers a standard
UNIX command line interface and allows users to specify multiple term
databases in the GMT format used by GSEA as well as in SaddleSum's own
ETD (Extended Term Database) format (see :ref:`etd-label` below).


Downloading
-----------

The source code for the standalone *SaddleSum* is available from the NCBI FTP site
(ftp://ftp.ncbi.nih.gov/pub/qmbpmn/SaddleSum/src) as a tar.gz archive. Releases
of *SaddleSum* share the version numbering with the corresponding
**qmbpmn-tools**.


Building and installing executables
-----------------------------------

After downloading a source code, extract it to a chosen directory and enter the
``saddlesum-<version>`` directory. The standard process of::

  ./configure
  make
  make install

should install the binaries into your executable directory. To clean the build,
type ``make clean``. This requires gcc and GNU make.

We have successfully built the source on Linux systems, Mac OS X with gcc and on
Windows XP using MinGW. It may be possible to build it on other platforms but we
have not attempted to do so.

.. _term-datasets-label:

Obtaining Term Datasets
-----------------------

Term datasets corresponding to the three Gene Ontology (GO) categories
and KEGG pathways are available for a number of species in both ETD and GMT
format (gzipped). ETD datasets are recommended but GMT versions of the
same data is provided for compatibility with the earlier version. The
FTP link is ftp://ftp.ncbi.nih.gov/pub/qmbpmn/SaddleSum/term_datasets . The
most recently built datasets are available under in the *current* subdirectory.

The naming of the GMT files follows the convention::

  <database>-<species>.gmt.gz

The <database> field is either ``GO-<category>`` or ``KEGG``, where ``category``
is one of ``mf`` (molecular function), ``bp`` (biological process) or ``cc``
(cellular component). The species field uses KEGG three-letter convention. The
currently available species are:

========== =========================
KEGG code  Species
========== =========================
ath        Arabidopsis thaliana
cel        Caenorhabditis elegans
dre        Danio rerio
dme        Drosophila melanogaster
hsa        Homo sapiens
mmu        Mus musculus
pfa        Plasmodium falciparum
rno        Rattus norvegicus
sce        Saccharomyces cerevisiae
========== =========================

For example ``GO-bp-pfa.gmt.gz`` contains the GO biological process database for
*Plasmodium falciparum*.

Each ETD file contains all provided term datasets for a single species
(the three GO namespaces plus KEGG). The naming of the ETD files
follows the convention::

  <species>.etd.gz

The species field again uses the KEGG three letter codes shown above.

Using SaddleSum
---------------

Standalone *SaddleSum* has a simple command line interface::

    saddlesum [options] <weights_file> [<namespace>:]<term_db> [<namespace>:<term_db> ...]

For full list of options, type::

    saddlesum -h

or see the :ref:`command-line-label` page in this manual.

The ``<weights_file>`` argument can either be a filename path or ``-``, which
means that the weights are to be taken from the command line. The weights are
supplied in the two-column tab-delimited format, where the first column gives
entity (gene) IDs, while the second gives their weights as floating point
numbers. No comment lines or multiple columns are allowed.

Entity IDs depend on the term database being searched. Each ETD database
contains, in addition to the term datasets, an abbreviated version of
the NCBI Gene database for its species. Thus, when using an ETD
database, SaddleSum is able to interpret entity IDs provided as Gene
IDs, as gene symbols and as gene aliases. If provided symbol is
ambiguous, SaddleSum reports a warning or an error, depending on
whether the ambiguity can be resolved.

The GMT term databases supplied on our FTP site use NCBI Gene IDs to
label genes and hence the weights supplied to ``saddlesum`` program
must also use NCBI Gene IDs to supply weights. We have decided on
using NCBI Gene IDs for this purpose because many genes have several
widely-used names, while GMT format does not allow specifying aliases.

The database arguments are in the format ``<namespace>:<term_db>``
where ``<namespace>`` is an arbitrary namespace label and
``<term_db>`` is the path to the term database. The first database can
be specified without the ``<namespace>`` label, which indicates that
the ``<term_db>`` is an ETD file. If ``<namespace>`` is present, a GMT
dataset is assumed. All subsequent arguments must specify the
namespace part and the term datasets must be in GMT format. In this
way, it is possible to combine an ETD database with several GMT
ones. The recommended use for GMT databases is to create customized
databases containing associations that cannot be found in GO or KEGG.

Namespace labels for GMT databases can be repeated. Significant terms
under the same namespace are sorted and printed together, with each
results table headed by the ``<namespace>`` label. Term datasets
within an ETD database already contain their namespace labels, which
can be inspected using the :ref:`saddlesum-show-etd-label` utility.


.. _standalone-examples-label:

Examples
--------

The ``examples/`` subdirectory contains several examples of weights
from microarrays (log2 ratios) and term databases. All the examples
were taken from the NCBI GEO database and have the extension
``.tab``. The database files are in GMT format and have extension
``.gmt``. Copy ``saddlesum`` executable to your path and try the
following examples::

    saddlesum GDS3184_GSM253560_up_geneid.tab test:GO-bp-rno.gmt
    saddlesum GDS3184_GSM253560_down_geneid.tab test:GO-bp-rno.gmt

    saddlesum GDS2338_GSM102668_up_geneid.tab test:KEGG-sce.gmt
    saddlesum GDS2338_GSM102668_down_geneid.tab test:KEGG-sce.gmt

    saddlesum GDS2352_GSM89756_up_geneid.tab test:GO-cc-hsa.gmt
    saddlesum GDS2352_GSM89756_down_geneid.tab test:GO-cc-hsa.gmt

One can also experiment with options, such as setting the cutoffs and
similar. Note that the down-regulated weights were obtained by flipping the sign
on the original log2 ratios and setting all weights smaller than zero to zero
(please see our paper for the explanation).


Manpages
--------

.. toctree::
   :maxdepth: 1

   saddlesum-cli.rst
   saddlesum-show-etd.rst

.. _etd-label:

Structure of Extended Term Databases
------------------------------------

Extended Term Databases (ETDs) are binary files that contain term
databases used by *SaddleSum*. They are created through Python scripts
from *qmbpmn-tools*. Here, we describe the structure of the binary
format.

Each ETD consists of a header, genes database and one or more
*namespaces*. The header contains the overall information about the
database. Genes database is derived from NCBI Gene and serves to
enable parsing of gene label aliases. Each namespace is a separate
term database, containing terms, their relationships and mapping to genes.

.. note:: In the listing that follow we use C-like pseudocode to
   indicate the types of the fields. These resemble C declarations with
   a difference that array sizes are usually not fixed but need to be
   read from the file. The type ``uint32`` here means 32-bit
   little-endian integer.


Header
^^^^^^

The header part is written as following (comments are on the right)::

    char start_separator[8]        - always 'EXTERMDB'
    uint32 version_magic           - set to 1644632861 for the current version
    uint32 db_name_buflen
    char db_name[db_name_buflen]   - Full database name (description)
    uint32 num_namespaces
    uint32 namespaces_buflen
    char namespaces_buf[namespaces_buflen] - Names of all namespaces


Genes Database
^^^^^^^^^^^^^^

The genes database is derived from NCBI Gene records for a given
species. It contains a list of genes, their descriptions, their
symbols and lists of conflicts. There are two types of conflicts. Type
1 conflict occurs when the canonical gene symbol for a gene is an
alias of another gene. These are always resolved to mean the gene that
has the conflicting label as the canonical symbol. Type 2 conflict
occurs when a label is an alias of multiple genes but is not set as a
canonical symbol for any gene. Such conflicts cannot be resolved. The
genes database structure is::

    char start_separator[8]        - always 'NCBIGENE'
    uint32 version_magic           - set to 1200900292 for the current version
    uint32 metadata_buflen
    char metadata[metadata_buflen] - gene_info_file, NCBIGENE_URL_FMT
    uint32 gene_info_checksum      - CRC32 checksum of gene_info_file
    uint32 tax_id                  - NCBI Taxonomy ID

    uint32 N                       - number of genes
    uint32 gene_ids[N]             - assumes every NCBI Gene ID is uint32
    uint32 offsets[N+1]            - offsets into gene_info_file
    uint32 symbols_counts[N]       - counts of symbol lists for each gene
    uint32 symbols_buflen
    char symbols_buf[symbols_len]  - N lists of symbols. Each symbol is
                                     terminated by '\0'
                                     Each list contains one primary
                                     identifier and 0 or more synonyms.
    uint32 desc1_buflen
    char desc1[desc1_buflen]       - gene descriptions

    uint32 num_conf1               - number of conflicts (type1)
    uint32 conf1_counts[num_conf1] - counts of conflict lists (type 1)
    uint32 conf1_buflen
    char conf1_buf[conf1_buflen]   - A list of conflict lists (of type 1).
                                     Each symbol terminated by '\0'.
                                     The first item in each list conflicts
                                     with all the others.

    uint32 num_conf2               - number of conflicts (type2)
    uint32 conf2_counts[num_conf2] - counts of conflict lists (type 2)
    uint32 conf2_buflen
    char conf2_buf[conf2_buflen]   - A list of conflict lists (of type 2).
                                     Each symbol terminated by '\0'.
                                     The first item in each list conflicts
                                     with all the others.

Term Namespaces
^^^^^^^^^^^^^^^

Each namespace is a separate term database. It contains terms (IDs and
descriptions) their relationships and maps to genes. KEGG and Gene
Ontology databases differ according to their metadata. The structure is::

    char start_separator[8]        - always 'TERMDBNS'
    uint32 version_magic           - changes for each different version
                                     Currently 2264738403 for KEGG,
                                               2187050528 for Gene Ontology.
    uint32 num_edgetypes           - number of term relationship types
    uint32 edgetype_buflen
    char edgetype_names_buf[edgetype_buflen] - names for term
                                     relationships. Each name terminated by '\0'.
    uint32 M                       - number of terms
    uint32 slim_flags[M]           - whether the terms are in reduced dataset
                                     (not currently used)
    uint32 num_hits[M]             - counts of hits for each term
    uint32 hits[M][num_hits]       - hits to genes. These map each term index k
                                     to an array of gene indices, of length
                                     num_hits[k]
    uint32 termid_buflen
    char termid_buf[termid_buflen] - term IDs. Each item terminated by '\0'.
    uint32 desc_buflen
    char desc_buf[desc_buflen]     - term descriptions, Each item terminated by '\0'.
    uint32 num_parents[M]          - numbers of parents of each term
    uint32 parents[M][num_parents] - term indices of parents
    uint32 edgetypes[M][num_parents]  - relationships to parents (edgetype indices)
    uint32 metadata_buflen
    char metadata[metadata_buflen] - metadata (provided by the caller).
                                     List of strings, each string terminated by '\0'.

KEGG namespaces use three items in metadata field: organism prefix, leaf URL format
(URL format for the page describing a leaf term) and higher URL format
(URL format for terms higher in the hierarchy). Gene Ontology
namespaces have a single item; a url format for Amigo website. All URL
formats must contain a single ``%s`` formatting specifier.



..
   Local Variables:
   mode: rst
   indent-tabs-mode: nil
   sentence-end-double-space: t
   fill-column: 70
   End:
