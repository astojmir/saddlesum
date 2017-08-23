.. _command-line-label:

saddlesum
=========

.. program:: saddlesum

SYNOPSIS
--------

The standalone saddlesum interface takes the following general form::

    saddlesum [options] <weights_file> [<namespace>:]<term_db> [<namespace>:<term_db> ...]

OPTIONS
-------

Arguments
^^^^^^^^^

.. cmdoption:: <weights_file>

   Tab-delimited file with entity ids and weights. If
   ``<weights_file>`` is specified as ``-``, use standard input.

.. cmdoption:: <namespace>:<term_db>

   The text up to the first column denotes the namespace label while
   the rest of the argument denotes a path to a term database. For the
   first database, the namespace component can be omitted, in which
   case the file is assumed to be in ETD format. All subsequent
   databases should be in GMT format and have namespace specified.
   This allows combining multiple databases to obtain joint
   results. Significant terms for each namespace are treated separately.

Generic options
^^^^^^^^^^^^^^^

.. cmdoption:: -h

   Print a description of all command line options.


.. cmdoption:: -V

   Print the version number and exit. The version number always
   matches that of the *qmbpmn-tools*.


Statistical options
^^^^^^^^^^^^^^^^^^^

.. cmdoption:: -m <min_term_size>

   Set the minimum number of entities for a term to be considered
   (*default* = 2). Only entities with supplied weights count towards
   the term size.

.. cmdoption:: -e <Evalue_cutoff>

   Set the largest E-value for a term to be considered *significant*
   (*default* = 1e-02).

   .. note::

      Since saddlesum uses an algorithm that quickly rejects those
      terms that cannot have an E-value smaller than ``<Evalue_cutoff>``,
      the choice of ``<Evalue_cutoff>`` can significantly affect the
      running time of the program.

.. cmdoption:: -n <effective_db_size>

   Set the effective term database size for applying Bonferroni
   correction (i.e. calculating E-values) to P-values output by the
   algorithm (*default* = the total number of terms
   considered). Setting ``<effective_db_size>`` to 1.0 will result in the
   original P-values being returned.

.. cmdoption:: -s <statistics_type>

   Set statistical method used to evaluate P-values. The argument must
   be one of the following:

   ``wsum`` (*default*)
     use SaddleSum statistics based on sum-of-weights score for each term
   ``hgem``
     use one-sided Fisher's Exact test.

   .. note::

      Option -s hgem implies option -d (whether specified or
      not). Also, one of -r or -w options must be used to specify the
      cutoff for selecting weights.

.. cmdoption:: -a

   When this option is specified, all weights that can be mapped to
   valid entities are used as statistical background. Otherwise, only
   those weights that both map to a valid entity and a vocabulary term
   in the term database are used.

.. cmdoption:: -x <namespace>

   Exclude ``<namespace>`` from ETD database. Each ETD database may
   contain multiple namespaces. This option allows specific namespaces
   to be totally excluded from consideration. It affects the effective
   database and hence the term E-values. More than one namespace can
   be excluded by using this option multiple times.

   .. note::

      Use saddlesum-show-etd program to discover the
      names of all namespaces present in an ETD file.

.. cmdoption:: -T <term_id>

   Compute statistics only for the term with ID ``<term_id>`` and
   display the list of entities associated with that term, together
   with their weights. All other statistical and weight processing
   options can be specified in conjunction with -T but -e has no
   effect.

   .. note::

      If specifying -m results in the term to be excluded from
      computation of P-value, no statistics will be computed and displayed.


Weight processing options
^^^^^^^^^^^^^^^^^^^^^^^^^
.. cmdoption:: -t <weight_transformation>

   Apply a transformation to each of the provided weights prior to
   other applying other processing options (see below) and
   calculating enrichment statistics. The argument must be one of the
   following:

   ``flip``
     flip the sign of each weight
   ``abs``
     take the absolute value of each weight.

   When this option is omitted, no transformation will be performed.

.. cmdoption:: -r <rank_cutoff>

   Set all weights ranked lower than ``<rank_cutoff>`` to 0. If there
   are several weights tied at ``<rank_cutoff>``, keep all of them.

.. cmdoption:: -w <weight_cutoff>

   Set all weights smaller than ``<weight_cutoff>`` to 0.

.. note::

   Only one of the ``-r`` and ``-w`` options can be set at the same time.

.. cmdoption:: -d

   Discretize weights. Set all weights greater than 0 to 1 and all
   those smaller than 0 to 0.

.. note::

   The weight processing options are applied in this order: ``-t``,
   then ``-r`` or ``-w`` and finally ``-d``.

Output options
^^^^^^^^^^^^^^

.. cmdoption:: -O <output_file>

   Output results to ``<output_file>`` instead of to the standard output.

.. cmdoption:: -F <output_format>

   Select output format. The argument must be one of the following:

   ``txt`` (*default*)
     print results as formatted (pretty) text.
   ``tab``
     print results as a  tab-delimited file. Different sections are
     separated by heading lines starting with ``#`` character.

.. cmdoption:: -W

   Print warnings about entity identifiers from the
   ``<weights_file>``.

.. cmdoption:: -U

   Print ids from ``<weights_file>`` that are not present in the term
   databases.

.. note::

   Options -W and -U apply only to text output. Tab-delimited output
   always contains sections containing warnings and unknown ids.



..
   Local Variables:
   mode: rst
   indent-tabs-mode: nil
   sentence-end-double-space: t
   fill-column: 70
   End:
