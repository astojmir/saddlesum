.. _saddlesum-show-etd-label:

saddlesum-show-etd
==================

.. program:: saddlesum-show-etd


SYNOPSIS
--------

::

   saddlesum-show-etd [options] <term_db>

OPTIONS
-------

Arguments
^^^^^^^^^

.. cmdoption:: <term_db>

   A database in ETD format.

Generic options
^^^^^^^^^^^^^^^

.. cmdoption:: -h

   Print a description of all command line options.


.. cmdoption:: -V

   Print the version number and exit. The version number always
   matches that of the *qmbpmn-tools*.

Output options
^^^^^^^^^^^^^^

.. cmdoption:: -N

   Only list the namespaces present in the database, one per
   line. Otherwise, a more detailed info table is output.

   .. note::

      Specify -Ftab to obtain only the list of namespaces, without any
      other text.

.. cmdoption:: -O <output_file>

   Output results to ``<output_file>`` instead of to the standard output.

.. cmdoption:: -F <output_format>

   Select output format. The argument must be one of the following:

   ``txt`` (*default*)
     print results as formatted (pretty) text.
   ``tab``
     print results as a  tab-delimited file. Different sections are
     separated by heading lines starting with ``#`` character.


..
   Local Variables:
   mode: rst
   indent-tabs-mode: nil
   sentence-end-double-space: t
   fill-column: 70
   End:
