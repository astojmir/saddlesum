Additional Information
======================

License
-------

All code for *SaddleSum* and *CytoSaddleSum* written at the NCBI is
released into Public Domain. The licenses of external components are
indicated in the source packages.

Reference
---------

If you found *SaddleSum* useful, please cite:

  * A. Stojmirovic and Y.-K. Yu. `Robust and accurate data enrichment
    statistics via distribution function of sum of
    weights. </pmc/articles/PMC2958744/>`_
    *Bioinformatics*, **26** (21):2752-2759, 2010.

Credits
-------

  * Aleksandar Stojmirovic and Yi-Kuo Yu for designed the study, conducted the research and write the paper;
  * Aleksandar Stojmirovic wrote most of the code for the web and the
    standalone version and a part of the Cytoscape plugin;
  * Alexander Bliskovsky wrote most of the Cytoscape plugin as well as
    the code to import KEGG datasets and to download and maintain all term databases.

.. _acknowledge-label:

Acknowledgments
---------------

Standalone SaddleSum uses C code from the Cephes library by Stephen
L. Moshier and a hashtable library by Christopher Clark.

*SaddleSum* web server is written in the
`Python <http://www.python.org/>`_ programming language and
relies on several open-source components:

  * `Numpy and Scipy <http://www.scipy.org/>`_ libraries of scientific tools for Python,
  * `Graphviz <http://www.graphviz.org/>`_ graph vizualization software,
  * `Sphinx <http://sphinx.pocoo.org/index.html>`_ Python documentation generator.
  * `Jinja2 <http://jinja.pocoo.org/2/>`_ template engine.

The web browser client code uses various Javascript routines from the
NCBI and elsewhere. ECMAScrpit code and widgets for the SVG 'Network Navigator'
were taken from `Carto:Net <http://www.carto.net/>`_

All discrete network image color schemes were taken from the
`www.ColorBrewer.org <http://www.ColorBrewer.org>`_ site by
Cynthia A. Brewer, Geography, Pennsylvania State University.

*SaddleSum* term datasets are based on the
`Gene Ontology <http://www.geneontology.org/>`_ (associations
downloaded from the
`NCBI FTP site <ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/>`_) and the
`KEGG Pathways Database <http://www.genome.jp/kegg/pathway.html>`_.

*CytoSaddleSum* uses
`Apache HttpComponents <http://hc.apache.org/>`_ library for HTTP
requests.

We are grateful to Zvezdana Stojmirovic for help with the graphic
design of the *SaddleSum* web pages and interfaces.


News and Updates
----------------

  * **02-Mar-2017.** The documentation was updated and corrected.

  * **14-Dec-2011.**  Updated *CytoSaddleSum* to version 1.5. The new
    version allows users to export results into files in plain text
    or TAB format and to import them from TAB files. The web form
    now alows results to be recieved in TAB format as well. Standalone
    *SaddleSum* 1.5.0 is identical to *SaddleSum* 1.4.0 except for
    the documentation.

  * **19-Sep-2011.**  *CytoSaddleSum* (version 1.4), a Cytoscape interface to *SaddleSum* released. JAR file and source code available at available at the `FTP site <ftp://ftp.ncbi.nlm.nih.gov/pub/qmbpmn/CytoSaddleSum/>`_.

  * **19-Sep-2011.**  *SaddleSum* 1.4.0 released (source code only). Standalone version source code available at the `FTP site <ftp://ftp.ncbi.nlm.nih.gov/pub/qmbpmn/SaddleSum/>`_.

  * **11-Jul-2011.**  *SaddleSum* 1.3.0 released. Standalone version source code and term databases available at the `FTP site <ftp://ftp.ncbi.nlm.nih.gov/pub/qmbpmn/SaddleSum/>`_.
  * **30-Jul-2010.**  Standalone *SaddleSum* source code available at the the `FTP site <ftp://ftp.ncbi.nlm.nih.gov/pub/qmbpmn/SaddleSum/>`_.


..
   Local Variables:
   mode: rst
   indent-tabs-mode: nil
   sentence-end-double-space: t
   fill-column: 70
   End:
