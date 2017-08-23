Overview
========

Term enrichment analysis facilitates biological interpretation by assigning to
experimentally/computationally obtained data annotation associated with terms
from controlled vocabularies. This process usually involves obtaining
statistical significance for each vocabulary term and using the most significant
terms to describe a given set of biological entities, often associated with
weights. Many existing enrichment methods require selections of (arbitrary
number of) the most significant entities and/or do not account for weights of
entities. Others either mandate extensive simulations to obtain statistics or
make assumptions about data distribution that need not be justifiable. In
addition, most methods have difficulty assigning correct statistical
significance to terms with few entities.

Implementing the well-known Lugananni-Rice formula, we have developed a novel
approach, called *SaddleSum*, that is free from all these undesirable
constraints. It approximates the distribution of sum of weights asymptotically
by saddlepoint method to arrive at analytically and computationally tractable
form.

We evaluated *SaddleSum* against several existing methods and demonstrated its
ability to adapt equally well to distributions with widely different
properties. With entity weights properly taken into account, *SaddleSum* is
internally consistent and stable with respect to the choice of number of most
significant entities selected. Making few assumptions on the input data, the
proposed method is universal and can thus be applied to areas beyond analysis of
microarrays, such as deep sequencing, quantitative proteomics and in silico
network simulations. *SaddleSum* also provides a term-size dependent score
distribution function that gives rise to accurate statistical significance even
for terms with few entities. As a consequence, *SaddleSum* enables researchers to
place confidence in its significance assignments to small terms that are often
biologically most specific.


.. _statistics-label:

Enrichment Statistics of Sum-of-weights Scores
----------------------------------------------

Suppose that you have obtained some measurements for a number of genes or
proteins in an organism. These can be gene expression log-ratios or
values computed through some simulation. The problem now is to
biologically interpret your experiment, that is, to link it to known
terms and concepts.

One way to do so is to obtain a *controlled vocabulary* that covers
the domain of your investigation and serves as basis for annotation of
the genes from your organism. This means that every term in that
vocabulary is associated with one or more genes; the term *describes*
its associated genes. Your original problem is now transformed to the
question of finding what terms best describe your overall measurements.

Let us further assume that gene measurements, which we will call *weights*,
indicate relative importance of genes. If that is not so, for example
if all highly expressed genes should be treated equally, you can
transform your values so that this assumption holds and continue from
there. A natural way to express the importance of each vocabulary
term is to compute its *score* as the sum of weights of all genes
mapping to that term.

Unfortunately, the score just calculated is not a good way to compare
terms. Different terms may contain different numbers of genes and
hence cover different parts of your measurement set. Thus, it may be
possible that difference in scores of two terms is only due to their
coverage. Normalizing by the number of genes may help but does not go
far enough; there is still a question of the scale of the scores and
their distribution.

The standard approach to this issue is statistical: we want to know
the probability of a given score occurring by chance with respect to
some null distribution. This probability is known as P-value and can
be used as the *statistical significance* of the term. Note that there
is a null distribution here associated with each term. In our case,
it is reasonable to assume that the scores are built by summing *m*
independently and identically distributed weights, where *m* is the
number of genes mapping to the term. If the term P-value is
sufficiently small, say 0.01, we can conclude that it is very unlikely
that term score has arisen by chance. Therefore, we are able to
interpret our dataset using that term.

*SaddleSum* is a tool that automates the above procedure. It requires
a term database and a collection of gene labels associated with
weights. It outputs the most statistically significant terms together
with their estimated statistical significance. To compute P-values,
*SaddleSum* first estimates the distribution of weights using all the
supplied weights. Then using this estimated distribution, it
approximates the P-value of a score using the asymptotic saddlepoint
approach developed by Lugananni and Rice. This estimate is extremely
accurate, even for terms with very small number of associated genes
or for very small P-values.

To correct for the fact that a term database usually has many terms
that enhance the chance of random matches,
it is necessary to perform *multiple hypothesis testing correction*.
*SaddleSum* uses Bonferroni correction for multiple hypothesis
testing, where raw term P-values are multiplied with the effective
database size to give E-values. We define the effective database size
as the number of terms in the term database that map to at least *k*
genes corresponding to selected nodes, where *k* is the minimum term
size. This parameter is calculated at the time of the query. By
default, the calculated value is used to compute E-values but you can
can override it by changing the corresponding
form entry. For example, raw P-values can be displayed by setting the
effective database size to 1.


SaddleSum software
------------------

*SaddleSum* can be used through three interfaces: command line
(standalone version), web and as a
`Cytoscape <http://www.cytoscape.org/>`_ plugin.

The standalone version is written completely in the C programming language and has
no dependencies beyond the gcc compiler and GNU make (although it is likely that
it can be adapted to compile using other development tools). It provides a
command line interface and uses term databases in a simple tab-delimited format
(GMT) as well as in SaddleSum's own ETD (Extended Term Database) format.

The web version provides a wrapper for the standalone version via a
web form. The server script that interacts with the standalone version
is written in Python programming language and is a part of *qmbpmn-tools*
package, which also includes *ITM Probe*, a framework for analysis of information
flow in interaction networks based on random walks. The *qmbpmn-tools* package
depends on a number of external Python libraries and on the Graphviz suite for
visualizing graphs.

Cytoscape is an open source platform for complex network analysis and
visualization written in Java programming language. Apart for a rich
set of graph visualization tools, it provides an interface for
externally written plugins that provide additional functionality such
as network analysis algorithms, database import and functional
enrichment analysis. Cytoscape users are therefore able to combine
algorithms and data from different sources to perform complex
network-based analyses.

*CytoSaddleSum* is a Cytoscape plugin that enables *SaddleSum* queries
from Cytoscape platforms. It can interact either with a locally
installed command-line program directly, or through our web
server. The results are presented as a Cytoscape network view showing
term relationships. Term statistics, written to node attributes in
the relationship network, can be manipulated by the user.




..
   Local Variables:
   mode: rst
   indent-tabs-mode: nil
   sentence-end-double-space: t
   fill-column: 70
   End:
