.. _PyVCF: https://github.com/jamescasbon/PyVCF
.. _hgvs: https://github.com/biocommons/hgvs
.. _seqrepo: https://github.com/biocommons/biocommons.seqrepo
.. _ExAC: http://exac.hms.harvard.edu/
.. _SnpEff: http://snpeff.sourceforge.net/
.. _AnnoVar: http://annovar.openbioinformatics.org/en/latest/

Overview
========
Several avenues for implementation were considered:

1. grab as much as possible from the vcf itself and run a batch request (``/rest/bulk/variant``) for variant data on ExAC_.
    this would have been fast and relatively easy to implement, but would not annotate most of the non-SNP variation,
    since ExAC_ is largely a view of dbsnp.

2. just run the vcf through snpeff_ or annovar_.
    cheating ?

3. try and annotate variant consequences myself by investigating effects on transcription and translation.
    avoid the weeds of genomic to tx projections and codon tables by using the hgvs_ package.


This solution implements option (3).


Installation
============
    ``python >= 3.7``

.. code-block:: sh

    $ python setup.py install

Performance of hgvs_ can be greatly improved by installing a local instance of seqrepo_ ( takes ~20 min )

.. code-block:: sh

    $ pip install seqrepo
    $ seqrepo init
    $ seqrepo pull -i 2019-06-20


Usage
=====
If seqrepo_ installed:

.. code-block:: sh

    $ export HGVS_SEQREPO_DIR=/usr/local/share/seqrepo/2019-06-20

Example call used to obtain ``annotations.csv``:

.. code-block:: sh

    $ python -m tempus annotate test/data/Challenge_data.vcf annotations.csv --workers 5 --chunk-size 1000


Challenge Notes
===============
Few of my questions and concerns regarding the challenge statement:

- table vs annotated VCF::

  [..] output a table annotating each variant in the file.  upload [..] along with the annotated VCF file [..]

  Assuming that CSV should suffice.


- type of variation::

  (Substitution, Insertion, Silent, Intergenic, etc.) [..] annotate with the most deleterious possibility

  This list is likely intentionally vague. Genetics are somewhat mixed with genomics here - variant type and
  function can be orthogonal.
  Will try to make some safe assumptions and implement just some of the basic functional nomenclature possible.

Implementation Notes
====================
* Parallelism via splitting the vcf into chunks and letting separate processes chug away is rather inelegant.
    Threads or coroutines would make more sense for the network-bound IO characteristic of this annotator.
    But seqrepo_ is not thread-safe and I didn't feel like digging too deep into the issues there.
    Coroutines could certainly be considered, but again I just wanted a quick and dirty solution and keep the library
    APIs maximally transparent.

* Some extra effort went into 5'-normalizing (left shuffle) variants to query ExAC_.
    But the gains from such normalization were never assessed and could well be negligible.

* PyVCF is a bit rough around the edges.
    For example, pickling _Record(s) (objects that represent loci) destroys some of the INFO and GT data inside.
    Writing is also no joke. Perhaps, this is why CVS table is good enough for now :)

* jupyter notebooks found in ``nb/`` are merely scratch paper.
