.. _PyVCF: https://github.com/jamescasbon/PyVCF
.. _hgvs: https://github.com/biocommons/hgvs
.. _seqrepo: https://github.com/biocommons/biocommons.seqrepo
.. _ExAC: http://exac.hms.harvard.edu/
.. _SnpEff: http://snpeff.sourceforge.net/
.. _AnnoVar: http://annovar.openbioinformatics.org/en/latest/

Overview
========



Installation
============
.. code-block:: sh

    $ python setup.py install


Usage
=====
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
  Will try to make some safe assumptions and implement just some of the basic functional stratification possible.

Visual inspection of given VCF
------------------------------
* somatic
    appears to be such from the `normal` and `vaf5` samples. but `SOMATIC` record tag is not explicitly used.

* caller
    freebayes


Implementation Notes
====================

