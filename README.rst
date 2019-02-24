.. _PyVCF: https://github.com/jamescasbon/PyVCF
.. _SnpEff: http://snpeff.sourceforge.net/
.. _AnnoVar: http://annovar.openbioinformatics.org/en/latest/
.. _1kgRef: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/
.. _Exac: http://exac.hms.harvard.edu/
.. _USCSTableBrowser: http://genome.ucsc.edu/cgi-bin/hgTables
.. _BioPython:


Challenge Notes
===============
Few of my questions and concerns regarding the challenge statement:

- table vs annotated VCF::

  [..] output a table annotating each variant in the file.  upload [..] along with the annotated VCF file [..]

  Assuming that both are required.


- type of variation::

  (Substitution, Insertion, Silent, Intergenic, etc.) [..] annotate with the most deleterious possibility

  This list is likely intentionally vague. Genetics are somewhat mixed with genomics here - variant type and function can be orthogonal.
  Will try to make some safe assumptions and implement just some of the basic functional stratification possible.

Visual inspection of given VCF
------------------------------
somatic
  appears to be such from the `normal` and `vaf5` samples. `SOMATIC` record tag is not explicitly used, however.

caller
  freebayes


Implementation Notes
====================

Exploration
-----------

Basic Stats in Jupyter
^^^^^^^^^^^^^^^^^^^^^^
Parsed sample genotypes using PyVCF_ into a pandas dataframe:

pass/filter
  all 6977 genotypes for both `normal` and `vaf5` samples are *passing*.

ploidy
  all 6977 genotypes for both `normal` and `vaf5` samples are *triploid*.


Reference Sources
-----------------
Downloaded reference data from established sources.

genome (fasta)
  1kgRef_ (per vcf header ``##reference=/data/human_g1k_v37.fasta``)

genes (bed)
  USCSTableBrowser_ (RefSeq)


