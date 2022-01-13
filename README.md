These are the source files for our paper on Celtis genetic diversity using
SSRs. Note that final formatting and additional edits were applied directly
to the submitted manuscript in Word doc format. Those changes aren't
reflected here.

# File List

- **celtis.Rmd**: The manuscript source, in RMarkdown (plain text) format
- **plantarum.json**: My bibliography database, required only to insert
    citations into the compiled manuscript.
- **data_prep.R**: The main script used to read and process the data
- **systematic-botany.csl**: The Systematic Botany citation format
    specification. Needed only to format citations and the literature cited
    section of the manuscript.
- **celtis.pdf**: The manuscript, compiled into a pdf for easier reading.
- **data/**:
  - **cssr2018-09-06.txt**: microsatellite raw data, in tab-delimited
    format. Sample.Name: code for each sample; Marker: locus; Allele.1-4:
    size of each allele at that locus.
  - **flow.csv**: raw flow cytometry data as generated from flowPloidy.
    channel: filename, including sample name. countsA, countsB: cell counts
    in each histogram peak. sizeA, sizeB: mean/location of each peak. cvA,
    cbB: coefficient of variation for each peak. rcs: residual Chi Square.
    linearity: linearity parameter for NLS model. pg: estimated nuclear DNA
    content in picograms. ploidy: individual ploidy level (blank =
    diploid). 
  - **popTable.csv**: metadata for each individual. sample: individual
    code. population: population code. species: field ID (not to be
    trusted!). ploidy (ploidy from flow cytometry analysis). Genotype (not
    used). Latitude/Longitude: coordinates of population.
  - **dipAssign.csv**: Admixture assignments from STRUCTURE, in
    tab-delimited format. sample: sample code. pop1, pop2: estimated
    admixture from each population. population: which of the two
    populations the sample was assigned to (i.e., had > 0.7 admixture
    proportion from that species). "mix" indicates hybrids.
  - **struc-out.csv**: STRUCTURE probabilities for identifying optimal K
    value. 
  - **celtis_morph.csv**: morphological data. Some of the variables are
    coded in differenet ways in different columns. See the code in
    `data_prep.R` for the correct way to import it
- **ssr-screening/**:
  - **ssr-screening.pdf**: Our manuscript describing primer development for
    this project. APPS no longer publishes primer notes, so we don't have
    somewhere to publish it.
- **structure-all**: directory for files associated with the STRUCTURE
  analysis of the complete (mixed ploidy) data.
  - **celtis-structure.org**: Details of the STRUCTURE analysis. This is an
    `orgmode` file. This is a plain text format, you can open and read it
    in any text editor. If that editor is `Emacs`, you'll have extra
    feature available to you, but this is not necessary to read the text
    and review the code.
  - **mainparams** and **extraparams**: run options for STRUCTURE
  - **structure.sh**: the script used to run STRUCTURE, looping over all
    the replicates for each K value.
  - **all2_1.csv**, **all3_1.csv** and **all4_1.csv**: The
    admixture/assignment results for one replicate of the STRUCTURE runs
    for K = 2, 3, and 4.
  - **ksummary.csv**: "Estimated Ln Prob of Data", extracted from STRUCTURE
    output and used to calculate Evanno's Delta values.
