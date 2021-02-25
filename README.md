These are the source files for our paper on Celtis genetic diversity using SSRs.

# File List

- **celtis.Rmd**: The manuscript source, in RMarkdown (plain text) format
- **plantarum.json**: My bibliography database, required only to insert
    citations into the compiled manuscript.
- **data_prep.R**: The main script used to read and process the data
- **american-journal-of-botany.csl**: The American Journal of Botany
  citation format specification. Needed only to format citations and the
  literature cited section of the manuscript.
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
