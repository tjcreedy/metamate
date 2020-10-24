## metaMATE: your metabarcoding friend!

### Overview

metaMATE (metabarcoding Multiple Abundance Threshold Evaluator) analyses a set of amplicons derived through metabarcoding of a mitochondrial coding locus to determine putative NUMT and other erroneous sequences based on relative read abundance thresholds within libraries, phylogenetic clades and/or taxonomic groupings. 

The paper for metaMATE is [available on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.06.17.157347v1) under the previous name of NUMTdumper. If you use metaMATE in your work, please cite this paper.

The development of this tool was supported by the iBioGen project, funded by the H2020 European Research Council, Grant/Award Number: 810729.

## Table of contents
* [Introduction](#Introduction)
  + [Input data](#input-data-required)
  + [`find`](#find-introduction)
  + [`dump`](#dump-introduction)
* [Installation](#installation)
* [Usage](#usage)
  + [Specifications](#specifications)
  + [Core arguments](#core-arguments)
  + [Clade delimitation and binning arguments](#clade-delimitation-and-binning-arguments)
  + [Taxon binning arguments](#taxon-binning-arguments)
  + [`find`-specific arguments](#find-specific-arguments)
  + [Reference-matching arguments](#reference-matching-arguments)
  + [Length-based arguments](#length-based-arguments)
  + [Translation-based arguments](#translation-based-arguments)
  + [`dump`-specific arguments](#dump-specific-arguments)
* [Examples](#Examples)
  + [`find`](#find-examples)
  + [`dump`](#dump-examples)
* [Outputs](#outputs)
* [Details](#details)
  + [Validation](#validation)
  + [Generating clades](#generating-clades)
  + [ASV assessment](#asv-assessment)
  + [Scoring and estimation](#scoring-and-estimation)
* [Development](#development)

## Introduction

metaMATE takes the guesswork and faff out of applying frequency-based filtering thresholds in NGS amplicon pipelines by utilising a modular threshold specification approach combined with detailed outputs to efficiently provide detailed insight into the effects of different filtering and binning strategies. This approach wase explicitly developed towards removing putative NUMT sequences from a population of ASVs, but the methodology is broad enough that it will work simultaneously on any other types of low-frequency erroneous sequences, such as sequencing errors. A principal benefit of metaMATE over closely related tools is its validation approach, whereby input ASVs are assessed for membership of control groups for known authenticity. The effect of frequency thresholds and binning strategies on retention or rejection of the members of these is used to assess the optimal methodology for NUMT and error removal.

The downside of this approach is that metaMATE is more data- and parameter- hungry than other similar tools. Rather than being content with being fed a set of ASVs, metaMATE requires inputs that a) enable the binning of ASV reads to provide pools upon which frequency thresholds can be applied and b) parameterise the determination of the two control groups. 

metaMATE operates using two submodules, `find` and `dump`, with the former being the main workhorse of the tool. The standard workflow would be as follows:

1. Run `find` on your ASVs with other data and output tabulated results and results cache
2. Independently analyse the tabulated results to decide on the best strategy
3. Run `dump` on the results cache for a given threshold set to extract the filtered ASVs

### Input data required

The ideal dataset for metaMATE is **a set of ASVs** arising from a **multi-sample metabarcoding dataset** accompanied by a solid set of **reference sequences** that are expected to be present in the dataset. Optionally, metaMATE can also utilise data assigning each ASV to a taxonomic group. 

The ASVs (amplicon sequence variants) should comprise all unique sequences found across the dataset, optionally denoised, singletons removed, and filtered to remove outlying length variants. Note that metaMATE currently cannot process more than 65,536 input ASVs if perfoming clade binning due to the exponential complexity of the UPGMA algorithm.

The multi-sample aspect is crucial, as the main novelty and power of metaMATE comes from assessing read counts of ASVs across descrete sampling units, as opposed to just using read counts of ASVs across the dataset as a whole. metaMATE refers to sampling units as 'libraries', under the assumption that a sample has likely been sequenced as a descrete pool of amplicons sharing the same identifying index. Depending on your pipeline, your descrete samples of amplicons may be stored in different files or in one file with sequences identified in the sequence header. metaMATE can deal with either of these ([see below](#-l--libraries-path-path)), assuming the format is correct.

Reference sequences allow metaMATE to identify some of the input ASVs as verified authentic, that is to say definitely *not* NUMTs. The reference sequences are not expected to be comprehensive, nor are all references necessarily expected to occur, but the more ASVs that can be designated as verified authentic, the better metaMATE is able to estimate the impact on filtering and the more accurate the selection of filtering thresholds can be. If necessary, references can be drawn from global databases, such as BOLD or GenBank, or curated versions of these such as MIDORI, but [we suggest](#reference-matching-arguments) that the hit thresholds be more stringent in this case.

### `find` introduction

The purpose of `find` mode is to comprehensively assess a range of frequency filtering specifications to analyse the impact of these on determining NUMTs. By default, `find` doesn't actually output filtered ASV sequences; instead, it outputs comprehensive information about the effect of each term and threshold set on the number of ASVs filtered and, crucially, the numbers of validated ASVs retained or rejected by each threshold set. This information can then be used to guide a `dump` run to actually output filtered ASV sequences without putative NUMTs.

Running in `find` mode requires the most input data, as metaMATE will be undertaking all of its core processes. This will generally be the first mode that is used, in order to generate a set of output statistics to analyse. 

A default `find` run carries out five main tasks:
1. Parse the input frequency filtering specification into a set of binning strategies and thresholds.
2. Assess all ASVs for potential membership of the authentic or non-authentic control groups, by
   1. Finding ASVs that match to the supplied reference set(s) or blast database(s) using BLAST
   2. Finding ASVs that fall outside acceptable length or translation parameters
3. Bin ASVs according to the specified binning strategies and generate counts of ASV reads within these bins
4. For each specified set of thresholds, assess all ASVs for retention or rejection according to their binned read frequencies
5. Output a report detailing counts of ASV rejection and retention overall and for the two control groups, over all thresholds.

This report can then be easily interrogated by the user according to project-specific requirements to balance rejection and retention. 

### `dump` introduction

The purpose of `dump` mode is to output a set of filtered ASVs without any NUMTs. It does this by enacting a single desired threshold set, either by providing the results from a `find` run and selecting the desired threshold set, or by providing an ASV set, other necessary inputs, and a single threshold specification. In this latter case, metaMATE runs a slimmed-down version of a `find` run, skipping step 2, running step 4 only once (rather than once for every combination of thresholds), and skipping step 5. This functionality is provided for enhanced versatility of the tool for differing applications, but it is recommended that for the most accuracy, `dump` is used on the analysed outputs from a `find` run.

## Installation


The best place to get metaMATE is to install from [the PyPI package](https://pypi.org/project/metaMATE/). The metaMATE source is available on [GitHub](https://github.com/tjcreedy/metamate).


metaMATE was developed and tested on Ubuntu Linux. It has not been tested anywhere else, but will probably work on most Linux systems, and likely Mac OS as well. No idea about Windows.


Note that just installing from the PyPI package is not sufficient to run metaMATE, it also requires some system dependencies and R packages.


### Dependencies


metaMATE requires python3 and the python3 libraries biopython and scipy. These should automatically be installed if using the pip installer above.


metaMATE requires the following executables to be available on the command line:
* Rscript (part of [R](https://cran.r-project.org/))
* blastn and makeblastdb (part of [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download))
* [mafft](https://mafft.cbrc.jp/alignment/software/)


The R packages getopt, ape and phangorn are also required.


### Quick install


Make sure you have all of the system dependencies: Python3 (3.6+), pip, [MAFFT](https://mafft.cbrc.jp/alignment/software/linux.html), [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [R](https://cran.r-project.org/). If you're on ubuntu linux, you can run:

```
sudo apt install python3 python3-pip mafft ncbi-blast+ r-base
```

Install metaMATE from the PyPI package:

```
python3 -m pip install metaMATE
```

or

```
sudo -H python3 -m pip install metaMATE
```

If the above fails, ensure that the python version that `python3` refers to in your current environment is >=3.6 (`python3 --version`).

Finally, ensure the necessary R libraries are installed:

```
Rscript <(echo "install.packages(c('getopt', 'ape', 'fastcluster'), repos = 'https://cloud.r-project.org')")
```

Done!

## Usage

metaMATE has summary help built in for the overall tool and the two run modes. These can be accessed by running the below commands.

```
metamate --help
metamate find --help
metamate dump --help
```

This documentation provides further explanation for input data and argument selection, and the [details](#details) section goes in depth into the way that the key parts of metaMATE work.

Note that all commandline arguments can be provided in a file, one per line, with the path to the file supplied on the commandline preceeded by `@`. For example, running `metamate find -A asvs.fasta -R references.fasta --realign` would be the same as `metamate find @args.txt` where the contents of `args.txt` is:

```
-A
asvs.fasta
-R
references.fasta
--realign
```

### Specifications

metaMATE enables considerable flexibility in the way that frequency filters can be applied in amplicon filtering. Unfortunately, this means it has a slightly complex way of specifying these filters. This is described in detail here, but for a quick start, you can simply use the [specifications.txt file available in the GitHub repository](https://github.com/tjcreedy/metamate/blob/master/specifications.txt). We recommend starting with this, looking at the results, and then modifying it as necessary for your data.

A filtering specification consists of one or more 'terms'. Each term comprises three parts, as follows:

#### 1. Categories
Read frequencies can be binned according to four categories or combinations of these categories

The four available categories are:
* "total" = the total number of reads for a haplotype across the entire dataset
* "library" = the number of reads of a haplotype per library (or sample, [see above](#input-data-required)
* "clade" = the clade assignment of a haplotype
* "taxon" = the taxon assignment of a haplotype


The binning strategy must start with either "total" or "library". Counts will be further subdivided by any further terms. For example: 
* "total|clade" will bin reads by clade over the whole dataset
* "library" will use only the per-library ASV reads
* "library|taxon+clade" will bin reads by unique clade and taxon within each library

#### 2. Metrics
Frequency can be assessed as one of two metrics, specified as follows:
* "n" = the absolute total of reads per ASV within a category or combination of categories
* "p" = the proportion of reads per ASV relative to the total number of reads of all haplotypes within a category or combination of categories

#### 3. Thresholds
Thresholds for designating NUMTs can be specified as a single value, a range of values, or mixture. 

Ranges are specified in the form "start-stop/nsteps", e.g. `1-2/5` will run with threshold values of `1, 1.25. 1.5, 1.75, 2`

Multiple values or ranges are specified in the form `a,b,c`, where `a`, `b`, and `c` can each be single values or ranges, for example `1,2,3,4-10/4` will expand to `1, 2, 3, 4, 6, 8, 10`

#### Building a term

A term comprising of one of each of the three specifications (category(/ies), metric, threshold(s)) is written within square brackets, with the three parts separated by semicolons, and any spaces are ignored. For example:

`[total; n; 2-10/5]`

This specification would run five iterations. In each iteration, ASVs would be designated as NUMTs if they have fewer than 2, 4, 6, 8 or 10 reads respectively in total across the entire dataset.

#### Multiple terms

When running metaMATE in `find` mode, multiple specification terms can be compared sequentially, or 'additively', in the same run. These are included with the `+` symbol between terms in the specifications text file, either on one line or split over multiple lines, e.g.:
`[total; p; 0.001,0.005,0.01] + [library|taxon; p; 0.4,0.6,0.8]`
can also be written as
```
[total; p; 0.001,0.005,0.01] 
+ [library|taxon; p; 0.4,0.6,0.8]`
```
These terms would run 6 total iterations: 
* the first 3 would designate ASVs as NUMTs if they had less than 0.1%, 0.5% or 1% of the total number of reads respectively
* the second 3 would designate ASVs as NUMTs if they were present as less than 40%, 60% or 80% of the total reads per taxon per
  library in all* of the taxon-library combinations in which they occur.

Multiple terms can also be compared simultaneously, or 'multiplicitively'. When running metaMATE in `find` mode, thes would be specified using the `*` betwene terms in the specifications text file, e.g.:
`[library; p; .001,0.005,0.01] * [total|clade; n; 2,5,10]`
These lines would run 9 total iterations, comprising all possible combinations of thresholds from each set. 
In the first iteration, ASVs would be designated as NUMTs if they had less than 1% of the reads in all* libraries in which they were present AND/OR if they had fewer than 2 reads within their clade across the entire dataset, and so on.

When running metaMATE in `dump` mode and not supplying a `-C/--resultcache`, simultaneous terms can be specified on the command line. Note that in this mode each term can only have one threshold. For example:
`-s '[library; p; .001]' '[total|clade; n; 5]'`
This line would run 1 iteration, designating any ASVs as NUMTs that appear as less than 0.1% of the reads in all libraries in which they occur or as fewer than 5 reads in the clade in which they occur across the entire dataset. Note the example uses `'` to avoid the shell parsing the `|` and `;` symbols.

Sequential and simultaneous terms can all be run together on the same run, for example:

```
[total; n; 2-10/5] + [library|clade; p; 0.15,0.3]
+ [library; p; 0.1] * [total|taxon; n; 3,5] * [library|clade+taxon; p; 0.03-0.07/3]
```

These lines would run 5+2+1\*2\*3 = 13 iterations:
* The first 5 iterations would designate as NUMTs any ASVs with fewer than 2, 4, 6, 8 or 10 reads across the entire dataset respectively.
* The next two iterations would designate as NUMTs any ASVs with less than 15% or 30% of the total reads for their clade within all* libraries in which they occur.
* The final 6 iterations would combine three sets. In all cases, ASVs would be designated as NUMTs if they occur in less than 10% of the total reads within all* libraries in which they occur. Depending on the iteration, ASVs would also be designated as NUMTs if they occur as fewer than 3 or 5 reads within their taxon over the entire dataset (binning by taxon is redundant here) or as less than 3%, 5% or 7% of the total reads for their taxon within their clade within all* libraries in which they occur

Remember, lines beginning with `#`, blank lines, spaces and line breaks are always ignored. So specifications can be written in any of the following ways:
```
A + B * C
```
```
A
+ B
* C
```
```
A +
# a comment
B * C
```
  
* Note that more strict filtering for library-based specifications can be applied by setting `--anyfail`, whereby ASVs will be designated as NUMTs if they fail to meet the threshold in *any* library in which they occur, as opposed to the default of *all* libraries in which they occur. This is currently a global setting, i.e. it applies to all terms involving library filtering. It could be applied on a per-term basis if there is demand.


### Core arguments

#### `-A/--asvs path`

`find`: *always required* | `dump`: *always required*

`path` should be the path to a fasta file of unique sequences with unique header names to filter. There are no header format requirements. For best operation of metaMATE, the sequences in this file should have had primers trimmed and low-quality sequences removed. It is recommended that sequences that are more than 100bp longer or shorter than the expected range of your target locus should have been removed, but other length variants remain. This file must contain some unwanted length variants ([see below](#identify-validated-non-target-ASVs) for why.).

The fasta file may be aligned or unaligned. If metaMATE requires an alignment (to build a UPGMA tree and delimit clades), an unaligned input will be aligned using [MAFFT FFT-NS-1](https://mafft.cbrc.jp/alignment/software/manual/manual.html). Supply an aligned set of ASVs if alignment is required but FFT-NS-1 is not expected to perform well with your data. 

#### `-L/--libraries path [path]`

`find`: *always required* | `dump`: *required* if not supplying `-C/--resultcache`

The read library or libraries supplied are used by metaMATE to assess the incidence of each ASV sequence per library. Each `path` should be the path to a fasta or fastq file containing ASV reads. Files may contain reads that are not in the `-A/--asvs` file, these will be ignored. If multiple paths are supplied, metaMATE assumes that each file is a separate library and uses the file names as library names. In this case, the headers in each file are ignored, only the sequences are relevant. If a single path is supplied, then metaMATE assumes that the library names are specified in the read headers, specified in the format `;barcodelabel=NAME;` or `;sample=NAME;` where `NAME` is the unique name of each library.

Note that the composition of libraries, i.e. the expected true richness and abundance of individuals within the sample, and whether they are complete samples, subsamples, or replicates, should be carefully considered when designing the specifications.

#### `-o/--outputdirectory path`

`find`: *always required* | `dump` *sometimes required*

If using a mode that outputs multiple files, `path` should be the path to a directory in which to place these files. If the directory already exists, metaMATE will exit with an error unless `--overwrite` is set. 

#### `-y/--anyfail` *flag*

If supplied, when comparing the per-bin frequencies of an ASV against a threshold, metaMATE will designate the ASV as a NUMT for that threshold set if it fails to meet the threshold in *any* of the bins in which it occurs (instead of the default, *all*). This is substantially more stringent and is not generally recommended, but may be suitable for small datasets.

#### `--realign` *flag*

If supplied, metaMATE will always run alignment on the supplied ASVs, whether already aligned or not. However this flag will be ignored if a tree is supplied, or if a `-C/--resultcache` is supplied in `dump` mode.

#### `--overwrite` *flag*

If supplied, metaMATE will overwrite any files in the destination output directory (including the current directory, if this is the specified `-o/-outputdirectory`) that have the same name as new outputs. This is provided as insurance against accidentally re-starting a run from scratch.

#### `-t/--threads n`

`n` should be a positive integer specifying the maximum number of parallel threads to use, where relevant. The main frequency threshold comparison section of metaMATE is run on multiple threads to improve speed. This argument is also passed to `MAFFT` and `blastn` if these are run. 

### Clade delimitation and binning arguments

`find`: *optional* | `dump`: *optional* 

These arguments apply if ASV binning is to be performed on a per-clade basis in mode `find` or in mode `dump` without a `-C/--resultcache`.

#### `-T/--tree path`

`path` should be the path to a newick-format ultrametric tree. Supplying a tree will skip alignment of unaligned ASV sequences supplied to `-A/--asvs`, and subsequent distance matrix computation and UPGMA tree building. This can speed up the run time substantially for datasets with large numbers of ASVs.

#### `-d/--divergence [0-1]`

A value between 0 and 1 specifiying the maximum percent divergence (1 = 100%) by which to delimit clades from a UPGMA tree. This is implemented by passing the supplied value to the `h` argument of the R function [`cutree`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cutree.html).

#### `--distancemodel MOD`

`MOD` should be the short name specifying the evolutionary model to be used for calculating pairwise distances between aligned ASV sequences for the purpose of building a UPGMA tree. The value of this argument is passed directly to the `model` argument of the R  function [`dist.dna`](https://www.rdocumentation.org/packages/ape/versions/5.4/topics/dist.dna) from the `ape` package, and should match exactly to one of the allowed functions for that model, namely "raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel" or "indelblock". By default, metaMATE uses "F84". Please read the [documentation](https://www.rdocumentation.org/packages/ape/versions/5.4/topics/dist.dna) for this function to select a model.

### Taxon binning arguments

#### `-G/--taxgroups path`

`find`: *optional* | `dump`: *optional* 

In mode `find` or mode `dump` without a `-C/--resultcache`, if ASV binning is to be performed on a per-taxon basis ([see specifications](#specifications)), `path` should be the path to a comma-delimited text file with as many rows as ASVs supplied to `-A/--asvs` and two columns. The first column should give the ASV read names, the second should give a grouping value for each ASV. This group is intended to allow for taxonomic bias in NUMT frequency, so is expected to be a taxonomic group, but metaMATE simply bins ASVs by unique group names so the actual content could be any string. The number of groups should be relatively low, less than 1% of the number of ASVs, and the ratio of ASVs to group number should be carefully considered when setting thresholds, remembering that when binning by library and taxonomic group, each library / taxon combination has only a small subset of the total ASV number.

### `find`-specific arguments

#### `-S/--specification [path]` *required*

`path` should be the path to a text file containing a specification for the terms apply when applying frequency filtering. These are explained in detail above, and an example specifications.txt file is [supplied in the GitHub repository](https://github.com/tjcreedy/metamate/blob/master/specifications.txt).

#### `-g/--generateASVresults [n]`

By default, `find` mode does not output filtered ASVs, instead providing detailed statistics for assessment to determine the optimal threshold set for the dataset in question. This option provides for the output of a fasta file for each threshold set, containing the sequences that were not rejected by that threshold set and/or were determined to be verified-authentic ASVs. If supplied without a value, a fasta file will be output for each threshold set that recieved the best score (see below). If a value is supplied, metaMATE will output a fasta file for each threshold set that scored within the top `n` proportion of threshold sets. If `n` is 1, a fasta will be output for every threshold set. This can generate a large number of fasta files! For more information about scoring, [see below](#scoring-and-estimation).

#### `-q/--scoremetric [metric]`

When using `-g/--generateASVresults`, you can select which scoring metric (accuracy, precision or recall) is used to select the best threshold set(s) by supplying the metric name to `-q/--scoremetric`. [See below](#scoring-and-estimation) for more details on these metrics.

### Reference-matching arguments

`find`: *required* | `dump` *not used*

These arguments control matching against a reference fasta or blast database for the purposes of determining verified-authentic ASVs. Put simply, an ASV is designated a verified-target-ASV if it matches against a reference sequence with a match length greater than specified and a percent identity greater than specified. Note that if multiple ASVs match to the same reference when percent identity is set to 100, only the ASV with the longest match length will be designated as a verified-target-ASV. 

#### `-R/--references path` and/or `-D/--blastdb path`

`path` should be the path to a fasta file (`-R/--references`) or blast database (`-D/--blastdb`) of sequences that represent known species that are likely to occur in the dataset. Both arguments are available for cases where multiple reference sources are desired, with different parameterisation of hits, but only one is required. For example, you may have a set of sanger-sequenced barcodes from your project's morphospecies, and which have been carefully curated to ensure accuracy and no NUMTs, against which you want to allow matches of 99% to allow for some minor sanger sequencing error and true haplotypes. You may also want to find 100% matches to a local copy of GenBank nt. You would supply the former to `-R` and the latter to `-D`. A BLAST search against these sequences will be used to designate a set of ASVs as verified-authentic. 

#### `--refmatchlength n` and/or `--dbmatchlength n`

`n` should be a positive integer specifying the minimum alignment length to consider a BLAST match when comparing ASVs against sequences in the file supplied to `-R/--references` or the database supplied to `-D/--blastdb`. The default value is calculated as 80% of the length below which ASVs will be designated as verified-non-authentic. 

#### `--refmatchpercent [0-100]` and/or `--dbmatchpercent [0-100]`

The supplied value is the minimum percent identity against to consider a BLAST match when comparing ASVs against sequences in the file supplied to `-R/--references` or the database supplied to `-D/--blastdb`. The default value is 99.9 for `--refmatchpercent` and 100 for `--dbmatchpercent`.

#### `--keeptemporaryfiles`

If this argument is supplied, the temporary blast database and/or blast result xml files generated during reference and/or blast database matching will not be deleted and can be found in a directory named 'blastdb' inside the output directory.

### Length-based arguments

`find`: *required* | `dump` *not used*

Length based parameters are used to calculate a set of lengths; an ASV with length outside this set will be designated a verified-non-authentic-ASV. A number of different arguments are available for flexible delimitation of this range. The user must specify sufficient information to compute the range, i.e.
1. the minimum and maximum length of the range, OR
2. the expected length of the centre of the range, AND one of
   1. the number of bases around this expectation,
   2. the percentage variation around this expectation,
   3. the number of codons around this expectation

Additionally, `--onlyvarybycodon` may be specified, which further restricts the range to only values differing from the expected length by a multiple of three. If used in conjunction with a minimum and maximum length specification, `--onlyvarybycodon` requires that `-l/--expectedlength` is also specified.

#### `-n/--minimumlength n` and `-x/--maximumlength x`

The values of `n` and `x` should be positive integers. These values would generally describe the distribution of lengths expected for real amplicons derived from the target locus. Any ASVs outside this range of lengths will be designated a verified-non-target ASV.

#### `-l/--expectedlength n`

`n` should be a positive integer that specifies the exact or centroid length expected for real amplicons derived from the target locus. This argument is required if specifying `--basesvariation`, `--percentvariation` or `--codonsvariation`, or if specifying `--onlyvarybycodon`

#### `--basesvariation b`, `--percentvariation [0-100]`, `--codonsvariation c`

One, and only one, of these arguments is required if specifying `-l/--expectedlength` but not specifying `-n/--minimumlength` and `-x/--maximumlength`. The supplied value is used to compute a range of values describing the distribution of lengths expected for real amplicons derived from the target locus. Any ASVs outside this range of lengths will be designated a verified-non-target ASV. `--percentvariation` refers to a percentage of the value supplied to `-l/--expectedlength`. `--codonsvariation` is simply a multiple of 3 bases, i.e. `--codonsvariation 2` == `--basesvariation 6`

#### `--onlyvarybycodon` *flag*

If specified, the computed range of values will be reduced to only include those values that differ from the value supplied to `-l/--expectedlength` by 3. The purpose of this is to allow realistic variation around the expected length, i.e. variation that does not involve a frameshift in the translation of the sequence. 

### Translation-based arguments

In most cases, these arguments can be left alone, aside from `-s/--table` which is alway required. By default, the reading frame is automatically detected by translating ASVs one-by-one in all frames and recording stop counts until a) a minimum sample size of total number of stops has been achieved and b) one frame achieves a significantly fewer stops than the other two frames. ASVs are translated in the order in which they are encountered in the ASVs file, under the assumption that ASVs are usually written in descending order of frequency after dereplication.

#### `-s/--table n`

`find`: *required* | `dump` *not used*

`n` should be a value corresponding to one of [the standard NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for metaMATE to use when translating ASVs. This is always required. 

#### `-r/--readingframe [1,2,3]`

If supplied, this value will be used as the reading frame for amino acid translation for all ASVs, where reading frame 1 translates from the 1st base of the ASV sequence, and so on. If not supplied, this will be automatically detected.

#### `--detectionconfidence [0-1]`

If `-r/--readingframe` is not specified, metaMATE will automatically detect the reading frame. Higher values of this argument increase certainty in this detection. The default is 0.95, i.e. 95% confidence.

#### `--detectionminstops n`

If `-r/--readingframe` is not specified, metaMATE will automatically detect the reading frame. Higher values of this argument increase the minimum number of stops that must be encountered in all reading frames before a reading frame can be selected. In cases of very small datasets, this value may need to be reduced from the default of 100.

### `dump`-specific arguments

#### `-f/--outfasta path`

`path` should be the path to a file to write retained ASVs to.  If the file already exists, metaMATE will exit with an error unless `--overwrite` is set. `-f/--outfasta path` is only required if the other arguments specify a single output file; otherwise `-o/--outputdirectory` would be required.

#### `-C/--resultcache path`

`path` to a '\_resultcache' file output from a previous metaMATE run in mode `find`. If specified, `-i/--resultindex` is required.

#### `-i/--resultindex n [n]`

Each value of `n` should be 0 or a positive integer referring to an result set index for which to extract retained filtered ASVs from a specified `-C/--resultcache`. Result set indices can be found in the statistics output from a `find` run.

#### `-S/--specification '[c; m; t]'`

`'[c; m; t]'` should be one or more text strings specifying terms to apply when frequency filtering. These are explained in detail below, but in brief, each term is in the format `[category(/ies); metric; threshold(s)]` where category is one or a combination of "total", "library", "clade" or "taxon", metric is one of "n" or "p", and threshold is a single value. In `dump` mode, terms are always considered simultaneously and only one threshold is permitted per term. For example, `[library|clade, p, 0.05]` would designate as NUMTs any ASVs that occurred as fewer than 5% of the reads within members of its clade within all library in which it occured.

## Examples

The following commands detail example runs of metaMATE using metabarcoding data provided in the [GitHub tests directory](https://github.com/tjcreedy/metamate/tree/master/tests/data). The README in that directory details how this data was generated. Where a specifications file is used, this refers to the [default specifications available on the GitHub](https://github.com/tjcreedy/metamate/blob/master/specifications.txt)

Note that the target amplicon size of this dataset is 418bp, so length specifications will be based on this; your data will vary! 

### `find` examples

This is probably the simplest `find` run that cound be run on this data:

```
metamate find -A 6_coleoptera.fasta -L 0_merge/*.fastq -S specifications.txt -G 6_coleoptera_taxon.csv -R dummy_references.fasta --expectedlength 418 --percentvar 0 --table 5 -o outputdir
```
The path to a fasta containing input ASVs to be filtered have been supplied to `-A`, and the paths to each of set of libraries to `-L` (note the bash parameter expansion: `0_merge/*.fastq` will expand to all files ending with `.fastq` in `0_merge/`). Each of the library files will be searched for each of the ASVs to find the count of ASVs per library. Note that these sequences have not been quality filtered - metaMATE is likely to be more accurate with quality filtered reads.
The specifications file path has been supplied to `-S`, this will be parsed to find all specification terms, thresholds and combinations and generate all iterations to run.
To find verified-authentic-ASVs, metaMATE will search the ASVs against the sequences in the `dummy_references.fasta` file path provided to `-R`. To find verified-non-authentic-ASVs, metaMATE will check the length of each ASV and translate it using NCBI translation table 5 (`--table`). Any ASV falling outside 418 +- 0% in length and/or containing any stops in translation will be designated as verified-non-authentic.
metaMATE will align and build a UPGMA tree from the input ASVs to group ASVs into clades. ASVs will also be grouped into taxa according to the tablular file path supplied to `-G`.
Once filtering has been performed for all specified threshold sets, the results will be written to `outputdir` (which will be created if not present). As the input ASVs are unaligned and no tree was provided, these results will contain a fasta with aligned ASVs and a newick format text file with a UPGMA tree of the ASVs.

If your reads are in a single file, with library information in the headers, rather than multiple files for each different library, you would supply the path to that file to `-L` instead:

```
metamate find -A 6_coleoptera.fasta -L 6_concat.fasta -S specifications.txt -G 6_coleoptera_taxon.csv -R dummy_references.fasta --expectedlength 418 --percentvar 0 --table 5 -o outputdir
```
Note that either `6_concat.fastq` or `6_concat.fasta` can be used here. The former has not been quality filtered, so read counts will be higher, which will affect the output of metaMATE. It is suggested that quality filtered reads be used as the input.


If you already have aligned ASVs and a tree file, but don't have any references, you could instead run:
```
metamate find -A 6_coleoptera_fftnsi.fasta  -T 6_coleoptera_UPGMA.nwk -L 0_merge/*.fastq -S specifications.txt -G 6_coleoptera_taxon.csv -D /blastdb/nt -expectedlength 418 --percentvar 0 --table 5 -o outputdir 
```
The path to aligned ASVs is supplied to the same argument as unaligned ASVs, metaMATE will detect whether the file is aligned or not. The path to a tree file has been supplied to `-T`. The `-L`ibraries and `-S`pecifications arguments are the same as in the first example. 
The `-D` argument is used to supply a path to a [blast-formatted database generated by `makeblastdb`](https://www.ncbi.nlm.nih.gov/books/NBK279688/). In this case, it looks like this is a local copy of GenBank nt (note this is not supplied in the test data). The remaining arguments are the same as the previous command.

This command will be faster than the previous command, as metaMATE does not need to generate the alignment and tree. It may also be more accurate, because you can build the optimal alignment for your specific data, rather than using the fast alignment algorithm metaMATE uses. You can also review the tree and decide on the most appropriate clade delimitation threshold for your data - by default this is 20%.

This final example overwrites many of these sorts of defaults:

```
metamate find -A 6_coleoptera_fftnsi.fasta -T 6_coleoptera_UPGMA.nwk -d 0.15 -L 0_merge/*.fastq -S specifications.txt -G 6_coleoptera_taxon.csv -R dummy_references.fasta --refmatchpercent 98 --refmatchlength 400 -D /blastdb/nt -expectedlength 418 --basesvariation 6 --onlyvarybycodon -s 5 -o outputdir 
```
The input ASV alignment and tree are the same, but clades will be delimited at 15% divergence rather than the default 20%. Verified-authentic ASVs will be determined by matches against two references, the references fasta file (against which passing hits must be 98% similar over at least 400 bp), and the nt dataset (which uses the default similarity and length settings). The expected length is unchanged, but we now specify a more complicated way of determining verified-non-authentic ASVs: any ASVs that are less than 6 bases of variation around 418 **and** do not vary by exactly 3 bases from 418. Thus any reads that are not 412, 415, 418, 421 or 424 bp will be designated as verified-non-authentic. Finally, the translation table has been supplied using the short option `-s` rather than the long option `--table`.

### `dump` examples

The most straightforward, and recommended, way to run `dump` is to do so simply to extract the results of a specified result set:
```
metamate dump -A 6_coleoptera.fasta -C outputdir/6_coleoptera_resultcache -i 35 -o 7_coleoptera_numtfiltered.fasta
```
Alongside the same ASV file as used for a `find` run (`-A`), this command specifies the path to a resultcache file output from that `find` run (`-C`). After analysing the results, the user has determined that the optimal threshold set (`-i`) is set 35 (an index given in the `find` output data). metaMATE will parse these inputs and write all ASVs that passed the frequency thresholds for this set *and/or* were determined to be verified-authentic-ASVs to the output file (`-o`). All of the specifications, threshold sets and other input files are not needed as all of this information is contained in the resultcache file.

Alternatively, `dump` mode can work very similarly to `find` mode. Rather than taking a resultcache and result index, you can supply any arguments needed for binning and threshold specification. For example:
```
metamate dump -A 6_coleoptera.fasta -L 2_concat.fasta -S '[library|clade; p; 0.05]' -o outputdir
```
The main differences between this and a `find` run are:
* filtering specifications are supplied directly to the commandline
* no arguments specifying how to delimit control groups are used
`dump` mode resolves any binning strategies, building a tree if needed for clade-based strategies, and then applies the specified thresholds, writing all ASVs that pass the thresholds to a file in the output directory. This method is provided for rapid filtering on a known dataset.

*Important note:* using `dump` by specifying thresholds is unlikely to return the same set ASVs for the same input data and threshold specification than using `find`, outputting a resultcache and running this in `dump`. [See below](#numtdumped-asv-fasta-format-file) for more details.

## Outputs

metaMATE returns different outputs depending on mode and arguments. Unless specified by using `-f/--outfasta`, all output files will use the same name as the input ASVs with a suffix appended.

### NUMTdumped ASV fasta-format file

As standard, for a given threshold set, metaMATE outputs those ASVs that fulfil the following conditions, in order of priority:
1. ASV is a verified-authentic ASV
2. ASV is not a verified-non-authentic ASV
3. ASV is present in a frequency equal to or exceeding any thresholds in any bins in which it occurs

The standard way to generate these files is by first running metaMATE in mode `find`, then selecting one or more result set indices from the output statistics file, then running metaMATE in mode `dump`. If one result set is selected, the resulting filtered ASV file will be output with the name specified by `-f/--outfasta`. Otherwise, the resulting filtered ASV files will be output to the directory specified by `-o/--outputdirectory`, with the suffix '_numtdumpresultsetN', where N is the index. This is also the effect of running `find` with the `--generateASVresults` argument.

This pipeline is designed to take account of the fact that threshold specifications and decisions on the optimal output should be dataset-dependent, and filtering with a single threshold set without validation is arbitrary.

There are two non-standard methods that generate outputs with different conditions. Firstly, using `--anyfail` in either `find` or `dump` mode will change condition 3 to "ASV is present in a frequency equal to or exceeding all thresholds in any bins in which it occurs", [see above](#-y--anyfail-flag) for details.

In the second case, users may optionally skip all validation and simply use mode `dump` to output a file of filtered ASVs based on a single set of term specifications and thresholds. In this case, the output ASVs are only those where the ASV is present in a frequency equal to or exceeding all thresholds in any bins in which it occurs. Thus the set of ASVs output from `find` and `dump` for the same set of term specifications and thresholds **may not be identical**, because `dump` does not perform any validation. 

### Results (`find` only)

The main output from metaMATE `find` mode is the *_results.csv file. This is a comma-delimited table that synthesises all of the results from applying all combinations of the specified terms and thresholds to the input ASVs, given the control groups of authentic- and non-authentic- ASVs. It is designed to be easily parseable and reformatable by downstream processes, in particular for analysis with R (a template for analysis is in development). Its columns comprise:

* *resultindex*: a unique identifier for each term specification and threshold combination, for ease of ASV retrieval using `dump` mode.
* *term*: each additive term, repeated for a number of rows equal to the number of thresholds supplied to this term. For example, if the first terms in the specifications are `[library; n; 10-100/20] + [library; p; 0.001-0.1/10]`, the first 20 rows will be for term `library_n` and the next 10 rows will be term `library_p`.
* *_threshold*: one or more columns giving threshold values for all bins/metric combinations in the specification. If a term does not include a bin/metric combination, this value will be 0. For the example above, the third column would be `library_n_threshold` and contain 20 values from 10 - 100, followed by 10 0s. The fourth column would be `library_p_threshold` and contain 20 0s followed by 10 values from 0.001 - 0.1.
* *accuracy_score, precision_score, recall_score*: the accuracy, [precision and recall](https://en.wikipedia.org/wiki/Precision_and_recall) scores for this threshold set based on the rejection/retention of va- and vna-ASVs, [see below](#scoring-and-estimation) for details
* *asvs_total*: the total number of ASVs in the dataset.
* *verifiedauthentic_total_observed*: the total number of verified authentic ASVs 
* *verifiednonauthentic_total_observed*: the total number of verified non-authentic ASVs
* *asvs_prelim_retained_n/p*: the number of ASVs and proportion of total ASVs that passed the given threshold(s), before rejecting any surviving verified-non-authentic ASVs and retaining any lost verified-authentics ASVs
* *asvs_prelim_rejected_n/p*: the number of ASVs and proportion of total ASVs that failed the given threshold(s), before rejecting any surviving verified-non-authentic ASVs and retaining any lost verified-authentics ASVs
* *verifiedauthentic_retained/rejected_n/p*: the number of and proportion of total verified authentic ASVs that passed/failed the given threshold(s)
* *verifiednonauthentic_retained/rejected_n/p*: the number of and proportion of total verified non-authentic ASVs that passed/failed the given threshold(s)
* *asvsactual_retained_n/p*: the total number of and proportion of total ASVs that were verified authentic, not verified non-authentic, and passed the given threshold(s)
* *asvsactual_rejected_n/p*: the total number of and proportion of total ASVs that were verified non-authentic and/or failed the given threshold(s)
* *verifiedauthentic/nonauthentic_total_estimate*: the estimated number of verified authentic and verified non-authentic ASVs in the initial dataset, prior to filtering, based on the calculations described below.
* *verifiedauthentic/nonauthentic_retained_estimate*: the estimated number of verified authentic and verified non-authentic ASVs within those preliminarily retained, based on the calculations described below.
* *rejects_hash*: the hashed alphabetically sorted list of actual rejected ASVs for this given set. Identical values on different rows denote those rows rejected an identical set of ASVs.

### ASV counts

The *_ASVcounts.csv file is a comma-separated table recording the number of reads of each input ASV found in each library.

### Clade groupings

The *_clades.csv file is a two-column comma-separated table recording the clade grouping for each input ASV, generated by the script `get_clades.R` supplied as part of metaMATE.

### Control list (`find` only)

The *_control.txt file is a two-column tab-separated table recording all ASVs determined to be validated-authentic or validated-non-authentic. The first column lists the reason for each determination: "lengthfail" means the ASV was too short, too long or otherwise did not fall into the range of acceptable lengths; "stopfail" means the ASV had stop codons in its amino acid translation; "refpass" means the ASV had a passing BLAST match to one of the supplied references, and did not fail either of the non-authentic tests.

### Result cache (`find` only)

The *_resultcache file is a compressed text file containing information on the ASVs rejected or retained for each of the supplied specification terms and threshold sets of a `find` run. It can be used in conjunction with the ASVs from the run and one or more result indices to extract the ASVs in a fasta file format using mode `dump`.

### Aligned/Unaligned fasta

If the input ASVs were aligned, metaMATE unaligns (degaps) them and outputs the sequences to the *_unaligned.fa file, or vice versa if the input ASVs were unaligned.

### UPGMA tree

Unless a tree is supplied, metaMATE uses the aligned ASVs to build a UPGMA tree using the script `make_tree.R` supplied as part of metaMATE. This is a newick-format tree that can be opened by any newick parser or tree viewing software.

## Details

### Validation

As part of a metaMATE run, the user must parameterise the determination of two control groups of ASVs. Control ASVs are those that can be determined _a priori_ to be either validly authentic or non-authentic: authentic ASVs are determined by a match against a reference set of sequences, non-authentic ASVs are determined by falling outside acceptable length and translation properties. The user must provide a reference set and specify acceptable sequence lengths, and may fine-tune reference matching and translation assessment. 

### Generating clades

Clades are generated as following:
1. If the input ASVs are not aligned, they are aligned using the mafft FFT-NSI algorithm. This algorithm is used because ASV sets may be very large and a more accurate algorithm could be very time consuming to run. It is strongly suggested that if this algorithm is not likely to be accurate for your ASVs, that you pre-align your ASVs using a suitable method
2. A distance matrix is computed for all sequences, using Felsenstein's distance implemented in Phylip (1984), unless an alternative metric is provided to `--distancemodel`.
3. A UPGMA tree is constructed based on the distance matrix
4. The UPGMA tree is cut into subtrees at the divergence level supplied to `-d/--divergence`.
5. The members of each subtree are considered to belong to the same clade

### Scoring and estimation

metaMATE scores the successful rejection of verified non-authentic ASVs and retention of verified authentic ASVs using three metrics: accuracy, precision and recall. These are standard metrics employed in cases where a model is tested on known data, for example when training machine learning models. Accuracy is simply the proportion of verified ASVs that were dealt with accurately, with a value of 1 denoting complete accuracy - all vna-ASVs were rejected by a given set of thresholds while all va-ASVs were retained. Precision and recall [are more specific metrics](https://en.wikipedia.org/wiki/Precision_and_recall).  Note that we do not suggest this as the only way to select the best threshold set - multiple threshold sets may recieve the same score, and it's important to carefully select the levels of va-ASV retention and vna-ASV rejection that you are happy with rather than using a reductive score to decide for you. However these scores may be useful in helping you decide.


## Development
