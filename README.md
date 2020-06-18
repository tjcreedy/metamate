## NUMTdumper: let's dump those NUMTs!

### Overview

NUMTdumper analyses a set of amplicons derived through metabarcoding of a mitochondrial coding locus to determine putative NUMT and other erroneous sequences based on relative read abundance thresholds within libraries, phylogenetic clades and/or taxonomic groupings. 

This documentation is a work in progress! Check back very soon.

The paper for NUMTdumper is in review:  Andújar, C., Creedy, T.J., Arribas, P., López, H., Salces-Castellano, A., Pérez-Delgado, A., Vogler, A.P. & B.C. Emerson (in review). NUMTdumper: a self-validating method for generating reliable haplotype data from mtDNA metabarcoding. Methods in Ecology and Evolution. 

If you use NUMTdumper in your work, please cite this paper.

The development of this tool was supported by the iBioGen project, funded by the H2020 European Research Council, Grant/Award Number: 810729.

## Table of contents
* [How NUMTdumper works](#how-numtdumper-works)
* [Installation](#installation)
* [Usage](#usage)
  + [Input data](#input-data)
  + [`find`](#find)
  + [Outputs](#outputs)
  + [`dump`](#dump)
* [Details](#details)
  + [Identifying non-targets](#identifying-validated-non-target-ASVs)
  + [Identifying targets](#identifying-validated-target-ASVs)
  + [Generating clades](#generating-clades)
  + [ASV assessment](#asv-assessment)
  + [Scoring and estimation](#scoring-and-estimation)
* [Development](#development)

## How NUMTdumper works

NUMTdumper takes the guesswork and faff out of applying frequency-based filtering thresholds in NGS amplicon pipelines by utilising a modular threshold specification approach combined with detailed outputs to efficiently provide detailed insight into the effects of different filtering and binning strategies. This approach is explicitly developed towards removing putative NUMT sequences from a population of ASVs, but the methodology is broad enough that it will work simultaneously on any other types of low-frequency erroneous sequences, such as sequencing errors. A principal benefit of NUMTdumper over closely related tools is its validation approach, whereby input ASVs are assessed for membership of control groups for known authenticity. The effect of frequency thresholds and binning strategies on retention or rejection of the members of these is used to assess the optimal methodology for NUMT and error removal.

The downside of this approach is that NUMTdumper is more data- and parameter- hungry than other similar tools. Rather than being content with being fed a set of ASVs, NUMTdumper requires inputs that a) enable the binning of ASV reads to provide pools upon which frequency thresholds can be applied and b) parameterise the determination of the two control groups. 

NUMTdumper operates using two submodules, `find` and `dump`, with the former being the main workhorse of the tool. A default `find` run carries out five main tasks:
1. Parse the input frequency filtering specification into a set of binning strategies and thresholds.
2. Assess all ASVs for potential membership of the authentic or non-authentic control groups, by
a) Finding ASVs that match to the supplied reference set(s) using BLAST
b) Finding ASVs that fall outside acceptable length or translation parameters
3. Bin ASVs according to the specified binning strategies and generate counts of ASV reads within these bins
4. For each specified set of thresholds, assess all ASVs for retention or rejection according to their binned read frequencies
5. Output a report detailing counts of ASV rejection and retention overall and for the two control groups, over all thresholds.
This report can then be easily interrogated by the user according to project-specific requirements to balance rejection and retention. 
A `dump` run is used to enact a single desired threshold set, either by providing the results from a `find` run and selecting the desired threshold set output, or by providing an ASV set and other inputs, and a single threshold specification.

## Installation

Currently NUMTdumper is supplied as a set of python and R scripts. You can install by downloading and unpacking this repository. You should ensure that you have all of the dependencies listed below installed and accessible on the PATH. 

Then you simply need to run `python3 ./numtdumper/numtdumper.py`. 

We plan to build an automatic installer and complete package soon.

### Dependencies

NUMTdumper requires python 3 and the python 3 library biopython

NUMTdumper also requires R and the R packages getopt, ape and phangorn

The following tools must be installed and accessible on the path:
* Rscript (R)
* blastn
* mafft

## Usage

### Input data

### `find`

### Outputs

### `dump`

## Details

As part of a NUMTdumper run, the user must parameterise the determination of two control groups of ASVs. Control ASVs are those that can be determined _a priori_ to be either validly authentic or non-authentic: authentic ASVs are determined by a match against a reference set of sequences, non-authentic ASVs are determined by falling outside acceptable length and translation properties. The user must provide a reference set and specify acceptable sequence lengths, and may fine-tune reference matching and translation assessment. 


### Identifying validated non-target ASVs


### Identifying validated target ASVs

### Generating clades

### ASV assessment

### Scoring and estimation

## Development
