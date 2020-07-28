# This directory contains dummy data for testing NUMTdumper
# 
#
# METABARCODING DATA
#
# The commands given here are a guideline for how data could be prepared for NUMTdumper,
# no metabarcoding pipeline is universal and each dataset should be processed in whatever
# way is most suitable for the amplicon type, project and research questions
#
# 0_merge/ - 12 metabarcoding sample read files, with primers trimmed and pairs merged
#
# 1_concat.fastq - the concatenation of the contents of 0_merge, the result of running:

samples=$(for f in 0_merge/*.fastq; do s1=${f##*/}; echo ${s1%.*}; done | sort | uniq)
for s in $samples; do sed -e "s/\(^@D00.*\) .*$/\1;sample=$s;/" 0_merge/$s.fastq ; done > 1_concat.fastq

# Note that for each file (T1.fastq, T2.fastq, etc), the sed command adds the sample name (T1, T2 etc) to the header of each read in that file before concatenation
#
# 2_concat.fasta - the result of quality filtering 1_concat.fastq, by running:

vsearch --fastx_filter 1_concat.fastq --fastq_maxee 1 --fastaout 2_concat.fasta

#
# 3_derep.fasta - the unique sequences from 2_concat.fasta, running:

vsearch --derep_fulllength 2_concat.fasta --output 3_derep.fasta --sizeout --relabel uniq

#
# 4_denoise.fasta - the denoised unique sequences, run with:

vsearch --cluster_unoise 3_derep.fasta --minsize 4 --unoise_alpha 2 --centroids 4_denoise.fasta

#
# 5_coarselengthfilter.fasta - denoised unique sequences, with outlying short and long reads removed

vsearch --fastx_filter 4_denoise.fasta --fastq_minlen 416 --fastq_maxlen 421 -fastaout 5_coarselengthfilter.fasta

#
# 6_coleoptera.fasta - denoised unique sequences that are putatively Coleoptera. These are the ASVs
#                      to be input into NUMTdumper
#
# OTHER DATA
#
# 6_coleoptera_taxon.csv - a two-column table specifying a (dummy, random) taxon for each ASV
#
#
# dummy_reference.fasta - a fasta of reference sequences
#

