# webster_thesis
This project contains 3 scripts that analyze genetic variation and markers that occur among genes. These scripts were originally developed for transcriptome parsing and filtering at a specific threshold within a subgroup of genes, such as a gene family or genes grouped by Gene Ontology, however the scripts can be easily altered for application to larger data sets. Arguments that pertain to input/output files and important parameters were included in each script.

Link to written thesis:

#chapter_1_preliminary.py:

This script was used as a preliminary analysis for the first chapter of my thesis, it confirmed that mutations occurred among the Phytochrome gene family in the Eriophorum vaginatum transcriptome at a given threshold.

#chapter_1_filter_SNPs_from_VCF.py:

This script was developed to parse genes that meet specific parameters for mutation quality from a transcriptome using a VCF and FASTA file. The script reads in the sequences from the FASTA file and mutation quality and frequency from the VCF file, if the mutations fall within the given parameters, the sequence data is transformed with the mutation and written to an output directory. This output directory contains only the gene names and sequences that qualify with the threshold parameters.

Script Arguments:
-f FASTA file input
-v VCF file input
-ad Allele depth (for mutation quality)
-sd Sample depth (for mutation frequency)
-o Output directory name

#chapter_2_isolate_SSRs_inside_coding.py:

This script was developed to be used with the program SciRoKo 3.4, this program output file has information on microsatellite markers that occur in a given FASTA file. The script reads in the FASTA file and translates these genes on each reading frame, then storing where the longest coding region occurs. The positions of the microsatellites are compared with the coding region positions, and an output file is created with only the microsatellites that occur inside the coding regions. The script can also be easily altered the look at microsatellites outside the coding regions.

Script Arguments:
-f FASTA file input name
-s SciRoKo output file
-l SSR (microsatellite) minimum repeat (ex. dinucleotide = 2)
-out_file Output file
-out_fasta FASTA file output

