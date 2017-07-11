# seqpy
Python package that provides library of tools and API for dealing with next generation sequencing data.

## Python/R tools under /bin

#### Aggregating result files

**aggregate_fastqc.py** Combine output files from running FASTQC.

**aggregate_flagstat.py** Combine output files from running samtools flagstat.

**aggregate_fpkm.py** Combine output files from running cufflinks.

**aggregate_fusioncatcher_results.py** Combine output files from running fusioncatcher.

**aggregate_gatk_coverage.py** Combine output files from running GATK coverage.

**aggregate_htseqcount.py** Combine output files from running HTSeq-count.

**aggregate_mpileup.py** Combine output files from running samtools mpileup.

**aggregate_star_QC.py** Combine final output log files from running STAR.

**aggregate_trinucleotide_counts.py**

### VCF and mutation calling tools

**filter_control_vcf.py**

**filter_mutect.py** Filter output VCF file form MuTect2 (in paired or single mode) using various filtering criterion.

**filter_pindel.py** Filter output VCF file form Pindel (in paired or single mode) using various filtering criterion.

**maf_to_matrix.py**

**modify_QSS_field.py**

**vcf_to_table.py**

### Misc

**baseq_stats.py**

**bed_to_kmer_frequency.py**

**biomart_convert_gene_identifiers.R**

**count_reads_for_junction.py**

**count_sequenced_contexts.py**

**count_sequenced_kmer.py**

**deseq2.R** Run DESeq2 on the user provided counts file and samples to be compared, along with other parameters.

**deseq2rnk.py** Convert the output from DESeq2 to an .rnk file using the specified ranking metric.

**expression2gct.py** Convert a CSV file containing expression estimates to a GCT file.

**filter_bad_cigar.py** Fix the bad CIGAR formatting that sometimes happens with RNA-seq.

**get_gsea_output.py** Get the table output from GSEA and save to excel file

**mouseSymbol2human.py** Convert mouse gene symbols to human gene symbols in a csv file containing rows as genes and samples as columns.