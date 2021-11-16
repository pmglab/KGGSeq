# The pipeline used to perform the equal weight version of RUNNER (RUNNER1) test in three real datasets of complex diseases

java -jar ./kggseq.jar
--runner1-gene-coding #the main parameter to perform the equal weight version of RUNNER (RUNNER1) test
--runner-coding-gene-cov mu_mis,mu_lof,oe_mis,oe_lof,ExonGC #Perform RUNNER1 by cosidering these covariables
--out /path/to/out_file(prefix)
--vcf-file /path/to/vcf_file
--ped-file /path/to/ped_file
--excel #Assign the output format as Excel
--nt 6 Set the number of parallel tasks for an analysis job
--buildver hg19 #Set the version of reference human genome build
--hwe-all 0.001 #QC paremeter: Exclude variants in all subjects with the Hardy-Weinberg test p value <= 0.01
--max-allele 4 #QC paremeter: Ignore variant with alleles over 4
--gty-qual 20.0 #QC paremeter: Exclude genotypes with the minimal genotyping quality (Phred Quality Score) per genotype < 20.0
--gty-sec-pl 20 #QC paremeter: Exclude genotypes with the second smallest normalized, Phred-scaled likelihoods for genotypes < 20
--gty-dp 8 #QC paremeter: Exclude genotypes with the minimal read depth per genotype < 4
--gty-af-ref 0.05 #QC paremeter: Exclude genotypes with the fraction of the reads carrying alternative allele >= 5% at a reference-allele homozygous genotype, usually indicating by AF in a VCF file
--gty-af-het 0.25 #QC paremeter: Exclude genotypes with the fraction of the reads carrying alternative allele <= 25% at a heterozygous genotype
--gty-af-alt 0.75 #QC paremeter: Exclude variants with the fraction of the reads carrying alternative allele <= 75% at a alternative-allele homozygous genotype
--min-obsu-rate 0.9 #QC paremeter: Exclude variants with proportion of observed number of non-missing genotypes in controls < 90% When there are control samples
--min-obsa-rate 0.9 #QC paremeter: Exclude variants with proportion of observed number of non-missing genotypes in cases < 90%
--filter-case-maf-oe 0.1 #Filter out variants with minor allele frequency over OR equal to 0.1 in the input cases to decrease the number of private variants
--db-gene refgene,gencode #Set database(s) to annotate and filter variants
--gene-feature-in 0,1,2,3,4,5,6 #Filtration according to gene features. 0-Frameshift,1-Nonframeshift,2-Startloss,3-Stoploss,4-Stopgain,5-Splicing,6-Missense
--db-filter gadexome.eas,gadgenome.eas #Set databases for allele frequency filtration. gadexome.eas is the East Asian panel of gnomAD exomes database; gadgenome.eas is the East Asian panel of gnomAD genomes database
--rare-allele-freq 0.01 #Variants with alternative allele frequency over 0.01 will be excluded
--ignore-cnv #Ignore copy number variation variants
--min-case-control-freq-ratio 3.0 #When there are control samples, this option can be used to filter out variants at which the alternative allele frequency in cases is less than that in controls multiplied by 3.0
--gene-freq-score eas #Assign the ancestrally matched population for the analyzed sample 
--qqplot #Generate QQ plot of p-values by this RUNNER1 test