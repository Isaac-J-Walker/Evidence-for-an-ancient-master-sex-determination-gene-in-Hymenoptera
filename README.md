# Evidence-for-an-ancient-master-sex-determination-gene-in-Hymenoptera
Ripository for methods and data related to "Evidence for an ancient master sex determination gene in Hymenoptera"

---
# Diversity Analysis Pipeline
### Bam Generation
```
nohup fasterq-dump [SRA] --stdout --split-spot | bowtie2 --local -p 8 --rg-id SRA --rg "SM:[SRA]" -x bowtie2DatabaseName --interleaved - | samtools view -bS -@ 8 - | samtools sort -@ 8 - > outputFile.sorted.bam &
samtools index outputFile.sorted.bam
```

### Determining Mean & Median Read Depth
```
# mean of whole genome
samtools depth -a file.bam | awk '{c++; s+=$3} END {print s/c}'

# mean of ANTSR region
samtools coverage -r <chrom>:<ANTSR_start>-<ANTSR_end> file.bam

# median of whole genome
samtools depth -a file.bam | awk '{print $3}' | sort -n | awk 'NF{a[i++]=$1} END {print a[int(i/2)]}'
```

### VCF Generation
```
# for single pooled BAMs
freebayes -f reference.fa --pooled-continuous pooled.bam > pooled.vcf

# for joint calling multiple BAMs
freebayes -f reference.fa --ploidy 2 male1.bam female1.bam female2.bam > joint_call.vcf

# for single organism BAM
freebayes -f reference.fa --ploidy 2 female.bam > female.vcf
```

### VCF Filtering
```
# for pooled VCFs
vcftools --vcf input.vcf --minQ 30 --minDP 5 --maxDP [mean read depth x 3] --recode --out filtered.vcf

# for joint called VCFs
vcftools --vcf input.vcf --minQ 30 --minDP 5 --maxDP [mean read depth of the sum of all samples x3] --recode --minGQ 20 --out filtered.vcf

# for single organism VCFs
vcftools --vcf input.vcf --minQ 30 --minDP 5 --maxDP [mean read deoth x 3] --minGQ 20 --recode --out filtered.vcf
```

### Diversity File Generation
```
vcftools --vcf filtered.vcf --window-pi 10000 --out filtered.pi
```

### Significance Calculation 
The PI file, the highest PI window, and the min number of windows needed to cover the entire ANTSR region, are input into the get-Perc-pi.R script. This calculates a percentile and corrected p-value for the window.
```
Rscript get-Perc-pi.R filtered.pi [window index] [number of windows covering region]
```

### Heterozygous Analysis
For heterozygous analysis, SNPs are filtered for heterozygous status using bcftools, and diversity and significance calculations were run on heterozygous VCF.
```
bcftools view -g het -v snps filtered.vcf > het_snps.vcf

vcftools --vcf het_snps.vcf --window-pi 10000 --out het_snps.pi

Rscript get-zScore-pi.R het_snps.pi [window index] [number of windows covering region]
```

---
# Phylogenetic Reconstruction
18SrRNA for each species were concatinated into a single fasta. This was then used as input for clusalo for multi-sequence alignment. The output of clustalo was then passed to iqtree2 for phylogenetic reconstruction.
```
clustalo -i 18s_concat.fasta -o 18s_msa.fasta --force --outfmt=fasta
iqtree2 -s 18s_aligned.fasta -redo -st DNA -m MFP -bb 1000 -alrt 1000
```
