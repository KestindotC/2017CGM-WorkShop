# 2017 CGM High-throughput genomic analysis WorkShop

National Taiwan University Center of Genomic Medicine   
Tuesday May 2, 2017   
Taipei, Taiwan   

## Workshop Schedule (permanent links to talks & handouts)

**Tuesday May 2th:**
- NGS Introduction
- DNA sequencing analysis
- RNA sequencing analysis

**Wednesday May 3th:**
- R for statistical genomic analysis
- Metagenomics
- anamiR


## Getting started

### Login using the username assigned to you

```
ssh username@172.1x.xx.x
```

### Environment and home directory
This setup is help to let you know the basic linux
- **You should be in your home directory**
- **Type `pwd` and get the path**
- **Change directory to demo_DNA `cd DNA/` and `ls` to see the contents**

```
cd ~/
pwd
cd DNA/
```


## DNA-seq

### Let's see the folder structure

```bash

vm@ip-172-16-XXXX:~$ ls

~/DNA-seq/
├── ERR687879_aln.bam
├── ERR687879_aln.sam
├── ERR687879_aln_sort.bam
├── ERR687879_aln_sort_rm.bam
├── ERR687879_aln_sort_rm.bam.bai
├── VarScan_ERR687879_chr17.log
├── VarScan_ERR687879_chr17.vcf
├── raw/
│	├── ERR_.fastq1
│	└── ERR_.fastq2
└── ref/
	├── chr17.fa
	├── chr17.fa.amb
	├── chr17.fa.ann
	├── chr17.fa.bwt
	├── chr17.fa.fai
	├── chr17.fa.pac
	└── chr17.fa.sa


```

### Analytic Pipeline
#### To use bwa-mem aligner, run the commands below

```bash
bwa mem -t 4 chr17.fa ERR687879_1.fastq ERR687879_2.fastq > ERR687879_aln.sam
```

#### Samtools Data Preprocessing
```bash
samtools view -Sb -o ERR687879_aln.bam ERR687879_aln.sam
samtools sort ERR687879_aln.bam -o ERR687879_aln_sort.bam
samtools rmdup ERR687879_aln_sort.bam ERR687879_aln_sort_rm.bam 
samtools index ERR687879_aln_sort.bam
```
#### Vairant Detection
```bash
samtools mpileup -B -f ~/DNA-seq/ref/chr17.fa \
ERR687879_aln_sort_rm.bam | \
varscan mpileup2snp \
--min-coverage 10 \
--min-var-freq 0.1 \
--min-reads2 16 \
--min-avg-qual 30 \
--output-vcf \
1>> VarScan_ERR687879_chr17.vcf \
2> VarScan_ERR687879_chr17.log

```
### Variant Annotation
We use SNPnexus to analyze SNPs in DNA sequence. And the resources are listed in the following seg


## RNA-seq

### Let's see the folder structure

```bash
vm@ip-172-16-XXXX:~$ ls

~/RNA-Seq/
├── align
│   ├── ERR188245_chrX.bam
│   ├── ERR188245_chrX.sam
│   ├── ERR188245_chrX_summarymetric.txt
│   ├── ERR188257_chrX.bam
│   ├── ERR188257_chrX.sam
│   ├── ERR188257_chrX_summarymetric.txt
│   ├── ERR188401_chrX.bam
│   ├── ERR188401_chrX.sam
│   ├── ERR188401_chrX_summarymetric.txt
│   ├── ERR188428_chrX.bam
│   ├── ERR188428_chrX.sam
│   └── ERR188428_chrX_summarymetric.txt
├── chrX_phenodata.csv
├── de
│   ├── boxplot.jpg
│   ├── chrX_gene_results.csv
│   ├── chrX_transcript_results.csv
│   ├── ERR188245
│   │   ├── e2t.ctab
│   │   ├── e_data.ctab
│   │   ├── ERR188245_chrX.gtf
│   │   ├── i2t.ctab
│   │   ├── i_data.ctab
│   │   └── t_data.ctab
│   ├── ERR188257
│   │   ├── e2t.ctab
│   │   ├── e_data.ctab
│   │   ├── ERR188257_chrX.gtf
│   │   ├── i2t.ctab
│   │   ├── i_data.ctab
│   │   └── t_data.ctab
│   ├── ERR188401
│   │   ├── e2t.ctab
│   │   ├── e_data.ctab
│   │   ├── ERR188401_chrX.gtf
│   │   ├── i2t.ctab
│   │   ├── i_data.ctab
│   │   └── t_data.ctab
│   ├── ERR188428
│   │   ├── e2t.ctab
│   │   ├── e_data.ctab
│   │   ├── ERR188428_chrX.gtf
│   │   ├── i2t.ctab
│   │   ├── i_data.ctab
│   │   └── t_data.ctab
│   └── FPKM_summary_table.csv
├── DEanalysis.R
├── quan
│   ├── ERR188245_chrX.gtf
│   ├── ERR188257_chrX.gtf
│   ├── ERR188401_chrX.gtf
│   ├── ERR188428_chrX.gtf
│   ├── mergelist.txt
│   └── stringtie_merged.gtf
├── raw
│   ├── ERR188245_chrX_1.fastq.gz
│   ├── ERR188245_chrX_2.fastq.gz
│   ├── ERR188257_chrX_1.fastq.gz
│   ├── ERR188257_chrX_2.fastq.gz
│   ├── ERR188401_chrX_1.fastq.gz
│   ├── ERR188401_chrX_2.fastq.gz
│   ├── ERR188428_chrX_1.fastq.gz
│   └── ERR188428_chrX_2.fastq.gz
└── ref
    ├── chrX.gtf
    ├── chrX_tran.1.ht2
    ├── chrX_tran.2.ht2
    ├── chrX_tran.3.ht2
    ├── chrX_tran.4.ht2
    ├── chrX_tran.5.ht2
    ├── chrX_tran.6.ht2
    ├── chrX_tran.7.ht2
    └── chrX_tran.8.ht2
```


#### Spliced alignment by HISAT2
```bash
cd ../align
hisat2 --dta -x ../ref/chrX_tran -1 ../raw/ERR188428_chrX_1.fastq.gz -2 ../raw/ERR188428_chrX_2.fastq.gz -S ERR188428_chrX.sam 2>ERR188428_chrX_summarymetric.txt
hisat2 --dta -x ../ref/chrX_tran -1 ../raw/ERR188401_chrX_1.fastq.gz -2 ../raw/ERR188401_chrX_2.fastq.gz -S ERR188401_chrX.sam 2>ERR188401_chrX_summarymetric.txt
hisat2 --dta -x ../ref/chrX_tran -1 ../raw/ERR188257_chrX_1.fastq.gz -2 ../raw/ERR188257_chrX_2.fastq.gz -S ERR188257_chrX.sam 2>ERR188257_chrX_summarymetric.txt
hisat2 --dta -x ../ref/chrX_tran -1 ../raw/ERR188245_chrX_1.fastq.gz -2 ../raw/ERR188245_chrX_2.fastq.gz -S ERR188245_chrX.sam 2>ERR188245_chrX_summarymetric.txt
```
#### Conversion to sorted BAM file by SAMtools
```bash
samtools sort -@ 1 -o ERR188428_chrX.bam ERR188428_chrX.sam
samtools sort -@ 1 -o ERR188401_chrX.bam ERR188401_chrX.sam
samtools sort -@ 1 -o ERR188257_chrX.bam ERR188257_chrX.sam
samtools sort -@ 1 -o ERR188245_chrX.bam ERR188245_chrX.sam
```
#### Transcript assembly and quantification by StringTie

**Step1: to assemble transcripts for each sample**
```bash
cd ../quan 
stringtie -G ../ref/chrX.gtf -o ERR188428_chrX.gtf -l ERR188428_chrX ../align/ERR188428_chrX.bam
stringtie -G ../ref/chrX.gtf -o ERR188401_chrX.gtf -l ERR188401_chrX ../align/ERR188401_chrX.bam
stringtie -G ../ref/chrX.gtf -o ERR188257_chrX.gtf -l ERR188257_chrX ../align/ERR188428_chrX.bam
stringtie -G ../ref/chrX.gtf -o ERR188245_chrX.gtf -l ERR188245_chrX ../align/ERR188245_chrX.bam
```
**Step2: to merge transcripts for all samples**
```bash
find ./*gtf >mergelist.txt
stringtie --merge -G ../ref/chrX.gtf -o stringtie_merged.gtf ./mergelist.txt
#Step3: to re-estimate the abundance
stringtie -e -B -G ../ref/chrX.gtf -o ../de/ERR188428/ERR188428_chrX.gtf ../align/ERR188428_chrX.bam
stringtie -e -B -G ../ref/chrX.gtf -o ../de/ERR188401/ERR188401_chrX.gtf ../align/ERR188401_chrX.bam
stringtie -e -B -G ../ref/chrX.gtf -o ../de/ERR188257/ERR188257_chrX.gtf ../align/ERR188257_chrX.bam
stringtie -e -B -G ../ref/chrX.gtf -o ../de/ERR188245/ERR188245_chrX.gtf ../align/ERR188245_chrX.bam
```

#### Differential Expression Analysis by 
```bash
cd ../
Rscript DEanalysis.R
```
## Resources

* [SNPnexus](http://www.snp-nexus.org)
* [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.ilmn)
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](http://multiqc.info)
* [BWA](http://bio-bwa.sourceforge.net/)
* [SAMtools](http://samtools.sourceforge.net/)
* [VarScan](http://varscan.sourceforge.net/)
* [SNPnexus](http://www.snp-nexus.org/)
* [Integrative Genomics Viewer (IGV)](http://www.broadinstitute.org/igv/)
* [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)
* [StringTie](https://ccb.jhu.edu/software/stringtie/)
* [R](https://cran.r-project.org)
* [R package anamiR](https://cran.r-project.org/web/packages/circlize/index.html)
* [R package ballgown](http://bioconductor.org/packages/release/bioc/html/ballgown.html)

