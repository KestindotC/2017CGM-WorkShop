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


## RNA-seq

### Let see what the data look like


```bash

vm@ip-172-16-XXXX:~$ ls

/home/DNA-seq/
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





### To use bwa-mem aligner, run the commands below

```
bwa mem -t 1 -R '@RG\tID:group1\tSM:test\tPL:illumina\tLB:lib1\tPU:unit1' /home/DNA-seq/ref/chr17.fa ERRRWE.fastq > exome_chr22.sam
```


### Samtools Data Preprocessing
```
samtools view -Sb -o ERR687879_aln.bam ERR687879_aln.sam
```



### Resources

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

