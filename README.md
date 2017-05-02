# 2017 CGM High-throughput genomic analysis WorkShop

National Taiwan University Center of Genomic Medicine   
Tuesday May 2,3 2017   
Taipei, Taiwan   


## Workshop Schedule (permanent links to talks & handouts)

**Tuesday May 2th:**
- NGS Introduction 
- DNA sequencing analysis [[Script]](#dna-seq) [[Slides]](http://slides.com/kestinchang/dnaseq-2/fullscreen)
- RNA sequencing analysis [[Script]](#rna-seq) [[Slides]](https://www.dropbox.com/s/m2qp2yq6ouh01yz/RNA-seq_workshop_2017_v4.pdf?dl=0)
- [CellExpress](http://cellexpress.cgm.ntu.edu.tw)

**Wednesday May 3th:**
- R for statistical genomic analysis [[Script]](#r-for-statistical-genomic-analysis)
- Metagenomics
- anamiR [[Script]](#anamir) [[Slides]](https://www.dropbox.com/s/s61oytmauewesvz/anamiR_training.pptx?dl=0)

> All tool links, websites, citations were list in [Resources](#resources).  

## Getting started

### Login using the username assigned to you
```
ssh username@140.112.129.20
```

### Environment and home directory
Basic linux command
- **You should be in your home directory after login**
- **Type `pwd` and get the current path**
- **Change directory to DNA-seq `cd DNA-seq/` and `ls` to list the contents in the folder**

```
[user1@localhost ~]$ ls
DNA-Seq  miniconda3  RNA-Seq

[user1@localhost ~]$ pwd
/training_home/user1

[user1@localhost ~]$ cd DNA-seq/
[user1@localhost ~]$ ls
ERR687879_1.fastq  ERR687879_2.fastq  raw  README  ref  results
```


## DNA-seq

### Let's see the folder structure

```bash

~/DNA-seq/
├── ERR687879_1.fastq
├── ERR687879_2.fastq
├── raw
│   ├── ERR687879_1.fastq
│   └── ERR687879_2.fastq
├── README
├── ref
│   ├── chr17.fa
│   ├── chr17.fa.amb
│   ├── chr17.fa.ann
│   ├── chr17.fa.bwt
│   ├── chr17.fa.fai
│   ├── chr17.fa.pac
│   └── chr17.fa.sa
└── results
    ├── ERR687879_aln.bam
    ├── ERR687879_aln.sam
    ├── ERR687879_aln_sort.bam
    ├── ERR687879_aln_sort_rm.bam
    ├── ERR687879_aln_sort_rm.bam.bai
    ├── VarScan_ERR687879_chr17.log
    └── VarScan_ERR687879_chr17.vcf

```

### Analytic Pipeline
#### To use bwa-mem aligner, run the commands below

```bash
bwa mem ref/chr17.fa ERR687879_1.fastq ERR687879_2.fastq > ERR687879_aln.sam
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
samtools mpileup -B -f ref/chr17.fa \
ERR687879_aln_sort_rm.bam | \
varscan mpileup2snp \
--min-coverage 40 \
--min-var-freq 0.4 \
--min-reads2 16 \
--min-avg-qual 30 \
--output-vcf \
1>> VarScan_ERR687879_chr17.vcf \
2> VarScan_ERR687879_chr17.log

```
#### Variant Annotation
We can use [wANNOVAR](http://wannovar.wglab.org) to annotate SNPs calling result (vcf files).
So please download the vcf files from the FTP or use demo_vcf [here](query.output.genome_summary.csv) to run the wANNOVAR.


## RNA-seq

### Let's see the folder structure

```bash

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
│   ├── FPKM_summary_table.csv 
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
│   └── ERR188428
│       ├── e2t.ctab
│       ├── e_data.ctab
│       ├── ERR188428_chrX.gtf
│       ├── i2t.ctab
│       ├── i_data.ctab
│       └── t_data.ctab
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

### Analytic Pipeline
#### link to http://140.112.129.20
```bash
[training@localhost ~]$ pwd
/training_home/training
[training@localhost ~]$ cd ../user#
[training@localhost ~]$ ls
RNA-Seq
[training@localhost ~]$ cd RNA-Seq
[training@localhost ~]$ ls
align chrX_phenodata.csv de DEanalysis.R quan raw ref
```
#### Spliced alignment using HISAT2
```bash
cd ./align
hisat2 -p 2 --dta -x ../ref/chrX_tran -1 ../raw/ERR188428_chrX_1.fastq.gz -2 ../raw/ERR188428_chrX_2.fastq.gz -S ERR188428_chrX.sam 2>ERR188428_chrX_summarymetric.txt
hisat2 -p 2 --dta -x ../ref/chrX_tran -1 ../raw/ERR188401_chrX_1.fastq.gz -2 ../raw/ERR188401_chrX_2.fastq.gz -S ERR188401_chrX.sam 2>ERR188401_chrX_summarymetric.txt
hisat2 -p 2 --dta -x ../ref/chrX_tran -1 ../raw/ERR188257_chrX_1.fastq.gz -2 ../raw/ERR188257_chrX_2.fastq.gz -S ERR188257_chrX.sam 2>ERR188257_chrX_summarymetric.txt
hisat2 -p 2 --dta -x ../ref/chrX_tran -1 ../raw/ERR188245_chrX_1.fastq.gz -2 ../raw/ERR188245_chrX_2.fastq.gz -S ERR188245_chrX.sam 2>ERR188245_chrX_summarymetric.txt
```
#### Conversion to sorted BAM files using SAMtools
```bash
samtools sort -@ 2 -o ERR188428_chrX.bam ERR188428_chrX.sam
samtools sort -@ 2 -o ERR188401_chrX.bam ERR188401_chrX.sam
samtools sort -@ 2 -o ERR188257_chrX.bam ERR188257_chrX.sam
samtools sort -@ 2 -o ERR188245_chrX.bam ERR188245_chrX.sam
```
#### Transcript assembly and quantification using StringTie

1. to assemble transcripts for each sample
```bash
cd ../quan 
stringtie -G ../ref/chrX.gtf -o ERR188428_chrX.gtf -l ERR188428_chrX ../align/ERR188428_chrX.bam
stringtie -G ../ref/chrX.gtf -o ERR188401_chrX.gtf -l ERR188401_chrX ../align/ERR188401_chrX.bam
stringtie -G ../ref/chrX.gtf -o ERR188257_chrX.gtf -l ERR188257_chrX ../align/ERR188428_chrX.bam
stringtie -G ../ref/chrX.gtf -o ERR188245_chrX.gtf -l ERR188245_chrX ../align/ERR188245_chrX.bam
```
2. to merge transcripts for all samples
```bash
find ./*.gtf >mergelist.txt
stringtie --merge -G ../ref/chrX.gtf -o stringtie_merged.gtf ./mergelist.txt
```
3. to re-estimate the abundance
```bash
stringtie -e -B -G ./stringtie_merged.gtf -o ../de/ERR188428/ERR188428_chrX.gtf ../align/ERR188428_chrX.bam
stringtie -e -B -G ./stringtie_merged.gtf -o ../de/ERR188401/ERR188401_chrX.gtf ../align/ERR188401_chrX.bam
stringtie -e -B -G ./stringtie_merged.gtf -o ../de/ERR188257/ERR188257_chrX.gtf ../align/ERR188257_chrX.bam
stringtie -e -B -G ./stringtie_merged.gtf -o ../de/ERR188245/ERR188245_chrX.gtf ../align/ERR188245_chrX.bam
```
4. view vcf file 
```bash
head –n 5 ERR188245_chrX.gtf 
```
#### Differential Expression Analysis using ballgown R package
```bash
cd ../
Rscript DEanalysis.R
```

---


## R for statistical genomic analysis

### General Workflow
1. Import data
```R
library(SNPassoc)
data(SNPs)
SNPs[1:10,1:10]
ncol(SNPs)
nrow(SNPs)
```
2. **Single** SNP analysis 
```R
mySNP <- snp(SNPs$snp10001,sep="")
SNPs$snp10001
mySNP
summary(mySNP)
plot(mySNP,label="snp10001",col="darkgreen")
```
3. Analysis for **many** SNPs
```R
myData <- setupSNP(data=SNPs,colSNPs=6:40,sep="")
myData[1:5,6:10]
plot(myData, which=20)
summary(myData)
```
4. Hardy-Weinberg Equilibrium (HWE)
```R
res <- tableHWE(myData)
res
```
5. Remove SNPs violating HWE and MAF<0.05
```R
sum <- summary(myData)
delete_index <- subset(c(1:nrow(sum)),sum[,2]>95|sum[,3]<0.05)
Data <- myData[,-(delete_index+5)]
summary(Data)
```
6. Single SNP association
```R
association(casco~snp10001, data=Data)
association(casco~snp10001, model=c("log-additive"), data=Data)
```
7. Adjusted association analysis
```R
association(casco~sex+snp10001,model=c("log-additive"), data=Data)
```
8. Stratified analysis
```R
association(casco~snp10001+strata(sex), model=c("log-additive"), data=Data)
```

9. Interaction analysis
```R
association(casco~snp10001*sex, model=c("log-additive"), data=Data)
```
10. GWAS
```R
result <- WGassociation(casco,data=Data)
result
plot(result)
summary(result)
result <- WGassociation(casco, data=Data, model="log-add")
result <- WGassociation(casco~sex, data=Data, model="log-add")
```

11. Bonferroni correction
```R
ans <- WGassociation(casco,data=Data, model="all")
Bonferroni.sig(ans, model="dominant", alpha=0.05, include.all.SNPs=FALSE)
```
12. GWAS
```R
data(HapMap)
summary(resHapMap)
plot(resHapMap)
plot(resHapMap, whole=FALSE)
```
---   

## anamiR
### General Workflow
1. Import data
```R
library(anamiR)
data(mrna)
data(mirna)
data(pheno.mirna)
data(pheno.mrna)
```
2. SummarizedExperimaent object
```R
mirna_se <- SummarizedExperiment(assays = SimpleList(counts=mirna), colData = pheno.mirna)
mrna_se <- SummarizedExperiment(assays = SimpleList(counts=mrna), colData = pheno.mrna)
```
3. Differential Expression analysis
```R
mirna_d <- differExp_discrete(se = mirna_se, class = "ER", method = "limma", log2 = FALSE, p_value.cutoff = 0.05, logratio = 0.5, p_adjust.method = "BH")
mrna_d <- differExp_discrete(se = mrna_se, class = "ER", method = "limma", log2 = FALSE, p_value.cutoff = 0.05, logratio = 0.5, p_adjust.method = "BH")
```
4. miRNA ID conversion
```R
mirna_21 <- miR_converter(data = mirna_d, remove_old = TRUE, original_version = 17, latest_version = 21)
```
5. Correlation abalysis
```R
cor <- negative_cor(mrna_data = mrna_d, mirna_data = mirna_21, cut.off = -0.5, method = "pearson")
```
6. Database Intersection
```R
sup <- database_support(cor_data = cor, Sum.cutoff = 2, org = "hsa")
```
7. Heatmap visualization
```R
heat_vis(cor, mrna_d, mirna_21)
```
8. Functional analysis
```R
enr <- enrichment(data_support = sup, per_time = 5000)
```
### GSEA Workflow
1. Load data
```R
aa <- system.file("extdata", "GSE19536_mrna.csv", package = "anamiR")
mrna <- data.table::fread(aa, fill = T, header = T)
bb <- system.file("extdata", "GSE19536_mirna.csv", package = "anamiR")
mirna <- data.table::fread(bb, fill = T, header = T)
cc <- system.file("extdata", "pheno_data.csv", package = "anamiR")
pheno.data <- data.table::fread(cc, fill = T, header = T)
```
2. transform data format
```R
mirna_name <- mirna[["miRNA"]]
mrna_name <- mrna[["Gene"]]
mirna <- mirna[, -1]
mrna <- mrna[, -1]
mirna <- data.matrix(mirna)
mrna <- data.matrix(mrna)
row.names(mirna) <- mirna_name
row.names(mrna) <- mrna_name
pheno_name <- pheno.data[["Sample"]]
pheno.data <- pheno.data[, -1]
pheno.data <- as.matrix(pheno.data)
row.names(pheno.data) <- pheno_name
```
#### Extra steps
```R
mirna <- mirna[, 50:70]
mrna <- mrna[, 50:70]
pheno.data <- as.matrix(pheno.data[50:70, ])
colnames(pheno.data) <- "ER"
```
3. SummarizedExperiment object
```R
mirna_21 <- miR_converter(mirna, original_version = 17)
mirna_se <- SummarizedExperiment(assays = SimpleList(counts=mirna_21), colData = pheno.data)
mrna_se <- SummarizedExperiment(assays = SimpleList(counts=mrna), colData = pheno.data)
```
4. GSEA_ana
```R
table <- GSEA_ana(mrna_se = mrna_se, mirna_se = mirna_se, class = "ER", pathway_num = 2)
names(table)
miRNA_path1 <- table[[1]]
Gene_path1 <- table[[2]]
```
5. GSEA_res
```R
result <- GSEA_res(table, pheno.data = pheno.data, class = "ER", DE_method = "limma", cor_cut = -0.3)
names(result)
result_path1 <- result[[1]]
```

## Resources

* [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.ilmn)
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](http://multiqc.info)
* [BWA](http://bio-bwa.sourceforge.net/)
* [SAMtools](http://samtools.sourceforge.net/)
* [VarScan](http://varscan.sourceforge.net/)
* [wANNOVAR](http://wannovar.wglab.org)
* [Integrative Genomics Viewer (IGV)](http://www.broadinstitute.org/igv/)
* [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)
* [StringTie](https://ccb.jhu.edu/software/stringtie/)
* [RNA-seq Data for Exercise](http://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/)
* [CellExpress](http://cellexpress.cgm.ntu.edu.tw)
* [R](https://cran.r-project.org)
* [R package anamiR](https://bioconductor.org/packages/release/bioc/html/anamiR.html)
* [R package ballgown](http://bioconductor.org/packages/release/bioc/html/ballgown.html)
* [R package SNPassoc](https://cran.r-project.org/web/packages/SNPassoc/index.html)
* [The Applications of Big Data Analysis in Genomic Research 2016](http://www.airitilibrary.com/Publication/alDetailedMesh?docid=10281916-201611-201612070002-201612070002-609-619)
