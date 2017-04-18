# 2017 CGM WorkShop DNA-seq


---
# Getting started
### Login using the username assigned to you

```
ssh username@ron.sr.unh.edu
```
### Setting up a temp directory on Ron
This setup is specific to Ron and may not be needed for other servers.

- **You should be in your home directory**
- **Type `pwd` and get the path**
- **Change directory to demo_DNA `cd demo_DNA/` and `ls` to see the contents**

```
export TMPDIR=/home/gomre/USERNAME/temp
```


### To use bwa aligner, run the commands below

```
bwa mem -t 1 -R '@RG\tID:group1\tSM:test\tPL:illumina\tLB:lib1\tPU:unit1' /home/vm/demo_DNA/ref/chr22_ref.fa exome_chr22.fastq > exome_chr22.sam
```


Samtools Preprocessing
```
samtools view -Sb -o exome_chr22.bam exome_chr22.sam
```

Now, you should be in the QIIME environment and ready to run through the pipeline. 

在做資料時做分析上，我們所使用的Pipeline大致上如下，在每一個Part都會有相對應的
- Raw Data Analysis
- Genome Mapping
- Variant Calling
- Annotation
