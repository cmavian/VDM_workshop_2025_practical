# VDM_workshop_2025_practical
Virus discovery pipeline created by Erin Harvey, Carla Mavian, Nokuzotha Nkiwane, Eduan Wilkinson. #added in alphabetical order

### Metagenomic Workflow
In this tutorial we will learn how to taxonomically classify and visualize our metagenomic reads obtained with Illumina using the following programs:

1. [Fastp](https://github.com/OpenGene/fastp)
2. [Trimmomatic](https://github.com/usadellab/Trimmomatic)
3. [MEGAHIT](https://www.metagenomics.wiki/tools/assembly/megahit)
4. [Diamond](https://github.com/bbuchfink/diamond?tab=readme-ov-file)

<figure>
    <img src="workflow.png" width="920" height="1200">
    <figcaption>Virus discovery pipeline Workflow by Nokuzotha Nkiwane </figcaption>
</figure>


### Connecting to the server to run the analysis

Using MobaXterm connect to Host: 

```
xxxxx
```

```
xxxxxxx
```



### Setting up our folder for the analysis

orientate yourself 

```
pwd
ls
```

Let's make a working folder
 
```
mkdir metagenomics
```

Go into the folder

```
cd metagenomics
```

4. let's make other folders 

```
mkdir data results scripts
ls
```
4. let's go into data and copy our files there

```
cd data
cp xxxxxxx/lib* .
ls
```

### 1. FastP PRE Trimmomatic

```
cd metagenomics/scripts
nano 
```

### 2. Trimmomatic

```
cd metagenomics/scripts
nano 02.trimmomatic.sh
```
inside you will write:

```
#!/bin/env bash

module trimmomatic

input="/analyses/users/nokuzothan/discovery_pipeline/init_tools/trimmo_test/input"
output="/analyses/users/nokuzothan/discovery_pipeline/init_tools/trimmo_test/output"
ADAPTERS="/analyses/software/programs/trimmomatic/0.39/adapters.fa"

#first get the SRA ID as a variable using %%
for file in $input/*_1.fastq;

do
file_name=$(basename "$file")
id=${file_name%%_1.fastq}

R1="$input/${id}_1.fastq"
R2="$input/${id}_2.fastq"

P1="${output}/${id}_P1.fastq"
U1="${output}/${id}_U1.fastq"
P2="${output}/${id}_P2.fastq"
U2="${output}/${id}_U2.fastq"


trimmomatic PE -threads 8 -phred33 -summary "${output}/${id}_statsSummary.txt" $R1 $R2 $P1 $U1 $P2 $U2 ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25 

2>>"${output}/${id}.log" 1>>"${output}/${id}.log"

echo "Trim run for $id complete"
done
```
and save it by "Write Out" in Nano:
#### Ctrl + O 

If you want to save the changes to the existing file, simply press 
#### Enter

After saving, if you wish to exit the editor, press:
#### Ctrl + X

### 4. FastP POST Trimmomatic


### 5. MEGAHIT


### 6. Diamond


### 7. Taxonomy


### 8. Lineage filter


### 9. BlastN


### 10. BlastX
