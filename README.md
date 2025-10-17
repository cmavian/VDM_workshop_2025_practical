<p><b><h1>VDM_workshop_2025_practical</h1></b></p>
<h3>Virus discovery pipeline</h3>
<i>created by</i> Erin Harvey, Carla Mavian, Nokuzotha Nkiwane, Eduan Wilkinson.
<i>(in alphabetical order)</i>

<b><h3>Metagenomic Workflow</h3></b>
In this tutorial we will learn how to taxonomically classify and visualize our metagenomic reads obtained with Illumina using the following programs:
<list>
1. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. [Trimmomatic](https://github.com/usadellab/Trimmomatic)
3. [MEGAHIT](https://www.metagenomics.wiki/tools/assembly/megahit)
4. [Diamond](https://github.com/bbuchfink/diamond?tab=readme-ov-file)
</list>
<figure>
    <img src="workflow.png" width="920" height="1200">
    <figcaption><b>Figure 1.</b> Virus discovery pipeline Workflow by Nokuzotha Nkiwane </figcaption>
</figure>

<b><h3>Connecting to the server to run the analysis</h3></b>
MacBook's and Linux computers come with terminal windows as part of the system.
They do not require any additional programs to connect to a server.
On Windows computers, you have to emulate a terminal.
Multiple free programs provide that function, e.g [MobaXterm](https://mobaxterm.mobatek.net/), [putty](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html), [GitBash](https://git-scm.com/downloads), [Terminus](https://termius.com/) etc.
Using one of the emulators, open a terminal and connect to host by typing:

<pre><code>ssh username@ceri.sarmc.ac.za</code></pre>

First time connecting you will be ask if you trust the connection. 
You have to type 'yes' and press enter.
Then you will be asked to provide password.
You will be provided one together with your username at the beginning of the workshop.
<b>Passwords are case sensitive and invisiable when typing !!!</b>
Therefore, it is advisable to copy & paste the password.
Successful connection will print the server's logo on your terminal.<br>
<i>On average on a server you only have 2-5 tries to type the password afterwhich you'll be blocked for at least 30 minutes.
It is better to close the terminal and open it again and try again after the second failed password try.</i><br>

<b><h3>Setting up our folder for the analysis</h3></b>
<ol start=1>
<li>Commands used for orentation on the server</li>
<pre><code></code>pwd</code></pre>
<pre><code>ls</code></pre>
<pre><code>ls -Fal *</code></pre>
<pre><code>ll *</code></pre>

<li>Moving around on the server</li>
<pre><code>cd</code></pre>
<pre><code>cd ..</code></pre>
<pre><code>cd ../..</code></pre>
<pre><code>cd ~</code></pre>
<pre><code>cd workspace</code></pre>

<li>Abssolute vs relative path</li>
<pre><code>cd /home/vdw00</code></pre>
<pre><code>cd ../</code></pre>
<pre><code>cd ./vdw00</code></pre>

<li>First step will be creating a working folder in workspace folder and moving into the freshly created directory</li>
<pre><code>cd workspace
mkdir -p metagenomics
cd metagenomics</code></pre>
To check if the location was changed, print working directory
<pre><code>pwd</code></pre>

<li>Next, we will create other folders which will, help keep the analysis sorted</li> 
<pre><code>mkdir -p data results scripts</code></pre>
Then change permissions on all of the directories and in them:
<pre><code>chmod -R a=rwx ./</code></pre>
To check if everything was created correctly, list the current directory.
<pre><code>ll *</code></pre>

<li>Let's go into data directory and copy our files</li><br>
<pre><code>cd data</code></pre>
Coping from website:
<pre><code>wget <a href="https://">https://</a></code></pre>
Alternatively, if the internet connection is slow, you can copy the data from backup location:
<pre><code> cp -R /analyses/vdworkshop/.backup/data/* ./ </code></pre>
Then change permission on the files and confirm the change on the data in the location:
<pre><code>chmod a=rwx *
ll *</code></pre>
</ol>

<b><h3>Analysis</h3></b>
<ol start=1>
<h4><li>FastQC pre-Trimmomatic</li></h4>
Let's create a script to execute this step
<pre><code>cd ../scripts</code></pre>
(<i>if you got lost, you can use the absolute path:</i> <code>/analyses/vdworkshop/${USER}/metagenomics/scripts</code>)<br><br>

We will use <i>nano</i> text editor. You can create empty file and open it to edit at the same time
<pre><code>nano 01.fastqc_pretrim.sh</code></pre>
<i>you can copy & paste the script text directly into the open document.</i><br>

<pre><code>#!/bin/env bash

THR=5
### enabling conda environment and fastqc program
ON="module miniconda && conda activate fastqc"
eval $ON

### input and output directories
workdir=`realpath $(pwd) 2>dev/null`'/../'
input=${workdir}'/data'
output=${workdir}'/results/01.fastqc_pretrim'

### make output directory if it doesn't exist
if ! [[ -d ${output} ]]; then mkdir -p -m a=rwx ${output}; fi

echo 'Running FastQC pre-trim'

### run fastqc
fastqc -t ${THR} -o ${output}/ ${input}/*.f*q*
echo "Pretrim FastQC complete"

### deactivating fastqc program
OFF='conda deactivate'
eval ${OFF}

exit 0;
</code></pre>

To save your script press ctrl+X then Y and ENTER<br>
This will "override" the 01.fastqc_pretrim.sh (which was empty on opening)<br>
You can change the name after pressing Y <br><br>
Before running the script you have to change permissions
<pre><code>chmod a=rwx 01.fastqc_pretrim.sh</code></pre>
To run the script type:
<pre><code>bash 01.fastqc_pretrim.sh</code></pre>

<h4><li>Trimmomatic</li></h4>
Let's create another script for trimmomatic
<pre><code>nano 02.trimmomatic.sh</code></pre>
In the script, write (or copy & paste) the fallowing commands:<br>

<pre><code>#!/bin/env bash

THR=5
### activating program
ON='module trimmomatic'
eval ${ON}
 
workdir=`realapth $(pwd) 2>/dev/null`'/../'
input=${workdir}'/data'
output=${workdir}'/results/02.trim_output'
ADAPTERS=${input}'/adapters.fa'

### getting the SRA IDs and creating output files
for FOW in (ls ${input}/*.f*q* | grep -Ei "_r?1"); do
  REV=`echo ${FOW} | sed -r 's/\_(r|R)?1/\_\12/'`;
  ID=`basename ${FOW} | cut -d '_' -f1`
  P1=${output}/${ID}'_1.P.fq.gz'
  U1=${output}/${ID}'_1.U.fq.gz'
  P2=${output}/${ID}'_2.P.fq.gz'
  U1=${output}/${ID}'_2.U.fq.gz'

  ### running program
  trimmomatic PE -threads ${THR} -phred33 -summary ${output}/${ID}'_statsSummary.txt' \
    ${FOW} ${REV} ${P1} ${U1} ${P2} ${U2} \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25 \
    2>>${output}/${ID}.log 1>>${output}/${ID}.log

  echo "Trimming for ${ID} complete"
done
exit 0;
</code></pre>

Save your script by pressing ctrl+X then Y and ENTER<br><br>
Before running the script you have to change permissions
<pre><code>chmod a=rwx 02.trimmomatic.sh</code></pre>
To run the script type:
<pre><code>bash 02.trimmomatic.sh</code></pre>

<h4><li>FastQC post-Trimmomatic</li></h4>
Open new script file
<pre><code>nano 03.fastqc_posttrim.sh</code></pre>

<pre><code>#!/bin/env bash

THR=5
### enabling conda environment and fastqc program
ON="module miniconda && conda activate fastqc"
eval $ON

### input and output directories
workdir=`realpath $(pwd) 2>dev/null`'/../'
input=${workdir}'/results/02.trim_output'
output=${workdir}'/results/03.fastqc_posttrim'

### make output directory if it doesn't exist
if ! [[ -d ${output} ]]; then mkdir -p -m a=rwx ${output}; fi

echo 'Running FastQC post-trim'

### run fastqc
fastqc -t ${THR} -o ${output}/ ${input}/*.P.fq.gz
echo "Post-trim FastQC complete"

### deactivating fastqc program
OFF='conda deactivate'
eval ${OFF}

exit 0;
</code></pre>

Save your script by pressing ctrl+X then Y and ENTER and change permissions
<pre><code>chmod a=rwx 03.fastqc_posttrim.sh</code></pre>
To run the script type:
<pre><code>bash 03.fastqc_posttrim.sh</code></pre>

<h4><li>MultiQC</li></h4>
Let's create a script for this step
<pre><code>nano 04.multiqc.sh</code></pre>
In the script, write (or copy & paste) the fallowing commands:<br>

<pre><code>#!/bin/env bash

THR=5
### enabling conda environment and fastqc program
ON="module miniconda && conda activate fastqc"
eval $ON

### input and output directories
workdir=`realpath $(pwd) 2>dev/null`'/../'
input_pre=${workdir}'/results/01.fastqc_pretrim'
input_post=${workdir}'/results/03.fastqc_posttrim'
output=${workdir}'/results/04.multiqc'

### make output directory if it doesn't exist
if ! [[ -d ${output} ]]; then mkdir -p -m a=rwx ${output}; fi

### runnig multiqc
### this will join both results, we would rather want to keep them separate to compare
#multiqc ${input_pre} ${input_post} -o ${output}

multiqc ${input_pre} -o ${output}/multiqc_pretrim
multiqc ${input_post} -o ${output}/multiqc_posttrim

### deactivating multiqc program
OFF='conda deactivate'
eval ${OFF}

exit 0;
</code></pre>

Save your script by pressing ctrl+X then Y and ENTER and change permissions
<pre><code>chmod a=rwx 04.multiqc.sh</code></pre>
To run the script type:
<pre><code>bash 04.multiqc.sh</code></pre>

#################################################################################################
<h4><li>MEGAHIT</li></h4>

```
cd metagenomics/scripts
nano 05.megahit.sh
```

```
#!/bin/env bash

#load modules
ON="module miniconda && conda activate megahit"
eval $ON

#directories
input="/analyses/users/nokuzothan/disc_pipe/init_tools/fastp_test/output"
output="/analyses/users/nokuzothan/disc_pipe/init_tools/megahit/output/default"

#remove output directory if already exists and make new one
if [[ -d ${output} ]]; then
	rm -rf ${output}
fi
mkdir -p ${output}

#run megahit
for file in ${input}/*__R1_001_out.fastq;do 

	file_name=$(basename "$file")
	id=${file_name%%__R1_001_out.fastq}
	
	R1=${input}/${id}__R1_001_out.fastq
	R2=${input}/${id}__R2_001_out.fastq

	megahit --verbose -t 10 -1 ${R1} -2 ${R2} -o ${output}/${id}

	#after megahit run, prepend sample name to contigs (allows easier tracking of which contigs came from which sample in downstream analysis)
	awk -v prefix="${id}_" '
        /^>/ {$0=">" prefix substr($0,2)} {print}
    ' ${output}/${id}/final.contigs.fa > ${output}/${id}/sample.contigs.fa

done

echo "Megahit assmebly completed successfully :-) and contigs named with sample name"
```

### 6. Diamond

```
cd metagenomics/scripts
nano 06.1.diamond_rvdb.sh
```

```
#!/bin/env bash

#load module
module diamond

# variables and directories
wdir="/analyses/users/nokuzothan/disc_pipe"
cdir="${wdir}/init_tools/diamond"
input_reads_dir="${wdir}/init_tools/megahit/output/default"
db="${wdir}/ncbidb/RVDB/v30.0/U-RVDBv30.0-prot.fasta"
output="${cdir}/output/RVDB"
threads=$((`/bin/nproc` -2))

#clear existing output directory if any, make new output directory 
if [[ -e $output ]]; then
  rm -rf ${output} 
fi
mkdir -p ${output}

#make diamond protein database
diamond makedb --in ${db} -d ${output}/nr

#loop through each of the files created in megahit output directory to find final.contigs.fa files and run diamond
for folder in ${input_reads_dir}/*; do

  if [[ -d ${folder} ]]; then
    sample=$(basename ${folder})
    contigs=${folder}/sample.contigs.fa


    #alignment using blastx (exclude --min-score because it overrides the evalue (acc. to manual))
    if [[ -f ${contigs} ]]; then
    sample_out=${output}/${sample}_rvdb.matches.m8
    diamond blastx -d ${output}/nr.dmnd \
    -q ${contigs} \
    --out ${sample_out} \
    --threads ${threads} \
    --evalue 1E-5 \
    --outfmt 6 qseqid qlen sseqid stitle pident length evalue bitscore \
    --id 80 \
    --strand both \
    --unal 0 \
    --mp-init 

    else 
      echo "Contigs file for ${sample} not found."
    fi
  fi
done

exit 0
```

```
cd metagenomics/scripts
nano 06.2.nr_viral_seqs.sh
```

```
#!/bin/env bash

# variables and directories
wdir="/analyses/users/nokuzothan/disc_pipe"
cdir="${wdir}/init_tools/diamond/input"
viruses_csv="${cdir}/virus_taxonomy_lvls.csv"
names="${cdir}/viral_namess.txt"
db="${wdir}/ncbidb/fasta/nr.faa"
db_fasta="${cdir}/ncbi_fasta.faa"
threads=$((`/bin/nproc` -2))


#filter csv for viral sequences
awk -F',' '{print $3}' ${viruses_csv} >> ${names}

#filter nr database for viral sequences 
while read -r virus; do
    awk -v name="${virus}" '
        BEGIN {IGNORECASE=1}
        /^>/ {ON = index($0, name) > 0}
        ON {print}
    ' ${db} >> ${db_fasta}
done < ${names}

```

```
cd metagenomics/scripts
nano 06.3.diamond_ncbi.sh
```

```
#!/bin/env bash

#load module
module diamond

# variables and directories
wdir="/analyses/users/nokuzothan/disc_pipe"
cdir="${wdir}/init_tools/diamond"
input_reads_dir="${wdir}/init_tools/megahit/output/default"
output="${cdir}/output/NCBI"
db_fasta="${cdir}/input/ncbi_fasta.faa"
threads=$((`/bin/nproc` -2))

#clear existing output directory if any, make new output directory 
if [[ -e $output ]]; then
  rm -rf ${output} 
fi
mkdir -p -m a=rwx ${output}


#make diamond protein database
diamond makedb --in ${db_fasta} -d ${output}/nr

#loop through each of the files created in megahit output directory to find final.congtigs.fa files and run diamond
for folder in ls ${input_reads_dir}/*; do

  if [[ -d ${folder} ]]; then
    sample=$(basename ${folder})
    contigs=${folder}/sample.contigs.fa


    #alignment using blastx
    if [[ -s ${contigs} ]]; then
    sample_out=${output}/${sample}.matches.m8
    diamond blastx -d ${output}/nr.dmnd \
    -q ${contigs} \
    --out ${sample_out} \
    --threads ${threads} \
    --evalue 1E-5 \
    --outfmt 6 qseqid qlen sseqid stitle pident length evalue bitscore \
    --id 80 \
    --strand both \
    --unal 0 \
    --mp-init 

    else 
      echo "Contigs file for ${sample} not found."
    fi
  fi
done
```


### 7. Taxonomy

```
cd metagenomics/scripts
nano 07.taxonomy.sh
```

```
#!/bin/env bash

#directories

xxxxxxxxx

#function to get metadata from eutils
function get_meta {
    contig=$1
    length=$2
    acc_id=$3
    columns=$4
    output=$5


    #print ncbi page of protein accession and parse taxonomic id for use in taxonkit for lineage
    url1="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc_id}&rettype=gb&retmode=text"
    info=$(curl -L --retry 3 --connect-timeout 5 -N -# ${url1})
    #taxonomic number
    tax=$(echo "${info}" | awk '/\/db_xref/ { match($0, /taxon:([0-9]+)/, tax_id); print tax_id[1] }')

    #put NA if no taxonomy id found
    if [[ -z "${tax}" ]]; then
        tax="NA"
    fi

    #print output
    echo -e "${contig}\t${length}\t${acc_id}\t${rest}\t${tax}" >>"${output}"
    sleep 0.34

} 

#clear output file before each run
> ${output_tsv}

while IFS=$'\t' read -r col1 col2 col3 rest;do
    echo "[${col3}]"
    tmpfile=$(mktemp)
        get_meta "${col1}" "${col2}" "${col3}" "${rest}" "${tmpfile}"
        cat "${tmpfile}" >> "${output_tsv}"
        rm "${tmpfile}"
done < ${out}/acc_ids.txt
```

### 8. Lineage filter

```
cd metagenomics/scripts
nano 08.lineage_filter.sh
```

```
#!/bin/env bash

#load modules
ON="module miniconda && conda activate viral_pipeline"
eval ${ON}

#directories used
wdir="/analyses/users/nokuzothan/disc_pipe/init_tools"
out="${wdir}/play"
acc_tax_id="${out}/acc_tax_id.tsv"
lineage_out="${out}/lineages.tsv"
u_match_out="${out}/unique_contig_ids.txt"
contig_matches="${out}/contig_matches.tsv"
output_fa="${out}/blast_fasta.fa"
mega_conts="${wdir}/megahit/output/default"
fin_fasta=${mega_conts}/*/sample.contigs.fa


#clear lineage output files
> ${lineage_out}

#get raw lineage information
echo "Retrieving lineage information"
taxonkit lineage -d $'\t' -i 9 ${acc_tax_id} > ${lineage_out}

#contig filtering according to kingdom viruses
#empty files before extraction
> ${contig_matches}
> ${u_match_out}
> ${output_fa}

#extract contig matches that are part of viruses
echo "Extracting contig matches that are part of viruses"
while IFS=$'\t' read -r -a fields; do
    if [[ ${fields[9]} == *Viruses* ]]; then
        echo "${fields[2]}"
        printf "%s\t" "${fields[@]}" "\n" >> ${contig_matches}
        printf "\n" >> ${contig_matches}
    fi
done < ${lineage_out}

#get unique contig matches
echo "Extract unique contig matches"
cat ${contig_matches} | awk '{print $1}' > ${u_match_out}
sort -u ${u_match_out} -o ${u_match_out}

#find the contig matches in the final.contigs.fa file
echo "Producing final blasta fasta with matches to blast against NT and NR"
while read -r hit; do
    if grep -qF ">${hit}" "${fin_fasta}"; then
        awk -v contig=">${hit}" '
            $0 ~ ("^"contig) {print; ON=1; next}
            ON && /^>/ {exit}
            ON {print}
        ' ${fin_fasta} >> ${output_fa}
    fi

    #progress check
    echo "Sequence for ${hit} found"
done < ${u_match_out}

exit 0
```


### 9. BlastN

```
cd metagenomics/scripts
nano 09.blastn.sh
```

```
#!/bin/env bash

#load modules
module ncbi
ON="module miniconda && conda activate viral_pipeline"

#directories
wdir="/analyses/users/nokuzothan/disc_pipe"
db_fa="${wdir}/ncbidb/nt/nt"
input_fa="/analyses/users/nokuzothan/Virolocate/work/c4/03b0a4833978dcf7aabea1ab5a3987/final_blast_contig.fasta"
blastn_out="${wdir}/init_tools/play"
output="${blastn_out}/blastn_output_3.tsv"
blastn_tax_tmp_1="${blastn_out}/contig_acc.txt"
blastn_tax_tmp_2="${blastn_out}/contig_acc_tax.txt"
blastn_tax="${blastn_out}/blastn_taxonomy.tsv"
threads=4

# #make blastn_results subdirectory in blastn output folder
# if [[ -e ${blastn_out} ]]; then
#     rm -rf ${blastn_out}
# fi
# mkdir -p -m a=rwx ${blastn_out}

#blastn run
if [[ -s ${input_fa} ]]; then
    echo "Running blastn"
    blastn -query ${input_fa} \
        -db ${db_fa} \
        -out ${output}\
        -strand both \
        -num_threads ${threads} \
        -evalue 1E-5 \
        -outfmt "6 qseqid qlen sseqid stitle pident length qstart qend evalue bitscore" \
        -perc_identity 80 \
        -max_target_seqs 5

else
    echo "Query fasta file of samples does not exist or is empty, skipping blastn."
fi


#extract contigs ids and accession numbers from blastn_output.tsv

#get taxonomic ids
function get_meta {
    contig=$1
    length=$2
    acc_id=$3
    columns=$4
    output=$5

    #progress check
    echo "Fetching metadata for "${acc_id}""
    #print ncbi page of protein accession and parse taxonomic id for use in taxonkit for lineage
    url1="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc_id}&rettype=gb&retmode=text"
    info=$(curl -N -# ${url1})

    #host source, gographical location name, collection date, gene, product, taxonomic number
    host=$(echo "${info}" | awk -F'"' '/\/host/ {print $2}')
    geo_loc_name=$(echo "${info}" | awk -F'"' '/\/geo_loc_name/ {print $2}')
    date=$(echo "${info}" | awk -F'"' '/\/collection_date/ {print $2}')
    gene=$(echo "${info}" | awk -F'"' '/\/coded_by/ {print $2}')
    product=$(echo "${info}" | awk -F'"' '/\\/product/ {print $2}')
    tax=$(echo "${info}" | awk '/\/db_xref/ { match($0, /taxon:([0-9]+)/, tax_id); print tax_id[1] }')

    #put NA if any of the fields are empty
    if [[ -z "${host}" ]]; then
        host="NA"
    fi

    if [[ -z "${geo_loc_name}" ]]; then
        geo_loc_name="NA"
    fi

    if [[ -z "${date}" ]]; then
        date="NA"
    fi

    if [[ -z "${gene}" ]]; then
        gene="NA"
    fi

    if [[ -z "${product}" ]]; then
        product="NA"
    fi

    if [[ -z "${tax}" ]]; then
        tax="NA"
    fi

    #print output
    echo -e "${contig}\t${length}\t${acc_id}\t${rest}\t${host}\t${gene}\t${product}\t${geo_loc_name}\t${date}\t${tax}" >>${output}
}

while IFS=$'\t' read -r col1 col2 col3 rest;do
    acc=$(echo "[\${col3}]" | cut -d '|' -f2)
    get_meta "${col1}" "${col2}" "${acc}" "${rest}" "blastn_metadata.tsv"
done < "${output}"

#get lineage information
eval ${ON}

#progress check
echo "Getting taxonomic lineage information"

taxonkit lineage -d $'\t' -i 10 ${blastn_tax_tmp_2} >> ${blastn_tax}
conda deactivate

#remove temp file
rm ${blastn_tax_tmp_1}
rm ${blastn_tax_tmp_2}

```

### 10. BlastX

```
cd metagenomics/scripts
nano 10.blastx.sh
```

```
#!/bin/env bash

#load module
module diamond

# variables and directories
wdir="/analyses/users/nokuzothan/disc_pipe"
input_reads_dir="${wdir}/init_tools/diamond/output/blastn.fasta"
db="${wdir}/ncbidb/fasta/nr.faa"
output="${wdir}/init_tools/blastx_nr/output"
threads=$((`/bin/nproc` -2))

#clear existing output directory if any, make new output directory 
if [[ -e $output ]]; then
  rm -rf ${output} 
fi
mkdir -p ${output}

#clear existing output file
out_file="${output}/blastx_output.tsv"
> ${out_file}


#make diamond protein database
diamond makedb --in ${db} -d ${output}/nr

#loop through each of the files created in megahit output directory to find final.congtigs.fa files and run diamond
if [[ -s ${input_reads_dir} ]]; then
  diamond blastx -d ${output}/nr.dmnd \
    -q ${input_reads_dir} \
    --out ${out_file} \
    --threads ${threads} \
    --evalue 1E-5 \
    --outfmt 6 \
    --id 80 \
    --strand both \
    --unal 0 \
    --mp-init 

else 
  echo "Contigs file for samples not found or is empty."
fi

exit 0
```
