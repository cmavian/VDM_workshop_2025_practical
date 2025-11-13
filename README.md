<p><b><h1>VDM_workshop_2025_practical</h1></b></p>
<h3>Virus discovery pipeline</h3>
<i>created by</i> Erin Harvey, Carla Mavian, Nokuzotha Nkiwane, TJ Sanko, Eduan Wilkinson.
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

<li>Absolute vs relative path</li>
<pre><code>cd /home/vdw00</code></pre>
<pre><code>cd ../</code></pre>
<pre><code>cd ./vdw00</code></pre>

<li>First step will be creating a working folder in workspace folder and moving into the freshly created directory</li>
<pre><code>cd ~/workspace
mkdir -p metagenomics
cd metagenomics</code></pre>
To check if the location was changed, print working directory
<pre><code>pwd</code></pre>

<li>Next, we will create other folders which will, help keep the analysis sorted</li> 
<pre><code>mkdir -p data results scripts</code></pre>
Then change permissions on all of the directories and in them:
<pre><code>chmod -R a=rwx *</code></pre>
To check if everything was created correctly, list the current directory.
<pre><code>ll *</code></pre>

<li>Let's go into data directory and copy our files</li><br>
<pre><code>cd data</code></pre>
Downloading files from website:
<pre><code>wget <a href="https://">https://</a></code></pre>
Alternatively, if the internet connection is slow, you can copy the data from backup location:
<pre><code> cp -R /analyses/vdworkshop/.backup/data/* ./ </code></pre>
Then change permission on the files and confirm the change on the data in the location:
<pre><code>chmod a=rwx *
ll *</code></pre>
</ol>

<b><h3>Analysis</h3></b>
In this section we will create scripts to execute each step. The scripts will be located in metagenomics/scripts directory.
<pre><code>cd ~/workspace/metagenomics</code></code></pre>
(<i>if you get lost, you can use the absolute path:</i> <code>/analyses/vdworkshop/${USER}/metagenomics</code>)<br><br>
We will use <i>nano</i> text editor to creat all scripts. You can create empty file and open it to edit at the same time
<pre><code>nano scripts/[script_name].sh</code></pre>
<i>you can copy & paste the script text directly into the open document.</i><br>
To save any of your scripts, press ctrl+X then Y and ENTER. This will "override" the [script_name].sh (which was empty on opening)<br>
You can change the name after pressing 'Y'.Before running the script you have to change permissions
<pre><code>chmod a=rwx scripts/[script_name].sh</code></pre>
To run the script type:
<pre><code>bash scripts/[script_name].sh</code></pre><br>

<ol start=1>
<h4><li>FastQC pre-Trimmomatic</li></h4>
<pre><code>nano scripts/01.fastqc_pretrim.sh</code></pre>
<pre><code>#!/bin/env bash

THR=5
###enabling conda environment and fastqc program
ON="module miniconda && conda activate fastqc"
eval $ON

###input and output directories
workdir=`realpath $(pwd) 2>/dev/null`
input=${workdir}'/data/fastq'
output=${workdir}'/results/01.fastqc_pretrim'

###make output directory if it doesn't exist
if ! [[ -d ${output} ]]; then mkdir -p -m a=rwx ${output}; fi

echo '[INFO] Running FastQC pre-trim'

###run fastqc
fastqc -t ${THR} -o ${output}/ ${input}/*.f*q*
echo '[DONE] Pretrim FastQC complete"

###deactivating fastqc program
OFF='conda deactivate'
eval ${OFF}

chmod -R a=rwx ${output}
exit 0; </code></pre>
<pre><code>chmod a=rwx scripts/01.fastqc_pretrim.sh</code></pre>
<pre><code>bash scripts/01.fastqc_pretrim.sh</code></pre>

<h4><li>Trimmomatic</li></h4>
<pre><code>nano scripts/02.trimmomatic.sh</code></pre><br>
<pre><code>#!/bin/env bash

THR=5
###activating program
ON='module trimmomatic 1>/dev/null 2>/dev/null'
eval ${ON}

workdir=`realpath $(pwd) 2>/dev/null`
input=${workdir}'/data/fastq'
output=${workdir}'/results/02.trim_output'
ADAPTERS=${workdir}'/data/adapters.fa'

###make output directory if it doesn't exist
if [[ ! -e ${output} ]]; then mkdir -p ${output}; fi

###getting the SRA IDs and creating output files
for FOW in $(ls ${input}/*.f*q* | grep -Ei "_r?1"); do
  REV=`echo ${FOW} | sed -r 's/\_(r|R)?1/\_\12/'`;
  ID=`basename ${FOW} | cut -d '_' -f1`
  P1=${output}/${ID}'_1.P.fq.gz'
  U1=${output}/${ID}'_1.U.fq.gz'
  P2=${output}/${ID}'_2.P.fq.gz'
  U2=${output}/${ID}'_2.U.fq.gz'

  ###running program
  echo "[INFO] trimming: ${ID}"
  trimmomatic PE -threads ${THR} -phred33 -summary ${output}/${ID}'_statsSummary.txt' \
    ${FOW} ${REV} ${P1} ${U1} ${P2} ${U2} \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25 \
    2>>${output}/${ID}.log 1>>${output}/${ID}.log

  echo "[DONE] Trimming for ${ID} complete"
done
chmod -R a=rwx ${output}
exit 0;</code></pre>
<pre><code>chmod a=rwx  scripts/02.trimmomatic.sh</code></pre>
<pre><code>bash scripts/02.trimmomatic.sh</code></pre>

<h4><li>FastQC post-Trimmomatic</li></h4>
<pre><code>nano scripts/03.fastqc_posttrim.sh</code></pre>
<pre><code>#!/bin/env bash

THR=5
###enabling conda environment and fastqc program
ON="module miniconda && conda activate fastqc 1>/dev/null 2>/dev/null"
eval $ON

###input and output directories
workdir=`realpath $(pwd) 2>/dev/null`
input=${workdir}'/results/02.trim_output'
output=${workdir}'/results/03.fastqc_posttrim'

###make output directory if it doesn't exist
if ! [[ -d ${output} ]]; then mkdir -p -m a=rwx ${output}; fi

echo '[INFO] Running FastQC post-trim'

###run fastqc
fastqc -t ${THR} -o ${output}/ ${input}/*.P.fq.gz
echo '[DONE] Post-trim FastQC complete'

###deactivating fastqc program
OFF='conda deactivate'
eval ${OFF}

chmod -R a=rwx ${output}
exit 0; </code></pre>
<pre><code>chmod a=rwx  scripts/03.fastqc_posttrim.sh</code></pre>
<pre><code>bash  scripts/03.fastqc_posttrim.sh</code></pre>

<h4><li>MultiQC</li></h4>
<pre><code>nano  scripts/04.multiqc.sh</code></pre>
<pre><code>#!/bin/env bash

THR=5
###enabling conda environment and fastqc program
ON="module miniconda && conda activate fastqc 1>/dev/null 2>/dev/null"
eval $ON

###input and output directories
workdir=`realpath $(pwd) 2>/dev/null`
input_pre=${workdir}'/results/01.fastqc_pretrim'
input_post=${workdir}'/results/03.fastqc_posttrim'
output=${workdir}'/results/04.multiqc'

###make output directory if it doesn't exist
if ! [[ -d ${output} ]]; then mkdir -p -m a=rwx ${output}; fi

echo '[INFO] Runnig multiqc'
###to join both results
#multiqc ${input_pre} ${input_post} -o ${output}

multiqc ${input_pre} -o ${output}/multiqc_pretrim
multiqc ${input_post} -o ${output}/multiqc_posttrim

###deactivating multiqc program
OFF='conda deactivate'
eval ${OFF}

zip -9r ${output}/multiqc_pretrim.zip ${output}/multiqc_pretrim
zip -9r ${output}/multiqc_posttrim.zip ${output}/multiqc_posttrim

chmod -R a=rwx ${output}
exit 0; </code></pre>
<pre><code>chmod a=rwx  scripts/04.multiqc.sh</code></pre>
<pre><code>bash  scripts/04.multiqc.sh</code></pre>

<h4><li>Megahit</li></h4>
<pre><code>nano  scripts/05.megahit.sh</code></pre>
<pre><code>#!/bin/env bash

THR=5
###enabling conda environment and megahit program
ON="module miniconda && conda activate megahit 1>/dev/null 2>/dev/null"
eval $ON

###input and output directories
workdir=`realpath $(pwd) 2>/dev/null`
input=${workdir}'/results/02.trim_output'
output=${workdir}'/results/05.megahit'

###make output directory if it doesn't exist
if ! [[ -d ${output} ]]; then mkdir -p -m a=rwx ${output}; fi

###run megahit
for FOW in $(ls ${input}/*_1.P.fq.gz); do
  REV=`echo ${FOW} | sed -r 's/\_(r|R)?1/\_\12/'`;
  ID=`basename ${FOW} | cut -d '_' -f1`
  LOG=${output}/${ID}.log

  echo "[INFO] assembling: ${ID}"
  megahit -t ${THR} -1 ${FOW} -2 ${REV} -o ${output}/${ID} 1>${LOG} 2>${LOG}

  ###after megahit run, prepend sample name to contigs
  ###(allows easier tracking of which contigs came from which sample in downstream analysis)
  awk -v prefix="${ID}_" '/^>/ {$0=">" prefix substr($0,2)} {print}
    ' ${output}/${ID}/final.contigs.fa > ${output}/${ID}/${ID}.contigs.fasta

  echo "[DONE] Megahit assmebly completed successfully and contigs named with sample name: ${ID}"
done

###deactivating megahit program
OFF='conda deactivate'
eval ${OFF}

chmod -R a=rwx ${output}
exit 0; </code></pre>
<pre><code>chmod a=rwx  scripts/05.megahit.sh</code></pre>
<pre><code>bash  scripts/05.megahit.sh</code></pre>

<h4><li>Diamond</li></h4>
<ol start=i>
<li>Filtering viral sequences from RVDB protein database</li>
<pre><code>nano scripts/06.1.diamond_rvdb.sh</code></pre>

<pre><code>#!/bin/env bash
THR=5
### enabling conda environment and diamond program
ON="module diamond 1>/dev/null 2>/dev/null"
eval $ON

### input and output directories
workdir=`realpath $(pwd) 2>/dev/null`
input=${workdir}'/results/05.megahit'
output=${workdir}'/results/06.1.diamond_rvdb'
input_db=${workdir}'/data/database'

### DB variables and directories
FASTA="${input_db}/U-RVDBv30.0-prot.fasta"
DB='RVDB'
if [[ ! -e ${FASTA} && -e ${FASTA}.xz ]]; then xz -dk -T ${THR} ${FASTA}.xz; fi

### make output directory if it doesn't exist
if ! [[ -d ${output} ]]; then mkdir -p -m a=rwx ${output}; fi

### make diamond protein database
echo '[INFO] Reading database file '{$DB}
if [[ ! -e ${output}/${DB}.dmnd ]]; then
  diamond makedb --in ${FASTA} --threads ${THR} -d ${output}/${DB}; fi
chmod a=rwx ${output}/${DB}.dmnd

### loop through each of the files created in megahit output directory
### to find final.contigs.fa files and run diamond
function DBLX {
 FOLDER=$1; OUT=$2; DB=$3; T=$4;
  ID=$(basename $1)
  CONTIGS=${FOLDER}/${ID}.contigs.fasta
  LOG=${OUT}/${ID}.log

  if [[ -f ${CONTIGS} ]]; then
	echo "[INFO] Analyzing sample ${ID}"
    sample_out=${OUT}/${ID}_rvdb.blastx
    diamond blastx -d ${OUT}/${DB}.dmnd \
    -q ${CONTIGS} \
    --out ${sample_out} \
    --threads ${T} \
    --evalue 1e-5 \
    --outfmt 6 qseqid qlen sseqid stitle pident length evalue bitscore \
    --id 80 \
    --strand both \
    --unal 0 \
    --mp-init 1>${LOG} 2>${LOG}
  else
    echo "[WARNING] Contigs file for ${ID} not found."
  fi
}
export -f DBLX

### check and run in parallel
PARALLEL_VER=(`parallel --version | grep 'GNU parallel'`)
if [[ ${PARALLEL_VER[2]} =~ [0-9]+ ]]; then
 ls -dl ${input}/* | grep ^d | awk '{print $9}' | parallel -j ${THR} -n1 -I% "DBLX %" ${output} ${DB} 1
else
 for FOLDER in $(ls -dl ${input}/* | grep ^d | awk '{print $9}'); do
  DBLX ${FOLDER} ${output} ${DB} 1
 done
fi

chmod -R a=rwx ${output}
exit 0; </code></pre>
<pre><code>chmod a=rwx scripts/06.1.diamond_rvdb.sh</code></pre>
<pre><code>bash scripts/06.1.diamond_rvdb.sh</code></pre>

<li>Filtering viral sequences from NCBI protein database</li><br>
<pre><code>nano 06.2.diamond_nrdb.sh</code></pre>

<pre><code>#!/bin/env bash

THR=5
### enabling conda environment and diamond program
ON="module diamond 1>/dev/null 2>/dev/null"
eval $ON

### input and output directories
workdir=`realpath $(pwd) 2>/dev/null`
input=${workdir}'/results/05.megahit'
output=${workdir}'/results/06.1.diamond_nrdb'
input_db=${workdir}'/data/database'

### DB variables and directories
FASTA="${input_db}/nr.faa"
DB='NRDB'
if [[ ! -e ${FASTA} && -e ${FASTA}.gz ]]; then gzip -dkf ${FASTA}.gzip; fi

### make output directory if it doesn't exist
if ! [[ -d ${output} ]]; then mkdir -p -m a=rwx ${output}; fi

### make diamond protein database
echo '[INFO] Reading database file '{$DB}
if [[ ! -e ${output}/${DB}.dmnd ]]; then
  diamond makedb --in ${FASTA} --threads ${THR} -d ${output}/${DB}; fi
chmod a=rwx ${output}/${DB}.dmnd

### loop through each of the files created in megahit output directory
### to find final.contigs.fa files and run diamond
function DBLX {
 FOLDER=$1; OUT=$2; DB=$3; T=$4;
  ID=$(basename $1)
  CONTIGS=${FOLDER}/${ID}.contigs.fasta
  LOG=${OUT}/${ID}.log

  if [[ -f ${CONTIGS} ]]; then
	echo "[INFO] Analyzing sample ${ID}"
    sample_out=${OUT}/${ID}_nrdb.blastx
    diamond blastx -d ${OUT}/${DB}.dmnd \
    -q ${CONTIGS} \
    --out ${sample_out} \
    --threads ${T} \
    --evalue 1e-5 \
    --outfmt 6 qseqid qlen sseqid stitle pident length evalue bitscore \
    --id 80 \
    --strand both \
    --unal 0 \
    --mp-init 1>${LOG} 2>${LOG}
  else
    echo "[WARNING] Contigs file for ${ID} not found."
  fi
}
export -f DBLX

### check and run in parallel
PARALLEL_VER=(`parallel --version | grep 'GNU parallel'`)
if [[ ${PARALLEL_VER[2]} =~ [0-9]+ ]]; then
 ls -dl ${input}/* | grep ^d | awk '{print $9}' | parallel -j ${THR} -n1 -I% "DBLX %" ${output} ${DB} 1
else
 for FOLDER in $(ls -dl ${input}/* | grep ^d | awk '{print $9}'); do
  DBLX ${FOLDER} ${output} ${DB} 1
 done
fi

chmod -R a=rwx ${output}
exit 0; </code></pre>
<pre><code>chmod a=rwx scripts/06.2.diamond_nrdb</code></pre>
<pre><code>bash scripts/06.2.diamond_nrdb</code></pre>
</ol>

<li>Combining results from both blastx seqrches</li><br>
<pre><code>nano scripts/07.combine.sh</code></pre>

<pre><code>#!/bin/env bash

THR=5
### input and output directories
workdir=`realpath $(pwd) 2>/dev/null`
input=${workdir}'/results/05.megahit'
output=${workdir}'/results/07.combine'
rvdb_in=${workdir}'/results/06.1.diamond_rvdb'
nrdb_in=${workdir}'/results/06.2.diamond_nrdb'
contigs_in=${workdir}'/results/05.megahit'

combined_dir="${workdir}/results/07.combined"
### make output directory if it doesn't exist
if ! [[ -d ${combined_dir} ]]; then mkdir -p -m a=rwx ${combined_dir}; fi

rvdb_out=${combined_dir}'/combined_rvdb.blastx'
nrdb_out=${combined_dir}'/combined_nr.blastx'
rvdb_virus=${combined_dir}'/combined_rvdb.virus_hits.blastx'
nrdb_virus=${combined_dir}'/combined_nrdb.virus_hits.blastx'
unique=${combined_dir}'/unique_viral_contigs.txt'
combined_contigs=${combined_dir}'/all_contigs.fasta'
viral_fasta=${combined_dir}'/viral_contigs.fasta'

echo "[INFO] Concatenating all contigs..."
cat ${contigs_in}/*.contigs.fasta > ${combined_contigs}

echo "[INFO] Concatenating RVDB blastx results..."
cat ${rvdb_in}/*_rvdb.blastx > ${rvdb_out}

echo "[INFO] Concatenating VP blastx results..."
cat ${nrdb_in}/*_nrdb.blastx > ${nrdb_out}

echo "[INFO] Filtering for viral hits..."

grep -Ei "vir[us|idae|oid]" ${rvdb_out} > ${rvdb_virus}
grep -Ei "vir[us|idae|oid]" ${nrdb_out} > ${nrdb_virus}

#grep -Ei "virus|viridae|viroid" ${rvdb_out} > ${rvdb_virus}
#grep -Ei "virus|viridae|viroid" ${nrdb_out} > ${nrdb_virus}

echo "[INFO] Extracting unique contig names..."
awk '{print $1}' ${rvdb_virus} ${nrdb_virus} | sort -u > ${unique}

# Activate seqkit environment
ON="module miniconda && conda activate seqkit 1>/dev/null 2>/dev/null"
eval ${ON}

echo "[INFO] Extracting viral contigs using seqkit..."
seqkit grep -f ${unique} ${combined_contigs} > ${viral_fasta}

OFF='conda deactivate'
eval ${OFF}

echo -en "
[DONE]
 Combined RVDB blastx:      ${rvdb_out}
 Combined VP blastx:        ${vp_out}
 RVDB viral hits:           ${rvdb_virus}
 VP viral hits:             ${vp_virus}
 Unique viral contigs list: ${unique}
 Viral contigs FASTA:       ${viral_fasta}
"

echo "[COUNT] Number of viral contigs extracted:"
grep -Ec "^>" ${viral_fasta}

exit 0; </code></pre>
<pre><code>chmod a=rwx scripts/07.combine.sh</code></pre>
<pre><code>bash scripts/07.combine.sh</code></pre>

<li>Filtering viral sequences from NCBI nucleotide database</li><br>
<pre><code>nano scripts/08.ncbi_ntdb.sh</code></pre>
<pre><code>#!/bin/env bash
THR=5
###### enabling BLAST program
ON='module ncbi'
eval ${ON}

###input and output directories
workdir=`realpath $(pwd) 2>/dev/null`
input=${workdir}'/results/07.combined/viral_contigs.fasta'
output=${workdir}'/results/07.combine/08.ncbi_ntdb'
input_db=${workdir}'/data/database'

###DB variables and directories
FASTA="${input_db}/nt.fasta"
DB='NTDB'
if [[ ! -e ${FASTA} && -e ${FASTA}.gz ]]; then gzip -dkf ${FASTA}.gzip; fi

###make output directory if it doesn't exist
if ! [[ -d ${output} ]]; then mkdir -p -m a=rwx ${output}; fi

echo -en "
Running BLASTN against nt database...
Query: ${input}
DB:    ${DB}
"

blastn \
    -db ${DB} \
    -query ${input} \
    -out ${output}'/viral_contigs_vs_nt.blastn' \
    -num_threads ${THR} \
    -evalue 1e-10 \
    -max_target_seqs 25 \
    -outfmt "6 qseqid qlen sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    > ${output}/blastn.log 2>&1

chmod -R a+rwx ${output}

echo -en "
BLASTN complete.
Results written to: ${output}/viral_contigs_vs_nt.blastn
Log:                ${output}/blastn.log
"
exit 0; </code></pre>
<pre><code>chmod a=rwx scripts/08.ncbi_ntdb.sh</code></pre>
<pre><code>bash scripts/08.ncbi_ntdb.sh</code></pre>




