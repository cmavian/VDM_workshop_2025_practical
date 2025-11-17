#!/bin/env bash


######################################################################
#                                                                    #
# Automated script for Viral discovery workflow samples download     # 
# from SRA on NCBI                                                   #
#                                                                    #
# CERI 17.11.2025                                                    #
#                                            created by TJ Sanko     #
#                                                                    #
# Usage:                                                             #
#  Downlaod the script. Copy it to location where you want to        #
#  download the databases. Then run commands:                        #
#                                                                    #
#     chmod +x sra_download.sh                                       #
#     bash sra_download.sh SRA_numbers.txt                           #
#                                                                    #
# it will take long time to download so to avoid disconnections      #
# put the downlaod in screen mode:                                   #
#     screen -S download                                             #
#                                                                    #
# this will open a new screen. type to start downlaod                #
#      bash sra_download.sh SRA_numbers.txt &                        #
#                                                                    #
# and then press ctrl+A and D to disconnect the screen and run       #
# the process in background.                                         #
#                                                                    #
#                                                                    #
# Requires sra-toolkit to be installed                               #
# https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit  #
#                                                                    #
# installing under conda                                             #
#    conda create -y -n sratoolkit -c bioconda sra-tools             #
#                                                                    #
#                                                                    #
######################################################################


THR=5
### Importing SRA numbers in list format
SRA_LST=$1;

if [[ ! -e ${SRA_LST} || ${SRA_LST} =~ ^$ ]]; then echo "[ERROR] You need to provide one column list of SRA numbers to download" && exit 2; fi;
if [[   -e ${SRA_LST} && `cat ${SRA_LST} | grep -vE "^\s?$" | wc -l` -eq 0 ]]; then echo "[ERROR] The SRA list: ${SRA_LST} is empty" && exit 2; fi;

### checking for software
PATTERN="[0-9]+.[0-9]+(.[0-9]+)?"
if ! [[ `fastq-dump -V | cut -d ':' -f2 | sed -r 's,\s+,,g'` =~ ${PATTERN} ]]; then echo "[ERROR] SRAtoolkit was not detected" && exit 2; fi;

#### downloading data
PARALLEL_VER=(`parallel --version | grep 'GNU parallel'`)
if [[ ${PARALLEL_VER[2]} =~ [0-9]+ ]]; then

 echo '[INFO] Starting parallel mode' && sleep 1s
 cat ${SRA_LST} | parallel -j ${THR} -n1 -I% "fastq-dump --split-files %"
 cat ${SRA_LST} | parallel -j ${THR} -n1 -I% "gzip -9 %*.fastq"

else

 echo '[INFO] Starting non-parallel mode' && sleep 1s
 for SRA in $(cat ${SRA_LST}); do fastq-dump --split-files ${SRA}; gzip -9 ${SRA}*.fastq; done
 
fi

echo 'All done'
exit 0;

