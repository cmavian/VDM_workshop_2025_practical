#!/bin/env bash


#####################################################################
#                                                                   #
# Automated script for Viral discovery workflow databseses download # 
#                                                                   #
#                                                                   #
# CERI 17.11.2025                                                   #
#                                            created by TJ Sanko    #
#                                                                   #
# Usage:                                                            #
#  Downlaod the script. Copy it to location where you want to       #
#  download the databases. Then run commands:                       #
#                                                                   #
#     chmod +x db_download.sh                                       #
#     bash db_download.sh                                           #
#                                                                   #
# it will take long time to download so to avoid disconnections     #
# put the downlaod in screen mode:                                  #
#     screen -S download                                            #
#                                                                   #
# this will open a new screen. type to start downlaod               #
#      bash db_download.sh &                                        #
#                                                                   #
# and then press ctrl+A and D to disconnect the screen and run      #
# the process in background.                                        #
#                                                                   #
# +--------------------------------------------------------------+  #
# !!! For all databases you need ~2.5 TB of space on your server    #
# +--------------------------------------------------------------+  #
#                                                                   #
#####################################################################


### Downloading NCBI databases:
########################
### 1. NR db (~1.05 TB)
########################
DB='nr'

for X in `seq -w 0 128`; do 
  ### downlading
  wget ftp.ncbi.nlm.nih.gov/blast/db/nr.${X}.tar.gz
  wget ftp.ncbi.nlm.nih.gov/blast/db/nr.${X}.tar.gz.md5
  
  ### checking for competness of download
  md5sum -c ftp.ncbi.nlm.nih.gov/blast/db/nr.${X}.tar.gz.md5
  
  ### decompressing db
  gzip -dkf nr.${X}.tar.gz
done

if [[ $? -eq 0 ]]; then echo "[INFO] Finished downloading $DB db succesfully";
else echo "[WARNING] The download of $DB db was not competely successful"; fi
exit 0;
########################
### 2. NT db (~ 1.3 TB)
########################
DB='nt'

for X in `seq -w 0 245`; do
  ### downlading
  wget ftp.ncbi.nlm.nih.gov/blast/db/nt.${X}.tar.gz
  wget ftp.ncbi.nlm.nih.gov/blast/db/nt.${X}.tar.gz.md5
  
  ### checking for competness of download
  md5sum -c ftp.ncbi.nlm.nih.gov/blast/db/nt.${X}.tar.gz.md5
  
  ### decompressing db
  gzip -dkf nr.${X}.tar.gz
done

if [[ $? -eq 0 ]]; then echo "[INFO] Finished downloading $DB db succesfully";
else echo "[WARNING] The download of $DB db was not competely successful"; fi
exit 0;

### Downloading RVDB proteing db:
#############
### 3. RVDB (~ 
#############
DB='RVDB'

### downloading the protein db 
wget https://rvdb-prot.pasteur.fr/files/U-RVDBv30.0-prot.fasta.xz

### decompressing db
xz -dk -T 5 U-RVDBv30.0-prot.fasta.xz

if [[ $? -eq 0 ]]; then echo "[INFO] Finished downloading $DB db succesfully";
else echo "[WARNING] The download of $DB db was not competely successful"; fi
exit 0;



