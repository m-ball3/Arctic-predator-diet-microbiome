# ------------------------------------------------------------------
# Useful Hyak Commands
# ------------------------------------------------------------------

# log in
ssh mball3@klone.hyak.uw.edu

# SHows storage avaialble in my groups
hyakstorage

# Prints working directory
pwd

# Makes new folder
mkdir

# Rename fils 
mv <old name> <new name>

# Change directory to WADE-003 Arctic Predator Diet Project
cd /gscratch/coenv/mball3/WADE003-arctic-pred
  # folders inside: 
    rawdata  scripts
    # folders inside rawdata
    12sp1 16sp2

# Converts DOS to Unix format
    dos2unix your_script.sh

# Transfer files from local computer to klone with server copy
scp data_to_transfer UWNetID@klone.hyak.uw.edu:/path/to/directory    


scp -r "C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic Predator Diets/12SP1" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/rawdata
scp -r C:/Users/Intern/Desktop/arctic-pred/scripts/ mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/scripts
  # MUST BE RUN FROM LOCAL TERMINAL! NOT HYAK SSH SESSION


scp -r "Z:/Arctic predator diet sequences/16S/*" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/rawdata/16SP2
# ------------------------------------------------------------------
# Uploads reference Db to HYAK
# ------------------------------------------------------------------

# 12S
scp -r "C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/Ref-DB/MURI_MFU_07_2025.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/
  scp -r "C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/Ref-DB/12S-AddSpecies_11-25.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/
  
  # 16S
  scp -r "C:/Users/Intern/arctic-predator-diet-microbiome/DADA2/Ref-DB/16S_Arctic_predator_reference_database_07_2025.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/
  scp -r "C:/Users/Intern/arctic-predator-diet-microbiome/DADA2/Ref-DB/16S-AddSpecies_11-25.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/
  
  scp -r "C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/Ref-DB/16S_Arctic_predator_reference_database_07_2025.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/
  scp -r "C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/Ref-DB/16S-AddSpecies_11-25.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/
  
  scp -r "C:/Users/MBall/OneDrive/Documents/UW-DOCS/WADE LAB/Arctic-predator-diet-microbiome/DADA2/Ref-DB/16S_Arctic_predator_reference_database_07_2025.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/
  scp -r "C:/Users/MBall/OneDrive/Documents/UW-DOCS/WADE LAB/Arctic-predator-diet-microbiome/DADA2/Ref-DB/16S-AddSpecies_11-25.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/
  
# ------------------------------------------------------------------
# Gets output from HYAK
# ------------------------------------------------------------------

#12S
scp -r "mball3@klone.hyak.uw.edu:/mmfs1/home/mball3/WADE003-arcticpred_dada2_QAQC_12SP1_output.Rdata" "C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs"
scp -r "mball3@klone.hyak.uw.edu:/mmfs1/home/mball3/SRKW-diet-16SP1.Rdata" "C:/Users/MBall/OneDrive/文档/WADE LAB/SRKW/DADA2/DADA2 Outputs"


#16S
scp -r "mball3@klone.hyak.uw.edu:/mmfs1/home/mball3/WADE003-arcticpred_dada2_QAQC_16SP1+2-2.Rdata" "C:/Users/MBall/OneDrive/Documents/UW-DOCS/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs"

scp -r "mball3@klone.hyak.uw.edu:/mmfs1/home/mball3/SRKW-diet-16SP1.Rdata" "C:/Users/Intern/SRKW/DADA2/DADA2 Outputs"

scp -r "C:\Users\Intern\Desktop\arctic-pred\16S_Arctic_predator_reference_database_05_2025.fasta" mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/WADE003-arctic-pred/scripts


scp -r * mball3@klone.hyak.uw.edu:/gscratch/coenv/mball3/SRKW/rawdata/16SP1/


# Read error output
cat file

# Prints the .err file
cat /gscratch/coenv/mball3/WADE003-16P2_dada2_QAQC.err

# Pulls the container .sif file
apptainer pull /gscratch/coenv/containerstidyverse_4.0.1.sif 
docker://rocker/tidyverse:4.0.1

# permissions for .sif file
chmod +rx /gscratch/coenv/containerstidyverse_4.0.1.sif

# Makes .R file executable
chmod +x /gscratch/coenv/mball3/WADE003-arctic-pred/scripts/dada2QAQC-arcticpred.R

# FOR SLURM
module load apptainer

apptainer run --bind /gscratch /gscratch/coenv/containerstidyverse_4.0.1.sif Rscript /gscratch/coenv/mball3/WADE003-arctic-pred/scripts/dada2QAQC-arcticpred.R
  #--bind is supposed to make it so that apptainer has access to the files in gscratch


# COPY OUTPUT FILE TO MY LOCAL COMPUTER!!
scp mball3@klone.hyak.uw.edu:/mmfs1/home/mball3/WADE003-arcticpred_dada2_QAQC_12SP1_output-130trunc-ADFGnotes-50BOOTTRUE-leftrim20.Rdata "C:/Users/MBall/OneDrive/Documents/UW-DOCS/WADE lab/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/boot-true"