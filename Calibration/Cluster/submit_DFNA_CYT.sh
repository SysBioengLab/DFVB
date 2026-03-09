#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=Calib_DFNA_5
# Archivo de salida
#SBATCH --output=salida_dfna_cyt_5.txt
# Partición 
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=1
#SBATCH --exclude=n14,n4
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=iptapia@uc.cl
#SBATCH --mail-type=ALL

module load MATLAB/R2021a
cd /home1/iptapia/Calibration
srun matlab -nodisplay -r "Bootstraping_DFNA_CYT_200_5; exit"