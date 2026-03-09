#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=DFBA_1
# Archivo de salida
#SBATCH --output=dfba_1.txt
# Partición 
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=1
#SBATCH --nodelist=n1
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=iptapia@uc.cl
#SBATCH --mail-type=ALL

module load MATLAB/R2021a
cd /home1/iptapia/Calibration
srun matlab -nodisplay -r "Bootstraping_DFBA_CYT_200_1; exit"