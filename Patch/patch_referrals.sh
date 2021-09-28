#!/bin/bash
# set -euo pipefail
#PBS -l walltime=40:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

# Description: Get FASTQ metrics for validation run .
# Authors: Christopher Medway and Laura McCluskey, Seemu Ali
# Date: 27th March 2020
# Usage: qsub run_star-fusion.sh [inside sample dir with .variables and \
# .fastq.gz files]

version=0.0.2

# source variables file
. *.variables
. SomaticFusion.config



source "$conda_bin_path"/activate SomaticFusion

python fusion_report_referrals.py /data/results/"$seqId"/RocheSTFusion/"$sampleId"/ $sampleId

source "$conda_bin_path"/deactivate
