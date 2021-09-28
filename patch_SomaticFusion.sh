
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

version=0.0.4

# source variables file
. *.variables

. $panel.variables
. SomaticFusion.config

# set conda env
source "$conda_bin_path"/activate SomaticFusion


###########################
# Generate Metrics Report #
###########################

numberSamplescomplete=$(cat /data/results/$seqId/$panel/Summary_Metrics/samples_list.txt | uniq | wc -l)
numberSamplesInProject=$(find ../ -maxdepth 2 -mindepth 2 | grep .variables | grep -v "$panel".variables | uniq | wc -l)

echo $numberSamplescomplete
echo $numberSamplesInProject 
find ../ -maxdepth 2 -mindepth 2 | grep .variables | grep -v "$panel".variables > proj_samp.txt
cat /data/results/$seqId/$panel/Summary_Metrics/samples_list.txt | uniq > smpcomp.txt




if [ $numberSamplescomplete -eq $numberSamplesInProject ]
then

    python /data/results/$seqId/$panel/Summary_Metrics/total_reads_list.py "$seqId" "/data/results/$seqId/$panel/Summary_Metrics"

else
    echo "not all samples have completed running. Finishing process for this sample."
fi



source "$conda_bin_path"/deactivate