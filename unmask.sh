sample="$(basename -- $(pwd))"
. SomaticFusion.config
. $sample.variables
. IlluminaSTFusion.variables

python /data/results/$seqId/IlluminaSTFusion/$sample/fusion_report_cases_unmasked.py /data/results/$seqId/IlluminaSTFusion/$sample /data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-0.0.4/IlluminaSTFusion $sampleId
