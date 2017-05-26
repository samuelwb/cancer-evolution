# Author: Ryan Layer

SAMPLES="B0
B1
B2
B3
B4"

PE_SR_OPS=""

for SAMPLE in $SAMPLES
do
  DISC_BAM=$SAMPLE/$SAMPLE/$SAMPLE.discordants.bam
  SPLIT_BAM=$SAMPLE/$SAMPLE/$SAMPLE.splitters.bam
  HISTO=$SAMPLE/$SAMPLE/$SAMPLE.psort.bam.histo
  MEAN=`cat $SAMPLE/$SAMPLE/$SAMPLE.psort.bam.stats | tail -n 1 | cut -f1 | cut -d":" -f2`
  STD=`cat $SAMPLE/$SAMPLE/$SAMPLE.psort.bam.stats | tail -n 1 | cut -f2 | cut -d":" -f2`
  RL=`cat $SAMPLE/$SAMPLE/$SAMPLE.psort.bam.stats | head -n 1`
  Z=4
  PE_BD=20
  SR_BD=10
  PE_PARAMS="-pe id:$SAMPLE,bam_file:$DISC_BAM,histo_file:$HISTO,mean:$MEAN,stdev:$STD,read_length:$RL,min_non_overlap:$RL,discordant_z:$Z,back_distance:$PE_BD,weight:1,min_mapping_threshold:20"
  SR_PARAMS="-sr id:$SAMPLE,bam_file:$SPLIT_BAM,back_distance:$SR_BD,weight:1,min_mapping_threshold:20"
  PE_SR_OPS="$PE_SR_OPS $PE_PARAMS $SR_PARAMS"
done

lumpy -mw 4 -tt 0 -x btu356_LCR-hs37d5.MT.hs37d5.phix.bed $PE_SR_OPS > B0.B1.B2.B3.B4.x.btu356_LCR-hs37d5.MT.hs37d5.bed.vcf

for SAMPLE in $SAMPLES
do
  BAM=$SAMPLE/$SAMPLE/$SAMPLE.psort.bam
  SPLIT_BAM=$SAMPLE/$SAMPLE/$SAMPLE.splitters.bam
  svtyper -B $BAM -S $SPLIT_BAM -i B0.B1.B2.B3.B4.x.btu356_LCR-hs37d5.MT.hs37d5.phix.bed.vcf -o B0.B1.B2.B3.B4.x.btu356_LCR-hs37d5.MT.hs37d5.phix.bed.vcf.$SAMPLE
done

#vcf_paste is distributed with svtyper
vcf_paste.py \
  -m B0.B1.B2.B3.B4.x.btu356_LCR-hs37d5.MT.hs37d5.phix.bed.vcf.B0 \
  B0.B1.B2.B3.B4.x.btu356_LCR-hs37d5.MT.hs37d5.phix.bed.vcf.B0 \
  B0.B1.B2.B3.B4.x.btu356_LCR-hs37d5.MT.hs37d5.phix.bed.vcf.B1 \
  B0.B1.B2.B3.B4.x.btu356_LCR-hs37d5.MT.hs37d5.phix.bed.vcf.B2 \
  B0.B1.B2.B3.B4.x.btu356_LCR-hs37d5.MT.hs37d5.phix.bed.vcf.B3 \
  B0.B1.B2.B3.B4.x.btu356_LCR-hs37d5.MT.hs37d5.phix.bed.vcf.B4  \
| cut -f1-10,16,22,28,34 \
> B0.B1.B2.B3.B4.x.btu356_LCR-hs37d5.MT.hs37d5.phix.bed.svt.B0.B1.B2.B3.B4.vcf
