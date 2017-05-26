# Author: Ryan Layer

BAM=
OUT_DIR=

PYTHON=/usr/local/bin/python
SAMBAMBA=$HOME/bin/sambamba
LUMPY_HOME=$HOME/src/lumpy-sv
PAIREND_DISTRO=$LUMPY_HOME/scripts/pairend_distro.py
SAMBLASTER=$HOME/bin/samblaster
SAMTOBAM="$SAMBAMBA view -S -f bam -l 0"
SAMSORT="$SAMBAMBA sort -m 1G --tmpdir "
BAMGROUPREADS=$LUMPY_HOME/scripts/bamkit/bamgroupreads.py
BAMFILTERRG=$LUMPY_HOME/scripts/bamkit/bamfilterrg.py
MAX_SPLIT_COUNT=2
MIN_NON_OVERLAP=20

TEMP_DIR=$OUT_DIR
mkdir -p $TEMP_DIR/spl $TEMP_DIR/disc
mkfifo $TEMP_DIR/spl_pipe
mkfifo $TEMP_DIR/disc_pipe

READ_LENGTH=`$SAMBAMBA view $BAM | head -n 10000 | awk 'BEGIN { MAX_LEN=0 } { LEN=length($10); if (LEN>MAX_LEN) MAX_LEN=LEN } END { print MAX_LEN }'`
$PYTHON $BAMGROUPREADS \
    -i $BAM \
| $SAMBLASTER \
    --acceptDupMarks \
    --excludeDups \
    --addMateTags \
    --maxSplitCount $MAX_SPLIT_COUNT \
    --minNonOverlap $MIN_NON_OVERLAP \
    --splitterFile $TEMP_DIR/spl_pipe \
    --discordantFile $TEMP_DIR/disc_pipe \
| $SAMTOBAM /dev/stdin \
| $SAMSORT $TEMP_DIR/spl \
    -o $TEMP_DIR/$BAM_BASE.psort.bam  \
    /dev/stdin &
$SAMTOBAM $TEMP_DIR/spl_pipe \
| $SAMSORT $TEMP_DIR/spl \
    -o $TEMP_DIR/$BAM_BASE.splitters.bam  \
    /dev/stdin &
$SAMTOBAM $TEMP_DIR/disc_pipe \
| $SAMSORT $TEMP_DIR/disc \
    -o $TEMP_DIR/$BAM_BASE.discordants.bam \
    /dev/stdin
wait
rm $TEMP_DIR/spl_pipe
rm $TEMP_DIR/disc_pipe
rmdir $TEMP_DIR/spl
rmdir $TEMP_DIR/disc

echo $READ_LENGTH > $BAM.stats
$SAMBAMBA view $BAM \
| awk '{if (NR<=1000000)
            print > "/dev/stdout";
        else
            print > "/dev/null" }' \
| $PYTHON $PAIREND_DISTRO \
    -r $READ_LENGTH \
    -X 4 \
    -N 100000 \
    -o $BAM.histo \
>> $BAM.stats
