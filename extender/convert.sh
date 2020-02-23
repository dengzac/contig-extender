bowtie2-build $1 ecolii
bowtie2 -x ecolii -1 $2 -2 $3 -L 5 --n-ceil L,1000,0 --np 0 -p 8 -i S,10,0 --dpad 250 --gbar 20 --rdg 100000,100000 --rfg 1000000,1000000 --no-unal> out3.sam
samtools view -b out3.sam -o out3.bam
samtools sort out3.bam -o out3.bam
samtools index out3.bam
