#对比参考基因组
#!/bin/bash
#$ -S    /bin/sh
#$ -N   bwa_job             # 任务名称
#$   -j     y                       #  任务运行的stdout和stderr信息的文件名相同
#$ -o   ./                          #  stdout文件存放目录
#$ -e   ./                          #  stderr文件存放目录
#$ -cwd                       #在运行任务节点上，先cd 到运行qsub所在的目录再运行
source ~/.bashrc
hash -r
export path=$TMPDIR:$path
vi sampleid #创建一个写明样品id,一行一行写
mkdir -p 01-qc
mkdir -p 02-read-align
mkdir -p 03-mkdup
mkdir -p 04-gatk_vcf
mkdir -p 05-snp_vcf
mkdir -p 06-vcf_snp_flt
path=`pwd`
NUMBER_THREADS=80
REFERENCE=ng #绝对路径(建议把参考基因组放目录下)
rn1=${REFERENCE##*/}
rn=${rg1%%.*}
#01 qc
echo qc
fastqc -t 2 -o 01-qc/ $path/*.fq.gz
multiqc 01-qc/*.zip
echo qc done
#02 bwa and sort sampleid 文件需自备
echo bwa
bwa index $REFERENCE -p $rn
cat sampleid | while read id 
do
    (bwa mem -R "@RG\tID:$id\tSM:$id\tLB:WGS\tPL:Illumina" -t $NUMBER_THREADS $rn ${id}_1_clean.fq.gz ${id}_2_clean.fq.gz || echo -n 'error' ) | samtools sort -@ 5 -o 02-read-align/${id}_sort.bam -
    samtools index 02-read-align/${id}_sort.bam
done
echo bwa done
#03 mkdup
echo mkdup
for i in `ls 02-read-align/*_sort.bam`; do
    sample1=${i##*/}
    sample=${sample1%%_*}
    sambamba markdup -r -t 10 02-read-align/${sample}_sort.bam 03-mkdup/${sample}_sort_mkdup.bam 

done
echo mkdup done
#04 gatk
echo gatk
samtools faidx $REFERENCE
gatk CreateSequenceDictionary -R $REFERENCE
for i in `ls 03-mkdup/*_sort_mkdup.bam`; do
    sample1=${i##*/}
    sample=${sample1%%_*}
    gatk HaplotypeCaller \
 	-R  $REFERENCE\
 	--emit-ref-confidence GVCF \
 	-I 03-mkdup/${sample}_sort_mkdup.bam \
 	-O 04-gatk_vcf/${sample}.gvcf && echo "** gvcf done **"
done
#这个是两个群体的脚本，如果以后分析多个群体，请用下面的
id1=`cat sampleid | sed -n '1p'`
id2=`cat sampleid | sed -n '2p'`
gatk CombineGVCFs -R $REFERENCE --variant 04-gatk_vcf/$id1.gvcf --variant 04-gatk_vcf/$id2.gvcf -O 04-gatk_vcf/allsample.gvcf
gatk GenotypeGVCFs -R $REFERENCE -V 04-gatk_vcf/allsample.gvcf  -O 04-gatk_vcf/allsample.vcf
echo gatk done
#echo 'gatk CombineGVCFs \' >combine,sh
#echo '-R '$REFERENCE' \'  >> combine,sh
#for i in `ls 04-gatk_vcf/*.gvcf`;
#do
#echo '--variant '${i}' \' >>combine.sh;
#done
#echo '-O  04-gatk_vcf/allsample.vcf' >>combine.sh;
#done
#05 selecting variants
echo selecting variants
gatk SelectVariants -R $REFERENCE -V 04-gatk_vcf/allsample.vcf -select-type  SNP -O 05-snp_vcf/allsample_snp.vcf
echo selecting variants done
#06 filtering variants
echo filtering variants
gatk -T VariantFiltration -R $REFERENCE -V 05-snp_vcf/allsample_snp.vcf --filterExpression "QD<2.0||FS>60.0||MQ<40.0||MQRankSum<-12.5||ReadPosRankSum<-8.0" --filterName "LowQual" -o 06-vcf_snp_flt/allsample_snp_flt.vcf
echo filtering variants done
echo finish
