wget https://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/analysis/IonTorrent_TVC_06302015/TSVC_variants_defaultlowsetting.vcf
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/ion_exome/IonXpress_020_rawlib.hg19.bam
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/ion_exome/IonXpress_020_rawlib.hg19.bam.bai

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
samtools faidx hg19.fa
picard CreateSequenceDictionary -R hg19.fa -O hg19.dict

bgzip -f tests/data/TSVC_variants_defaultlowsetting.vcf
bcftools index tests/data/TSVC_variants_defaultlowsetting.vcf.gz
bcftools view -r chr1:800000-1000000 tests/data/TSVC_variants_defaultlowsetting.vcf.gz -Oz -o tests/data/region.vcf.gz
gzip -dc tests/data/region.vcf.gz > tests/data/region.vcf