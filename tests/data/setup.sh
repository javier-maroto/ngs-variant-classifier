# Downloads the data needed for the tests
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/ion_exome/HG002_NA24385_SRR1767406_IonXpress_020_rawlib_24028.bam

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
samtools faidx hg19.fa
picard CreateSequenceDictionary -R hg19.fa -O hg19.dict