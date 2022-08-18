

wget http://ftp.ensembl.org/pub/release-107/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.107.gtf.gz

unpigz Rattus_norvegicus.mRatBN7.2.107.gtf.gz

wget http://ftp.ensembl.org/pub/release-107/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz

unpigz Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz

#build mRatBN7.2 release 107 ENSEMBL STAR reference

Rscript ../../Rscripts/remove_small_RNA_annotations_ensembl.R Rattus_norvegicus.mRatBN7.2.107.gtf

STAR --runMode genomeGenerate --genomeDir mRatBN7.2_107_nosmRNA/ \
--genomeFastaFiles Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
--sjdbGTFfile Rattus_norvegicus.mRatBN7.2.107_filteredsmRNA.gtf --sjdbOverhang 150 --runThreadN 64

pigz Rattus_norvegicus.mRatBN7.2.107.gtf
pigz Rattus_norvegicus.mRatBN7.2.107_filteredsmRNA.gtf
pigz Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa

mv ./mRatBN7.2_107_nosmRNA/ /datadrive/Rhapsody_analysis/index/mRatBN7.2_107_nosmRNA/



