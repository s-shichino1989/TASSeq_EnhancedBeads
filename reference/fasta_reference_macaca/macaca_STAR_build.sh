

wget http://ftp.ensembl.org/pub/release-107/gtf/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_6.0.107.gtf.gz

unpigz Macaca_fascicularis.Macaca_fascicularis_6.0.107.gtf.gz

wget http://ftp.ensembl.org/pub/release-107/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa.gz

unpigz Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa.gz

#build macaca6 release 107 ENSEMBL STAR reference

Rscript ../../Rscripts/remove_small_RNA_annotations_ensembl.R Macaca_fascicularis.Macaca_fascicularis_6.0.107.gtf

STAR --runMode genomeGenerate --genomeDir macaca6_107_nosmRNA/ \
--genomeFastaFiles Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa \
--sjdbGTFfile Macaca_fascicularis.Macaca_fascicularis_6.0.107_filteredsmRNA.gtf --sjdbOverhang 150 --runThreadN 64

pigz Macaca_fascicularis.Macaca_fascicularis_6.0.107.gtf
pigz Macaca_fascicularis.Macaca_fascicularis_6.0.107_filteredsmRNA.gtf
pigz Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa

mv ./macaca6_107_nosmRNA/ /datadrive/Rhapsody_analysis/index/macaca6_107_nosmRNA/



