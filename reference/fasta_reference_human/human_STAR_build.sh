

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.primary_assembly.annotation.gtf.gz

unpigz gencode.v41.primary_assembly.annotation.gtf.gz

wget http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

unpigz Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#build GRCh38 release 107 / v41 GENCODE/ENSEMBL STAR reference

Rscript ../../remove_small_RNA_annotations.R gencode.v41.primary_assembly.annotation.gtf

STAR --runMode genomeGenerate --genomeDir GRCh38_107_v41_nosmRNA/ \
--genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile gencode.v41.primary_assembly.annotation_filteredsmRNA.gtf --sjdbOverhang 150 --runThreadN 64

pigz -f gencode.v41.primary_assembly.annotation.gtf
pigz -f gencode.v41.primary_assembly.annotation_filteredsmRNA.gtf
pigz -f Homo_sapiens.GRCh38.dna.primary_assembly.fa

mv ./GRCh38_107_v41_nosmRNA/ /datadrive/Rhapsody_analysis/index/GRCh38_107_v41_nosmRNA/



