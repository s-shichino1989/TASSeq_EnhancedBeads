

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.primary_assembly.annotation.gtf.gz

unpigz gencode.vM30.primary_assembly.annotation.gtf.gz

wget http://ftp.ensembl.org/pub/release-107/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
unpigz Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

#build GRCm39 release 107 v30 GENCODE/ENSEMBL STAR reference

Rscript ../../remove_small_RNA_annotations.R gencode.vM30.primary_assembly.annotation.gtf

STAR --runMode genomeGenerate --genomeDir GRCm39_107_v30_nosmRNA/ \
--genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa \
--sjdbGTFfile gencode.vM30.primary_assembly.annotation_filteredsmRNA.gtf --sjdbOverhang 150 --runThreadN 64

pigz gencode.vM30.primary_assembly.annotation.gtf
pigz gencode.vM30.primary_assembly.annotation_filteredsmRNA.gtf
pigz Mus_musculus.GRCm39.dna.primary_assembly.fa

mv ./GRCm39_107_v30_nosmRNA/ /datadrive/Rhapsody_analysis/index/GRCm39_107_v30_nosmRNA/



