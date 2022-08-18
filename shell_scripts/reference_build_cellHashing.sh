#build bowtie2 index for cell hasing
#build date 20210127

threads=`grep -ce '^processor\s\+:' /proc/cpuinfo`
tmp=`echo $((threads/2))`
bowtie2-build --threads ${tmp} -f ./reference/hashtags_fasta/humanHashtag_Biolegend.fa humanHashtag14
bowtie2-build --threads ${tmp} -f ./reference/hashtags_fasta/humanSampleTag_BD.fa humanSampleTag12
bowtie2-build --threads ${tmp} -f ./reference/hashtags_fasta/mouseHashtag_Biolegend.fa mouseHashtag15
bowtie2-build --threads ${tmp} -f ./reference/hashtags_fasta/mouseSampleTag_BD.fa mouseSampleTag_BD
bowtie2-build --threads ${tmp} -f ./reference/hashtags_fasta/StreptavidinHashtag.fa Streptavidin_Hashtag
bowtie2-build --threads ${tmp} -f ./reference/hashtags_fasta/Hashtag_snRNAseq.fa Hashtag_snRNAseq


sleep 1

mv humanHashtag14.1.bt2 ./index/humanHashtag14.1.bt2
mv humanHashtag14.2.bt2 ./index/humanHashtag14.2.bt2
mv humanHashtag14.3.bt2 ./index/humanHashtag14.3.bt2
mv humanHashtag14.4.bt2 ./index/humanHashtag14.4.bt2
mv humanHashtag14.rev.1.bt2 ./index/humanHashtag14.rev.1.bt2
mv humanHashtag14.rev.2.bt2 ./index/humanHashtag14.rev.2.bt2

mv humanSampleTag12.1.bt2 ./index/humanSampleTag12.1.bt2
mv humanSampleTag12.2.bt2 ./index/humanSampleTag12.2.bt2
mv humanSampleTag12.3.bt2 ./index/humanSampleTag12.3.bt2
mv humanSampleTag12.4.bt2 ./index/humanSampleTag12.4.bt2
mv humanSampleTag12.rev.1.bt2 ./index/humanSampleTag12.rev.1.bt2
mv humanSampleTag12.rev.2.bt2 ./index/humanSampleTag12.rev.2.bt2

mv mouseHashtag15.1.bt2 ./index/mouseHashtag15.1.bt2
mv mouseHashtag15.2.bt2 ./index/mouseHashtag15.2.bt2
mv mouseHashtag15.3.bt2 ./index/mouseHashtag15.3.bt2
mv mouseHashtag15.4.bt2 ./index/mouseHashtag15.4.bt2
mv mouseHashtag15.rev.1.bt2 ./index/mouseHashtag15.rev.1.bt2
mv mouseHashtag15.rev.2.bt2 ./index/mouseHashtag15.rev.2.bt2

mv mouseSampleTag_BD.1.bt2 ./index/mouseSampleTag_BD.1.bt2
mv mouseSampleTag_BD.2.bt2 ./index/mouseSampleTag_BD.2.bt2
mv mouseSampleTag_BD.3.bt2 ./index/mouseSampleTag_BD.3.bt2
mv mouseSampleTag_BD.4.bt2 ./index/mouseSampleTag_BD.4.bt2
mv mouseSampleTag_BD.rev.1.bt2 ./index/mouseSampleTag_BD.rev.1.bt2
mv mouseSampleTag_BD.rev.2.bt2 ./index/mouseSampleTag_BD.rev.2.bt2

mv Streptavidin_Hashtag.1.bt2 ./index/Streptavidin_Hashtag.1.bt2
mv Streptavidin_Hashtag.2.bt2 ./index/Streptavidin_Hashtag.2.bt2
mv Streptavidin_Hashtag.3.bt2 ./index/Streptavidin_Hashtag.3.bt2
mv Streptavidin_Hashtag.4.bt2 ./index/Streptavidin_Hashtag.4.bt2
mv Streptavidin_Hashtag.rev.1.bt2 ./index/Streptavidin_Hashtag.rev.1.bt2
mv Streptavidin_Hashtag.rev.2.bt2 ./index/Streptavidin_Hashtag.rev.2.bt2

mv Hashtag_snRNAseq.1.bt2 ./index/Hashtag_snRNAseq.1.bt2
mv Hashtag_snRNAseq.2.bt2 ./index/Hashtag_snRNAseq.2.bt2
mv Hashtag_snRNAseq.3.bt2 ./index/Hashtag_snRNAseq.3.bt2
mv Hashtag_snRNAseq.4.bt2 ./index/Hashtag_snRNAseq.4.bt2
mv Hashtag_snRNAseq.rev.1.bt2 ./index/Hashtag_snRNAseq.rev.1.bt2
mv Hashtag_snRNAseq.rev.2.bt2 ./index/Hashtag_snRNAseq.rev.2.bt2


sleep 1

