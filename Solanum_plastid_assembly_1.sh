#!/bin/sh
#### This is step 1 in the plastid assembly process. The results from this script are manually assembled into a draft genome.
# chmod 755 ./Solanum_plastid_assembly_1.sh
# nohup ./Solanum_plastid_assembly_1.sh &

###List of accessions here
accessions="bre545968 aba458403 aba458404 ach558032 acr365313 acr365314 acr498204 alb498206 amb365317 amb365362 amb498209 amb498210 amb498212 amb498213 and320345 and561648 and561658 avi498091 avi498092 avi498093 ber218225 ber498105 ber527886 ber545850 bla498214 bra320265 bra558460 bre310931 bre473378 bre498111 bre498218 bre545970 bre545971 bre545981 buk266385 buk365353 buk414155 buk473492 buk473493 buk473494 buk568933 buk568954 bul545751 bul604074 caj230522 can210035 can246533 can265864 can265865 can283084 can442696 can473355 can498226 can498227 can545972 can568969 car283062 car283063 car347759 cha275138 cha320294 cha472816 cha472830 cha500020 cho365328 cho365339 etu498311 gon195186 gon195188 gon195214 gon458393 gou472911 gou472991 gou472995 gou473019 gou473077 gou473106 gou500022 gou537026 gou545865 gou545975 gou545978 gou558067 haw320364 hon473365 hon498067 hon498071 hon545879 hua320327 hua498255 hyp473477 inc473060 inc473067 inc473069 inc473070 inc500048 jam641944 jam664024 kur472924 kur472936 kur472948 kur472952 kur558185 kur558208 lax283088 lax498252 lax607887 lep458378 lep473342 lep473446 lep473451 lep545985 lep545987 lim473468 mar210040 mar310944 med210045 med230507 med320260 med458402 med473496 meg210034 meg473158 meg500029 meg546000 mul210044 mul210052 mul210055 mul275272 mul365336 mul365337 mul365338 mul473349 mul473352 mul498266 pal245763 pam275274 pam275275 pam442697 pam458381 pau473489 phu195191 phu195198 phu225665 phu225693 phu225703 phu243467 phu243468 phu243469 phu258855 pin253214 pin537023 pol161728 pol347770 sog230510 sog365360 spa246536 spa473375 spa473385 spa498134 spa498284 spa498285 spe320299 spe458335 spe458337 spe472966 spe472988 spe472990 ste195204 ste205527 ste230512 ste230513 ste234011 ste255527 ste283141 ste365344 tar217457 tar414152 tar458366 tar473217 tar473218 ver195170 ver275256 ver275260 ver320332 ver320333 ver458370 ver473303 ver473309 ver498010 ver498061 ver500070 ver545745 ver545747 ver545884 ver558150 ver558463 ver558488 vio473396 vio473398 vio498296"
###Reference sequence
plastid_ref="/path/S.tuberosum_chloroplast.fa"

###############################START PROGRAM###############################

for PI in $accessions ; do ###Loop through the list of accessions

##########OTHER PARAMETERS##########
###Reads
fastq1="/path/${PI}_1_clean.fq.gz"
fastq2="/path/${PI}_2_clean.fq.gz"
###RG header
SM=`echo $PI | cut -c 4-9`
ID=`echo $PI | cut -c 1-9`
Read_Group="@RG\tID:${ID}\tLB:barcode\tPL:illumina\tSM:${SM}\tPI:300"
##########END PARAMETERS##########

###Create index for reference (only do once)
if [ ! -s ${plastid_ref::(-2)}dict ]; then #check if already done
    samtools faidx ${plastid_ref}
    bwa index ${plastid_ref}
    java -jar -Xmx4g /usr/local/bin/picard-tools/picard.jar CreateSequenceDictionary REFERENCE=${plastid_ref} OUTPUT=${plastid_ref::(-2)}dict
else
    echo "Reference has been indexed"
fi

mkdir ${PI}_plastid

###Map reads to the reference and only keep those that map. Then convert to bam, and sort.
bwa mem -a -M -R ${Read_Group} -t 20 ${plastid_ref} ${fastq1} ${fastq2} | awk '$3 != "*"' >${PI}_plastid/${PI}_plastid.sam
samtools view -bt ${plastid_ref}.fai -o ${PI}_plastid/${PI}_plastid.bam ${PI}_plastid/${PI}_plastid.sam
samtools sort ${PI}_plastid/${PI}_plastid.bam ${PI}_plastid/${PI}_plastid.s
samtools view -b -F8 ${PI}_plastid/${PI}_plastid.s.bam >${PI}_plastid/${PI}_plastid.sF30.bam
rm ${PI}_plastid/*.sam ${PI}_plastid/${PI}_plastid.bam ${PI}_plastid/${PI}_plastid.s.bam

###Pull reads that map to the plastid 
java -Xmx2g -jar /usr/local/bin/picard-tools/picard.jar SamToFastq INPUT=${PI}_plastid/${PI}_plastid.sF30.bam FASTQ=${PI}_plastid/${PI}.mapped_plastid.1.fastq SECOND_END_FASTQ=${PI}_plastid/${PI}.mapped_plastid.2.fastq TMP_DIR=/dauc2/tmp 

###Denovo assemble reads that map to the plastid and take the *-8.fa file 
abyss-pe k=64 name=${PI}_plastid/${PI}.plastid in="${PI}_plastid/${PI}.mapped_plastid.1.fastq ${PI}_plastid/${PI}.mapped_plastid.2.fastq"
if [ ! -e "${PI}_plastid/${PI}.plastid-8.fa" ]; then
  filename=$(ls ${PI}_plastid/${PI}.plastid-[1-7].fa | sort -n | tail -n1)
  mv $filename  ${PI}_plastid/${PI}.plastid-8.fa
fi
rm ${PI}_plastid/*-1* ${PI}_plastid/*-2* ${PI}_plastid/*-3* ${PI}_plastid/*-4* ${PI}_plastid/*-5* ${PI}_plastid/*-6* ${PI}_plastid/*-7* ${PI}_plastid/*untig* ${PI}_plastid/*contig* ${PI}_plastid/*scaff* ${PI}_plastid/*.dot ${PI}_plastid/*bubb* ${PI}_plastid/*inde* ${PI}_plastid/*unitig* ${PI}_plastid/*stats* ${PI}_plastid/coverage.hist ${PI}_plastid/error.log
bb.revcomp ${PI}_plastid/${PI}.plastid-8.fa >${PI}_plastid/${PI}.plastid-8rc.fa

##run graph to identify how much of reference is there, put plastid contigs in order based on reference
nucmer -maxmatch -p ${PI}_plastid/${PI}.plastidmatch ${plastid_ref} ${PI}_plastid/${PI}.plastid-8.fa
mummerplot -y "[0,175000]" --layout --postscript -p ${PI}_plastid/${PI}.PLASTIDMATCH ${PI}_plastid/${PI}.plastidmatch.delta
show-coords -rd -L 1000 ${PI}_plastid/${PI}.plastidmatch.delta >${PI}_plastid/${PI}.plastidmatch.coords #-L is the minimal alignement length
rm ${PI}_plastid/*.delta ${PI}_plastid/*.fliter ${PI}_plastid/*.rplot ${PI}_plastid/*.gp ${PI}_plastid/*.fplot
###Use the *.coords and *.ps files to manually assemble the chloroplast

done ###for PI in $accessions ; do
###############################END PROGRAM###############################
