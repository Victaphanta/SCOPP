#!/bin/bash

# THINGS TO DO
# 1) update SUBSETTRIM=subsettrimpre/protrespex such that it placed subtaxa and subgene in the folder name
# 2) update SUBSETANCESTORTRIM with an additional switch for path to guiding tree for fastML
# 3) Maybe only apply BMGE100 for subsettrim-PRO-trespex? Note sure, discuss with Liz

TEMP=`getopt -o m:i:s:a:r:p:g:x:t:b:f: --long argm:,argi:,args:,arga:,argr:,argp:,argg:,argx:,argt:,argb:,argf: -n 'orthology.sh' -- "$@"`
eval set -- "$TEMP"
# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -m|--argm)
            case "$2" in
                "") ARG_M='all' ; shift 2 ;;
                *) ARG_M=$2 ; shift 2 ;;
            esac ;;
        -s|--argi)
            case "$2" in
                "") shift 2 ;;
                *) ARG_S=$2 ; shift 2 ;;
            esac ;;
        -a|--arge)
            case "$2" in
                "") shift 2 ;;
                *) ARG_A=$2 ; shift 2 ;;
            esac ;;
        -r|--argr)
            case "$2" in
                "") shift 2 ;;
                *) ARG_R=$2 ; shift 2 ;;
            esac ;;
        -p|--argp)
            case "$2" in
                "") ARG_P='1' ; shift 2 ;;
                *) ARG_P=$2 ; shift 2 ;;
            esac ;;
        -g|--argg)
            case "$2" in
                 "") shift 2 ;;
                *) ARG_G=$2 ; shift 2 ;;
            esac ;;
        -x|--argx)
            case "$2" in
                 "") shift 2 ;;
                *) ARG_X=$2 ; shift 2 ;;
            esac ;;
        -t|--argt)
            case "$2" in
                 "") shift 2 ;;
                *) ARG_T=$2 ; shift 2 ;;
            esac ;;
	-b|--argb)
            case "$2" in
                 "") shift 2 ;;
                *) ARG_B=$2 ; shift 2 ;;
            esac ;;
	-f|--argf)
            case "$2" in
                 "") shift 2 ;;
                *) ARG_F=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

###################################################################################################
#                                                                                                 #
# MODULE 1 (ASSEMBLE) - Trim and assemble                                                         #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "cleanassemble" ]; then

mkdir -p $PWD/0_assemblies &&
mkdir -p $PWD/0_assemblies/TRIM &&

# Read PE sample names into array
cat $ARG_S | while read i
do
#	java -jar /usr/local/bin/trimmomatic.jar PE -phred64 -threads $ARG_P \
#		$ARG_A/${i}_R1.fastq.gz \
#		$ARG_A/${i}_R2.fastq.gz \
#		$ARG_A/TRIM/${i}_R1.TRIMP.fastq.gz \
#		$ARG_A/TRIM/${i}_R1.TRIMU.fastq.gz \
#		$ARG_A/TRIM/${i}_R2.TRIMP.fastq.gz \
#		$ARG_A/TRIM/${i}_R2.TRIMU.fastq.gz \
#		ILLUMINACLIP:/usr/local/share/trimmomatic/TruSeq3-PE-2.fa:1:35:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:40 2>&1 >/dev/null |\
#	awk '/Input\ Read\ Pairs:/ {print $4,$7,$8,$12,$13,$17,$18,$20,$21}' OFS="\t" > $ARG_A/TRIM/$i.trimsumPE.log &&
#	cat $ARG_A/TRIM/${i}_R1.TRIMP.fastq.gz $ARG_A/TRIM/${i}_R1.TRIMU.fastq.gz $ARG_A/TRIM/${i}_R2.TRIMU.fastq.gz > $ARG_A/TRIM/${i}_R1.TRIMPU.fastq.gz
#	Trinity --seqType fq --max_memory 40G --min_contig_length 100 \
#		--left $ARG_A/TRIM/${i}_R1.TRIMPU.fastq.gz --right $ARG_A/TRIM/${i}_R2.TRIMP.fastq.gz  \
#		--output $ARG_A/${i}.trinity --CPU $ARG_P --full_cleanup --min_kmer_cov 2 &&
	makesomethingNotInterleaved.pl $ARG_A/${i}.trinity.Trinity.fasta > $PWD/0_assemblies/${i}.trinity.NI.fasta
done

cat $ARG_A/TRIM/*.trimsumPE.log > $ARG_A/TRIM/trimsumPE.log
fi
###################################################################################################
#                                                                                                 #
# MODULE 2 (LAST) - Lastal search                                                                 #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "last" ]; then
mkdir -p $PWD/1_last &&
lastdb -p $PWD/1_last/REF $PWD/$ARG_R
# Loop through the trinity assembly files
# First, read sample names into array
readarray -t S < $ARG_S 
for i in "${S[@]}"
do
	pigz -d $ARG_A/$i.trinity.fasta.gz
	makesomethingNotInterleaved.pl $ARG_A/$i.trinity.fasta > $ARG_A/$i.trinity.NI.fasta &&
	lastal -j3 -F15 -fTab -P0 -D1e+10 $PWD/1_last/REF $ARG_A/$i.trinity.NI.fasta > $PWD/1_last/$i.trinity.NI.LASTAL.TAB.txt &&
	lastal -j3 -F15 -fMAF -P0 -D1e+10 $PWD/1_last/REF $ARG_A/$i.trinity.NI.fasta > $PWD/1_last/$i.trinity.NI.LASTAL.MAF.txt &&
	lastal -j3 -F15 -fBlastTab -P0 -D1e+10 $PWD/1_last/REF $ARG_A/$i.trinity.NI.fasta > $PWD/1_last/$i.trinity.NI.LASTAL.BLASTTAB.txt &&
	pigz -f $ARG_A/$i.trinity.fasta

	grep -v '#' $PWD/1_last/$i.trinity.NI.LASTAL.TAB.txt | awk '!x[$7]++' > $PWD/1_last/$i.trinity.NI.LASTAL.TAB.tophits.txt &&
	grep -v '#' $PWD/1_last/$i.trinity.NI.LASTAL.MAF.txt | grep . | sed -n '0~3p' | sed -E 's/\s+/\t/g' | awk '{print $7}' > $PWD/1_last/$i.trinity.NI.LASTAL.MAF.Qlocal.txt &&
	grep -v '#' $PWD/1_last/$i.trinity.NI.LASTAL.TAB.txt > $PWD/1_last/$i.trinity.NI.LASTAL.TAB.noheader.txt &&
	paste $PWD/1_last/$i.trinity.NI.LASTAL.TAB.noheader.txt $PWD/1_last/$i.trinity.NI.LASTAL.MAF.Qlocal.txt > $PWD/1_last/$i.trinity.NI.LASTAL.TAB.Qlocal.txt &&
	awk '!x[$7]++' $PWD/1_last/$i.trinity.NI.LASTAL.TAB.Qlocal.txt > $PWD/1_last/$i.trinity.NI.LASTAL.TAB.Qlocal.tophit.txt &&
	grep -v '#' $PWD/1_last/$i.trinity.NI.LASTAL.BLASTTAB.txt > $PWD/1_last/$i.trinity.NI.LASTAL.BLASTTAB.noheader.txt &&
	paste $PWD/1_last/$i.trinity.NI.LASTAL.TAB.Qlocal.txt $PWD/1_last/$i.trinity.NI.LASTAL.BLASTTAB.noheader.txt > $PWD/1_last/$i.trinity.NI.LASTAL.TAB.Qlocal.plusblast.txt &&
	awk '!x[$7]++' $PWD/1_last/$i.trinity.NI.LASTAL.TAB.Qlocal.plusblast.txt > $PWD/1_last/$i.trinity.NI.LASTAL.TAB.Qlocal.plusblast.tophit.txt
#	pigz -f $ARG_A/$i.trinity.NI.fasta
done
echo "LASTAL SEARCH COMPLETE"
fi
###################################################################################################
#                                                                                                 #
# MODULE 3 (HOMOLOG) - Extract homologous CDS sequences from assemblies                           #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "homolog" ]; then
mkdir -p $PWD/2_homolog &&

pullexons_lastal_20170822.py $ARG_G $ARG_S &&

count=0 &&
cat $ARG_G | while read g
do
	alignAM.py \
		-prot_outfile $PWD/2_homolog/$g.homologs.mafft_aa.fasta \
		-nuc_outfile $PWD/2_homolog/$g.homologs.mafft_nt.fasta \
		-aligner mafft -options " --localpair --maxiterate 1000 --reorder --anysymbol --thread -1 " \
		$PWD/2_homolog/$g.homologs.fasta 
	let count=count+1;
	echo "$count $PWD/2_homolog/$g.homologs.fasta"
done
echo "HOMOLOG ALIGNMENTS COMPLETE"
fi
###################################################################################################
#                                                                                                 #
# MODULE 4 (CONSENSUS) - Collapse contigs by making consenus given mismatch threshold             #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "consensus" ]; then
taxa=$(wc -l < "$ARG_S") 
genes=$(wc -l < "$ARG_G")
mkdir -p $PWD/3_consensus.T${taxa}.G${genes} 
count=0 &&
>ConsensusSummary.T$taxa.G$genes.03.50.15.1.0.txt
cat $ARG_G | while read g
do
	consensus_maker_v2.7.py \
		-mismatches 0.03 \
		-minlength 0.5 \
		-dw_overlap 15 \
		-dw_factor 1 \
		-debug 0 \
		-taxacount $taxa \
		$PWD/2_homolog/$g.homologs.mafft_nt.fasta \
		$PWD/3_consensus.T${taxa}.G${genes}/$g.T$taxa.G$genes.03.50.15.1.0.fasta \
		ConsensusSummary.T$taxa.G$genes.03.50.15.1.0.txt \
		ConsensusStats.T$taxa.G$genes.03.50.15.1.0.txt  &&

	let count=count+1;
	echo "$count CONSENSUS DONE FOR $g"
done

echo "CONSENSUS COMPLETE"
fi

###################################################################################################
#                                                                                                 #
# MODULE 5a (SUBSETTRIM) - subset genes and/or taxa via the -x and -t switch respectively,        #
#                          and trim using BMGE                                                    #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "subsettrimpretrespex" ]||[ "$ARG_M" = "subsettrimposttrespex" ]; then
taxa=$(wc -l < "$ARG_S") 
genes=$(wc -l < "$ARG_G")
subsetgenes=$(wc -l < "$ARG_X")

if [ "$ARG_M" = "subsettrimpretrespex" ]; then
	mkdir -p $PWD/4_SubsetTrim.T${taxa}.G${genes}.S${subsetgenes}
	OUTDIR=$PWD/4_SubsetTrim.T${taxa}.G${genes}.S${subsetgenes}
else 
	mkdir -p $PWD/6_SubsetTrim.T${taxa}.G${genes}.S${subsetgenes}
	OUTDIR=$PWD/6_SubsetTrim.T${taxa}.G${genes}.S${subsetgenes}
echo $OUTDIR
fi

esweep_linux_v1 \
	-i $PWD/3_consensus.T${taxa}.G${genes} \
	-o $OUTDIR \
	-n $PWD/$ARG_T \
	-k $PWD/$ARG_X \
	-d %s.T$taxa.G$genes.03.50.15.1.0.fasta \
	-x \
	-v

cat $ARG_X | while read genesubset
do
	# The following will apply the changes in line numbers 2,4,6,8,......
	perl -nle "s/N|~|n/-/g unless $. & 1; print " $OUTDIR/${genesubset}_unaligned.fasta \
		>$OUTDIR/${genesubset}_aligned.fasta &&
	rm $OUTDIR/${genesubset}_unaligned.fasta &&
	perl -i -p -e "s{^}{TTT} if $. %2 ==0" $OUTDIR/${genesubset}_aligned.fasta &&
	selectSites.pl -x 3 \
		$OUTDIR/${genesubset}_aligned.fasta \
		> $OUTDIR/${genesubset}_aligned.trim.fasta &&
	makesomethingNotInterleaved.pl $OUTDIR/${genesubset}_aligned.trim.fasta \
		> $OUTDIR/${genesubset}_aligned.trim.NI.fasta &&
	degapcodon2.pl $OUTDIR/${genesubset}_aligned.trim.NI.fasta \
		> $OUTDIR/${genesubset}_aligned.DG.fasta &&
	java -jar /usr/local/bin/BMGE100.jar \
		-i $OUTDIR/${genesubset}_aligned.DG.fasta \
		-t CODON \
		-m BLOSUM62 \
		-h 0.5 \
		-w 1 \
		-g  0.5 \
		-oco123f $OUTDIR/${genesubset}_aligned.trim.NI.BMGE.fasta &&
	selectSites.pl -s '4-' $OUTDIR/${genesubset}_aligned.trim.NI.BMGE.fasta \
		> $OUTDIR/${genesubset}_aligned.trim.NI.BMGET.fasta && 
	rm $OUTDIR/${genesubset}_aligned.trim.NI.BMGE.fasta &&
	makesomethingNotInterleaved.pl $OUTDIR/${genesubset}_aligned.trim.NI.BMGET.fasta \
		> $OUTDIR/${genesubset}_aligned.trim.NI.BMGE.NI.fasta &&
	perl -i -nle "s/X|N|~|n/-/g unless $. & 1; print " $OUTDIR/${genesubset}_aligned.trim.NI.BMGE.NI.fasta
done

esweep_linux_v1 \
	-i $OUTDIR \
	-o $OUTDIR \
	-n $PWD/$ARG_T \
	-k $PWD/$ARG_X \
	-d %s_aligned.trim.NI.BMGE.NI.fasta \
	-x \
	-v

fi

############################# MODULE5b SUBSET & RE CONSTRUCT ANCESTOR FOR PREBAIT #################
#                                                                                                 #
# MODULE 5b (SUBSETANCESTORTRIM) -subset genes and/or taxa via the -x and -t switch respectively, #
#                       and construct ancestral sequences                                         #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "subsetancestorbait" ]; then
taxa=$(wc -l < "$ARG_S") 
subsettaxa=$(wc -l < "$ARG_T") 
genes=$(wc -l < "$ARG_G")
subsetgenes=$(wc -l < "$ARG_X")

mkdir -p $PWD/4_SubsetFilterBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes} &&
mkdir -p $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes &&
mkdir -p $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.exons &&

#perl -i -nle "s/N|~|n/-/g unless $. & 1; print " $PWD/3_consensus.T${taxa}.G${genes}/${genesubset}.T$taxa.G$genes.03.50.15.1.0.fasta &&
#selectSites.pl -x 3 \
#	$PWD/4_SubsetFilterBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${genesubset}_unaligned.fasta \
#	> $PWD/4_SubsetFilterBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${genesubset}_unaligned.trim.fasta &&


esweep_linux_v1 \
	-i $PWD/3_consensus.T${taxa}.G${genes} \
	-o $PWD/4_SubsetFilterBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes} \
	-n $PWD/$ARG_T \
	-k $PWD/$ARG_X \
	-d %s.T$taxa.G$genes.03.50.15.1.0.fasta \
	-x \
	-v

count=0 &&
cat $ARG_X | while read genesubset
do
#	wc -L $PWD/4_SubsetFilterBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${genesubset}_unaligned.fasta > $PWD/4_SubsetFilterBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${genesubset}.len1.txt &&
	# The following will apply the changes in line numbers 2,4,6,8,......
	perl -i -nle "s/N|~|n/-/g unless $. & 1; print " $PWD/4_SubsetFilterBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${genesubset}_unaligned.fasta &&
	selectSites.pl -x 3 \
		$PWD/4_SubsetFilterBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${genesubset}_unaligned.fasta \
		> $PWD/4_SubsetFilterBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${genesubset}_unaligned.trim.fasta &&

	# Ancestral reconstruction using FASTML
	FastML_Wrapper.pl \
		--MSA_File $PWD/4_SubsetFilterBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${genesubset}_unaligned.trim.fasta \
		--Tree RAxML_bestTree.barcoca.PUNCTOIDEA.nwk \
		--seqType NUC \
		--outDir $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset} &&
	cp $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}/seq.joint.txt $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.joint.fas &&
	cp $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}/seq.marginal.txt $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.fas &&
	rm -rf $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset} &&
	
	cat $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.fas > $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.fas &&
	extractSequence_beginingMatch.pl \
		$PWD/3_consensus.T${taxa}.G${genes}/${genesubset}.T$taxa.G$genes.03.50.15.1.0.fasta \
		Lottia_gigantea.GCA_000327385.1.cds.all >> $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.fas &&
	perl -i -nle "s/-|~|n/N/g unless $. & 1; print " $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.fas &&
	degapcodon2.pl $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.fas \
		> $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.DG.fas &&

	alignAM.py \
		-prot_outfile $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.mafft_aa.fasta \
		-nuc_outfile $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.mafft_nt.fasta \
		-aligner mafft -options " --localpair --maxiterate 1000 --reorder --anysymbol --thread -1 " \
		$PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.DG.fas &&

	perl -i -nle "s/N|~|n/-/g unless $. & 1; print " $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.mafft_nt.fasta &&
	selectSites.pl -x 3 $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.mafft_nt.fasta \
		> $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.mafft_nt.trim.fasta
	if [ $count -eq 0 ] 
	then 
		grep '^>' $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes/${genesubset}.seq.marginal.plusog.mafft_nt.trim.fasta | awk 'sub(/^>/, "")' | awk -F'~' '{print $1}' > T.${subsettaxa}.plusOG.plusAncestral.txt
	fi

	let count=count+1;
done
echo "Lottia_gigantea.GCA_000327385.1.cds.all" >> T.${subsettaxa}.plusOG.plusAncestral.txt &&
head -n 1 $ARG_F > ${subsetgenes}.$ARG_F && 
awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' $ARG_X $ARG_F >> ${subsetgenes}.$ARG_F

ulimit -n 25000 &&
esweep_linux_v2 \
	-i $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.genes \
	-o $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.exons \
	-s $PWD/T.${subsettaxa}.plusOG.plusAncestral.txt \
	-p $PWD/${subsetgenes}.$ARG_F \
	-g %s.seq.marginal.plusog.mafft_nt.trim.fasta \
	-x \
	-v 
fi


############################# MODULE6 MAKE preBAIT ################################################
#                                                                                                 #
# The module collates all exons into a single fasta file per species (or ancestral node)          #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "prebait" ]; then
taxa=$(wc -l < "$ARG_S")
subsettaxa=$(wc -l < "$ARG_T")
genes=$(wc -l < "$ARG_G")
subsetgenes=$(wc -l < "$ARG_X")

head -n 1 $ARG_F > ${subsetgenes}.$ARG_F && 
awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' $ARG_X $ARG_F >> ${subsetgenes}.$ARG_F


# subtitute ambiguous sites in tip sequences with corresponding site from most basal ancestral reconstruction.
# UPDATE - need to automate which sequence is reference 

# tail -n +2 $PWD/${subsetgenes}.$ARG_F | cut -f 2 > $PWD/${subsetgenes}.$ARG_F.exonnames &&
cat $PWD/${subsetgenes}.$ARG_F.exonnames | while read EXON
do
	paste -sd'\t\n' $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.exons/${EXON}_unaligned.fasta  | tr -d '>' | tr -d '\r' > $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.exons/${EXON}_unaligned.txt
	awk '{t[NR]=$1;s[NR]=$2} END{r=1; while(r<15){printf ">%s\n",t[r]; for(i=1;i<=length(s[1]);i++) {if(substr(s[r],i,1)~/[ACGT-]/) printf substr(s[r],i,1); else printf substr(s[8],i,1)} printf "\n";r++}}' $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.exons/${EXON}_unaligned.txt > $PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.exons/${EXON}_substit.fasta
done

	if [ "$ARG_B" = "max" ]; then
		mkdir -p $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.MAX &&
		tail -n +2 $PWD/${subsetgenes}.$ARG_F  |sort -k1,1 -k3,3gr | sort -u -k1,1 --merge |cut -f 2 > $PWD/${subsetgenes}.$ARG_F.biggestexonname &&

		cat taxa.OG.Ancestor.txt | while read TAXA
		do
			> $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.MAX/${TAXA}.preBAIT.fasta
			cat $PWD/${subsetgenes}.$ARG_F.biggestexonname | while read EXON
			do
				filterbyname.sh \
					in=$PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.exons/${EXON}_substit.fasta \
					out=stdout.fasta \
					include=t \
					names=${TAXA} >> $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.MAX/${TAXA}.preBAIT.fasta &&
				sed -i "s/>$TAXA/>${EXON}_${TAXA}/g" $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.MAX/${TAXA}.preBAIT.fasta
			done
			makesomethingNotInterleaved.pl $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.MAX/${TAXA}.preBAIT.fasta \
				> $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.MAX/${TAXA}.preBAIT.NI.fasta &&
			perl -i -nle "s/N//g unless $. & 1; print " $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.MAX/${TAXA}.preBAIT.NI.fasta &&
			reformat.sh \
				in=$PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.MAX/${TAXA}.preBAIT.NI.fasta \
				out=$PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.MAX/${TAXA}.preBAIT.MAX.MIN100.fasta \
				minlen=121 \
				ow=t \
				fastawrap=1000000
		done
		echo "PREBAIT COMPLETE"
	else
		mkdir -p $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes} &&
#		tail -n +2 $PWD/${subsetgenes}.$ARG_F | cut -f 2 > $PWD/${subsetgenes}.$ARG_F.exonnames &&
		cat taxa.OG.Ancestor.txt | while read TAXA
#		cut -f 2 $PWD/${subsetgenes}.$ARG_F >$PWD/${subsetgenes}.$ARG_F.exonnames
		do 
			> $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${TAXA}.preBAIT.fasta
			cat $PWD/${subsetgenes}.$ARG_F.exonnames | while read EXON
			do
				selectSeqs.pl -m ${TAXA} \
					$PWD/5_ancestral.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}.exons/${EXON}_substit.fasta \
					>> $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${TAXA}.preBAIT.fasta &&
				sed -i "s/>$TAXA/>${EXON}_${TAXA}/g" $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${TAXA}.preBAIT.fasta
			done
			makesomethingNotInterleaved.pl $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${TAXA}.preBAIT.fasta \
				> $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${TAXA}.preBAIT.NI.fasta &&
			perl -i -nle "s/N|-//g unless $. & 1; print " $PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${TAXA}.preBAIT.NI.fasta &&
			reformat.sh \
				in=$PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${TAXA}.preBAIT.NI.fasta \
				out=$PWD/6_preBait.T${taxa}.G${genes}.ST${subsettaxa}.SG${subsetgenes}/${TAXA}.preBAIT.MIN100.fasta \
				minlen=100 \
				ow=t \
				fastawrap=1000000
		done
		echo "PREBAIT COMPLETE"
	fi
fi
