#!/bin/bash

# THINGS TO DO
# 1) update SUBSETTRIM=subsettrimpre/protrespex such that it placed subtaxa and subgene in the folder name
# 3) Maybe only apply BMGE100 for subsettrim-PRO-trespex?

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

mkdir -p $ARG_A/TRIM &&
mkdir -p $PWD/1_assemblies &&

# Read PE sample names into array
cat $ARG_S | while read i
do
	java -jar /usr/local/bin/trimmomatic.jar PE -phred33 -threads $ARG_P \
		$ARG_A/${i}_R1.fastq.gz \
		$ARG_A/${i}_R2.fastq.gz \
		$ARG_A/TRIM/${i}_R1.TRIMP.fastq.gz \
		$ARG_A/TRIM/${i}_R1.TRIMU.fastq.gz \
		$ARG_A/TRIM/${i}_R2.TRIMP.fastq.gz \
		$ARG_A/TRIM/${i}_R2.TRIMU.fastq.gz \
		ILLUMINACLIP:/usr/local/share/trimmomatic/TruSeq3-PE-2.fa:1:35:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:40  |\
	awk '/Input\ Read\ Pairs:/ {print $4,$7,$8,$12,$13,$17,$18,$20,$21}' OFS="\t" > $ARG_A/TRIM/$i.trimsumPE.log &&
	cat $ARG_A/TRIM/${i}_R1.TRIMP.fastq.gz $ARG_A/TRIM/${i}_R1.TRIMU.fastq.gz $ARG_A/TRIM/${i}_R2.TRIMU.fastq.gz > $ARG_A/TRIM/${i}_R1.TRIMPU.fastq.gz &&
	Trinity --seqType fq --max_memory 40G --min_contig_length 100 \
		--left $ARG_A/TRIM/${i}_R1.TRIMPU.fastq.gz --right $ARG_A/TRIM/${i}_R2.TRIMP.fastq.gz  \
		--output $PWD/1_assemblies/${i}.trinity --CPU $ARG_P --full_cleanup --min_kmer_cov 2 &&
	makesomethingNotInterleaved.pl $PWD/1_assemblies/${i}.trinity.Trinity.fasta > $PWD/1_assemblies/${i}.trinity.NI.fasta
done

cat $ARG_A/TRIM/*.trimsumPE.log > $ARG_A/TRIM/trimsumPE.log
fi
###################################################################################################
#                                                                                                 #
# MODULE 2 (LAST) - Lastal search                                                                 #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "last" ]; then
mkdir -p $PWD/2_last &&
lastdb -p $PWD/2_last/REF $PWD/$ARG_R
# Loop through the trinity assembly files
# First, read sample names into array
readarray -t S < $ARG_S 
for i in "${S[@]}"
do
	lastal -j3 -F15 -fTab -P0 -D1e+10 $PWD/2_last/REF $PWD/1_assemblies/$i.trinity.NI.fasta > $PWD/2_last/$i.trinity.NI.LASTAL.TAB.txt &&
	lastal -j3 -F15 -fMAF -P0 -D1e+10 $PWD/2_last/REF $PWD/1_assemblies/$i.trinity.NI.fasta > $PWD/2_last/$i.trinity.NI.LASTAL.MAF.txt &&
	lastal -j3 -F15 -fBlastTab -P0 -D1e+10 $PWD/2_last/REF $PWD/1_assemblies/$i.trinity.NI.fasta > $PWD/2_last/$i.trinity.NI.LASTAL.BLASTTAB.txt &&
#	pigz -f $PWD/1_assemblies/$i.trinity.fasta

	grep -v '#' $PWD/2_last/$i.trinity.NI.LASTAL.TAB.txt | awk '!x[$7]++' > $PWD/2_last/$i.trinity.NI.LASTAL.TAB.tophits.txt &&
	grep -v '#' $PWD/2_last/$i.trinity.NI.LASTAL.MAF.txt | grep . | sed -n '0~3p' | sed -E 's/\s+/\t/g' | awk '{print $7}' > $PWD/2_last/$i.trinity.NI.LASTAL.MAF.Qlocal.txt &&
	grep -v '#' $PWD/2_last/$i.trinity.NI.LASTAL.TAB.txt > $PWD/2_last/$i.trinity.NI.LASTAL.TAB.noheader.txt &&
	paste $PWD/2_last/$i.trinity.NI.LASTAL.TAB.noheader.txt $PWD/2_last/$i.trinity.NI.LASTAL.MAF.Qlocal.txt > $PWD/2_last/$i.trinity.NI.LASTAL.TAB.Qlocal.txt &&
	awk '!x[$7]++' $PWD/2_last/$i.trinity.NI.LASTAL.TAB.Qlocal.txt > $PWD/2_last/$i.trinity.NI.LASTAL.TAB.Qlocal.tophit.txt &&
	grep -v '#' $PWD/2_last/$i.trinity.NI.LASTAL.BLASTTAB.txt > $PWD/2_last/$i.trinity.NI.LASTAL.BLASTTAB.noheader.txt &&
	paste $PWD/2_last/$i.trinity.NI.LASTAL.TAB.Qlocal.txt $PWD/2_last/$i.trinity.NI.LASTAL.BLASTTAB.noheader.txt > $PWD/2_last/$i.trinity.NI.LASTAL.TAB.Qlocal.plusblast.txt &&
	awk '!x[$7]++' $PWD/2_last/$i.trinity.NI.LASTAL.TAB.Qlocal.plusblast.txt > $PWD/2_last/$i.trinity.NI.LASTAL.TAB.Qlocal.plusblast.tophit.txt
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
mkdir -p $PWD/3_homolog &&

pullexons_lastal_20170822.py $ARG_G $ARG_S &&

count=0 &&
cat $ARG_G | while read g
do
	alignAM.py \
		-prot_outfile $PWD/3_homolog/$g.homologs.mafft_aa.fasta \
		-nuc_outfile $PWD/3_homolog/$g.homologs.mafft_nt.fasta \
		-aligner mafft -options " --localpair --maxiterate 1000 --reorder --anysymbol --thread -1 " \
		$PWD/3_homolog/$g.homologs.fasta 
	let count=count+1;
	echo "$count $PWD/3_homolog/$g.homologs.fasta"
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
mkdir -p $PWD/4_consensus.T${taxa}.G${genes} 
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
		$PWD/3_homolog/$g.homologs.mafft_nt.fasta \
		$PWD/4_consensus.T${taxa}.G${genes}/$g.T$taxa.G$genes.03.50.15.1.0.fasta \
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
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "subsettrimpre" ]||[ "$ARG_M" = "subsettrimpost" ]; then
taxa=$(wc -l < "$ARG_S") 
genes=$(wc -l < "$ARG_G")
subsetgenes=$(wc -l < "$ARG_X")

if [ "$ARG_M" = "subsettrimpre" ]; then
	mkdir -p $PWD/5_SubsetTrim.T${taxa}.G${genes}.S${subsetgenes}
	OUTDIR=$PWD/5_SubsetTrim.T${taxa}.G${genes}.S${subsetgenes}
else 
	mkdir -p $PWD/6_SubsetTrim.T${taxa}.G${genes}.S${subsetgenes}
	OUTDIR=$PWD/6_SubsetTrim.T${taxa}.G${genes}.S${subsetgenes}
echo $OUTDIR
fi

esweep_linux_v1 \
	-i $PWD/4_consensus.T${taxa}.G${genes} \
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
