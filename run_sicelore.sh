#!/bin/bash -x

#check if all software is present 
# need Java.1.8 in JAVA_HOME
#java=`which $JAVA_HOME`
java=`which java`
spoa=`which spoa`
minimap2=`which minimap2`
samtools=`which samtools`
SICEPATH="/usr/local/bin/sicelore-master/Jar/"
#output_dir="${PWD}/output_dir"
#tmp_dir="${output_dir}/tmp/"

if [ -z "$java" ] || [ -z "$spoa" ] || [ -z "$samtools" ] || [ -z "$minimap2" ]
then
    echo -e "\nMissing path to required softwares:"
    echo -e "\tjava=$java"
    echo -e "\tspoa=$spoa"
    echo -e "\tsamtools=$samtools"
    echo -e "\tminimap2=$minimap2"
    echo -e "\nPlease update your \$PATH and rerun.\n\n"
    exit
fi


function usage() #shows how to use this script
{
	echo 
	echo "Run the whole sicelore pipeline once"
	
	echo
	echo "Usage: $0 [-n name of the samples (comma separated)] [-f fastq file paths (comma separated)]"
	echo
	echo "OPTIONS:"
	echo "  [-f | -fastqs Fastq-File]    - Path to Long Read FastQ file (comma separated)"
	echo "  [-r | -ref RefFlat-File]    - Path that contains the prepared refFlat file"
	echo "  [-c | -cel Cellranger outs Folder]    - Path to the cellranger output folder"
	echo "  [-m | -mmi Minimap2 Index files]    - Path to the Minimap2 Index file (.mmi)"
	echo "  [-j | -junc Junction files]    - Path to the Junction file (generated with 'paftools.js gff2bed')"
	echo "  [-v | -variations SNP|Variation files]    - Path to the SNP file (generated as described here: https://github.com/ucagenomix/sicelore#10-calling-single-nucleotide-polymorphisms-cell-by-cell	')"
	echo "  [-o | -outFolder Output folder]	- Path were the sicelore resultsss are written to"
	exit -1
}


#proceed parameters
while test $# -gt 0
	do
		case $1 in
			-h | -help) usage;;
			-f | -fastq) FASTQ=$2;;
			-r | -ref) REFFLAT=$2;;
			-c | -cel) CELLOUT=$2;;
			-m | -mmi) MINIMAPINDEX=$2;;
			-j | -junc) JUNCTIONFILE=$2;;
			-v | -variations) VARFILE=$2;;
			-o | -outFolder) OUTFOLDER=$2;;
			*) echo "ERROR: unknown argument $1" ; exit -1
	esac
	shift 2
	
done


#check input parameters

if [ -z $FASTQ ] ; then
  echo "Parameter: Fastq (-f) required!" && exit 2;
fi
if [ -z $REFFLAT ] ; then
  echo "Parameter: REFFLAT (-r) required!" && exit 2;
fi
if [ -z $OUTFOLDER ] ; then
  echo "Parameter: OUTFOLDER (-o) required!" && exit 2;
fi
if [ -z $CELLOUT ] ; then
  echo "Parameter: CELLOUT (-c) required!" && exit 2;
fi
if [ -z $MINIMAPINDEX ] ; then
  echo "Parameter: MINIMAPINDEX (-m) required!" && exit 2;
fi
if [ -z $JUNCTIONFILE ] ; then
  echo "Parameter: JUNCTIONFILE (-j) required!" && exit 2;
fi
if [ -z $VARFILE ] ; then
  echo "Parameter: SNPFILE (-v) required!" && exit 2;
fi


OUTPREFIX=$(basename ${FASTQ/.fastq/})

SAMPLEFOLDER=$OUTFOLDER$OUTPREFIX"/"
TMPDIR=$SAMPLEFOLDER"tmp"
LOGFILE=$SAMPLEFOLDER$OUTPREFIX".log"

# create output directory


printf "\nCreate Output Directories\n" 
[ -d $OUTFOLDER ] && echo $OUTFOLDER " Exists" || mkdir $OUTFOLDER 
[ -d $SAMPLEFOLDER ] && echo $SAMPLEFOLDER" Exists" || mkdir $SAMPLEFOLDER
[ -d $TMPDIR ] && echo $TMPDIR" Exists" || mkdir $TMPDIR

echo "" | tee -a $LOGFILE
echo "Start Sicelore Pipeline for " $FASTQ | tee $LOGFILE


#Print Input Parameters
echo "" | tee -a $LOGFILE
echo "### Provided Input Parameters ###" | tee -a $LOGFILE
echo "OUTPREFIX = " $OUTPREFIX  | tee -a $LOGFILE
echo "OUTFOLDER = " $OUTFOLDER  | tee -a $LOGFILE
echo "SAMPLEFOLDER = " $SAMPLEFOLDER  | tee -a $LOGFILE
echo ""  && echo ""  | tee -a $LOGFILE
echo "INPUT FILES:"  | tee -a $LOGFILE
echo "REFFLAT = " $REFFLAT  | tee -a $LOGFILE
echo "CELLOUT = " $CELLOUT  | tee -a $LOGFILE
echo "MINIMAPINDEX = " $MINIMAPINDEX  | tee -a $LOGFILE
echo "JUNCTIONFILE = " $JUNCTIONFILE  | tee -a $LOGFILE
echo "VARFILE = " $VARFILE  | tee -a $LOGFILE




# parse illumina bam file
CURRENTOUTFILE=$SAMPLEFOLDER$OUTPREFIX".obj"
echo "Running Illuminaparser [Step 1 of 10]"  | tee -a $LOGFILE
if [ -f "$CURRENTOUTFILE" ]; then
    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
else 
    $java -Xms100G -Xmx170G -jar /usr/local/bin/sicelore-master/Jar/IlluminaParser-1.0.jar -i $CELLOUT/"possorted_genome_bam.bam" -o $SAMPLEFOLDER$OUTPREFIX".obj" -t $CELLOUT"filtered_feature_bc_matrix/barcodes.tsv" -b CB -g GN -u UB 2>&1 |tee -a $LOGFILE
fi
if [ ! -f "$CURRENTOUTFILE" ]; then
echo "Error: Previous Step Failed." | tee -a $LOGFILE
	exit -1
fi

# scan nanopore reads
CURRENTOUTFILE=$SAMPLEFOLDER"/passed/"$OUTPREFIX"FWD.fastq"
echo "Running NanoporeReadScanner [Step 2 of 10]" | tee -a $LOGFILE
if [ -f "$CURRENTOUTFILE" ]; then
    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
else 
	$java -Xms100G -Xmx170G -jar /usr/local/bin/sicelore-master/Jar/NanoporeReadScanner-0.5.jar -i $FASTQ -o $SAMPLEFOLDER 2>&1 |tee -a $LOGFILE
fi
if [ ! -f "$CURRENTOUTFILE" ]; then
echo "Error: Previous Step Failed." | tee -a $LOGFILE
	exit -1
fi

# map reads to genome
CURRENTOUTFILE=$SAMPLEFOLDER"minimap.bam"
echo "Running Minima2 [Step 3 of 10]" | tee -a $LOGFILE
if [ -f "$CURRENTOUTFILE" ]; then
    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
else 
	$minimap2 -ax splice -uf --MD --sam-hit-only -t 20 --junc-bed $JUNCTIONFILE $MINIMAPINDEX $SAMPLEFOLDER"/passed/"$OUTPREFIX"FWD.fastq" -o $SAMPLEFOLDER/minimap.sam | tee -a $LOGFILE
	$samtools view -Sb $SAMPLEFOLDER/minimap.sam -o $SAMPLEFOLDER/minimap.unsorted.bam | tee -a $LOGFILE
	$samtools sort $SAMPLEFOLDER/minimap.unsorted.bam -o $SAMPLEFOLDER/minimap.bam | tee -a $LOGFILE
	$samtools index $SAMPLEFOLDER/minimap.bam | tee -a $LOGFILE
fi
if [ ! -f "$CURRENTOUTFILE" ]; then
echo "Error: Previous Step Failed." | tee -a $LOGFILE
	exit -1
fi

#minimap2 -ax splice -uf --MD --sam-hit-only -t 20 --junc-bed junctions.bed /media/Storage/Wesley/LongReadSequencing/HumanRefGenome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.mmi \
#	/media/Helios_scStorage/Wesley/LongReadSequencing/combined_1st_2nd/CombinedLongRead_01_HF1.fastq > minimap2/CombinedLongRead_01_HF1.sam

# tag reads with gene name
CURRENTOUTFILE=$SAMPLEFOLDER/GE.bam
echo "Running AddGeneNameTag [Step 4 of 10]" | tee -a $LOGFILE
if [ -f "$CURRENTOUTFILE" ]; then
    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
else 
	$java -Xms100G -Xmx170G -jar /usr/local/bin/sicelore-master/Jar/Sicelore-2.0.jar AddGeneNameTag I=$SAMPLEFOLDER/minimap.bam O=$SAMPLEFOLDER/GE.bam REFFLAT=$REFFLAT GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT | tee -a $LOGFILE
	$samtools index $SAMPLEFOLDER/GE.bam | tee -a $LOGFILE
fi
if [ ! -f "$CURRENTOUTFILE" ]; then
echo "Error: Previous Step Failed." | tee -a $LOGFILE
	exit -1
fi



# tag reads with fastq sequence
CURRENTOUTFILE=$SAMPLEFOLDER/GEUS.bam
echo "Running AddBamReadSequenceTag [Step 5 of 10]" | tee -a $LOGFILE
if [ -f "$CURRENTOUTFILE" ]; then
    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
else 
	$java -Xms100G -Xmx170G -jar /usr/local/bin/sicelore-master/Jar/Sicelore-2.0.jar AddBamReadSequenceTag I=$SAMPLEFOLDER/GE.bam O=$SAMPLEFOLDER/GEUS.bam FASTQDIR=$SAMPLEFOLDER/passed/ VALIDATION_STRINGENCY=SILENT | tee -a $LOGFILE
	$samtools index $SAMPLEFOLDER/GEUS.bam | tee -a $LOGFILE
fi
if [ ! -f "$CURRENTOUTFILE" ]; then
echo "Error: Previous Step Failed." | tee -a $LOGFILE
	exit -1
fi



# tag reads with cellBC/UMI barcodes
CURRENTOUTFILE=$SAMPLEFOLDER/GEUS10xAttributes_umifound_.bam
echo "Running NanoporeBC_UMI_finder [Step 6 of 10]" | tee -a $LOGFILE
if [ -f "$CURRENTOUTFILE" ]; then
    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
else 
	$java -Xms170G -Xmx170G -jar /usr/local/bin/sicelore-master/Jar/NanoporeBC_UMI_finder-1.0.jar -i $SAMPLEFOLDER/GEUS.bam -o $SAMPLEFOLDER/GEUS10xAttributes.bam -k $SAMPLEFOLDER$OUTPREFIX".obj" --ncpu 24 --maxUMIfalseMatchPercent 2 --maxBCfalseMatchPercent 5 --logFile $SAMPLEFOLDER/UMI_founder_out.log 2>&1 | tee -a $LOGFILE
	$samtools index $SAMPLEFOLDER/GEUS10xAttributes.bam | tee -a $LOGFILE
	$samtools index $SAMPLEFOLDER/GEUS10xAttributes_umifound_.bam | tee -a $LOGFILE
fi
if [ ! -f "$CURRENTOUTFILE" ]; then
echo "Error: Previous Step Failed." | tee -a $LOGFILE
	exit -1
fi




# generate isoform matrix
#$java -jar -Xmx4g /usr/local/bin/sicelore-master/Jar/Sicelore-2.0.jar IsoformMatrix DELTA=2 METHOD=STRICT GENETAG=GE I=$output_dir/GEUS10xAttributes_umifound_.bam REFFLAT=Gencode/gencode.v18.mm10.refFlat.txt CSV=Barcodes/cellBC.190.tsv OUTDIR=$output_dir PREFIX=sicread VALIDATION_STRINGENCY=SILENT

# compute consensus sequence
CURRENTOUTFILE=$SAMPLEFOLDER/consensus.fq
echo "Running ComputeConsensus [Step 7 of 10]" | tee -a $LOGFILE
if [ -f "$CURRENTOUTFILE" ]; then
    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
else 
	$java -Xms100G -Xmx170G -jar /usr/local/bin/sicelore-master/Jar/Sicelore-2.0.jar ComputeConsensus T=10 I=$SAMPLEFOLDER/GEUS10xAttributes_umifound_.bam O=$SAMPLEFOLDER/consensus.fq TMPDIR=$TMPDIR | tee -a $LOGFILE
fi
if [ ! -f "$CURRENTOUTFILE" ]; then
echo "Error: Previous Step Failed." | tee -a $LOGFILE
	exit -1
fi


# map molecules to genome
CURRENTOUTFILE=$SAMPLEFOLDER/molecule.bam
echo "Running Minimap on Consensus [Step 8 of 10]" | tee -a $LOGFILE
if [ -f "$CURRENTOUTFILE" ]; then
    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
else 
	$minimap2 -ax splice -uf --MD --sam-hit-only -t 4 --junc-bed $JUNCTIONFILE $MINIMAPINDEX $SAMPLEFOLDER/consensus.fq -o $SAMPLEFOLDER/molecule.sam | tee -a $LOGFILE
	$samtools view -Sb $SAMPLEFOLDER/molecule.sam -o $SAMPLEFOLDER/molecule.unsorted.bam | tee -a $LOGFILE
	$samtools sort $SAMPLEFOLDER/molecule.unsorted.bam -o $SAMPLEFOLDER/molecule.bam | tee -a $LOGFILE
	$samtools index $SAMPLEFOLDER/molecule.bam | tee -a $LOGFILE
fi
if [ ! -f "$CURRENTOUTFILE" ]; then
echo "Error: Previous Step Failed." | tee -a $LOGFILE
	exit -1
fi

# add cellBC/UMI tags
CURRENTOUTFILE=$SAMPLEFOLDER/molecule.tags.bam
echo "Running AddBamMoleculeTags on Consensus [Step 9 of 10]" | tee -a $LOGFILE
if [ -f "$CURRENTOUTFILE" ]; then
    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
else 
	$java -Xms100G -Xmx170G -jar /usr/local/bin/sicelore-master/Jar/Sicelore-2.0.jar AddBamMoleculeTags I=$SAMPLEFOLDER/molecule.bam O=$SAMPLEFOLDER/molecule.tags.bam | tee -a $LOGFILE
	$samtools index $SAMPLEFOLDER/molecule.tags.bam | tee -a $LOGFILE
fi
if [ ! -f "$CURRENTOUTFILE" ]; then
echo "Error: Previous Step Failed." | tee -a $LOGFILE
	exit -1
fi


	
# add gene name tag
CURRENTOUTFILE=$SAMPLEFOLDER/molecule.tags.GE.bam
echo "Running AddGeneNameTag on Consensus [Step 10 of 10]" | tee -a $LOGFILE
if [ -f "$CURRENTOUTFILE" ]; then
    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
else 
	$java -Xms100G -Xmx170G -jar /usr/local/bin/sicelore-master/Jar/Sicelore-2.0.jar AddGeneNameTag I=$SAMPLEFOLDER/molecule.tags.bam O=$SAMPLEFOLDER/molecule.tags.GE.bam REFFLAT=$REFFLAT GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT | tee -a $LOGFILE
	$samtools index $SAMPLEFOLDER/molecule.tags.GE.bam | tee -a $LOGFILE
fi
if [ ! -f "$CURRENTOUTFILE" ]; then
echo "Error: Previous Step Failed." | tee -a $LOGFILE
	exit -1
fi



	
# generate molecule isoform matrix
#$java -jar -Xmx4g /usr/local/bin/sicelore-master/Jar/Sicelore-2.0.jar IsoformMatrix DELTA=2 METHOD=STRICT ISOBAM=true GENETAG=GE I=$output_dir/molecule.tags.GE.bam REFFLAT=Gencode/gencode.v18.mm10.refFlat.txt CSV=Barcodes/cellBC.190.tsv OUTDIR=$output_dir PREFIX=sicmol VALIDATION_STRINGENCY=SILENT
#$samtools index $output_dir/sicmol_isobam.bam

# cleaning
#cd $output_dir
#rm -fr failed 190c.clta.illumina.bam.obj consensus.fq GEUS10xAttributes_umifound_.bam GEUS10xAttributes_umifound_.bam.bai molecule.tags.GE.bam molecule.tags.GE.bam.bai GE.bam GE.bam.bai GEUS.bam GEUS.bam.bai GEUS10xAttributes.bam GEUS10xAttributes.bam.bai minimap.bam minimap.bam.bai minimap.sam minimap.unsorted.bam molecule.bam molecule.bam.bai molecule.sam molecule.tags.bam molecule.tags.bam.bai molecule.unsorted.bam out.log passed tmp


# Call SNPs
CURRENTOUTFILE=$SAMPLEFOLDER/snps/snp_snpmatrix.txt
echo "Running SNPMatrix on Consensus [Step 11 of 10]" | tee -a $LOGFILE
if [ -f "$CURRENTOUTFILE" ]; then
    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
else 
	mkdir $SAMPLEFOLDER/snps
	java -jar -Xmx170g -Xms170g /usr/local/bin/sicelore-master/Jar/Sicelore-2.0.jar SNPMatrix I=$SAMPLEFOLDER/molecule.tags.GE.bam CSV=$CELLOUT"filtered_feature_bc_matrix/barcodes.tsv" SNP=$VARFILE O=$SAMPLEFOLDER/snps/
fi
if [ ! -f "$CURRENTOUTFILE" ]; then
echo "Error: Previous Step Failed." | tee -a $LOGFILE
	exit -1
fi

# Call SNPs
#CURRENTOUTFILE=$SAMPLEFOLDER/snps_combined/snp_snpmatrix.txt
#echo "Running SNPMatrix on combined Reads [Step 12 of 10]" | tee -a $LOGFILE
#if [ -f "$CURRENTOUTFILE" ]; then
#    echo "$CURRENTOUTFILE exists." | tee -a $LOGFILE
#else 
#	mkdir $SAMPLEFOLDER/snps_combined
#	samtools merge -c -p $SAMPLEFOLDER/GEUS10xAttributes.bam $SAMPLEFOLDER/GEUS10xAttributes_umifound_.bam -o $SAMPLEFOLDER/GEUS10xAttributes.combined.bam --threads 20  --write-index
#	java -jar -Xmx170g -Xms170g /usr/local/bin/sicelore-master/Jar/Sicelore-2.0.jar SNPMatrix I=$SAMPLEFOLDER/GEUS10xAttributes.combined.bam CSV=$CELLOUT"filtered_feature_bc_matrix/barcodes.tsv" SNP=$VARFILE O=$SAMPLEFOLDER/snps_combined/
#fi
#if [ ! -f "$CURRENTOUTFILE" ]; then
#echo "Error: Previous Step Failed." | tee -a $LOGFILE
#	exit -1
#trfi

