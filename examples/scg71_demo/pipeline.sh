if [ "$#" -ne 4 ]; then
    echo "usage: $0 dropest_directory config_file star_index_folder gtf_with_genes"
    echo "example: $0 ~/dropEst/build ~/dropEst/configs/indrop_v3.xml ~/star/mm10/index/ ~/star/mm10/genes.gtf"
    exit 1
fi

dropest_dir=$1
config_file=$2
star_index=$3
gtf_file=$4
cd 01_dropTag
$dropest_dir/droptag -c $config_file ../SCG_71_C2_S2_R2_001.fastq.gz ../SCG_71_C2_S2_R3_001.fastq.gz ../SCG_71_C2_S2_R1_001.fastq.gz
cd ../02_alignment
STAR --genomeDir $star_index --readFilesCommand zcat --outSAMtype BAM Unsorted --readFilesIn ../01_dropTag/SCG_71_C2_S2_R1_001.fastq.gz.tagged.1.fastq.gz
cd ../03_dropEst
$dropest_dir/dropest -w -M -u -G 20 -g $gtf_file -c $config_file ../02_alignment/Aligned.out.bam
