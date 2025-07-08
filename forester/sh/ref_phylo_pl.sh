VERSION="0.0.1"
CLEAN_FASTA="python3 /Users/czmasek/Dropbox/PROG/PYTHON/PYCHARM_PROJECTS/TWO/clean_fasta.py"
CDHIT="/Users/czmasek/anaconda3/bin/cd-hit"
MK_TREE="bash /Users/czmasek/IdeaProjects/forester/forester/sh/mktree.sh"
VIPR_X="java -Xmx8048m -cp /Users/czmasek/IdeaProjects/forester/forester/java/forester.jar org.forester.application.vipr_x6"

# Example options file:
# SEQ_TYPE=na
# MINIMAL_LENGTH_A=450
# MINIMAL_LENGTH_B=12000
# CHAR_RATIO=0.999
# CLUSTERING_CUTOFF_A=0.999
# CLUSTERING_CUTOFF_B=0.99
# MAFFT_OPTS_A=--auto --thread 8
# MAFFT_OPTS_B=--auto --thread 8
# MSA_PRO_OPTS_A=-rr=0.2 -rsl=450
# MSA_PRO_OPTS_B=-rr=0.2 -rsl=450
# TREE_INFERENCE_OPTS_A=-B100Zxm
# TREE_INFERENCE_OPTS_B=-B100Zxm
# TREE_A_DESC=this is tree a
# TREE_B_DESC=this is tree b

if [ "$#" -ne 7 ]; then
  echo "Usage: " >&2
  echo "    Example 1: ref_phylo_pl.sh options BVBRC_1900_2025.fasta ../JUNE_10_2025/BVBRC_genome.txt 2025-06-6 tree_a tree_b workdir"
  exit 1
fi

options_file=$1
genomes_file=$2
annotation_file=$3
date_str=$4
tree_a_output_name=$5
tree_b_output_name=$6
workdir=$7

if [ ! -f "$options_file" ]; then
  echo "$options_file does not exist"
  exit -1
fi

seq_type="?"
minimal_length_a=-1
minimal_length_b=-1
char_ratio=-1
clustering_cutoff_a=-1
clustering_cutoff_b=-1
mafft_opts_a="?"
mafft_opts_b="?"
msa_pro_opts_a="?"
msa_pro_opts_b="?"
tree_inference_opts_a="?"
tree_inference_opts_b="?"

while IFS= read -r line; do
  if [[ "$line" =~ SEQ_TYPE=(.+) ]]; then
    seq_type="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MINIMAL_LENGTH_A=(.+) ]]; then
    minimal_length_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MINIMAL_LENGTH_B=(.+) ]]; then
    minimal_length_b="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ CHAR_RATIO=(.+) ]]; then
    char_ratio="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ CLUSTERING_CUTOFF_A=(.+) ]]; then
    clustering_cutoff_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ CLUSTERING_CUTOFF_B=(.+) ]]; then
    clustering_cutoff_b="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MAFFT_OPTS_A=(.+) ]]; then
    mafft_opts_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MAFFT_OPTS_B=(.+) ]]; then
    mafft_opts_b="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MSA_PRO_OPTS_A=(.+) ]]; then
    msa_pro_opts_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MSA_PRO_OPTS_B=(.+) ]]; then
    msa_pro_opts_b="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ TREE_INFERENCE_OPTS_A=(.+) ]]; then
    tree_inference_opts_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ TREE_INFERENCE_OPTS_B=(.+) ]]; then
    tree_inference_opts_b="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ TREE_A_DESC=(.+) ]]; then
    tree_a_desc="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ TREE_B_DESC=(.+) ]]; then
    tree_b_desc="${BASH_REMATCH[1]}"
  fi
done <$options_file

echo "Version             : $VERSION"
echo "Options file        : $options_file"
echo "Genomes file A      : $genomes_file"
echo "Annotations file    : $annotation_file"
echo "Sequence type       : $seq_type"
echo "Min length A        : $minimal_length_a"
echo "Min length B        : $minimal_length_b"
echo "Ratio               : $char_ratio"
echo "Clustering cutoff A : $clustering_cutoff_a"
echo "Clustering cutoff B : $clustering_cutoff_b"
echo "MAFFT options A     : $mafft_opts_a"
echo "MAFFT options B     : $mafft_opts_b"
echo "MSA pro options A   : $msa_pro_opts_a"
echo "MSA pro options B   : $msa_pro_opts_b"
echo "Tree inf options A  : $tree_inference_opts_a"
echo "Tree inf options B  : $tree_inference_opts_b"
echo "Date string         : $date_str"
echo "Tree desc A         : $tree_a_desc"
echo "Tree desc B         : $tree_b_desc"
echo "Tree A output name  : $tree_a_output_name"
echo "Tree B output name  : $tree_b_output_name"
echo "Work dir            : $workdir"
echo ""
echo ""

if [ ! -d "$workdir" ]; then
  mkdir $workdir
fi

if [ ! -z "$(ls -A $workdir)" ]; then
  echo "Working directory is not empty"
  exit -1
fi

if [ -f "$tree_a_output_name" ]; then
  echo "Tree $tree_a_output_name already exists"
  exit -1
fi

if [ -f "$tree_b_output_name" ]; then
  echo "Tree $tree_b_output_name already exists"
  exit -1
fi

if [ ! -f "$genomes_file" ]; then
  echo "$genomes_file does not exist"
  exit -1
fi

if [ ! -f "$annotation_file" ]; then
  echo "$annotation_file does not exist"
  exit -1
fi

$CLEAN_FASTA -ml $minimal_length_a -r $char_ratio -t $seq_type $genomes_file $workdir/genomes_clean_a.fasta
$CLEAN_FASTA -ml $minimal_length_b -r $char_ratio -t $seq_type $genomes_file $workdir/genomes_clean_b.fasta

#only cluster if cutoff between 0 and 1:

$CDHIT -c $clustering_cutoff_a -d 0 -g 1 -i $workdir/genomes_clean_a.fasta -o $workdir/genomes_clean_a_cdhit.fasta
$CDHIT -c $clustering_cutoff_b -d 0 -g 1 -i $workdir/genomes_clean_b.fasta -o $workdir/genomes_clean_b_cdhit.fasta

mkdir $workdir/trees_a
mkdir $workdir/trees_b

cp $workdir/genomes_clean_a_cdhit.fasta $workdir/trees_a/
cp $workdir/genomes_clean_b_cdhit.fasta $workdir/trees_b/

$MK_TREE "$mafft_opts_a" "$msa_pro_opts_a" "-$tree_inference_opts_a" .fasta $workdir/trees_a/
$MK_TREE "$mafft_opts_b" "$msa_pro_opts_b" "-$tree_inference_opts_b" .fasta $workdir/trees_b/

$VIPR_X $workdir/trees_a/deco/genomes_clean_a_cdhit_ni_mafft_tree_raxml_d.xml $annotation_file $tree_a_output_name "$tree_a_desc"
$VIPR_X $workdir/trees_b/deco/genomes_clean_b_cdhit_ni_mafft_tree_raxml_d.xml $annotation_file $tree_b_output_name "$tree_b_desc"
