VERSION="0.0.2"
CLEAN_FASTA="python3 /Users/czmasek/Dropbox/PROG/PYTHON/PYCHARM_PROJECTS/TWO/clean_fasta.py"
CDHIT="/Users/czmasek/anaconda3/bin/cd-hit"
MAKE_MULTI_SEQ="perl /Users/czmasek/SOFT/VARIA/make_multi_seq.pl"
MK_TREE="bash /Users/czmasek/IdeaProjects/forester/forester/sh/mktree.sh"
VIPR_X="java -Xmx8048m -cp /Users/czmasek/IdeaProjects/forester/forester/java/forester.jar org.forester.application.vipr_x6"
PHYLOXML_TO_NH="java -Xmx8048m -cp /Users/czmasek/IdeaProjects/forester/forester/java/forester.jar org.forester.application.phyloxml2nh"

# Example options file:
# SEQ_TYPE=na
# MINIMAL_LENGTH_A=450
# MINIMAL_LENGTH_B=12000
# CHAR_RATIO=0.999
# CLUSTERING_CUTOFF_A=0.999
# CLUSTERING_CUTOFF_B=0.99
# ADDITIONAL_SEQS_A=
# ADDITIONAL_SEQS_B=
# MAFFT_OPTS_A=--auto --thread 8
# MAFFT_OPTS_B=--auto --thread 8
# MSA_PRO_OPTS_A=-rr=0.2 -rsl=450
# MSA_PRO_OPTS_B=-rr=0.2 -rsl=450
# TREE_INFERENCE_OPTS_A=-B100Zxm
# TREE_INFERENCE_OPTS_B=-B100Zxm
# TREE_A_DESC=this is tree a
# TREE_B_DESC=this is tree b

if [ "$#" -ne 8 ]; then
  echo "Usage: " >&2
  echo "    Example 1: ref_phylo_pl.sh options seqs.fasta seqs.fasta BVBRC_genome.txt 2025-06-6 tree_a tree_b workdir"
  exit 1
fi

options_file=$1
genomes_file_a=$2
genomes_file_b=$3
annotation_file=$4
date_str=$5
tree_a_output_name=$6
tree_b_output_name=$7
workdir=$8

if [ ! -f "$options_file" ]; then
  echo "$options_file does not exist"
  exit -1
fi

seq_type=""
minimal_length_a=100
minimal_length_b=100
char_ratio=0.99
clustering_cutoff_a=""
clustering_cutoff_b=""
mafft_opts_a="--auto"
mafft_opts_b="--auto"
msa_pro_opts_a=""
msa_pro_opts_b=""
tree_inference_opts_a=""
tree_inference_opts_b=""
add_a=""
add_b=""

while IFS= read -r line; do
  if [[ "$line" =~ SEQ_TYPE=[[:space:]]*(.+) ]]; then
    seq_type="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MINIMAL_LENGTH_A=[[:space:]]*(.+) ]]; then
    minimal_length_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MINIMAL_LENGTH_B=[[:space:]]*(.+) ]]; then
    minimal_length_b="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ CHAR_RATIO=[[:space:]]*(.+) ]]; then
    char_ratio="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ CLUSTERING_CUTOFF_A=[[:space:]]*(.+) ]]; then
    clustering_cutoff_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ CLUSTERING_CUTOFF_B=[[:space:]]*(.+) ]]; then
    clustering_cutoff_b="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MAFFT_OPTS_A=[[:space:]]*(.+) ]]; then
    mafft_opts_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MAFFT_OPTS_B=[[:space:]]*(.+) ]]; then
    mafft_opts_b="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MSA_PRO_OPTS_A=[[:space:]]*(.+) ]]; then
    msa_pro_opts_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MSA_PRO_OPTS_B=[[:space:]]*(.+) ]]; then
    msa_pro_opts_b="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ TREE_INFERENCE_OPTS_A=[[:space:]]*(.+) ]]; then
    tree_inference_opts_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ TREE_INFERENCE_OPTS_B=[[:space:]]*(.+) ]]; then
    tree_inference_opts_b="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ TREE_A_DESC=[[:space:]]*(.+) ]]; then
    tree_a_desc="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ TREE_B_DESC=[[:space:]]*(.+) ]]; then
    tree_b_desc="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ ADDITIONAL_SEQS_A=[[:space:]]*(.+) ]]; then
    add_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ ADDITIONAL_SEQS_B=[[:space:]]*(.+) ]]; then
    add_b="${BASH_REMATCH[1]}"
  fi
done <$options_file

echo ""
echo "ref_phylo_pl settings:"
echo ""
echo "Version             : $VERSION"
echo "Options file        : $options_file"
echo "Genomes file A      : $genomes_file_a"
echo "Genomes file B      : $genomes_file_b"
echo "Add to A            : $add_a"
echo "Add to B            : $add_b"
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
echo "Tree A output base  : $tree_a_output_name"
echo "Tree B output base  : $tree_b_output_name"
echo "Work dir            : $workdir"
echo ""
echo ""
echo ""
echo ""

if [ "$seq_type" != "na" ] && [ "$seq_type" != "aa" ]; then
  echo "SEQ_TYPE must be either 'na' or 'aa'"
  exit -1
fi

if [ "${#tree_inference_opts_a}" -lt 3 ]; then
  echo "TREE_INFERENCE_OPTS_A must be set"
  exit -1
fi

if [ "${#tree_inference_opts_b}" -lt 3 ]; then
  echo "TREE_INFERENCE_OPTS_B must be set"
  exit -1
fi

if [ -f "${tree_a_output_name}.xml" ]; then
  echo "Tree ${tree_a_output_name}.xml already exists"
  exit -1
fi

if [ -f "${tree_b_output_name}.xml" ]; then
  echo "Tree ${tree_b_output_name}.xml already exists"
  exit -1
fi

if [ ! -f "$genomes_file_a" ]; then
  echo "$genomes_file_a does not exist"
  exit -1
fi

if [ ! -f "$genomes_file_b" ]; then
  echo "$genomes_file_b does not exist"
  exit -1
fi

if [ ! -f "$annotation_file" ]; then
  echo "$annotation_file does not exist"
  exit -1
fi

if [ "${#clustering_cutoff_a}" -gt 0 ]; then
  if [ -d "clusters_a" ]; then
    echo "Directory for cluster output 'clusters_a' already exist"
    exit -1
  fi
fi

if [ "${#clustering_cutoff_b}" -gt 0 ]; then
  if [ -d "clusters_b" ]; then
    echo "Directory for cluster output 'clusters_b' already exist"
    exit -1
  fi
fi

if [ ! -d "$workdir" ]; then
  mkdir $workdir
fi

# Sequence cleanup
# ----------------
if [ ! -f $workdir/genomes_clean_a.fasta ]; then
  $CLEAN_FASTA -u t -ml $minimal_length_a -r $char_ratio -t $seq_type $genomes_file_a $workdir/genomes_clean_a.fasta
  rc=$?
  if [[ $rc != 0 ]]; then
    exit $rc
  fi
fi

if [ ! -f $workdir/genomes_clean_b.fasta ]; then
  $CLEAN_FASTA -u t -ml $minimal_length_b -r $char_ratio -t $seq_type $genomes_file_b $workdir/genomes_clean_b.fasta
  rc=$?
  if [[ $rc != 0 ]]; then
    exit $rc
  fi
fi

# Clustering
# ----------
if [ ! -f $workdir/genomes_clean_a_cdhit.fasta ]; then
  if [ "${#clustering_cutoff_a}" -gt 0 ]; then
    $CDHIT -c $clustering_cutoff_a -d 0 -g 1 -i $workdir/genomes_clean_a.fasta -o $workdir/genomes_clean_a_cdhit.fasta
    rc=$?
    if [[ $rc != 0 ]]; then
      exit $rc
    fi
    $MAKE_MULTI_SEQ $workdir/genomes_clean_a.fasta $workdir/genomes_clean_a_cdhit.fasta.clstr clusters_a 2
  else
    cp $workdir/genomes_clean_a.fasta $workdir/genomes_clean_a_cdhit.fasta
  fi
fi

if [ ! -f $workdir/genomes_clean_b_cdhit.fasta ]; then
  if [ "${#clustering_cutoff_b}" -gt 0 ]; then
    $CDHIT -c $clustering_cutoff_b -d 0 -g 1 -i $workdir/genomes_clean_b.fasta -o $workdir/genomes_clean_b_cdhit.fasta
    rc=$?
    if [[ $rc != 0 ]]; then
      exit $rc
    fi
    $MAKE_MULTI_SEQ $workdir/genomes_clean_b.fasta $workdir/genomes_clean_b_cdhit.fasta.clstr clusters_b 2
  else
    cp $workdir/genomes_clean_b.fasta $workdir/genomes_clean_b_cdhit.fasta
  fi
fi

# Addition of base sequences (optional)
# -------------------------------------
if [ "${#add_a}" -gt 0 ]; then
  cat $workdir/genomes_clean_a_cdhit.fasta $add_a >$workdir/genomes_clean_a_cdhit_added.fasta
  $CLEAN_FASTA -u t -ml 5 -r 0.1 -t $seq_type $workdir/genomes_clean_a_cdhit_added.fasta $workdir/genomes_clean_a_cdhit_added_u.fasta
  rc=$?
  if [[ $rc != 0 ]]; then
    exit $rc
  fi
  rm $workdir/genomes_clean_a_cdhit_added.fasta
  rm $workdir/genomes_clean_a_cdhit.fasta
  mv $workdir/genomes_clean_a_cdhit_added_u.fasta $workdir/genomes_clean_a_cdhit.fasta
fi

if [ "${#add_b}" -gt 0 ]; then
  cat $workdir/genomes_clean_b_cdhit.fasta $add_b >$workdir/genomes_clean_b_cdhit_added.fasta
  $CLEAN_FASTA -u t -ml 5 -r 0.1 -t $seq_type $workdir/genomes_clean_b_cdhit_added.fasta $workdir/genomes_clean_b_cdhit_added_u.fasta
  rc=$?
  if [[ $rc != 0 ]]; then
    exit $rc
  fi
  rm $workdir/genomes_clean_b_cdhit_added.fasta
  rm $workdir/genomes_clean_b_cdhit.fasta
  mv $workdir/genomes_clean_b_cdhit_added_u.fasta $workdir/genomes_clean_b_cdhit.fasta
fi

# MSA and tree inference
# ----------------------
mkdir $workdir/trees_a
mkdir $workdir/trees_b

cp $workdir/genomes_clean_a_cdhit.fasta $workdir/trees_a/
cp $workdir/genomes_clean_b_cdhit.fasta $workdir/trees_b/

$MK_TREE "$mafft_opts_a" "$msa_pro_opts_a" "-$tree_inference_opts_a" .fasta $workdir/trees_a/
rc=$?
if [[ $rc != 0 ]]; then
  exit $rc
fi
$MK_TREE "$mafft_opts_b" "$msa_pro_opts_b" "-$tree_inference_opts_b" .fasta $workdir/trees_b/
rc=$?
if [[ $rc != 0 ]]; then
  exit $rc
fi

# Annotation
# ----------
$VIPR_X $workdir/trees_a/deco/genomes_clean_a_cdhit_ni_mafft_tree_raxml_d.xml $annotation_file ${tree_a_output_name}.xml "${tree_a_desc} ${date_str}"
rc=$?
if [[ $rc != 0 ]]; then
  exit $rc
fi
$VIPR_X $workdir/trees_b/deco/genomes_clean_b_cdhit_ni_mafft_tree_raxml_d.xml $annotation_file ${tree_b_output_name}.xml "${tree_b_desc} ${date_str}"
rc=$?
if [[ $rc != 0 ]]; then
  exit $rc
fi

# Finish up
# ---------
$PHYLOXML_TO_NH ${tree_a_output_name}.xml ${tree_a_output_name}.nh
rc=$?
if [[ $rc != 0 ]]; then
  exit $rc
fi
$PHYLOXML_TO_NH ${tree_b_output_name}.xml ${tree_b_output_name}.nh
rc=$?
if [[ $rc != 0 ]]; then
  exit $rc
fi
cp $workdir/trees_a/genomes_clean_a_cdhit_mafft.fasta ${tree_a_output_name}_mafft.fasta
cp $workdir/trees_b/genomes_clean_b_cdhit_mafft.fasta ${tree_b_output_name}_mafft.fasta
