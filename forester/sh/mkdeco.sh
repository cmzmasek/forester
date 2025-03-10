VERSION="1.0.0"
DECORATOR="java -Xms256m -Xmx2048m -cp /Users/czmasek/IdeaProjects/forester/forester/java/forester.jar org.forester.application.decorator"

if [ "$#" -ne 5 ]; then
  echo "Usage: mkdeco.sh intree nim-file fasta-file dff-file outtree"
  echo ""
  echo "    Example: mkdeco.sh g_ni.xml g.nim g_ni.fasta g_ni_hmmscan.dff g_d.xml
"
  exit 1
fi

intree=$1
nim_file=$2
fasta_file=$3
dff_file=$4
outfile=$5

echo "Version   : $VERSION"
echo "Intree    : $intree"
echo "nim-file  : $nim_file"
echo "Fasta-file: $fasta_file"
echo "dff-file  : $dff_file"
echo "Outtree   : $outfile"
echo ""
echo ""

$DECORATOR -f=m $intree $fasta_file ___deco0___.xml
rc=$?
if [[ $rc != 0 ]]; then
  exit $rc
fi

$DECORATOR -f=d ___deco0___.xml $dff_file ___deco1___.xml
rc=$?
if [[ $rc != 0 ]]; then
  exit $rc
fi
rm ___deco0___.xml

$DECORATOR -mp -or -f=n ___deco1___.xml $nim_file $outfile
rc=$?
if [[ $rc != 0 ]]; then
  exit $rc
fi
rm ___deco1___.xml

echo "mkdeco successfully completed"
echo ""
