indir=$1

CD_HIT="/usr/bin/cd-hit"
CUTOFF_START=0.999
CUTOFF_STEP=0.001
max=30
fasta_re="/.+\/(.+)\.fasta$"


for i in $indir/* 
do
   
    if test -d "$i" 
    then
        echo "Dir $i"
        if [ ! -d "$i/hmmscan_seq_extract/cdhit" ]; then
            mkdir $i/hmmscan_seq_extract/cdhit
        else
            rm -rf $i/hmmscan_seq_extract/cdhit/*    
        fi
        outdir=$i/hmmscan_seq_extract/cdhit
        for j in $i/hmmscan_seq_extract/*
        do
            echo "j: $j"
            if test -f "$j"
            then
                echo "    file $j"
                
                if [[ $j =~ $fasta_re ]]
                then
                    da_name=${BASH_REMATCH[1]}
                    infile=$j
                    echo "da_name: $da_name"
                    
                    count=$(grep -c '>' $infile)
                    if [ $count -gt $max ]
                    then
                        cutoff=$CUTOFF_START
                        while [ $count -gt $max ]
                        do
                            rm $outdir/$da_name-cdhit.fasta*
                            echo "        Executing: $CD_HIT -d 0 -g 1 -T 0 -c $cutoff -i $infile -o $outdir/${da_name}_cdhit.fasta:"
                            $CD_HIT -d 0 -g 1 -T 0 -c $cutoff -i $infile -o $outdir/${da_name}_cdhit.fasta
                            rc=$?
                            if [[ $rc != 0 ]]
                            then
                                exit $rc
                            fi
                            count=$(grep -c '>' $outdir/${da_name}_cdhit.fasta)
                            cutoff=$(echo "$cutoff - $CUTOFF_STEP" | bc)
                        done
                        # make_multi_seq.pl
                        #------------------
                        perl /home/lambda/SOFT/cdhit-master/make_multi_seq.pl $infile $outdir/${da_name}_cdhit.fasta.clstr $outdir/${da_name}_multi_seq 0
                        rc=$?  
                        if [[ $rc != 0 ]]
                        then
                            exit $rc
                        fi
                        # mafft
                        # -----
                        for k in $outdir/${da_name}_multi_seq/* 
                        do
                            if test -f "$k" 
                            then  
                                /usr/local/bin/mafft --thread 12 $k > ${k}_mafft.fasta
                                rc=$?
                                if [[ $rc != 0 ]]
                                then
                                    exit $rc
                                fi
                            fi
                        done
                    else
                        echo "        Executing: cp $infile $outdir/:"
                        cp $infile $outdir/  
                        rc=$?  
                        if [[ $rc != 0 ]]
                        then
                            exit $rc
                        fi
                    fi    
                fi
                
            fi    
        done
    fi
done