indir=$1
main_outdir=$2

CD_HIT="/usr/bin/cd-hit"
CUTOFF_START=0.999
CUTOFF_STEP=0.001
max=30
fasta_re="/.+\/(.+)\.fasta$"


for i in $indir/* 
do
    if test -d "$i" 
    then
        species_name="$(basename $i)"
        echo "Dir $i"
        echo "Species $pecies_name"
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
                    echo "DA name: $da_name"
                    
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
                        # mafft & consensus
                        # -----------------
                        for cluster in $outdir/${da_name}_multi_seq/* 
                        do
                            cluster_number="$(basename $cluster)"
                            if test -f "$cluster" 
                            then 
                                seqcount=$(grep -c '>' $cluster)
                                if [ $seqcount -gt 1 ]
                                then
                                    /usr/local/bin/mafft --thread 12 $cluster > ${cluster}_mafft.fasta
                                    rc=$?
                                    if [[ $rc != 0 ]]
                                    then
                                        exit $rc
                                    fi
                                    java -Xmx8048m -cp /home/lambda/git/forester/forester/java/forester.jar org.forester.application.msa_consensus ${cluster}_mafft.fasta ${cluster}_mafft_cons.fasta ${da_name}_cons
                                    rc=$?
                                    if [[ $rc != 0 ]]
                                    then
                                        exit $rc
                                    fi
                                    if [ -f "${cluster}_mafft_cons.fasta" ];
                                    then
                                        cp ${cluster}_mafft_cons.fasta $main_outdir/${species_name}_${da_name}_cluster_${cluster_number}_mafft_cons.fasta
                                        rc=$?  
                                        if [[ $rc != 0 ]]
                                        then
                                            exit $rc
                                        fi
                                    else    
                                        cp ${cluster}_mafft_cons_TO_BLAST.fasta $main_outdir/${species_name}_${da_name}_cluster_${cluster_number}_mafft_cons_TO_BLAST.fasta
                                        rc=$?  
                                        if [[ $rc != 0 ]]
                                        then
                                            exit $rc
                                        fi
                                    fi
                                else
                                    cp $cluster ${cluster}_cluster_single.fasta
                                    rc=$?  
                                    if [[ $rc != 0 ]]
                                    then
                                        exit $rc
                                    fi
                                    cp $cluster $main_outdir/${species_name}_${da_name}_cluster_single.fasta
                                    rc=$?  
                                    if [[ $rc != 0 ]]
                                    then
                                        exit $rc
                                    fi
                                fi    
                            fi
                        done
                    else
                        # DAs with less than 30 seqs
                        cp $infile $main_outdir/${species_name}_${da_name}_limited.fasta
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