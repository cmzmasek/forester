for i in * 
do
    if test -f "$i" 
    then
       echo "Doing something to $i"
       cd-hit -c 0.90 -i $i -o cdhit090/$i
    fi
done