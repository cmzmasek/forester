for i in * 
do
    if test -f "$i" 
    then
       echo "Doing something to $i"
       cd-hit -c 0.95 -i $i -o cdhit095/$i
    fi
done