from sys import argv
import re

infile = argv[ 1 ]

print( 'Infile: ', infile )

match = 0
lessspecific = 0
qu = 0
na = 0
no_match = 0

with open( infile ) as f:
    for line in f:
        line = line.strip()
        if not line.startswith('#'):
            elements = line.split('\t')
            if elements[ 1 ] == 'Matching Clades':
                if (elements[ 0 ].endswith(elements[ 2 ])):
                    match+=1
                else:
                    my_regex = r".+\|" + re.escape(elements[ 2 ]) + r"\."
                    if  re.search(my_regex, elements[ 0 ]): 
                        lessspecific+=1
                    elif elements[ 2 ] == '?':
                        qu+=1
                        print('?       :  ', line)    
                    elif elements[ 0 ].endswith( 'NA' ):
                        na+=1
                        print('NA      :  ', line)
                    else:
                        no_match+=1
                        print('no match:  ', line)
                        
print()
print( 'Match                            :', match )
print( 'Less specific match              :', lessspecific )
print( 'No match: result undeceided ("?"):',  qu )
print( 'No match: target is "NA"         :',  na )
print( 'No match                         :',  no_match )
