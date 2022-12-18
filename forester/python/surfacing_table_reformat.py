import csv
import sys

csv.field_size_limit(sys.maxsize)

name_col = 1
da_col = 0
acc_col = 3

acc_to_da = dict()
acc_to_acc_name_da = dict()

infile = 'DA_SPECIES_IDS_MAP.txt'
outfile = 'ID_SOGNAME_DA_MAP.txt'

with open(infile, 'r') as in_table, open(outfile, 'w') as out_table:
    reader = csv.reader(in_table, delimiter='\t')

    out_table.write('#ACC\tSOG NAME\tDOMAIN ARCHITECTURE\n')

    for row in reader:
        name = row[name_col]
        da = row[da_col]
        accs = row[acc_col]

        name = name.strip('\"')
        da = da.strip('\"')
        accs = accs.strip('\"')

        accessors = accs.split(',')
      
        if len(accessors) > 0:
            for acc in accessors:
                if acc[-1:] == '_':
                    acc = acc[:-1]
                new_acc_na_da_str =  acc + '\t' + name + '\t' + da
                if acc in acc_to_da:
                    print('multiple: ' + acc_to_acc_name_da[ acc ])
                    print('        : ' + new_acc_na_da_str)
                    if da != acc_to_da[ acc ]:
                          if len(da) > len( acc_to_da[ acc ] ):
                               acc_to_da[ acc ] = da
                               acc_to_acc_name_da[ acc ] =  new_acc_na_da_str
                    print('keeping : ' + acc_to_acc_name_da[ acc ])     
                    print()
                else:
                    acc_to_da[ acc ] = da
                    acc_to_acc_name_da[ acc ] =  new_acc_na_da_str
            
                
        else:
            print(f"Error: no mapping for {name}")
            quit()

    for acc in acc_to_acc_name_da:
        out_table.write(f'{acc_to_acc_name_da[acc]}\n')
out_table.close()
print('wrote: ' + outfile )


