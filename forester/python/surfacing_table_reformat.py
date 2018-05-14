import csv

name_col = 1
da_col = 0
ids_col = 3

infile = 'DA_SPECIES_IDS_MAP_mod.csv'
outfile = 'ID_SOGNAME_DA_MAP.txt'

with open(infile, 'r') as in_table, open(outfile, 'w') as out_table:
    reader = csv.reader(in_table, delimiter='\t')

    out_table.write('#ID(GENBANK OR VIPR)\tSOG NAME\tDOMAIN ARCHITECTURE\n')

    for row in reader:
        name = row[name_col]
        da = row[da_col]
        ids = row[ids_col]

        name = name.strip('\"')
        da = da.strip('\"')
        ids = ids.strip('\"')

        idss = ids.split(',')
        print(idss)
        #   if len(idss) > 0:
        for identifier in idss:
            out_table.write(f'{identifier}\t{name}\t{da}\n')
        else:
            print(f"warning: no mapping for {name}")

out_table.close()