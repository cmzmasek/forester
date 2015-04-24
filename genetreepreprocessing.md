# Gene tree preprocessing #

## Purpose ##

Use EBI database to add species information to a gene tree which has sequence identifiers as labels.

## Usage ##
```
java -Xmx1024m -cp
path/to/forester.jar org.forester.application.gene_tree_preprocess <input tree in NH, NHX, Nexus, or phyloXML>
```



## Details ##

Output consists of three files:
  * input-name\_preprocessed\_gene\_tree.phylo.xml
  * input-name\_species\_present.txt
  * input-name\_removed\_nodes.txt


## Download ##

Download forester.jar here: http://code.google.com/p/forester/downloads/list