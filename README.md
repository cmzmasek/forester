Using Programs From Forester
============================



hmmscan_seq_extract
-------------------
```
java -cp path/to/forester.jar org.forester.application.hmmscan_seq_extract

hmmscan_seq_extract [options] <hmmscan output> <fasta file (used as input for hmmscan)> <output dir>

 options:
 ie: max (inclusive) iE-value
 mrel: min (inclusive) relative envelope length ratio
```


Cladinator3
-----------
```
java -Xmx8048m path/to/forester.jar org.forester.application.cladinator3

Usage:

cladinator3 [options] <input tree(s) file> [output table file]

 options:
  -s=<separator>     : the annotation-separator to be used (default: ".")
  -m=<mapping table> : to map node names to appropriate annotations (tab-separated, two columns) (default: no mapping)
  -x                 : to enable extra processing of annotations (e.g. "Q16611|A.1.1" becomes "A.1.1")
  -xs=<separator>    : the separator for extra annotations (default: "|")
  -xk                : to keep extra annotations (e.g. "Q16611|A.1.1" becomes "A.1.1.Q16611")
  -S=<pattern>       : special processing with pattern (e.g. "(\d+)([a-z]+)_.+" for changing "6q_EF42" to "6.q")
  -rs                : to remove the annotation-separator in the output (e.g. the ".")
  --q=<pattern>      : expert option: the regular expression pattern for the query (default: "_#\d+_M=(.+)" for pplacer output)

Examples:

 cladinator3 pp_out_tree.sing.tre result.tsv
 cladinator3 -s=. pp_out_tree.sing.tre result.tsv
 cladinator3 -s=_ -m=map.tsv pp_out_trees.sing.tre result.tsv
 cladinator3 -x -xs=& -xk pp_out_trees.sing.tre result.tsv
 cladinator3 -x -xs="|" pp_out_trees.sing.tre result.tsv
 cladinator3 -x -xk -m=map.tsv pp_out_trees.sing.tre result.tsv
 cladinator3 -m=map.tsv -S='(\d+)([a-z?]*)_.+' pp_out_trees.sing.tre result.tsv
```


Archaeoptryx
------------
```
java -Xmx2048m -cp path/to/forester.jar org.forester.archaeopteryx.Archaeopteryx

Usage:

Archaeopteryx [-c configuration file] [treefile]

```





