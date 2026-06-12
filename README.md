Archaeopteryx
=============

Archaeopteryx is an interactive viewer and editor for phylogenetic trees, and is
part of the **forester** toolkit. It reads phyloXML, Newick/New Hampshire
(NH/NHX), and Nexus trees, and supports rich annotation, on-the-fly coloring, and
export to PDF, PNG, and other graphics formats.


For Users
---------

### 1. Install Java

Archaeopteryx needs **Java 21 or newer** (a Java runtime is enough to run it; a
full JDK also works). Check what you already have:

```
java -version
```

If it reports version `21` or higher, you are ready. Otherwise install Eclipse
Temurin, a free, well-supported OpenJDK build:

- **macOS** — with [Homebrew](https://brew.sh):
  ```
  brew install --cask temurin@21
  ```
  or download the installer from
  <https://adoptium.net/temurin/releases/?version=21>

- **Windows** — with [winget](https://learn.microsoft.com/windows/package-manager/):
  ```
  winget install EclipseAdoptium.Temurin.21.JDK
  ```
  or download the `.msi` installer from
  <https://adoptium.net/temurin/releases/?version=21>

- **Linux** — from your distribution's package manager, for example:
  ```
  sudo apt install openjdk-21-jre     # Debian / Ubuntu
  sudo dnf install java-21-openjdk    # Fedora / RHEL
  ```
  or download a build from
  <https://adoptium.net/temurin/releases/?version=21>

### 2. Download Archaeopteryx

Download the ready-to-run `forester.jar`:

```
curl -L -o forester.jar https://github.com/cmzmasek/forester/raw/master/forester/java/forester.jar
```

(or, in the GitHub repository, open `forester/java/forester.jar` and click
**Download raw file**). The jar is self-contained — every required library is
bundled inside it, so there is nothing else to install.

### 3. Launch it

```
java -jar forester.jar
```

You can also open a tree directly, and/or load a configuration file:

```
java -jar forester.jar mytree.xml
java -jar forester.jar -c my_configuration_file mytree.xml
```

For very large trees, give the JVM more memory with `-Xmx`:

```
java -Xmx4g -jar forester.jar mytree.xml
```


For Developers
--------------

### Clone

```
git clone https://github.com/cmzmasek/forester.git
cd forester
```

### Prerequisites

- **JDK 21 or newer** — the build targets `release 21`; an older JDK fails with
  *"release version 21 not supported"*.
- **Apache Ant**.

All Java library dependencies (FlatLaf, iText, Apache Commons Codec, and
OpenChart) are vendored in `forester/java/resources/` and unpacked into the
output jar by the build, so there is no separate dependency-management step.

### Build

From the `forester/java` directory:

```
cd forester/java
ant all
```

This compiles the sources and produces the self-contained
`forester/java/forester.jar`, whose manifest main class is
`org.forester.archaeopteryx.Archaeopteryx`.

> If your default `java` is older than 21, point Ant at a newer JDK, e.g.
> `JAVA_HOME=/path/to/jdk-21 ant all`.

### Run the test suite

The suite resolves its test data relative to the repository root, so run it from
there with the compiled classes on the classpath:

```
java -Duser.dir="$(pwd)" -cp forester/java/classes org.forester.test.Test
```


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

Output example:
```
Input tree                 : clade_analysis_test_1_2_A.xml
Annotation-separator       : .
Query pattern              : _#\d+_M=(.+)
Output table               : test
Number of input trees      : 1
Ext. nodes in input tree   : 28

Results:

#Tree # Query	Assignment Confidence	Brackets           Conclusion              Placement count
1       Q	    A.3.1.2    1.0	       [A.3.1.2, A.3.1.2]	member of clade A.3.1.2	1
```
