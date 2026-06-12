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

### 4. Rooting a tree

The **Tools** menu offers two ways to (re)root a tree with branch lengths:

- **Midpoint-Root** — places the root at the midpoint of the longest tip-to-tip
  path. Fast, and works on any tree with branch lengths.
- **MAD-Root** — Minimal Ancestor Deviation rooting (Tria, Landan & Dagan,
  *Nature Ecology & Evolution* 2017,
  [doi:10.1038/s41559-017-0193](https://doi.org/10.1038/s41559-017-0193)). It
  chooses the root that minimizes how far every pair of tips deviates from a
  clock-like (equidistant-from-their-ancestor) expectation, and annotates each
  internal branch with its MAD value — smaller means a better-supported root
  location. Turn on **Confidence Values** in the left control panel to see the
  values. Requires branch lengths and at least three tips.

### 5. Coloring leaves by a property

The **Color by:** dropdown in the left control panel colors the leaves on the
fly by the value of a chosen phyloXML property (e.g. host, country, year). It
recomputes for whatever (sub)tree is currently displayed, and the legend is
included in PDF/PNG exports.

When a property has more distinct values than the palette has colors, the colors
are assigned **most-frequent-first**, so the most common values get the most
distinct colors (cycling only affects the rarest values). The legend then lists
the **most frequent** values — re-sorted alphabetically for readability — with a
`… +N more` footer for the remainder. (A few properties are special-cased:
`year` is shown as a continuous gradient, and `country`/`host` are grouped by the
part before a `:` / `;` qualifier.)


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

The headless test suite (`org.forester.test.Test`) resolves its test data
relative to the repository root, so run it from there with the compiled classes
on the classpath:

```
java -Duser.dir="$(pwd)" -cp forester/java/classes org.forester.test.Test
```

It prints `Failed tests: 0` on success. Numerically involved code is tested
thoroughly: for example, MAD rooting is validated against an independent
brute-force implementation — agreeing on the chosen root **and** every
per-branch support value — across a wide variety of random tree shapes (binary,
multifurcating, caterpillar, star) and sizes, so an accidental regression in the
algorithm is caught.

A handful of GUI integration tests need a display (and FlatLaf), so they are kept
out of the headless suite (and out of the shipped jar) and run individually, for
example:

```
java -cp forester/java/classes org.forester.archaeopteryx.SubSuperTreeButtonsTest
```

New or changed code should come with tests — see the existing `*Test` classes for
the established patterns.


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
