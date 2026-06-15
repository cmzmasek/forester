Archaeopteryx
=============

Archaeopteryx is an interactive viewer and editor for phylogenetic trees, and is
part of the **forester** toolkit. It reads phyloXML, Newick/New Hampshire
(NH/NHX), and Nexus trees, and supports rich annotation, on-the-fly coloring, and
export to PDF, SVG, EPS, PNG, and other graphics formats.


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

You can also open a tree directly:

```
java -jar forester.jar mytree.xml
```

For very large trees, give the JVM more memory with `-Xmx`:

```
java -Xmx4g -jar forester.jar mytree.xml
```

### 4. The Tools menu

The **Tools** menu collects operations that root, prune, colorize, collapse, and
annotate the current tree. In menu order:

**Rooting**

- **MAD-Root** — Minimal Ancestor Deviation rooting (Tria, Landan & Dagan,
  *Nature Ecology & Evolution* 2017,
  [doi:10.1038/s41559-017-0193](https://doi.org/10.1038/s41559-017-0193)). Chooses
  the root that minimizes how far every pair of tips deviates from a clock-like
  (equidistant-from-their-ancestor) expectation, and annotates each internal
  branch with its MAD value — smaller means a better-supported root location. Turn
  on **Confidence Values** in the left control panel to see them. Requires branch
  lengths and at least three tips.
- **Midpoint-Root** — places the root at the midpoint of the longest tip-to-tip
  path. Fast, and works on any tree with branch lengths.

**Pruning**

- **Delete Selected Nodes** — removes the selected external nodes from the tree.
  Select nodes first by clicking them, or via the **Search** field.
- **Retain Selected Nodes** — the inverse: keeps the selected external nodes and
  deletes all the others.

**Coloring**

- **Colorize Subtrees via Taxonomic Rank** — colors clades by their taxonomy at a
  rank you choose (e.g. order, family, genus). The chooser lists only the ranks
  actually present, with how many nodes carry each and the share of leaves it
  would color; a movable legend mapping colors to taxa is drawn on the tree and
  included in PDF/PNG exports. (This is distinct from **Color by:** in the control
  panel — see below — which colors *leaves* by an arbitrary property.)

**Clearing styles and colors**

- **Delete All Visual Styles From Nodes** — removes every per-node visual style
  (fonts, colors).
- **Delete All Colors From Branches** — removes every branch color.

**Collapsing branches** (the **Collapse Branches** submenu)

- **Collapse Weakly-Supported Branches…** — permanently collapses internal
  branches whose confidence is below a threshold you enter into polytomies
  (multifurcations). Cannot be undone.
- **Collapse Short Branches…** — permanently collapses internal branches shorter
  than a branch-length threshold you enter. Cannot be undone.

**Data retrieval**

- **Fetch Sequence & Taxonomic Data** — looks up additional sequence information
  and detailed taxonomy for the tree's nodes from UniProt / EMBL-GenBank and the
  NCBI taxonomy. Looked-up taxonomy is cached on disk, so subsequent runs and
  rank-coloring are fast.

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

### 6. Other recent additions

- **Settings dialog** — display, node/branch, font, export, and file options in
  one live-apply dialog (replaces the old Options and Type menus).
- **Bundled fonts** — ships three publication-quality fonts (**Source Sans 3**
  default, plus **Liberation Sans** and **Noto Sans**), so figures render
  identically on every machine. Set the tip-label size with the **font-size
  slider** in the left control panel.
- **Taxonomy cache** — taxonomy looked up from NCBI is cached on disk (30 days),
  so re-opening trees of organisms you've already seen is instant. Manage it
  under **Settings → Taxonomy Cache** (toggle, size, clear).
- **Vector export** — true **SVG** and **EPS** output alongside PDF/PNG/TIFF/JPG,
  for publication figures.
- **Adaptive control panel** — "Display Data" checkboxes appear only for data the
  tree actually has; node-symbol **support visualization** (threshold marks or
  size-scaled).


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

All Java library dependencies (FlatLaf, OpenPDF, Apache Commons Codec, OpenChart,
and VectorGraphics2D) and the bundled fonts are vendored in
`forester/java/resources/` and unpacked into the output jar by the build, so there
is no separate dependency-management step.

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
