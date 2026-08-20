Archaeopteryx
=============

Archaeopteryx is an interactive viewer and editor for phylogenetic trees, and is
part of the **forester** toolkit. It reads phyloXML, Newick/New Hampshire
(NH/NHX), and Nexus trees — including annotated
[**BEAST / BEAST X** output](#beast-and-beast-x-output) — and supports rich
annotation, on-the-fly coloring, and export to PDF, SVG, EPS, PNG, and other
graphics formats.


Download & Install
------------------

The easiest way to run Archaeopteryx is a **native installer**. Each one bundles
its own Java 21 runtime, so there is **nothing else to install** — no separate
Java needed.

Download the installer for your platform from the latest release:

**<https://github.com/cmzmasek/forester/releases>**

> The apps are not code-signed or notarized yet, so each platform shows a
> first-launch security prompt (steps below). Getting past it once is enough.

### macOS (Apple Silicon)

1. Download `Archaeopteryx-<version>.dmg` and open it.
2. Drag **Archaeopteryx** into **Applications**.
3. On the first launch, **right-click** (or Control-click) the app and choose
   **Open**, then **Open** again — this clears macOS's "unidentified developer"
   warning. After that, launch it normally.

> **Intel Macs:** the `.dmg` is Apple-Silicon (arm64) only. On an Intel Mac,
> use the [jar](#run-from-the-jar) below.

### Windows

1. Download `Archaeopteryx-<version>.msi` and run it. It adds a Start-menu entry
   (and, optionally, a desktop shortcut) and lets you choose the install folder.
2. Because the installer is unsigned, SmartScreen may warn — click
   **More info → Run anyway**.

### Linux (Debian / Ubuntu)

```
sudo apt install ./archaeopteryx_<version>_amd64.deb
```

(or `sudo dpkg -i archaeopteryx_<version>_amd64.deb`). Launch it from your
applications menu, or by running `archaeopteryx`.


Run From the Jar
----------------

If there is no installer for your platform (for example, an **Intel Mac**), or
you prefer a single self-contained file, run the jar with your own Java.

### 1. Install Java

Running the jar needs **Java 21 or newer** (a Java runtime is enough; a full JDK
also works). Check what you already have:

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

### 2. Download the jar

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


Try It: Demo Trees
------------------

Don't have a tree handy? Archaeopteryx ships with a small gallery of example
trees, each **pre-configured to show off a capability** — open one and it comes
up already colored, banded, or laid out, no setup required.

Just launch Archaeopteryx and pick a tree from **File → Demo Trees**:

- **Color Tips by Metadata** — a tree colored by a categorical property
- **Annotation Columns** — tip-aligned color strips and a numeric heat-map
- **Protein Domain Architectures** — multi-domain proteins drawn to scale
- **Ancestral State Pies** — a discrete geographic trait as posterior pies
- **Node Age Spindles** — divergence-time uncertainty (point age + 95% HPD) as
  tapered spindles
- **SARS-CoV-2 Time Tree** — a tip-dated viral tree on a calendar-year axis
- **Dinosaur Time Tree** — a dated archosaur tree (with *Archaeopteryx*!) on the
  geologic time scale
- **Ammonite Time Tree** — an all-extinct fossil clade with FAD/LAD range bars
- **Tree of Life (Deep Time)** — a time-calibrated tree back to LUCA (~3.8 Ga)
- **Tanglegram** — two trees (gophers vs. their lice) compared side by side

The demos are bundled inside the jar, so they are always available. (The same
trees, plus more, live in [`forester/demo/`](forester/demo/) for browsing on
GitHub.)


Searching Trees
---------------

Two independent **search boxes** on the left control panel (**A** and **B**) find
and highlight matching nodes; a node matched by A (**red**), by B, or by **both**
(**teal**) is shown in a distinct colour, so two searches can be compared at a
glance. Search B’s colour is chosen under **Settings → Labels & Colors →
Found/Selected Colors** — **Electric Violet** (default), **Neon Magenta**, or **Emerald Green** —
each picked to stay legible on a white background, and the choice is remembered
across restarts.

When both boxes carry a query, a **Combine:** control appears below them: keep the
two highlights **independent** (the default), or fold them into one result set —
**A AND B** (matches both) or **A OR B** (matches either) — which then drives the
highlight, step-through, counter, and export.

Each box has two dropdowns — **what** to search and **how** to match — above its
query field:

- **Field** — the node data to search. The list is tailored to the loaded tree, so
  you only see fields it actually has, and the labels match the **Display Data**
  checkboxes. **Any Text** (the default) searches every text field — and your custom
  annotation properties — at once; or pick a specific one — **Node Name**, a taxonomy
  field (**Taxonomy Scientific**, **Taxonomy Common**, **Taxonomy Code**, **Taxonomy
  Identifier**, **Taxonomy Synonym**, **Taxonomy Lineage**), a sequence field (**Seq
  Name**, **Gene Name**, **Gene Symbol**, **Seq Accession**), **Annotation**,
  **Domain**, or **any custom phyloXML property** by its reference (e.g. `data:host`).
  Numeric fields — **Branch Length**, **Support / Confidence**, and numeric
  properties — are offered as well, as are the tree's **structure** fields (prefixed
  `Structure:`): **Clade Size (tips)**, **Number of Children**, **Depth from Root
  (edges)**, **Distance from Root** (when the tree has branch lengths), and **Node
  Type** (leaf / internal / root) — so you can, for example, find every clade with
  more than 50 tips, or every unresolved node (children > 2).
- **Match** — how the query is compared. For a text field: **contains** (the
  default), **starts with**, **ends with**, **whole word**, or **regular
  expression**. For a numeric field the operators switch to plain-language
  comparisons — **equals**, **not equal**, **less than**, **at most**, **greater
  than**, **at least** — and **range** (which reveals a second box for the upper
  bound).

When you search a specific text field, the query box **suggests the values that
field actually has in the tree**, filtered as you type — pick one to match it exactly.
This makes categorical fields (**Node Type**, **Taxonomy Code**, an annotation column
like `data:host`) point-and-click, and saves you from mistyping a value.

Two shared options sit above the boxes: **Match Case** and **Inverse** (select the
nodes that do *not* match). Within a text query, `,` is a logical **OR** and `+` a
logical **AND** (e.g. `kinase, phosphatase`, or `human + receptor`); both are treated
literally in a regular-expression search. Your field and match choices are remembered
as you work, so switching fields or navigating between trees doesn't reset them.

Step through the hits with the **◀ / ▶** buttons beside the boxes, or **View → Find
Next / Find Previous** (**⌘G / ⌘⇧G**) — each jump centres the next match in the view.
Under **Settings → Labels & Colors**, **Bold Found Labels**, **Dim Non-Matches**, and
**Pulse Found Nodes** make the matches stand out further.

Whenever any nodes are highlighted, a **Found / Selected: N** counter appears at the
right of the menu bar. Search hits and manual selection are one and the same in
Archaeopteryx, so this is a single running total of the distinct highlighted nodes; it
hides itself when nothing is highlighted. Hover it for the breakdown by search box
(**A** / **B**) and manual **Selected** nodes.


BEAST and BEAST X Output
------------------------

Archaeopteryx reads the annotated trees produced by **BEAST**, **BEAST 2**, and
**BEAST X** (the current BEAST 2 release line), including **TreeAnnotator**
maximum-clade-credibility (MCC) summaries. Both output shapes are supported:

- annotated **Nexus** (TreeAnnotator's `.tree` / `.trees` output), and
- annotated **Newick / NHX** with FigTree-style `[&key=value, ...]` comment
  blocks on nodes and branches.

Just open the file — parsing of these tags is **on by default** (toggle under
**Settings → File Reading → "Parse BEAST-style extended Newick/Nexus tags"**).
Each annotation is mapped onto the viewer's existing display features:

| BEAST annotation | Becomes | Turn it on with |
| --- | --- | --- |
| `posterior` | Branch support (confidence) | **Confidence Values**; support coloring / symbols |
| node age `height` / `height_median` / `height_mean` + `height_95%_HPD={lo,hi}` (or `height_range`) | Node age with a 95% HPD interval | **Node Age Bars (HPD)** — auto-enabled on load for a dated tree with HPD intervals |
| discrete / geographic traits (e.g. a phylogeographic `location`) with posterior state sets | **Ancestral-state pie charts** | the **"Ancestral pie:"** dropdown (appears automatically when the tree carries such a trait) |
| any other field (`rate`, `length_*`, custom traits, …) | A node property `beast:<key>` | **Color by**, **Size by**, and **Annotation Columns** (numeric traits render as gradients / bars) |

Nothing is discarded: recognized fields become native structures (support, node
dates, pies), and every remaining field is preserved as a `beast:*` property you
can color, size, or tabulate. A malformed field is skipped rather than aborting
the load, so real-world TreeAnnotator files open cleanly.

To turn a dated MCC tree into a time tree, display it as a **phylogram** (the
`P`/`A` buttons) with **Node Age Bars (HPD)** on; add the **Scale Axis** for a
labeled time axis.

The node-age overlay has two shapes (**Settings → Overlays → Data Overlays → Node
age shape**): a flat **Bar** across the 95% HPD interval (the FigTree convention),
or a **Spindle** — a tapered lens that peaks at the point estimate and narrows to
the HPD bounds, so you can see *where* the estimate sits within its interval. The
spindle is a schematic of the *summarized* uncertainty (the point estimate + 95%
HPD), not the raw posterior density — a summary (MCC) tree doesn't carry the
per-node posterior sample.


Geologic Time Axis
------------------

For a dated, time-calibrated tree — where the branch lengths are geologic time
(millions of years) — Archaeopteryx can draw the **international geologic time
scale** beneath the tree instead of a plain numeric axis. Turn it on under
**Settings → Overlays → Time Axis → Geologic (ICS)**. The axis appears when the
tree is shown as a **phylogram** (branch lengths = time), as two coloured, named
bands — **System/Period** over **Series/Epoch** — so a clade's position along the
time axis reads directly against the named geologic intervals (Cretaceous,
Jurassic, Triassic, …).

The Time Axis is **per tree**: Archaeopteryx reads the appropriate axis from each
tree's own `<date>` values (their unit and magnitude), so a geologic Dinosaur tree
in one tab and a calendar-dated SARS-CoV-2 tree in another each show the right axis
at the same time — no global switch to flip. The Settings dropdown lets you
override the axis for the current tab (or turn it off), and when you **save** the
tree, a deliberate choice travels with it (restored on reload).

The axis follows the layout: it runs along the bottom in the **root-left**
orientation, down the breadth side in the **root-on-top / root-on-bottom**
orientations, and becomes concentric coloured **rings** (period bands from the
centre outward) in the **circular** layout — the iTOL-style geologic disc. (It is
not shown in the unrooted layout, which has no single time axis to band.) In the
rectangular orientations the axis stays pinned to the edge as you zoom and scroll,
so it is always in view.

Beneath the coloured bands a **numeric age axis** is drawn — a ruler in *millions of
years before present* (Ma), with tick marks and labels at round intervals that
increase toward the root, so you can read any node's age directly off the axis.
(In the circular layout the named, coloured annuli themselves are the age scale.)

The two bands **adapt to the tree's depth** so both always fully cover its range:
**System/Period** over **Series/Epoch** for a Phanerozoic tree, **Erathem/Era**
over **System/Period** once the tree reaches into the Proterozoic, and
**Eonothem/Eon** over **Erathem/Era** for a deep Archean tree — so even a
billions-of-years "tree of life" is fully banded (the Precambrian is never blank).

The tree is anchored in time by its **root age**: Archaeopteryx uses the oldest
`<date>` value in the tree, or you can set it explicitly with **"Set root age…"**
next to the Time Axis selector.

The axis is aligned to the tree's **own branches**, so it works on a **fossil-only**
clade — one with no living (age-0) tip. Each Ma spans exactly one branch-length
unit, anchored at the root age, so the coloured bands line up with the branch
nodes and the youngest tips sit at their true age (e.g. an all-extinct ammonite
tree ends at the end-Cretaceous, 66 Ma, rather than being stretched to the
present). A tree that *does* reach the present is the ordinary special case.

Two optional refinements (both **off by default**, and — like the axis itself —
per tree, saved with the tree) are available under **Settings → Overlays** when the
geologic axis is on: **Time-Axis Grid Lines** draws faint reference lines across the tree at the
finer band's interval boundaries (e.g. the Early/Middle/Late Triassic and the
Triassic/Permian boundaries), and **Geologic Boundary Ages** labels the coarser
band's boundaries with their age (e.g. *201.4* at the base of the Jurassic).

The interval names, boundaries, and official colours are those of the
**International Chronostratigraphic Chart** of the **International Commission on
Stratigraphy (ICS / IUGS)** ([stratigraphy.org](https://stratigraphy.org)).
Reference:

- Cohen, K.M., Harper, D.A.T., Gibbard, P.L. & Car, N. (2025, updated): "The ICS
  International Chronostratigraphic Chart this decade", *Episodes* 48:105–115.


Fossil Range Bars (FAD/LAD)
---------------------------

For a tree of **fossil taxa**, a tip is not a single point in time — each taxon is
known from a *stratigraphic range*, from its **First Appearance Datum** (FAD, its
oldest occurrence) to its **Last Appearance Datum** (LAD, its youngest). Turn on
**Settings → Overlays → Data Overlays → Fossil Range Bars (FAD/LAD)** and, on a
dated phylogram, each tip that carries a `<date>` **min/max** gets a capped
stratigraphic-range bar spanning its known duration, drawn back over the terminal
branch so the tip label stays clear. Read against the **Geologic Time Axis**, this
turns a phylogeny into a proper stratigraphic-range figure — the duration and
overlap of taxa laid out against the named geologic intervals, no hand-drawing in
Illustrator required.

Like the Node Age Bars, it is **auto-enabled** on load when the tree has fossil tip
ranges, and it renders in every rectangular orientation and as radial segments in
the **circular** layout. The range is read from the tip's native phyloXML `<date>`
(value/min/max — the same model shown in the node popup), so it works directly on a
tree time-scaled by any of the usual tools. Reference:

- Bell, M.A. & Lloyd, G.T. (2015): "strap: an R package for plotting phylogenies
  against stratigraphy and assessing their stratigraphic congruence",
  *Palaeontology* 58(2):379–389.


Calendar (Absolute-Date) Axis
-----------------------------

For a **tip-dated** tree — a time-scaled phylogeny whose branch lengths are
calendar time and whose tips carry sampling dates, as in molecular epidemiology
(e.g. a SARS-CoV-2 phylodynamic tree) — Archaeopteryx can draw a **calendar-year
axis** instead of the geologic one. Turn it on under **Settings → Overlays →
Time Axis → Calendar (dates)**. It draws a labelled year/decade ruler (like the
distance scale axis, but in calendar time): along the bottom in the root-left
orientation, down the side in root-on-top / root-on-bottom, and as concentric
labelled **year rings** in the circular layout.

The axis is anchored by the **most-recent tip** (the present): Archaeopteryx uses
the largest tip `<date>` (a calendar-year value) automatically, or you can set it
explicitly with **"Set most-recent-tip date…"** next to the Time Axis selector.
Each node's calendar date is then its distance-from-root back from that present.
The **Time-Axis Grid Lines** toggle (Settings → Overlays) also works here, drawing
faint reference lines across the tree at each labelled year tick.


Tanglegrams
-----------

A **tanglegram** compares two trees side by side with connectors linking their
matching tips — the standard way to see how congruent two phylogenies are (a gene
tree vs. a species tree, host vs. parasite, two reconstruction methods, …). With
two or more trees open, choose **Analysis → Create Tanglegram…**, pick the two
trees and the field to link their tips on (node name, taxonomy, or sequence), and
the tanglegram opens in its own window (the second tree mirrored). It reports the
number of **crossing** connectors and a size-normalised **entanglement** score
(described below) — measures of how incongruent the two trees are.

To untangle it, **click a clade's vertical bar** to flip it (a topology-preserving
rotation), or press **Auto-untangle** to have Archaeopteryx do it for both trees;
both are undoable. The **Colour** selector recolours the connectors — uniform,
**Crossings** (the crossing connectors highlighted in red), or by a tip attribute
(taxonomy or an imported category, with a legend). **Export…** saves the figure
as **PDF**, **SVG**, **EPS**, or **PNG** (the vector formats are publication-ready,
document-white with the colouring preserved).

**Auto-untangle heuristic.** Auto-untangle reorders the two trees by rotating
clades (reversing the child order at internal nodes — a topology-preserving
operation) to reduce the number of crossing connectors between their matched tips.
Archaeopteryx uses a barycentre heuristic — each clade's children are ordered by
the mean vertical position of the tips they link to in the other tree — applied
alternately to both trees and iterated to convergence, combined with random
restarts, keeping the arrangement with the fewest crossings (it never increases
the crossing count). Minimising the crossings of a tanglegram is NP-hard, so this
is a heuristic. Foundational references:

- Barycentre crossing-reduction heuristic: Sugiyama K, Tagawa S, Toda M (1981):
  "Methods for visual understanding of hierarchical system structures", *IEEE
  Transactions on Systems, Man, and Cybernetics* 11(2):109–125.
- Tanglegram layout and crossing minimisation: Scornavacca C, Zickmann F, Huson
  DH (2011): "Tanglegrams for rooted phylogenetic trees and networks",
  *Bioinformatics* 27(13):i248–i256.

**Entanglement (congruence score).** Alongside the raw crossing count, the window
reports a size-normalised **entanglement** in the range [0, 1]: **0** when the two
leaf orderings agree perfectly (no crossings) and **1** when they are fully
reversed. Two connectors cross exactly when their matched tips sit in the opposite
vertical order in the two trees, so the total number of crossings is the number of
**discordant** tip pairs — an inversion count, which is the Kendall rank-correlation
(τ) distance between the two leaf orderings. Archaeopteryx divides that count by the
maximum possible number of pairs, *n*(*n*−1)/2, so the score is comparable across
trees of different sizes, and computes it in *O*(*n* log *n*) with a counting merge
sort. This is a crossing-based (Kendall-τ) entanglement; it is a different, simpler
definition from the leaf-position-based *entanglement* of the dendextend package.
References:

- Rank-correlation (concordant/discordant pairs) basis: Kendall MG (1938): "A New
  Measure of Rank Correlation", *Biometrika* 30(1–2):81–93.
- Related, leaf-position-based entanglement: Galili T (2015): "dendextend: an R
  package for visualizing, adjusting and comparing trees of hierarchical
  clustering", *Bioinformatics* 31(22):3718–3720.


For Developers
--------------

### Clone

```
git clone https://github.com/cmzmasek/forester.git
cd forester
```

### Prerequisites

- **JDK 21 or newer** — the build targets `release 21`; an older JDK fails with
  *"release version 21 not supported"*. (Building the installers needs a JDK, not
  just a JRE, because it uses the JDK's `jpackage`.)
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

It prints `Failed tests: 0` on success. A handful of GUI integration tests need a
display (and FlatLaf), so they are kept out of the headless suite (and out of the
shipped jar) and run individually, for example:

```
java -cp forester/java/classes org.forester.archaeopteryx.SubSuperTreeButtonsTest
```

New or changed code should come with tests — see the existing `*Test` classes for
the established patterns.

### Packaging and Continuous Integration

The native installers are produced with the JDK's **`jpackage`**, which bundles a
trimmed Java runtime with the app.

#### Building an installer locally

`jpackage` only builds an installer for the **OS it runs on**, so on macOS you can
build the mac artifacts locally, from `forester/java`:

```
ant jpackage-dmg          # -> dist/Archaeopteryx-<version>.dmg (via a self-contained .app)
ant jpackage-app-image    # -> dist/Archaeopteryx.app only (no .dmg)
```

Output goes to `forester/java/dist/` (git-ignored). The Ant targets:

- parse the version from `AptxConstants.java` (the single source of truth);
- **map the version** for the bundle: macOS and Windows reject a leading-zero
  major, so the bundle version strips a leading `0.` (`0.9.82` → `9.82`, still
  correctly ordered across the `0.9.x` / `0.10.x` line). The real, user-visible
  version is restored into the app's `Info.plist` (so Finder shows `0.9.82`) and
  used for the `.dmg` / asset filenames;
- embed the app icon (see `archaeopteryx_icon_assets/`).

#### The CI workflow

`.github/workflows/installers.yml` builds all three platforms' installers on
GitHub-hosted runners.

- **Triggers:** pushing a **version tag** (`*.*.*`, e.g. `0.9.83`), or a manual
  run from the **Actions** tab (`workflow_dispatch`).
- **Build matrix** (`fail-fast: false`, so each OS reports independently):

  | Runner | Output | How |
  |--------|--------|-----|
  | `macos-latest` (Apple Silicon) | `.dmg` | reuses `ant jpackage-dmg` (the exact local target) |
  | `windows-latest` | `.msi` | `ant all` + `jpackage --type msi`; a step puts **WiX 3** on `PATH` (jpackage needs it for `.msi`); the msi is renamed to the real version; Start-menu/desktop shortcuts + install-dir chooser |
  | `ubuntu-latest` | `.deb` | `ant all` + `jpackage --type deb` (installs `fakeroot`); PNG icon |

  Each leg uploads its installer as a **run artifact**.
- **Release job** (`release`): runs **only on a tag push** (`needs: build`),
  downloads the three artifacts, and creates/updates a **GitHub Release** for the
  tag with them attached as assets. It is idempotent (`gh release upload
  --clobber`), marked **prerelease** while the version is `0.x`, and writes notes
  listing each file plus the unsigned-app caveat. Only this job is granted
  `contents: write`; the build matrix stays read-only.
- **App icons** live in `archaeopteryx_icon_assets/` — `.icns` (macOS), `.ico`
  (Windows), `.png` (Linux); see that folder's README for regenerating the
  `.icns` from the source SVG.

Notes / current limitations:

- Nothing is **code-signed or notarized** yet, hence the first-launch prompts
  documented above.
- There is **no Intel-Mac (`x86_64`) dmg**: GitHub has retired its Intel
  (`macos-13`) hosted runners, and `jpackage` cannot cleanly cross-build an
  x86_64 app on an arm64 runner. Intel users run the [jar](#run-from-the-jar).

#### Cutting a release

```
# 1. bump VERSION in forester/java/src/org/forester/archaeopteryx/AptxConstants.java, then commit
# 2. tag and push:
git tag 0.9.83 && git push origin 0.9.83
```

CI builds the macOS `.dmg`, Windows `.msi`, and Linux `.deb`, and publishes them
to a GitHub Release for the tag.
