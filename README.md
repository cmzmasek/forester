forester
========

**forester** is an open-source Java library and command-line toolkit for
phylogenetics and evolutionary genomics: reading and writing phylogenetic trees,
gene-tree / species-tree reconciliation, ortholog inference, protein-domain and
sequence analysis, and the data structures behind them. It is the reference
implementation of the **phyloXML** format.

Its best-known component is **Archaeopteryx**, an interactive viewer and editor
for phylogenetic trees.

> ### 🌳 Looking for Archaeopteryx?
>
> Archaeopteryx — the interactive tree viewer and editor — has its own home:
>
> - **Home page & downloads → <https://cmzmasek.github.io/archaeopteryx/>**
> - **Documentation & releases → <https://github.com/cmzmasek/archaeopteryx>**
>
> Native installers (macOS `.dmg`, Windows `.msi`, Linux `.deb`) bundle a Java
> runtime, so there is nothing else to install. The Archaeopteryx **source code**
> lives in *this* repository (under
> [`forester/java/src/org/forester/archaeopteryx/`](forester/java/src/org/forester/archaeopteryx));
> the links above are its user-facing front door and full documentation.


What's in forester
------------------

- **Archaeopteryx** (`org.forester.archaeopteryx`) — the interactive tree viewer
  and editor (see the links above).
- **phyloXML & tree I/O** (`org.forester.io`, `org.forester.phylogeny`) — the
  reference implementation of the [phyloXML](http://www.phyloxml.org) format
  (reader, writer, object model), plus Newick / New Hampshire (NH/NHX) and Nexus
  parsers, and the phylogeny data structures they build.
- **SDI / GSDI / GSDIR** (`org.forester.sdi`) — speciation–duplication inference:
  reconcile a gene tree against a species tree to infer gene duplications and
  speciations.
- **RIO** (`org.forester.rio`) — Resampled Inference of Orthologs.
- **cladinator** (`org.forester.clade_analysis`) — classify / place query
  sequences against annotated reference clades.
- **surfacing** (`org.forester.surfacing`) — genome-wide protein-domain-architecture
  analysis and comparison.
- **MSA tools** (`org.forester.msa`, `org.forester.msa_compactor`) —
  multiple-sequence-alignment reading/writing, entropy, consensus, and compaction.
- **evoinference** (`org.forester.evoinference`) — distance-matrix methods
  (e.g. neighbor-joining).
- Supporting packages: sequence, species, protein, Gene Ontology (`go`),
  sequence/taxonomy web-service clients (`ws`), and general utilities.

**Command-line tools** live in `org.forester.application` — among them `gsdi`,
`rio`, `cladinator`, `decorator`, `count_support`, `confadd`, `nj`,
`phyloxml_converter`, and the `msa_*` alignment utilities. Each is a `main`-class
you can run from the built jar, e.g.:

```
java -cp forester.jar org.forester.application.gsdi
```

(run with no arguments for usage). Building the jar is described under
[For Developers](#for-developers).


Citing
------

If you use **phyloXML** or the forester libraries, please cite:

- Han, M.V. & Zmasek, C.M. (2009): "phyloXML: XML for evolutionary biology and
  comparative genomics", *BMC Bioinformatics* 10:356.

For **Archaeopteryx**, see the citation information at
<https://github.com/cmzmasek/archaeopteryx> (a dedicated publication is in
preparation; a Zenodo DOI provides a stable, versioned citation).


License
-------

GPL-3.0 — free and open source. See [`LICENSE`](LICENSE). © Christian M. Zmasek.


For Developers
--------------

### Clone

```
git clone https://github.com/cmzmasek/forester.git
cd forester
```

### Prerequisites

- **JDK 21 or newer** — the build targets `release 21`; an older JDK fails with
  *"release version 21 not supported"*. (Building the Archaeopteryx installers
  needs a JDK, not just a JRE, because it uses the JDK's `jpackage`.)
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
`org.forester.archaeopteryx.Archaeopteryx` (so `java -jar forester.jar` launches
Archaeopteryx; the command-line tools are run via `-cp`, as shown above).

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

### Packaging the Archaeopteryx installers, and Continuous Integration

The native Archaeopteryx installers are produced with the JDK's **`jpackage`**,
which bundles a trimmed Java runtime with the app.

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
  documented on the Archaeopteryx site.
- There is **no Intel-Mac (`x86_64`) dmg**: GitHub has retired its Intel
  (`macos-13`) hosted runners, and `jpackage` cannot cleanly cross-build an
  x86_64 app on an arm64 runner. Intel users run the plain `forester.jar` with
  their own Java 21+.

#### Cutting an Archaeopteryx release

```
# 1. bump VERSION in forester/java/src/org/forester/archaeopteryx/AptxConstants.java, then commit
# 2. tag and push:
git tag 0.9.83 && git push origin 0.9.83
```

CI builds the macOS `.dmg`, Windows `.msi`, and Linux `.deb`, and publishes them
to a GitHub Release for the tag.

### Publishing the Archaeopteryx home page (cmzmasek.github.io/archaeopteryx)

The home page at **<https://cmzmasek.github.io/archaeopteryx/>** is **not** in
this `forester` repository — it lives in a small, separate repository,
**[`cmzmasek/archaeopteryx`](https://github.com/cmzmasek/archaeopteryx)**, and is
served by **GitHub Pages** (classic project site) from that repo's **`main`
branch, `docs/` folder**. The whole site is a single self-contained
`docs/index.html`; the Archaeopteryx **user documentation** lives in that repo's
`README.md`.

To update the site:

```
git clone https://github.com/cmzmasek/archaeopteryx.git
cd archaeopteryx
# edit docs/index.html
git add docs/index.html
git commit -m "Site: <what changed>"
git push origin main
```

Pushing to `main` triggers a Pages rebuild automatically; the change is live at
<https://cmzmasek.github.io/archaeopteryx/> within about a minute. There is no
build step — GitHub Pages serves `docs/index.html` as-is.

> **URL case-sensitivity:** GitHub project-site paths are case-sensitive, so the
> home URL is the lowercase `…/archaeopteryx/`. Keep every `github.com/cmzmasek`
> repo reference and the Pages URL lowercase to match.

The site's **Download** button points at that repo's own
[releases](https://github.com/cmzmasek/archaeopteryx/releases/latest). Those
installers are produced by a **`release.yml`** workflow **in the
`cmzmasek/archaeopteryx` repo** (run manually from its **Actions** tab via
`workflow_dispatch`): it checks out `forester` at a given ref, builds the three
installers with `jpackage` (the same recipe as this repo's `installers.yml`), and
publishes them as a GitHub Release there. So a normal `forester` release (tag →
`installers.yml`) is independent of the home page; to also refresh the site's
Download button, dispatch that `release.yml` (pass `prerelease=false` so the
"latest" link updates).
