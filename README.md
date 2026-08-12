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
