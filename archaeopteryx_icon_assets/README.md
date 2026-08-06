# Archaeopteryx app icon assets

Branding assets for the native installers (built via the `jpackage-*` targets in
`forester/java/build.xml`).

| File | Use |
|------|-----|
| `archaeopteryx-anime-badge.svg` | **source** for the macOS icon (rounded-tile "badge" style) |
| `Archaeopteryx.icns` | **macOS** app icon — wired into `jpackage-app-image` via `--icon` |
| `archaeopteryx-anime.ico` | **Windows** app icon (for the future Windows `jpackage` build) |
| `archaeopteryx-anime.svg` / `-512.png` | the plain (transparent, no tile) bird — not used by the mac build |
| `archaeopteryx-badge-512.png` | raster preview of the badge tile |

## Regenerating `Archaeopteryx.icns`

Rendered crisply from the badge SVG (needs `librsvg` — `brew install librsvg` — plus the
built-in `iconutil`/`sips`):

```sh
SVG=archaeopteryx-anime-badge.svg; SET=Archaeopteryx.iconset
mkdir -p "$SET"
r() { rsvg-convert -w "$1" -h "$1" "$SVG" -o "$SET/$2"; }
r 16 icon_16x16.png;    r 32 icon_16x16@2x.png
r 32 icon_32x32.png;    r 64 icon_32x32@2x.png
r 128 icon_128x128.png; r 256 icon_128x128@2x.png
r 256 icon_256x256.png; r 512 icon_256x256@2x.png
r 512 icon_512x512.png; r 1024 icon_512x512@2x.png
iconutil -c icns "$SET" -o Archaeopteryx.icns && rm -rf "$SET"
```

The committed `Archaeopteryx.icns` lets the build run with only a JDK (no `librsvg` needed).
