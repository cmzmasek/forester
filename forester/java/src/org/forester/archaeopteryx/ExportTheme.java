// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.archaeopteryx;

import java.awt.Color;

/**
 * Transiently switches a tree panel to the LIGHT theme -- color scheme AND panel background, exactly as an
 * on-screen UI theme switch does (see {@link MainFrame#updateTreeCanvasColors}) -- so an export or clipboard
 * copy is document-ready (white background, dark legible ink) regardless of the light/dark theme the user is
 * working in, then restores it. Shared by the raster render ({@link AptxUtil#renderPhylogenyToImage}) and the
 * vector exporters ({@link VectorGraphicsExporter}, {@link PdfExporter}) so all export formats behave alike.
 *
 * <p>Use in a try/finally around the paint:
 * <pre>
 *   final ExportTheme theme = ExportTheme.applyIf( tree_panel, whiteBackground );
 *   try { ...paintPhylogeny... } finally { theme.restore(); }
 * </pre>
 * {@link #restore()} is a no-op when the option is off or the panel is already in the light theme. All calls run
 * synchronously on the EDT during export, so no repaint observes the transient state.
 */
final class ExportTheme {

    private final TreePanel _tree_panel;
    private final boolean   _relit;
    private final int       _prev_scheme;
    private final Color     _prev_background;

    private ExportTheme( final TreePanel tree_panel, final boolean relit, final int prev_scheme,
                         final Color prev_background ) {
        _tree_panel = tree_panel;
        _relit = relit;
        _prev_scheme = prev_scheme;
        _prev_background = prev_background;
    }

    /** Switches {@code tree_panel} to the light theme when {@code white_background} is on and it is not already
     *  light; returns a handle whose {@link #restore()} undoes exactly that. */
    static ExportTheme applyIf( final TreePanel tree_panel, final boolean white_background ) {
        if ( !white_background ) {
            return new ExportTheme( tree_panel, false, 0, null ); // no-op handle; restore() does nothing
        }
        final TreeColorSet color_set = tree_panel.getTreeColorSet();
        final int prev_scheme = color_set.getCurrentColorScheme();
        if ( prev_scheme == TreeColorSet.LIGHT_COLOR_SCHEME ) {
            return new ExportTheme( tree_panel, false, prev_scheme, null ); // already light -> nothing to switch
        }
        final Color prev_background = tree_panel.getBackground();
        color_set.setColorSchema( TreeColorSet.LIGHT_COLOR_SCHEME );
        tree_panel.setBackground( color_set.getBackgroundColor() );
        return new ExportTheme( tree_panel, true, prev_scheme, prev_background );
    }

    void restore() {
        if ( _relit ) {
            _tree_panel.getTreeColorSet().setColorSchema( _prev_scheme );
            _tree_panel.setBackground( _prev_background );
        }
    }
}
