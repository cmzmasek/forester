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

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * A parsed external association table for a tanglegram: a two-column CSV/TSV mapping a value in the FIRST (left) tree
 * to the value it corresponds to in the SECOND (right) tree, so two trees whose linking tips have DIFFERENT names can
 * still be paired (the parasite-vs-host / gene-vs-species-with-different-ids case a plain value join cannot handle).
 * See {@link TanglegramLinker#linkByAssociation}.
 *
 * <p>Format (deliberately minimal): one association per line, two columns separated by a TAB or a comma (auto-detected
 * per line), no header row, and a leading UTF-8 byte-order mark tolerated. Blank lines and lines with fewer than two
 * non-empty cells are skipped. A left value may map to several right values (repeat it on several lines), and several
 * left values may map to one right value -- the mapping is many:many. Cells are trimmed and a surrounding pair of
 * double quotes is stripped. This is NOT a full CSV parser: a value that itself contains the delimiter (a quoted comma)
 * is not supported -- use the TAB-delimited form for values that may contain commas (tip identifiers rarely do).
 */
final class TanglegramAssociation {

    private final Map<String, List<String>> _left_to_right;
    private final int                       _pair_count;

    private TanglegramAssociation( final Map<String, List<String>> left_to_right, final int pair_count ) {
        _left_to_right = left_to_right;
        _pair_count = pair_count;
    }

    /** left value -> the right value(s) it is associated with (insertion order, deduplicated per left key). */
    Map<String, List<String>> leftToRight() {
        return _left_to_right;
    }

    /** The number of distinct (left, right) association pairs read. */
    int pairCount() {
        return _pair_count;
    }

    /** The number of distinct left keys. */
    int leftKeyCount() {
        return _left_to_right.size();
    }

    /**
     * Parses a two-column CSV/TSV association table. Never returns null; a table with no usable rows yields an empty
     * association (which links nothing). Robust to blank lines, a trailing newline, CR/LF line endings, and stray
     * short rows.
     */
    static TanglegramAssociation parse( final String content ) {
        final Map<String, List<String>> map = new LinkedHashMap<>();
        int pairs = 0;
        if ( content != null ) {
            // strip a leading UTF-8 byte-order mark (Excel / Google Sheets CSV exports prepend one), else the first
            // left key would carry a hidden U+FEFF and silently fail to match its tip
            final String text = ( !content.isEmpty() && ( content.charAt( 0 ) == '\uFEFF' ) ) ? content.substring( 1 )
                    : content;
            for( final String raw : text.split( "\r\n|\r|\n", -1 ) ) {
                if ( raw.trim().isEmpty() ) {
                    continue;
                }
                final char delimiter = ( raw.indexOf( '\t' ) >= 0 ) ? '\t' : ',';
                final String[] cells = raw.split( String.valueOf( delimiter ), -1 );
                if ( cells.length < 2 ) {
                    continue;
                }
                final String left = clean( cells[ 0 ] );
                final String right = clean( cells[ 1 ] );
                if ( left.isEmpty() || right.isEmpty() ) {
                    continue;
                }
                final List<String> rights = map.computeIfAbsent( left, k -> new ArrayList<>() );
                if ( !rights.contains( right ) ) {
                    rights.add( right );
                    pairs++;
                }
            }
        }
        return new TanglegramAssociation( map, pairs );
    }

    /** Trims a cell and strips a single pair of surrounding double quotes (as a spreadsheet CSV export may add). */
    private static String clean( final String cell ) {
        String s = cell.trim();
        if ( ( s.length() >= 2 ) && ( s.charAt( 0 ) == '"' ) && ( s.charAt( s.length() - 1 ) == '"' ) ) {
            s = s.substring( 1, s.length() - 1 ).trim();
        }
        return s;
    }
}
