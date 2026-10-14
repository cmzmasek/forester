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

import java.util.List;

/** A report a {@link ReportWindow} can show: a summary, column names and rows of cells, and its own TSV. */
interface TableReport {

    /** The lines above the table. */
    String getSummary();

    String[] getColumns();

    /** One array per row, one cell per column. */
    List<String[]> getCells();

    /** Tab-separated, one header line then one line per row; a cell's tabs and line breaks become spaces. */
    String toTsv();

    static String tsvCell( final String s ) {
        return ( s == null ) ? "" : s.replace( '\t', ' ' ).replace( '\r', ' ' ).replace( '\n', ' ' );
    }

    /** The standard TSV of {@code columns} + {@code cells}. */
    static String tsvOf( final String[] columns, final List<String[]> cells ) {
        final StringBuilder sb = new StringBuilder( String.join( "\t", columns ) ).append( '\n' );
        for( final String[] row : cells ) {
            for( int i = 0; i < row.length; ++i ) {
                sb.append( ( i == 0 ) ? "" : "\t" ).append( tsvCell( row[ i ] ) );
            }
            sb.append( '\n' );
        }
        return sb.toString();
    }
}
