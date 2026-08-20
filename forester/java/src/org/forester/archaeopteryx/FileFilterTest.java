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

import java.io.File;

import javax.swing.filechooser.FileFilter;

/**
 * The Read-Tree file-chooser filters accept Auspice/Nextstrain JSON: the dedicated {@code JsonFilter} matches
 * {@code .json} / {@code .auspice.json} (case-insensitively, plus directories), and the "All supported files" default
 * filter now also lists {@code .json} so a JSON dataset is visible/selectable in the open dialog. Headless.
 */
public final class FileFilterTest {

    public static void main( final String[] args ) {
        System.out.println( "File filters: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        final FileFilter json = new JsonFilter();
        // the dedicated JSON filter accepts .json and .auspice.json, case-insensitively
        if ( !json.accept( new File( "ncov.json" ) ) || !json.accept( new File( "ncov.auspice.json" ) )
                || !json.accept( new File( "NCOV.JSON" ) ) ) {
            return fail( "JsonFilter must accept .json / .auspice.json (case-insensitive)" );
        }
        // ... but not a non-JSON tree file, and it accepts directories (so the user can navigate)
        if ( json.accept( new File( "tree.xml" ) ) || json.accept( new File( "tree.nwk" ) ) ) {
            return fail( "JsonFilter must not accept non-JSON files" );
        }
        if ( !json.accept( new File( System.getProperty( "user.dir" ) ) ) ) {
            return fail( "JsonFilter must accept directories" );
        }
        if ( json.getDescription() == null || !json.getDescription().toLowerCase().contains( "json" ) ) {
            return fail( "JsonFilter description must mention JSON" );
        }
        // the "All supported files" default filter now lists .json too (so it shows in the open dialog)...
        final FileFilter def = new DefaultFilter();
        if ( !def.accept( new File( "ncov.json" ) ) || !def.accept( new File( "ncov.auspice.json" ) ) ) {
            return fail( "the default 'All supported' filter must accept .json" );
        }
        // ... without regressing the other supported types, and still rejecting an unsupported one
        if ( !def.accept( new File( "tree.xml" ) ) || !def.accept( new File( "tree.nwk" ) )
                || !def.accept( new File( "tree.nex" ) ) ) {
            return fail( "the default filter must still accept the existing tree types" );
        }
        if ( def.accept( new File( "photo.png" ) ) ) {
            return fail( "the default filter must still reject an unsupported type" );
        }
        return true;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [FileFilterTest] " + msg );
        return false;
    }

    private FileFilterTest() {
    }
}
