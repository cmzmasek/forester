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
import java.awt.GraphicsEnvironment;
import java.nio.file.Files;
import java.util.EnumMap;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.DirectoryPreferences.Category;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Verifies the per-purpose last-used directories in the live app: OPEN / SAVE / EXPORT are tracked <em>separately</em>
 * (setting one does not move the others), an unset purpose falls back to a readable home/Desktop directory, and the
 * directories round-trip through {@link DirectoryPreferences} using the real resolution. Uses an isolated state
 * directory so the user's real {@code ~/.archaeopteryx} is never read or written. Headful; a green no-op when headless.
 */
public final class PersistentDirectoriesTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "PersistentDirectories: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final String prev = System.getProperty( GuiPreferences.DIR_PROPERTY );
        final MainFrame[] mf = new MainFrame[ 1 ];
        final boolean[] ok = { true };
        try {
            // isolate ALL Archaeopteryx on-disk state to a throwaway directory for this test
            final File state = Files.createTempDirectory( "aptx-state" ).toFile();
            System.setProperty( GuiPreferences.DIR_PROPERTY, state.getAbsolutePath() );

            final File open_dir = Files.createTempDirectory( "aptx-open" ).toFile();
            final File export_dir = Files.createTempDirectory( "aptx-export" ).toFile();

            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( "((a,b),c);", new NHXParser() )[ 0 ];
            final Configuration conf = new Configuration();
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "dirs" ) );
            final MainFrame frame = mf[ 0 ];

            SwingUtilities.invokeAndWait( () -> {
                // set only OPEN and EXPORT; leave SAVE unset
                frame.setCurrentDir( Category.OPEN, open_dir );
                frame.setCurrentDir( Category.EXPORT, export_dir );

                // (1) separation: each purpose returns its own directory, and SAVE (unset) is neither of them
                if ( !samePath( frame.getCurrentDir( Category.OPEN ), open_dir ) ) {
                    fail( ok, "OPEN should return the open dir, got " + frame.getCurrentDir( Category.OPEN ) );
                }
                if ( !samePath( frame.getCurrentDir( Category.EXPORT ), export_dir ) ) {
                    fail( ok, "EXPORT should return the export dir, got " + frame.getCurrentDir( Category.EXPORT ) );
                }
                final File save = frame.getCurrentDir( Category.SAVE );
                if ( samePath( save, open_dir ) || samePath( save, export_dir ) ) {
                    fail( ok, "SAVE must not inherit OPEN/EXPORT -- the purposes are separate; got " + save );
                }
                // (2) an unset purpose falls back to a readable directory (never null/unreadable)
                if ( ( save == null ) || !save.canRead() ) {
                    fail( ok, "an unset purpose must fall back to a readable directory, got " + save );
                }

                // (3) end-to-end: the live map persists through DirectoryPreferences (isolated store) and reloads
                new DirectoryPreferences().saveFrom( frame._current_dirs );
                final Map<Category, File> reloaded = new EnumMap<>( Category.class );
                new DirectoryPreferences().applyTo( reloaded );
                if ( !samePath( reloaded.get( Category.OPEN ), open_dir ) ) {
                    fail( ok, "OPEN should persist+reload, got " + reloaded.get( Category.OPEN ) );
                }
                if ( !samePath( reloaded.get( Category.EXPORT ), export_dir ) ) {
                    fail( ok, "EXPORT should persist+reload, got " + reloaded.get( Category.EXPORT ) );
                }
                if ( reloaded.containsKey( Category.SAVE ) ) {
                    fail( ok, "an unset SAVE (a fallback, never stored) must not be persisted" );
                }

                // (4) lazy validation: a remembered directory that no longer exists is NOT returned -- getCurrentDir
                //     falls back (this check moved here from load-time so a stale network path can't hang startup)
                try {
                    final File doomed = Files.createTempDirectory( "aptx-doomed" ).toFile();
                    frame.setCurrentDir( Category.SAVE, doomed );
                    if ( !samePath( frame.getCurrentDir( Category.SAVE ), doomed ) ) {
                        fail( ok, "a freshly set, existing dir should be returned" );
                    }
                    if ( !doomed.delete() ) {
                        fail( ok, "could not delete the temp dir for the lazy-validation check" );
                    }
                    final File after = frame.getCurrentDir( Category.SAVE );
                    if ( samePath( after, doomed ) ) {
                        fail( ok, "getCurrentDir must not return a directory that no longer exists; got " + after );
                    }
                    if ( ( after == null ) || !after.canRead() ) {
                        fail( ok, "after a vanished dir, getCurrentDir must fall back to a readable dir, got " + after );
                    }
                }
                catch ( final Exception e ) {
                    fail( ok, "lazy-validation check threw: " + e );
                }
            } );
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            fail( ok, "unexpected: " + t );
        }
        finally {
            if ( mf[ 0 ] != null ) {
                try {
                    SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() );
                }
                catch ( final Exception ignore ) {
                    // ignore teardown failure
                }
            }
            if ( prev == null ) {
                System.clearProperty( GuiPreferences.DIR_PROPERTY );
            }
            else {
                System.setProperty( GuiPreferences.DIR_PROPERTY, prev );
            }
        }
        return ok[ 0 ];
    }

    private static boolean samePath( final File a, final File b ) {
        return ( a != null ) && ( b != null ) && a.getAbsolutePath().equals( b.getAbsolutePath() );
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [PersistentDirectoriesTest] " + msg );
        ok[ 0 ] = false;
    }

    private PersistentDirectoriesTest() {
    }
}
