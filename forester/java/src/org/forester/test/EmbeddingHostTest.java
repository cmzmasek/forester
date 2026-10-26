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

package org.forester.test;

import java.awt.GraphicsEnvironment;
import java.awt.Window;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Stream;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Configuration;
import org.forester.archaeopteryx.MainFrameApplication;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/** Hosts an embedded Archaeopteryx viewer the way another Swing program would, through public members only. */
public final class EmbeddingHostTest {

    private static final String CACHE_DIR_PROPERTY = "archaeopteryx.cache.dir";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "EmbeddingHost: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final String prev_dir = System.getProperty( CACHE_DIR_PROPERTY );
        Path dir = null;
        try {
            dir = Files.createTempDirectory( "aptx-host" );
            System.setProperty( CACHE_DIR_PROPERTY, dir.toString() );
            return hostOk();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        finally {
            if ( prev_dir == null ) {
                System.clearProperty( CACHE_DIR_PROPERTY );
            }
            else {
                System.setProperty( CACHE_DIR_PROPERTY, prev_dir );
            }
            deleteQuietly( dir );
        }
    }

    private static boolean hostOk() throws Exception {
        final Set<Window> before = new HashSet<>( Arrays.asList( Window.getWindows() ) );
        final MainFrameApplication viewer = MainFrameApplication
                .createEmbeddedInstance( new Phylogeny[] { tree( "((A,B),(C,D),(E,F),(G,H));" ),
                        tree( "(((P,Q),R),((S,T),(U,(V,W))));" ) }, new Configuration(), "host" );
        final JFrame[] host = { null };
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> {
                host[ 0 ] = new JFrame( "Host" );
                host[ 0 ].setJMenuBar( viewer.getJMenuBar() );
                final JPanel panel = viewer; // compiles only while the viewer is a JPanel
                host[ 0 ].getContentPane().add( panel );
                host[ 0 ].setSize( 900, 700 );
                host[ 0 ].setVisible( true );
            } );
            SwingUtilities.invokeAndWait( () -> {} ); // the first layout's resize events
            SwingUtilities.invokeAndWait( () -> {} );
            SwingUtilities.invokeAndWait( () -> {
                final JMenu file = menu( viewer.getJMenuBar(), "File" );
                if ( ( file == null ) || hasItem( file, "Exit" ) ) {
                    ok[ 0 ] = fail( "an embedded viewer's File menu must have no Exit: the host closes the viewer" );
                }
                viewer.end();
            } );
            SwingUtilities.invokeAndWait( () -> {} );
            SwingUtilities.invokeAndWait( () -> {
                if ( viewer.getComponentCount() != 0 ) {
                    ok[ 0 ] = fail( "end() must empty the viewer" );
                }
                if ( !host[ 0 ].isDisplayable() ) {
                    ok[ 0 ] = fail( "end() must leave the host's window open" );
                }
                for( final Window w : Window.getWindows() ) {
                    if ( !before.contains( w ) && ( w != host[ 0 ] ) && w.isDisplayable() ) {
                        ok[ 0 ] = fail( "end() must leave no window of the viewer open: " + w );
                        w.dispose();
                    }
                }
                host[ 0 ].dispose();
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> {
                if ( ( host[ 0 ] != null ) && host[ 0 ].isDisplayable() ) {
                    viewer.end();
                    host[ 0 ].dispose();
                }
            } );
        }
        return ok[ 0 ];
    }

    private static Phylogeny tree( final String newick ) throws Exception {
        return ParserBasedPhylogenyFactory.getInstance().create( newick, new NHXParser() )[ 0 ];
    }

    private static JMenu menu( final JMenuBar bar, final String title ) {
        for( int i = 0; i < bar.getMenuCount(); ++i ) {
            final JMenu m = bar.getMenu( i );
            if ( ( m != null ) && title.equals( m.getText() ) ) {
                return m;
            }
        }
        return null;
    }

    private static boolean hasItem( final JMenu menu, final String text ) {
        for( int i = 0; i < menu.getItemCount(); ++i ) {
            if ( ( menu.getItem( i ) != null ) && text.equals( menu.getItem( i ).getText() ) ) {
                return true;
            }
        }
        return false;
    }

    private static void deleteQuietly( final Path dir ) {
        if ( dir == null ) {
            return;
        }
        try (Stream<Path> paths = Files.walk( dir )) {
            paths.sorted( Comparator.reverseOrder() ).forEach( p -> p.toFile().delete() );
        }
        catch ( final IOException e ) {
            // best effort: a leftover temp directory is harmless
        }
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [EmbeddingHostTest] " + message );
        return false;
    }
}
