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

import java.awt.GraphicsEnvironment;
import java.io.File;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.TIME_AXIS_TYPE;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * The Time Axis is PER-TREE: the Dinosaur demo (geologic) and the SARS-CoV-2 demo (calendar) auto-derive DIFFERENT
 * axes in two tabs at once, an override on one tab does not affect the other, a deviating config round-trips through a
 * phyloXML save/reload (the {@code aptx:time_axis} property), and a non-deviating (auto-derive) config writes NO
 * property. Headful; a green no-op when headless. Dogfoods the two demos.
 */
public final class TimeAxisPerTreeTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TimeAxisPerTree: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny dino = parse( "dinosaur-time-tree.xml" );
            final Phylogeny sars = parse( "sars-cov-2-time-tree.xml" );
            if ( ( dino == null ) || ( sars == null ) ) {
                return fail( "demo trees missing" );
            }
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { dino, sars }, conf, "ta" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final java.util.List<TreePanel> panels = frame.getMainPanel().getTreePanels();
                    if ( panels.size() != 2 ) {
                        fail( ok, "expected two tabs, got " + panels.size() );
                        return;
                    }
                    final TreePanel tp_dino = panels.get( 0 );
                    final TreePanel tp_sars = panels.get( 1 );

                    // 1) auto-derive is per-tree and INDEPENDENT: geologic vs calendar at the same time
                    if ( tp_dino.effectiveTimeAxisType() != TIME_AXIS_TYPE.GEOLOGIC ) {
                        fail( ok, "the Dinosaur tab must auto-derive GEOLOGIC, got " + tp_dino.effectiveTimeAxisType() );
                    }
                    if ( tp_sars.effectiveTimeAxisType() != TIME_AXIS_TYPE.CALENDAR ) {
                        fail( ok, "the SARS-CoV-2 tab must auto-derive CALENDAR, got " + tp_sars.effectiveTimeAxisType() );
                    }

                    // 2) an override on one tab does NOT affect the other
                    tp_dino.setTimeAxisType( TIME_AXIS_TYPE.NONE );
                    tp_dino.setTimeAxisGrid( true );
                    if ( ( tp_dino.effectiveTimeAxisType() != TIME_AXIS_TYPE.NONE ) || !tp_dino.isTimeAxisGrid() ) {
                        fail( ok, "the Dinosaur override did not take effect" );
                    }
                    if ( ( tp_sars.effectiveTimeAxisType() != TIME_AXIS_TYPE.CALENDAR ) || tp_sars.isTimeAxisGrid() ) {
                        fail( ok, "overriding the Dinosaur tab must NOT affect the SARS tab (still CALENDAR, grid off)" );
                    }

                    // 3) a DEVIATING config round-trips through a phyloXML save/reload
                    tp_sars.setTimeAxisType( TIME_AXIS_TYPE.NONE ); // deviate from the derived CALENDAR (turn axis off)
                    tp_sars.setTimeAxisAges( true );
                    final Phylogeny reloaded = saveAndReread( tp_sars );
                    if ( reloaded == null ) {
                        fail( ok, "save/reread failed" );
                    }
                    else {
                        final TimeAxisConfig cfg = TimeAxisConfig.readFromTree( reloaded );
                        if ( ( cfg == null ) || ( cfg.getType() != TIME_AXIS_TYPE.NONE ) || !cfg.isAges() ) {
                            fail( ok, "a deviating Time-Axis config must survive save/reload, got "
                                    + ( cfg == null ? "null" : cfg.serialize() ) );
                        }
                        // and the reloaded tree, opened in a new tab, restores that config (auto-derive is overridden)
                        frame.getMainPanel().addPhylogenyInNewTab( reloaded, conf, "reloaded", "" );
                        final TreePanel tp_re = frame.getMainPanel().getTreePanels()
                                .get( frame.getMainPanel().getTreePanels().size() - 1 );
                        if ( ( tp_re.effectiveTimeAxisType() != TIME_AXIS_TYPE.NONE ) || !tp_re.isTimeAxisAges() ) {
                            fail( ok, "a reloaded tree must restore its saved axis override (NONE + ages), got "
                                    + tp_re.effectiveTimeAxisType() + " ages=" + tp_re.isTimeAxisAges() );
                        }
                    }

                    // 4) a GRID-ONLY deviation must NOT freeze the type: it round-trips as {type=auto, grid} and the
                    //    reloaded tree STILL auto-derives its type (a fresh dino tab: derives GEOLOGIC, type override null)
                    final Phylogeny dino2 = parse( "dinosaur-time-tree.xml" );
                    frame.getMainPanel().addPhylogenyInNewTab( dino2, conf, "dino2", "" );
                    final TreePanel tp_d2 = last( frame );
                    tp_d2.setTimeAxisGrid( true ); // deviate via grid ONLY -- the type stays on auto-derive
                    final Phylogeny reloaded2 = saveAndReread( tp_d2 );
                    final TimeAxisConfig cfg2 = ( reloaded2 == null ) ? null : TimeAxisConfig.readFromTree( reloaded2 );
                    if ( ( cfg2 == null ) || ( cfg2.getType() != null ) || !cfg2.isGrid() ) {
                        fail( ok, "a grid-only deviation must persist as {type=auto (null), grid=true}, got "
                                + ( cfg2 == null ? "null" : cfg2.serialize() ) );
                    }
                    else {
                        frame.getMainPanel().addPhylogenyInNewTab( reloaded2, conf, "reloaded2", "" );
                        final TreePanel tp_r2 = last( frame );
                        if ( ( tp_r2.effectiveTimeAxisType() != TIME_AXIS_TYPE.GEOLOGIC ) || !tp_r2.isTimeAxisGrid() ) {
                            fail( ok, "a reloaded grid-only tree must STILL auto-derive its type (GEOLOGIC) + keep grid, "
                                    + "got " + tp_r2.effectiveTimeAxisType() + " grid=" + tp_r2.isTimeAxisGrid() );
                        }
                    }

                    // 5) a NON-deviating (auto-derive) config writes NO property (clean files)
                    tp_dino.resetTimeAxisToAutoDerive(); // back to the derived GEOLOGIC, no flags/overrides
                    tp_dino.syncTimeAxisConfigToTree();
                    if ( TimeAxisConfig.readFromTree( tp_dino.getPhylogeny() ) != null ) {
                        fail( ok, "an auto-derive (non-deviating) config must NOT write an aptx:time_axis property" );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Sync the panel's Time-Axis config into its tree, write the tree to a temp phyloXML file, re-read it. */
    private static Phylogeny saveAndReread( final TreePanel tp ) throws Exception {
        tp.syncTimeAxisConfigToTree();
        final File tmp = File.createTempFile( "aptx_timeaxis", ".xml" );
        tmp.deleteOnExit();
        try {
            new PhylogenyWriter().toPhyloXML( tmp, tp.getPhylogeny(), 0 );
            final Phylogeny[] read = ParserBasedPhylogenyFactory.getInstance().create( tmp,
                    PhyloXmlParser.createPhyloXmlParser() );
            return ( read.length > 0 ) ? read[ 0 ] : null;
        }
        finally {
            tmp.delete();
        }
    }

    private static Phylogeny parse( final String demo ) {
        final File f = new File( System.getProperty( "user.dir" ), "forester/demo/" + demo );
        if ( !f.exists() ) {
            return null;
        }
        try {
            return ParserBasedPhylogenyFactory.getInstance().create( f, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    private static TreePanel last( final MainFrame frame ) {
        final java.util.List<TreePanel> ps = frame.getMainPanel().getTreePanels();
        return ps.get( ps.size() - 1 );
    }

    private static boolean fail( final String m ) {
        System.out.println( "  TimeAxisPerTreeTest: " + m );
        return false;
    }

    private static void fail( final boolean[] ok, final String m ) {
        System.out.println( "  [TimeAxisPerTreeTest] " + m );
        ok[ 0 ] = false;
    }
}
