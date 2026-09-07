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
import java.awt.image.BufferedImage;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.TreeFacts.Group;
import org.forester.phylogeny.Phylogeny;

/**
 * Tests for {@link TreePropertiesForm} (needs a display): the sections it builds (editable ones first, then one
 * per fact group, the identity section collapsed only when empty), dirtiness against the normalized baseline,
 * validation outlines, a write without a tree panel (no undo, still applied), the fact refresh after a write, the
 * header texts, and {@link TreePropertiesForm#rebind()} keeping unwritten edits while reloading clean fields
 * from a replaced tree. Also the {@link HistogramPanel}: bin hit-testing, tooltips, and that it paints.
 */
public final class TreePropertiesFormTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TreePropertiesForm: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                sections( ok );
                editing( ok );
                rebind( ok );
                histogramPanel( ok );
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void sections( final boolean[] ok ) {
        final Phylogeny phy = TreeFactsTest.fixture();
        phy.setName( "Fixture" );
        final TreePropertiesForm form = new TreePropertiesForm( phy, null );
        check( ok, "name section", form.hasSectionForTest( TreePropertiesForm.SEC_NAME ) );
        check( ok, "name section expanded", form.isSectionExpandedForTest( TreePropertiesForm.SEC_NAME ) );
        check( ok, "identity section", form.hasSectionForTest( TreePropertiesForm.SEC_IDENTITY ) );
        check( ok, "identity collapsed when empty", !form.isSectionExpandedForTest( TreePropertiesForm.SEC_IDENTITY ) );
        for( final String t : new String[] { TreeFacts.FILE, TreeFacts.STRUCTURE, TreeFacts.BRANCH_LENGTHS,
                "Support values (bootstrap)", "Support values (probability)", TreeFacts.COVERAGE } ) {
            check( ok, "fact section: " + t, form.hasSectionForTest( t ) && form.isSectionExpandedForTest( t ) );
        }
        check( ok, "no time axis without a panel", !form.hasSectionForTest( TreeFacts.TIME_AXIS ) );
        check( ok, "groups exposed", form.groups().size() == 6 );
        check( ok, "title is the name", "Fixture".equals( form.headerTitleForTest() ) );
        check( ok, "subtitle without a file: " + form.headerSubtitleForTest(),
               "not saved to a file yet · 5 tips".equals( form.headerSubtitleForTest() ) );
        form.toggleSectionForTest( TreeFacts.STRUCTURE );
        check( ok, "toggle collapses", !form.isSectionExpandedForTest( TreeFacts.STRUCTURE ) );
        // a refresh rebuilds the fact sections but keeps each one's expanded state
        form.refresh();
        check( ok, "collapsed state survives a refresh", !form.isSectionExpandedForTest( TreeFacts.STRUCTURE ) );
        check( ok, "other sections still expanded", form.isSectionExpandedForTest( TreeFacts.COVERAGE ) );
        // identity expanded when the tree has any of its fields
        phy.setType( "gene tree" );
        final TreePropertiesForm typed = new TreePropertiesForm( phy, null );
        check( ok, "identity expanded with data", typed.isSectionExpandedForTest( TreePropertiesForm.SEC_IDENTITY ) );
        // an unnamed tree
        phy.setName( "" );
        check( ok, "untitled", TreePropertiesForm.UNTITLED.equals( new TreePropertiesForm( phy, null ).titleText() ) );
        // every editable field is registered
        for( final String key : new String[] { TreePropertiesDraft.NAME, TreePropertiesDraft.DESCRIPTION,
                TreePropertiesDraft.ID_VALUE, TreePropertiesDraft.ID_PROVIDER, TreePropertiesDraft.TYPE,
                TreePropertiesDraft.DISTANCE_UNIT } ) {
            check( ok, "field " + key, form.fieldForTest( key ) != null );
        }
        check( ok, "editable", form.isEditable() && ( form.component() == form ) );
    }

    private static void editing( final boolean[] ok ) {
        final Phylogeny phy = TreeFactsTest.fixture();
        phy.setName( "Fixture" );
        final TreePropertiesForm form = new TreePropertiesForm( phy, null );
        final int[] changes = { 0 };
        form.addChangeListener( () -> changes[ 0 ]++ );
        check( ok, "clean on open", !form.isDirty() && form.problems().isEmpty() );
        // whitespace-only edits are not dirty
        form.setTextForTest( TreePropertiesDraft.NAME, "  Fixture " );
        check( ok, "listener fired", changes[ 0 ] > 0 );
        check( ok, "whitespace is not a change", !form.isDirty() );
        form.setTextForTest( TreePropertiesDraft.NAME, "" );
        check( ok, "blank name is dirty and invalid", form.isDirty() && !form.problems().isEmpty() );
        check( ok, "blank name outlined", form.isOutlinedForTest( TreePropertiesDraft.NAME ) );
        check( ok, "write refused", !form.write() && "Fixture".equals( phy.getName() ) );
        form.setTextForTest( TreePropertiesDraft.NAME, "Renamed" );
        check( ok, "outline cleared", !form.isOutlinedForTest( TreePropertiesDraft.NAME ) );
        form.setTextForTest( TreePropertiesDraft.ID_PROVIDER, "ncbi" );
        check( ok, "provider without value outlined", form.isOutlinedForTest( TreePropertiesDraft.ID_PROVIDER ) );
        form.setTextForTest( TreePropertiesDraft.ID_VALUE, "123" );
        form.setTextForTest( TreePropertiesDraft.DESCRIPTION, "desc" );
        form.setTextForTest( TreePropertiesDraft.DISTANCE_UNIT, "Ma" );
        check( ok, "valid again", form.problems().isEmpty() );
        final List<Group> before = form.groups();
        check( ok, "write without a panel", form.write() );
        check( ok, "applied", "Renamed".equals( phy.getName() ) && "desc".equals( phy.getDescription() )
                && ( phy.getIdentifier() != null ) && "ncbi".equals( phy.getIdentifier().getProvider() )
                && "Ma".equals( phy.getDistanceUnit() ) );
        check( ok, "clean after write", !form.isDirty() );
        check( ok, "baseline updated", "Renamed".equals( form.baseline().name ) );
        check( ok, "header follows the name", "Renamed".equals( form.headerTitleForTest() ) );
        check( ok, "facts recomputed on write", form.groups() != before );
        // a no-change write is a harmless true
        check( ok, "no-op write", form.write() );
    }

    private static void rebind( final boolean[] ok ) {
        final Phylogeny phy = TreeFactsTest.fixture();
        phy.setName( "Fixture" );
        // a form over a panel-less tree: rebind re-reads the SAME tree object (no panel to swap it)
        final TreePropertiesForm clean = new TreePropertiesForm( phy, null );
        phy.setName( "Changed elsewhere" );
        phy.setType( "species tree" );
        clean.rebind();
        check( ok, "clean form reloads the name", "Changed elsewhere".equals( clean.collect().name ) );
        check( ok, "clean form reloads the type", "species tree".equals( clean.collect().type ) );
        check( ok, "still clean", !clean.isDirty() );
        check( ok, "header re-read", "Changed elsewhere".equals( clean.headerTitleForTest() ) );
        // a DIRTY form keeps its edits and measures them against the new baseline
        final TreePropertiesForm dirty = new TreePropertiesForm( phy, null );
        dirty.setTextForTest( TreePropertiesDraft.DESCRIPTION, "my unwritten text" );
        phy.setName( "Changed again" );
        dirty.rebind();
        check( ok, "dirty edit kept", "my unwritten text".equals( dirty.collect().description ) );
        check( ok, "still dirty", dirty.isDirty() );
        check( ok, "baseline is the new tree", "Changed again".equals( dirty.baseline().name ) );
        check( ok, "unedited field NOT reloaded while dirty (the widgets are the user's)",
               "Changed elsewhere".equals( dirty.collect().name ) );
        check( ok, "writing after a rebind applies to the tree", dirty.write()
                && "my unwritten text".equals( phy.getDescription() ) && "Changed elsewhere".equals( phy.getName() ) );
    }

    private static void histogramPanel( final boolean[] ok ) {
        final org.forester.util.DescriptiveStatistics st = new org.forester.util.BasicDescriptiveStatistics();
        for( int i = 0; i <= 12; ++i ) {
            st.addValue( i );
        }
        final HistogramPanel hp = new HistogramPanel( TreeFacts.histogram( st ) );
        final JFrame f = new JFrame();
        f.getContentPane().add( hp );
        f.setSize( 400, 120 );
        f.validate(); // lays out the panel at the frame's width (not shown: no OS window)
        hp.setSize( 300, hp.getPreferredSize().height );
        check( ok, "preferred height covers bars + labels", hp.getPreferredSize().height > HistogramPanel.BAR_AREA_HEIGHT );
        check( ok, "max width capped", hp.getMaximumSize().width == HistogramPanel.MAX_WIDTH );
        check( ok, "first bin at x=0", hp.binAt( 0 ) == 0 );
        check( ok, "last bin near the right edge", hp.binAt( 299 ) == TreeFacts.HISTOGRAM_BINS - 1 );
        check( ok, "outside", ( hp.binAt( -1 ) == -1 ) && ( hp.binAt( 5000 ) == -1 ) );
        final java.awt.event.MouseEvent me = new java.awt.event.MouseEvent( hp, 0, 0, 0, 2, 5, 0, false );
        final String tip = hp.getToolTipText( me );
        check( ok, "tooltip names the range and count: " + tip, "0 – 1: 1 branch".equals( tip ) );
        final java.awt.event.MouseEvent last = new java.awt.event.MouseEvent( hp, 0, 0, 0, 298, 5, 0, false );
        check( ok, "last bin tooltip: " + hp.getToolTipText( last ), "11 – 12: 2 branches".equals( hp.getToolTipText( last ) ) );
        // it paints: accent pixels in the bar area, nothing thrown
        final BufferedImage img = new BufferedImage( 300, hp.getHeight(), BufferedImage.TYPE_INT_ARGB );
        final java.awt.Graphics2D g = img.createGraphics();
        hp.paint( g );
        g.dispose();
        final int accent = FormWidgets.accentColor().getRGB() & 0xFFFFFF;
        boolean found = false;
        for( int x = 0; ( x < 300 ) && !found; ++x ) {
            for( int y = 0; ( y < HistogramPanel.BAR_AREA_HEIGHT ) && !found; ++y ) {
                if ( ( img.getRGB( x, y ) & 0xFFFFFF ) == accent ) {
                    found = true;
                }
            }
        }
        check( ok, "accent-coloured bars painted", found );
        f.dispose();
    }

    private static void check( final boolean[] ok, final String what, final boolean condition ) {
        if ( !condition ) {
            System.out.println( "  [TreePropertiesFormTest] " + what );
            ok[ 0 ] = false;
        }
    }

    private TreePropertiesFormTest() {
    }
}
