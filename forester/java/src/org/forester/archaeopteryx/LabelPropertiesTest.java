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
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * Node properties in the tip label: the one-line "values only" text, the field inventory the Annotation Fields
 * chooser offers, and the per-tree selection + order that drives them.
 * <p>
 * The behaviour under test is a deliberate replacement of the old one, which appended EVERY property to the label
 * as "{@code <full ref>: <value>}", newline-joined into a string drawn as a single line -- so a ten-property tree
 * got an unreadable label with embedded newlines, in backwards ref order. The assertions below pin each of those
 * four things (values only, comma-joined, one line, ascending order) so none of them can quietly come back.
 */
public final class LabelPropertiesTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LabelProperties: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return testLabelText() && testPropertyListOrder() && testFieldInventory() && testFieldRoles()
                && testLabelSpecIsNotAColumn() && testFieldReorder() && testLabelInTreePanel();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [LabelPropertiesTest] " + msg );
        return false;
    }

    // ---- the one-line label text -------------------------------------------------------------------------------
    private static boolean testLabelText() {
        final PropertiesList props = new PropertiesList();
        props.addProperty( new Property( "data:host", "cat", "", "xsd:string", AppliesTo.NODE ) );
        props.addProperty( new Property( "data:country", "Brazil", "", "xsd:string", AppliesTo.NODE ) );
        props.addProperty( new Property( "data:reads", "42", "METRIC:kb", "xsd:decimal", AppliesTo.NODE ) );
        props.addProperty( new Property( "data:blank", "", "", "xsd:string", AppliesTo.NODE ) );
        props.addProperty( new Property( "data:spaces", "   ", "", "xsd:string", AppliesTo.NODE ) );
        props.addProperty( new Property( "aptx:import_profile", "v1;/data/x.csv", "", "xsd:string", AppliesTo.NODE ) );
        props.addProperty( new Property( "style:node_color", "0x112233", "", "xsd:string", AppliesTo.PHYLOGENY ) );
        // default (no explicit selection): every user-visible property, values only, ascending ref order, one line
        final String all = TreePanelUtil.labelPropertiesText( props, null );
        if ( !"Brazil, cat, 42 kb".equals( all ) ) {
            return fail( "default label text should be \"Brazil, cat, 42 kb\", got [" + all + "]" );
        }
        if ( all.indexOf( '\n' ) >= 0 ) {
            return fail( "the label is drawn as ONE line -- it must never contain a newline: [" + all + "]" );
        }
        if ( all.contains( ":" ) || all.contains( "data" ) || all.contains( "METRIC" ) ) {
            return fail( "the label shows values only -- no \"ref:\" prefix, and an unprefixed unit: [" + all + "]" );
        }
        if ( all.contains( "v1;" ) || all.contains( "0x112233" ) ) {
            return fail( "internal aptx:/style: metadata must never reach the label: [" + all + "]" );
        }
        // an explicit selection picks the fields AND their order (the reverse of the ref-sorted order, to prove
        // the chosen order wins rather than coinciding with it)
        final String chosen = TreePanelUtil.labelPropertiesText( props, Arrays.asList( "data:reads", "data:host" ) );
        if ( !"42 kb, cat".equals( chosen ) ) {
            return fail( "the chosen field order should win: expected \"42 kb, cat\", got [" + chosen + "]" );
        }
        // a ref the node does not carry is skipped, not rendered as a gap
        final String missing = TreePanelUtil.labelPropertiesText( props,
                                                                 Arrays.asList( "data:nope", "data:host" ) );
        if ( !"cat".equals( missing ) ) {
            return fail( "a ref the node lacks should be skipped: [" + missing + "]" );
        }
        // naming an internal ref explicitly still must not show it
        final String internal = TreePanelUtil.labelPropertiesText( props,
                                                                   Arrays.asList( "aptx:import_profile" ) );
        if ( internal.length() != 0 ) {
            return fail( "an explicitly named aptx: ref must still be hidden: [" + internal + "]" );
        }
        if ( TreePanelUtil.labelPropertiesText( props, new ArrayList<String>() ).length() != 0 ) {
            return fail( "an EMPTY selection means no properties in the label" );
        }
        if ( TreePanelUtil.labelPropertiesText( null, null ).length() != 0 ) {
            return fail( "a null property list should yield empty text" );
        }
        // a whitespace-only value contributes no text, so it must not emit a separator either
        if ( all.contains( ", ," ) || all.endsWith( "," ) || all.endsWith( " " ) ) {
            return fail( "a whitespace-only value must not leave a doubled or trailing separator: [" + all + "]" );
        }
        // the rollover popup / node panel text hides the same internal namespaces the label does
        final String popup = TreePanelUtil.userVisiblePropertiesText( props ).toString();
        if ( popup.contains( "0x112233" ) || popup.contains( "style:" ) || popup.contains( "aptx:" ) ) {
            return fail( "internal aptx:/style: metadata must not surface in the node-data display: " + popup );
        }
        if ( !TreePanelUtil.isVisualStylePropertyRef( "style:node_color" )
                || TreePanelUtil.isVisualStylePropertyRef( "data:host" )
                || TreePanelUtil.isVisualStylePropertyRef( null ) ) {
            return fail( "isVisualStylePropertyRef must match only the style: namespace" );
        }
        return true;
    }

    // ---- the property list sorts ASCENDING by ref ---------------------------------------------------------------
    private static boolean testPropertyListOrder() {
        final PropertiesList props = new PropertiesList();
        props.addProperty( new Property( "data:zebra", "z", "", "xsd:string", AppliesTo.NODE ) );
        props.addProperty( new Property( "data:aardvark", "a", "", "xsd:string", AppliesTo.NODE ) );
        if ( !"data:aardvark".equals( props.getProperties().get( 0 ).getRef() ) ) {
            return fail( "properties must sort ASCENDING by ref (the comparator's arguments used to be swapped), got "
                    + props.getProperties().get( 0 ).getRef() + " first" );
        }
        return true;
    }

    // ---- the chooser's field inventory is BROADER than the colorable refs ----------------------------------------
    private static boolean testFieldInventory() {
        final Phylogeny phy = inventoryTree();
        final List<String> refs = TreePanelUtil.userVisiblePropertyRefs( phy );
        // the three kinds of field colorableRefs drops -- constant, per-tip-unique categorical, internal-only --
        // are exactly the ones worth putting in a label, so the inventory must keep all three
        for( final String must : new String[] { "data:host", "data:study", "data:accession", "data:support" } ) {
            if ( !refs.contains( must ) ) {
                return fail( "the field inventory is missing " + must + ": " + refs );
            }
        }
        if ( refs.contains( "aptx:import_profile" ) || refs.contains( "style:node_color" ) ) {
            return fail( "internal aptx:/style: metadata must not be offered as a field: " + refs );
        }
        // ... and this is WHY the inventory cannot just be colorableRefs: that list drops three of the four
        final List<String> colorable = PropertyColorScheme.colorableRefs( phy );
        if ( !colorable.contains( "data:host" ) ) {
            return fail( "a varying categorical field should be colorable: " + colorable );
        }
        if ( colorable.contains( "data:study" ) || colorable.contains( "data:accession" )
                || colorable.contains( "data:support" ) ) {
            return fail( "colorableRefs is expected to drop the constant / per-tip-unique / internal-only fields; "
                    + "if it no longer does, the inventory rationale needs revisiting: " + colorable );
        }
        // sorted by the name the user sees: accession, host, study, support
        if ( !refs.equals( Arrays.asList( "data:accession", "data:host", "data:study", "data:support" ) ) ) {
            return fail( "the inventory should be sorted by display name, got " + refs );
        }
        return true;
    }

    // ---- which roles each kind of field is offered ---------------------------------------------------------------
    private static boolean testFieldRoles() {
        final Phylogeny phy = inventoryTree();
        // an INTERNAL-node-only field: a tip-aligned column has nothing to draw, so the label is the only role
        final List<AnnotationColumns.Type> support = AnnotationColumns.allowedTypes( phy, "data:support" );
        if ( !support.equals( Arrays.asList( AnnotationColumns.Type.LABEL ) ) ) {
            return fail( "an internal-node-only field should offer the label only, got " + support );
        }
        // a per-tip-UNIQUE categorical field (an accession): no useful colour, but a fine text column
        final List<AnnotationColumns.Type> acc = AnnotationColumns.allowedTypes( phy, "data:accession" );
        if ( !acc.equals( Arrays.asList( AnnotationColumns.Type.LABEL, AnnotationColumns.Type.TEXT ) ) ) {
            return fail( "a per-tip-unique field should offer label + text only, got " + acc );
        }
        // a CONSTANT field: likewise -- one colour for every tip tells the reader nothing
        final List<AnnotationColumns.Type> study = AnnotationColumns.allowedTypes( phy, "data:study" );
        if ( !study.equals( Arrays.asList( AnnotationColumns.Type.LABEL, AnnotationColumns.Type.TEXT ) ) ) {
            return fail( "a constant field should offer label + text only, got " + study );
        }
        // a normal categorical field keeps its full set of column roles, with the label added
        final List<AnnotationColumns.Type> host = AnnotationColumns.allowedTypes( phy, "data:host" );
        if ( !host.contains( AnnotationColumns.Type.COLOR_STRIP ) || !host.contains( AnnotationColumns.Type.SYMBOL )
                || !host.contains( AnnotationColumns.Type.LABEL ) ) {
            return fail( "a colorable categorical field should keep its column roles plus the label, got " + host );
        }
        return true;
    }

    // ---- a LABEL field is never built as a column ---------------------------------------------------------------
    private static boolean testLabelSpecIsNotAColumn() {
        final Phylogeny phy = inventoryTree();
        final AnnotationColumns cols = new AnnotationColumns( phy,
                Arrays.asList( new AnnotationColumns.ColumnSpec( "data:study", AnnotationColumns.Type.LABEL ),
                               new AnnotationColumns.ColumnSpec( "data:host",
                                       AnnotationColumns.Type.COLOR_STRIP ) ) );
        if ( cols.size() != 1 ) {
            return fail( "a LABEL field must not become a column: expected 1 column, got " + cols.size() );
        }
        if ( cols.getColumn( 0 ).getType() != AnnotationColumns.Type.COLOR_STRIP ) {
            return fail( "the surviving column should be the color strip, got " + cols.getColumn( 0 ).getType() );
        }
        return true;
    }

    // ---- the chooser's up/down reordering -----------------------------------------------------------------------
    private static boolean testFieldReorder() {
        final List<String> l = new ArrayList<String>( Arrays.asList( "a", "b", "c" ) );
        if ( MainFrame.moveInList( l, "a", -1 ) || !l.equals( Arrays.asList( "a", "b", "c" ) ) ) {
            return fail( "moving the FIRST item up must be a no-op, got " + l );
        }
        if ( MainFrame.moveInList( l, "c", 1 ) || !l.equals( Arrays.asList( "a", "b", "c" ) ) ) {
            return fail( "moving the LAST item down must be a no-op, got " + l );
        }
        if ( MainFrame.moveInList( l, "zzz", 1 ) || !l.equals( Arrays.asList( "a", "b", "c" ) ) ) {
            return fail( "moving an item that is not in the list must be a no-op, got " + l );
        }
        if ( !MainFrame.moveInList( l, "a", 1 ) || !l.equals( Arrays.asList( "b", "a", "c" ) ) ) {
            return fail( "moving \"a\" down should give [b, a, c], got " + l );
        }
        if ( !MainFrame.moveInList( l, "c", -1 ) || !l.equals( Arrays.asList( "b", "c", "a" ) ) ) {
            return fail( "moving \"c\" up should give [b, c, a], got " + l );
        }
        return true;
    }

    // ---- the real label path on a real TreePanel -----------------------------------------------------------------
    private static boolean testLabelInTreePanel() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration part; needs a display toolkit
        }
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { inventoryTree() }, new Configuration(), "labelprops" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final PhylogenyNode tip = tp.getPhylogeny().getFirstExternalNode();
                mf[ 0 ].getControlPanel().setCheckbox( DisplayOption.SHOW_PROPERTIES, true );
                // the default: every user-visible property, values only
                final String all = tp.nodeLabelTextForTest( tip );
                if ( !all.contains( "cat" ) || !all.contains( "ACC1" ) || !all.contains( "S2026" ) ) {
                    ok[ 0 ] = fail( "by default the label should show every property value, got [" + all + "]" );
                }
                else if ( all.contains( "data:" ) || ( all.indexOf( '\n' ) >= 0 ) ) {
                    ok[ 0 ] = fail( "the label must be one line of values, got [" + all + "]" );
                }
                // narrowed to one field, and only that field
                tp.setLabelPropertyRefs( Arrays.asList( "data:accession" ) );
                final String one = tp.nodeLabelTextForTest( tip );
                if ( !one.contains( "ACC1" ) || one.contains( "cat" ) || one.contains( "S2026" ) ) {
                    ok[ 0 ] = fail( "only the chosen field should show, got [" + one + "]" );
                }
                // two fields, in the order given (accession AFTER host, the reverse of the ref-sorted order)
                tp.setLabelPropertyRefs( Arrays.asList( "data:host", "data:accession" ) );
                final String two = tp.nodeLabelTextForTest( tip );
                if ( two.indexOf( "cat" ) > two.indexOf( "ACC1" ) ) {
                    ok[ 0 ] = fail( "the label should honour the chosen field order, got [" + two + "]" );
                }
                // none: the node keeps its name, with no dangling separator
                tp.setLabelPropertyRefs( new ArrayList<String>() );
                final String none = tp.nodeLabelTextForTest( tip );
                if ( none.contains( "cat" ) || none.contains( "ACC1" ) ) {
                    ok[ 0 ] = fail( "an empty selection should put no property in the label, got [" + none + "]" );
                }
                else if ( none.endsWith( " " ) || none.startsWith( " " ) ) {
                    ok[ 0 ] = fail( "an empty selection must not leave a dangling space, got [" + none + "]" );
                }
                else if ( !none.contains( tip.getName() ) ) {
                    ok[ 0 ] = fail( "the node's own name should survive, got [" + none + "]" );
                }
                // Reset to Defaults path: back to showing all of them
                tp.clearAnnotationColumns();
                if ( tp.getLabelPropertyRefs() != null ) {
                    ok[ 0 ] = fail( "a reset should clear the label-field selection" );
                }
                else if ( !tp.nodeLabelTextForTest( tip ).contains( "cat" ) ) {
                    ok[ 0 ] = fail( "after a reset the label should show every property again" );
                }
                // A field drawn as a COLUMN leaves the default label set: one role per field has to hold from the
                // moment a tree is opened, not only once the chooser has been visited.
                tp.setAnnotationColumns( java.util.Arrays
                        .asList( new AnnotationColumns.ColumnSpec( "data:host",
                                                                   AnnotationColumns.Type.COLOR_STRIP ) ) );
                final String with_col = tp.nodeLabelTextForTest( tip );
                if ( with_col.contains( "cat" ) ) {
                    ok[ 0 ] = fail( "a field shown as a column must not ALSO be in the default label set, got ["
                            + with_col + "]" );
                }
                else if ( !with_col.contains( "ACC1" ) || !with_col.contains( "S2026" ) ) {
                    ok[ 0 ] = fail( "the other fields should still be in the label, got [" + with_col + "]" );
                }
                tp.setAnnotationColumns( null );
                // the cached longest-label width must be RECOMPUTED when the label text changes, not zeroed --
                // it is what reserves room for the labels, the columns, the clade bands and the radial radius
                tp.calcParametersForPainting( 900, 620 );
                final int w_all = tp.getLongestExtNodeInfo();
                tp.setLabelPropertyRefs( new ArrayList<String>() );
                final int w_none = tp.getLongestExtNodeInfo();
                if ( ( w_none <= 0 ) || ( w_none >= w_all ) ) {
                    ok[ 0 ] = fail( "dropping the label properties should SHRINK the reserved label width to a "
                            + "positive value, got " + w_none + " (was " + w_all + ")" );
                }
                tp.setLabelPropertyRefs( null );
                if ( tp.getLongestExtNodeInfo() <= w_none ) {
                    ok[ 0 ] = fail( "restoring the label properties should grow the reserved width again, got "
                            + tp.getLongestExtNodeInfo() );
                }
                // PARITY: the radial layouts build their labels in their own painter, which used to stop at the
                // sequence data -- so the "Properties" option and everything the chooser does to the label were
                // silently dead in Circular and Unrooted. These tips are named "A"/"B"/"C", far short of the radial
                // label-width cap, so the extra ink can only be the property values.
                for ( final Options.PHYLOGENY_GRAPHICS_TYPE type : new Options.PHYLOGENY_GRAPHICS_TYPE[] {
                        Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                    tp.setPhylogenyGraphicsType( type );
                    mf[ 0 ].getControlPanel().setCheckbox( DisplayOption.SHOW_PROPERTIES, false );
                    final int without = inkOffscreen( mf[ 0 ], tp );
                    mf[ 0 ].getControlPanel().setCheckbox( DisplayOption.SHOW_PROPERTIES, true );
                    final int with = inkOffscreen( mf[ 0 ], tp );
                    if ( with <= without ) {
                        ok[ 0 ] = fail( "the " + type + " layout must draw the label properties too (ink " + with
                                + " vs " + without + ")" );
                    }
                }
            } );
            SwingUtilities.invokeAndWait( () -> mf[ 0 ].dispose() );
            return ok[ 0 ];
        }
        catch ( final Exception e ) {
            return fail( "exception: " + e );
        }
    }

    /** Paints the panel into an offscreen image and counts the non-white pixels. */
    private static int inkOffscreen( final MainFrame mf, final TreePanel tp ) {
        final int w = 900, h = 620;
        mf.showWhole();
        tp.setSize( w, h );
        tp.validate();
        tp.doLayout();
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.printAll( g );
        g.dispose();
        int n = 0;
        for( int x = 0; x < w; ++x ) {
            for( int y = 0; y < h; ++y ) {
                if ( ( img.getRGB( x, y ) & 0x00FFFFFF ) != 0x00FFFFFF ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /**
     * A tree carrying one field of each kind the chooser has to cope with: {@code data:host} varies across the tips
     * (colorable), {@code data:study} is the same on every tip (constant), {@code data:accession} differs on every
     * tip (per-tip-unique categorical), and {@code data:support} sits on an INTERNAL node only. Plus the two
     * internal namespaces, which must never be offered or drawn.
     */
    private static Phylogeny inventoryTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        root.setName( "root" );
        final PropertiesList rp = new PropertiesList();
        rp.addProperty( new Property( "data:support", "0.98", "", "xsd:decimal", AppliesTo.NODE ) );
        rp.addProperty( new Property( "aptx:import_profile", "v1;/data/x.csv", "", "xsd:string", AppliesTo.NODE ) );
        rp.addProperty( new Property( "style:node_color", "0x112233", "", "xsd:string", AppliesTo.PHYLOGENY ) );
        root.getNodeData().setProperties( rp );
        root.addAsChild( tip( "A", "cat", "ACC1" ) );
        root.addAsChild( tip( "B", "dog", "ACC2" ) );
        root.addAsChild( tip( "C", "cat", "ACC3" ) );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode tip( final String name, final String host, final String accession ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( "data:host", host, "", "xsd:string", AppliesTo.NODE ) );
        pl.addProperty( new Property( "data:accession", accession, "", "xsd:string", AppliesTo.NODE ) );
        pl.addProperty( new Property( "data:study", "S2026", "", "xsd:string", AppliesTo.NODE ) );
        n.getNodeData().setProperties( pl );
        return n;
    }
}
