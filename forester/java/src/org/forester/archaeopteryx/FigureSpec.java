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

import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterUtil;

/**
 * A tree's FIGURE: what is drawn, and how, saved with the tree so reopening it gives you back the figure you
 * composed rather than the defaults.
 * <p>
 * Until this existed, only three things travelled with a tree ({@code aptx:branch_event},
 * {@code aptx:import_profile}, {@code aptx:time_axis}). Everything that actually makes a figure -- which labels
 * are drawn, the annotation columns, the clade bands, colour-by, the layout -- was lost the moment the file was
 * saved and reopened. It is also the missing half of "regenerate from a config": a recipe that can be stored can
 * be re-applied, whether by the GUI or by a renderer.
 * <p>
 * <b>Deliberately per-TREE display state only.</b> The theme, fonts and colour scheme stay global preferences, so
 * a figure sent to a colleague looks right in whatever theme they use rather than overriding their setup.
 * <p>
 * <b>Format:</b> {@code v1;key=value;key=value...}, values escaped by {@link PropertyTextCodec}. A map rather
 * than fixed fields, and UNKNOWN KEYS ARE IGNORED on load, so a figure written by a newer Archaeopteryx opens
 * (minus what the older one cannot draw) instead of failing -- the property is not a schema, it is a recipe.
 */
final class FigureSpec {

    /** phyloXML property ref under which a figure is stored, at the PHYLOGENY level. In the {@code aptx:} namespace
     *  so it stays out of the user-facing property displays. */
    static final String FIGURE_REF = "aptx:figure";
    private static final String VERSION = "v1";

    // --- keys -------------------------------------------------------------------------------------------------
    private static final String K_GRAPHICS     = "layout";        // PHYLOGENY_GRAPHICS_TYPE
    private static final String K_ORIENTATION  = "orientation";   // TREE_ORIENTATION
    private static final String K_DISPLAY_TYPE = "displaytype";   // PHYLOGENY_DISPLAY_TYPE (the phylogram tri-state)
    private static final String K_COLOR_BY     = "colorby";
    private static final String K_SIZE_BY      = "sizeby";
    private static final String K_PIE_TRAIT    = "pietrait";
    private static final String K_CLADE_MODE   = "clademode";
    private static final String K_CLADE_SKIP   = "cladeskip";
    private static final String K_CLADE_LEVELS = "cladelevels";   // rank:angle|rank:angle...
    private static final String K_LABEL_PROPS  = "labelprops";    // ref|ref...
    private static final String K_COLUMNS      = "columns";       // ref:type:shape:normalized|...
    private static final String K_SHOW_PREFIX  = "show.";         // + DisplayOption.name()

    private final Map<String, String> _values = new LinkedHashMap<String, String>();

    // --- capture / apply --------------------------------------------------------------------------------------

    /** The figure the panel is currently drawing. */
    static FigureSpec capture( final TreePanel tp ) {
        final FigureSpec f = new FigureSpec();
        if ( tp == null ) {
            return f;
        }
        f.put( K_GRAPHICS, tp.getPhylogenyGraphicsType().name() );
        f.put( K_ORIENTATION, tp.getTreeOrientation().name() );
        // the phylogram / aligned-phylogram / cladogram choice is a TRI-state and already per-tab: it has its own
        // key rather than riding a boolean, which is what silently turned "aligned" into "cladogram" once before
        final Options.PHYLOGENY_DISPLAY_TYPE dt = tp.getTreeDisplayTypeForThisTab();
        if ( dt != null ) {
            f.put( K_DISPLAY_TYPE, dt.name() );
        }
        for( final DisplayOption opt : DisplayOption.values() ) {
            if ( opt != DisplayOption.DISPLAY_AS_PHYLOGRAM ) { // carried by K_DISPLAY_TYPE, above
                f.put( K_SHOW_PREFIX + opt.name(), Boolean.toString( tp.shows( opt ) ) );
            }
        }
        f.put( K_COLOR_BY, tp.getColorByPropertyRef() );
        f.put( K_SIZE_BY, tp.getSizeByPropertyRef() );
        f.put( K_PIE_TRAIT, tp.getAncestralPieTrait() );
        f.put( K_CLADE_MODE, ( tp.getCladeBandsMode() == null ) ? null : tp.getCladeBandsMode().name() );
        f.put( K_CLADE_SKIP, Boolean.toString( tp.isCladeBandsSkipSingletons() ) );
        f.put( K_CLADE_LEVELS, joinCladeLevels( tp.getCladeLevelSpecs() ) );
        f.put( K_LABEL_PROPS, joinList( tp.getLabelPropertyRefs() ) );
        f.put( K_COLUMNS, joinColumns( tp.getAnnotationColumnSpecs() ) );
        return f;
    }

    /**
     * Applies this figure to the panel. Anything the spec does not mention is left alone, and a value that no
     * longer makes sense (a rank the tree no longer has, an enum from a newer version) is skipped rather than
     * thrown -- a figure is a request, not a contract.
     */
    void applyTo( final TreePanel tp ) {
        if ( tp == null ) {
            return;
        }
        final Options.PHYLOGENY_GRAPHICS_TYPE g = enumOf( Options.PHYLOGENY_GRAPHICS_TYPE.class, get( K_GRAPHICS ) );
        if ( g != null ) {
            tp.setPhylogenyGraphicsType( g );
        }
        final Options.TREE_ORIENTATION o = enumOf( Options.TREE_ORIENTATION.class, get( K_ORIENTATION ) );
        if ( o != null ) {
            tp.setTreeOrientation( o );
        }
        final Options.PHYLOGENY_DISPLAY_TYPE dt = enumOf( Options.PHYLOGENY_DISPLAY_TYPE.class,
                                                          get( K_DISPLAY_TYPE ) );
        if ( ( dt != null ) && ( tp.getControlPanel() != null ) ) {
            tp.getControlPanel().setTreeDisplayType( dt );
        }
        for( final DisplayOption opt : DisplayOption.values() ) {
            final String v = get( K_SHOW_PREFIX + opt.name() );
            if ( v != null ) {
                tp.setShows( opt, Boolean.parseBoolean( v ) );
            }
        }
        tp.setColorByPropertyRef( emptyToNull( get( K_COLOR_BY ) ) );
        tp.setSizeByPropertyRef( emptyToNull( get( K_SIZE_BY ) ) );
        tp.setAncestralPieTrait( emptyToNull( get( K_PIE_TRAIT ) ) );
        final List<String> labels = splitList( get( K_LABEL_PROPS ) );
        tp.setLabelPropertyRefs( ( get( K_LABEL_PROPS ) == null ) ? null : labels );
        tp.setAnnotationColumns( parseColumns( get( K_COLUMNS ) ) );
        final List<CladeLevel.Spec> levels = parseCladeLevels( get( K_CLADE_LEVELS ) );
        final TreePanel.CLADE_VIS mode = enumOf( TreePanel.CLADE_VIS.class, get( K_CLADE_MODE ) );
        if ( !levels.isEmpty() && ( mode != null ) ) {
            tp.setCladeLevels( levels, mode, Boolean.parseBoolean( get( K_CLADE_SKIP ) ) );
        }
        else {
            tp.clearCladeBands();
        }
    }

    /** The figure a tree gets when everything is switched off -- what "Clear All Overlays" applies. Only the
     *  overlays: the layout and which labels are drawn are left as they are, because clearing overlays should not
     *  also re-orient the tree or strip its labels. */
    static FigureSpec overlaysOff() {
        final FigureSpec f = new FigureSpec();
        f.put( K_COLOR_BY, "" );
        f.put( K_SIZE_BY, "" );
        f.put( K_PIE_TRAIT, "" );
        f.put( K_CLADE_LEVELS, "" );
        f.put( K_LABEL_PROPS, "" );
        f.put( K_COLUMNS, "" );
        return f;
    }

    boolean isEmpty() {
        return _values.isEmpty();
    }

    // --- serialization ----------------------------------------------------------------------------------------

    String toPropertyValue() {
        final StringBuilder sb = new StringBuilder( VERSION );
        for( final Map.Entry<String, String> e : _values.entrySet() ) {
            sb.append( ';' ).append( e.getKey() ).append( '=' ).append( PropertyTextCodec.esc( e.getValue() ) );
        }
        return sb.toString();
    }

    /** Never throws: a corrupt or unrecognised value yields null (no figure) rather than a failed tree load. */
    static FigureSpec parse( final String s ) {
        if ( ForesterUtil.isEmpty( s ) ) {
            return null;
        }
        final String[] fields = s.split( ";", -1 );
        if ( ( fields.length < 1 ) || !VERSION.equals( fields[ 0 ] ) ) {
            return null; // an unknown version: better no figure than a misread one
        }
        final FigureSpec f = new FigureSpec();
        for( int i = 1; i < fields.length; ++i ) {
            final int eq = fields[ i ].indexOf( '=' );
            if ( eq > 0 ) {
                f._values.put( fields[ i ].substring( 0, eq ),
                               PropertyTextCodec.unesc( fields[ i ].substring( eq + 1 ) ) );
            }
        }
        return f.isEmpty() ? null : f;
    }

    // --- helpers ----------------------------------------------------------------------------------------------

    private void put( final String k, final String v ) {
        _values.put( k, ( v == null ) ? "" : v );
    }

    String get( final String k ) {
        return _values.get( k );
    }

    private static String emptyToNull( final String s ) {
        return ForesterUtil.isEmpty( s ) ? null : s;
    }

    private static <E extends Enum<E>> E enumOf( final Class<E> type, final String name ) {
        if ( ForesterUtil.isEmpty( name ) ) {
            return null;
        }
        try {
            return Enum.valueOf( type, name );
        }
        catch ( final IllegalArgumentException e ) {
            return null; // a value this version does not know -- skip it rather than fail the load
        }
    }

    private static String joinList( final List<String> items ) {
        if ( items == null ) {
            return null;
        }
        final StringBuilder sb = new StringBuilder();
        for( final String s : items ) {
            if ( sb.length() > 0 ) {
                sb.append( '|' );
            }
            sb.append( s );
        }
        return sb.toString();
    }

    private static List<String> splitList( final String s ) {
        final List<String> out = new ArrayList<String>();
        if ( ForesterUtil.isEmpty( s ) ) {
            return out;
        }
        for( final String p : s.split( "\\|", -1 ) ) {
            if ( !ForesterUtil.isEmpty( p ) ) {
                out.add( p );
            }
        }
        return out;
    }

    private static String joinColumns( final List<AnnotationColumns.ColumnSpec> specs ) {
        if ( specs == null ) {
            return null;
        }
        final StringBuilder sb = new StringBuilder();
        for( final AnnotationColumns.ColumnSpec c : specs ) {
            if ( sb.length() > 0 ) {
                sb.append( '|' );
            }
            sb.append( c._ref ).append( '~' ).append( c._type.name() ).append( '~' ).append( c._shape.name() )
                    .append( '~' ).append( c._normalized );
        }
        return sb.toString();
    }

    private static List<AnnotationColumns.ColumnSpec> parseColumns( final String s ) {
        final List<AnnotationColumns.ColumnSpec> out = new ArrayList<AnnotationColumns.ColumnSpec>();
        for( final String item : splitList( s ) ) {
            final String[] p = item.split( "~", -1 );
            if ( p.length < 2 ) {
                continue;
            }
            final AnnotationColumns.Type t = enumOf( AnnotationColumns.Type.class, p[ 1 ] );
            if ( t == null ) {
                continue;
            }
            final AnnotationColumns.SymbolShape shape = ( p.length > 2 )
                    ? enumOf( AnnotationColumns.SymbolShape.class, p[ 2 ] ) : null;
            final boolean norm = ( p.length > 3 ) && Boolean.parseBoolean( p[ 3 ] );
            out.add( new AnnotationColumns.ColumnSpec( p[ 0 ], t,
                    ( shape == null ) ? AnnotationColumns.SymbolShape.CIRCLE : shape, norm ) );
        }
        return out;
    }

    private static String joinCladeLevels( final List<CladeLevel.Spec> specs ) {
        if ( specs == null ) {
            return null;
        }
        final StringBuilder sb = new StringBuilder();
        for( final CladeLevel.Spec c : specs ) {
            if ( sb.length() > 0 ) {
                sb.append( '|' );
            }
            sb.append( c.getRank() ).append( '~' ).append( c.getLabelAngle().name() );
        }
        return sb.toString();
    }

    private static List<CladeLevel.Spec> parseCladeLevels( final String s ) {
        final List<CladeLevel.Spec> out = new ArrayList<CladeLevel.Spec>();
        for( final String item : splitList( s ) ) {
            final String[] p = item.split( "~", -1 );
            if ( ForesterUtil.isEmpty( p[ 0 ] ) ) {
                continue;
            }
            final TreePanel.CLADE_LABEL_ANGLE a = ( p.length > 1 )
                    ? enumOf( TreePanel.CLADE_LABEL_ANGLE.class, p[ 1 ] ) : null;
            out.add( new CladeLevel.Spec( p[ 0 ], a ) );
        }
        return out;
    }

    // --- storage on the tree ----------------------------------------------------------------------------------

    /** Reads the figure stored on {@code phy}'s root, or null. */
    static FigureSpec readFrom( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() || ( phy.getRoot() == null )
                || !phy.getRoot().getNodeData().isHasProperties() ) {
            return null;
        }
        for( final org.forester.phylogeny.data.Property p : phy.getRoot().getNodeData().getProperties()
                .getProperties() ) {
            if ( FIGURE_REF.equals( p.getRef() ) ) {
                return parse( p.getValue() );
            }
        }
        return null;
    }

    /**
     * Stores {@code spec} on {@code phy}'s root, replacing any figure already there (or removing it when
     * {@code spec} is null). Called just before the tree is written, so a saved file carries the figure that was
     * on screen.
     */
    static void writeToTree( final Phylogeny phy, final FigureSpec spec ) {
        if ( ( phy == null ) || phy.isEmpty() || ( phy.getRoot() == null ) ) {
            return;
        }
        final org.forester.phylogeny.PhylogenyNode root = phy.getRoot();
        if ( root.getNodeData().getProperties() != null ) {
            final List<org.forester.phylogeny.data.Property> keep =
                    new ArrayList<org.forester.phylogeny.data.Property>();
            for( final org.forester.phylogeny.data.Property p : root.getNodeData().getProperties()
                    .getProperties() ) {
                if ( !FIGURE_REF.equals( p.getRef() ) ) {
                    keep.add( p );
                }
            }
            final org.forester.phylogeny.data.PropertiesList fresh =
                    new org.forester.phylogeny.data.PropertiesList();
            for( final org.forester.phylogeny.data.Property p : keep ) {
                fresh.addProperty( p );
            }
            root.getNodeData().setProperties( fresh );
        }
        if ( ( spec == null ) || spec.isEmpty() ) {
            return;
        }
        if ( root.getNodeData().getProperties() == null ) {
            root.getNodeData().setProperties( new org.forester.phylogeny.data.PropertiesList() );
        }
        root.getNodeData().getProperties()
                .addProperty( new org.forester.phylogeny.data.Property( FIGURE_REF, spec.toPropertyValue(), "",
                        "xsd:string", org.forester.phylogeny.data.Property.AppliesTo.PHYLOGENY ) );
    }

    private FigureSpec() {
    }
}
