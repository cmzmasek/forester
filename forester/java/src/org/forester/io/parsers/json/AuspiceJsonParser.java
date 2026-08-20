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

package org.forester.io.parsers.json;

import java.io.BufferedReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.List;
import java.util.Map;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.parsers.util.PhylogenyParserException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.util.ForesterUtil;

/**
 * Reads an <b>Auspice / Nextstrain v2</b> {@code dataset.json} (the modern single-file format) into a {@link Phylogeny},
 * mapping its per-node data onto forester's native model so Archaeopteryx's existing overlays light it up:
 * <ul>
 * <li>{@code node_attrs.num_date.value} &rarr; a {@code <date>} value (decimal year) &rarr; the Calendar time axis; its
 *     {@code .confidence [lo,hi]} &rarr; the date min/max &rarr; the Node Age bars/spindles.</li>
 * <li>{@code node_attrs.div} &rarr; a {@code nextstrain:div} property (the divergence measure, kept for a future
 *     time&harr;divergence view).</li>
 * <li>every discrete trait (country, region, clade_membership, host, ...) &rarr; a {@code nextstrain:<key>} node
 *     property (colour-by / annotation columns / search); its {@code .confidence {state:prob}} &rarr; a
 *     {@code nextstrain:<key>_set} + {@code _set_prob} pair &rarr; the Ancestral-State Pies.</li>
 * <li>{@code branch_attrs.labels.clade} &rarr; a {@code nextstrain:clade_label} property.</li>
 * </ul>
 * The default branch lengths are TIME (successive {@code num_date} differences). Deliberately NOT ingested: the map
 * ({@code geo_resolutions}), entropy ({@code genome_annotations}), and frequencies panels -- Archaeopteryx is a tree
 * viewer, not a phylodynamics dashboard. The generic JSON tokenising is done by {@link JsonParser}.
 */
public final class AuspiceJsonParser implements PhylogenyParser {

    /** Namespace for the node properties this parser writes (so Auspice traits are colour-able/searchable, and clearly
     *  distinguished from BEAST's {@code beast:} namespace). */
    public static final String  PREFIX = "nextstrain:";

    private Object              _source;

    @Override
    public void setSource( final Object source ) throws PhylogenyParserException, IOException {
        _source = source;
    }

    @Override
    public String getName() {
        return "Auspice/Nextstrain JSON parser";
    }

    @Override
    public Phylogeny[] parse() throws IOException {
        if ( _source == null ) {
            throw new IOException( "no source has been set" );
        }
        final Object root_value = JsonParser.parse( readSource( _source ) );
        final Map<String, Object> doc = asObject( root_value );
        if ( doc == null ) {
            throw new IOException( "not an Auspice dataset (the JSON root is not an object)" );
        }
        if ( !"v2".equals( doc.get( "version" ) ) || !( doc.get( "tree" ) instanceof Map ) ) {
            throw new IOException( "not an Auspice v2 dataset (expected \"version\":\"v2\" and a \"tree\" object)" );
        }
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = buildNode( asObject( doc.get( "tree" ) ) );
        phy.setRoot( root );
        phy.setRooted( true );
        final Map<String, Object> meta = asObject( doc.get( "meta" ) );
        if ( meta != null ) {
            final String title = asString( meta.get( "title" ) );
            if ( !ForesterUtil.isEmpty( title ) ) {
                phy.setName( title );
            }
        }
        final boolean time_resolved = hasAnyDate( root );
        if ( time_resolved ) {
            setTimeBranchLengths( root, null ); // default view = time (num_date deltas); divergence kept as a property
        }
        else {
            // a divergence-only Auspice build carries no num_date on any node (only the cumulative "div"); fall back to
            // div deltas as the branch lengths so the tree still lays out with meaningful lengths instead of a cladogram
            setDivBranchLengths( root, null );
        }
        // a TIP is a dated sample -- keep its point date (for the calendar axis) but drop the date INTERVAL: the
        // divergence-time UNCERTAINTY (the confidence -> Node Age spindles) belongs to the INTERNAL nodes, and a tip
        // interval would otherwise auto-enable the Fossil Range Bars (a geologic concept) on a viral tree.
        for ( final PhylogenyNode ext : phy.getExternalNodes() ) {
            if ( ext.getNodeData().isHasDate() ) {
                final Date d = ext.getNodeData().getDate();
                if ( ( d.getMin() != null ) || ( d.getMax() != null ) ) {
                    ext.getNodeData().setDate( new Date( d.getDesc(), d.getValue(), null, null, d.getUnit() ) );
                }
            }
        }
        phy.setDistanceUnit( time_resolved ? "year" : "subs/site" );
        phy.externalNodesHaveChanged();
        return new Phylogeny[] { phy };
    }

    /** Recursively build a node (and its subtree) from an Auspice tree-node object. */
    private static PhylogenyNode buildNode( final Map<String, Object> jn ) {
        final PhylogenyNode node = new PhylogenyNode();
        final String name = asString( jn.get( "name" ) );
        if ( !ForesterUtil.isEmpty( name ) ) {
            node.setName( name );
        }
        final Map<String, Object> attrs = asObject( jn.get( "node_attrs" ) );
        if ( attrs != null ) {
            applyNodeAttrs( node, attrs );
        }
        final Map<String, Object> battrs = asObject( jn.get( "branch_attrs" ) );
        if ( battrs != null ) {
            final Map<String, Object> labels = asObject( battrs.get( "labels" ) );
            if ( labels != null ) {
                final String clade = asString( labels.get( "clade" ) );
                if ( !ForesterUtil.isEmpty( clade ) ) {
                    addProperty( node, PREFIX + "clade_label", clade );
                }
            }
        }
        final List<Object> children = asArray( jn.get( "children" ) );
        if ( children != null ) {
            for ( final Object c : children ) {
                final Map<String, Object> cj = asObject( c );
                if ( cj != null ) {
                    node.addAsChild( buildNode( cj ) );
                }
            }
        }
        return node;
    }

    private static void applyNodeAttrs( final PhylogenyNode node, final Map<String, Object> attrs ) {
        for ( final Map.Entry<String, Object> e : attrs.entrySet() ) {
            final String key = e.getKey();
            final Object val = e.getValue();
            if ( "num_date".equals( key ) ) {
                applyNumDate( node, asObject( val ) );
            }
            else if ( "div".equals( key ) ) {
                final Double d = asDouble( val ); // div is a bare number (cumulative divergence from the root)
                if ( d != null ) {
                    addProperty( node, PREFIX + "div", numberToString( d ) );
                }
            }
            else {
                final Map<String, Object> obj = asObject( val );
                if ( obj != null ) {
                    // a discrete trait: {value, confidence{state:prob}, entropy}
                    final Object v = obj.get( "value" );
                    if ( isScalar( v ) ) {
                        addProperty( node, PREFIX + key, scalarToString( v ) );
                    }
                    final Map<String, Object> conf = asObject( obj.get( "confidence" ) );
                    if ( conf != null ) {
                        applyTraitConfidence( node, key, conf );
                    }
                }
                else if ( isScalar( val ) ) {
                    addProperty( node, PREFIX + key, scalarToString( val ) ); // a bare attr (accession, url, ...)
                }
                // a nested object without a "value" (e.g. author) -> not ingested in this cut
            }
        }
    }

    private static void applyNumDate( final PhylogenyNode node, final Map<String, Object> nd ) {
        if ( nd == null ) {
            return;
        }
        final Double value = asDouble( nd.get( "value" ) );
        if ( value == null ) {
            return;
        }
        BigDecimal min = null;
        BigDecimal max = null;
        final List<Object> conf = asArray( nd.get( "confidence" ) ); // [earlier, later] decimal years
        if ( ( conf != null ) && ( conf.size() == 2 ) ) {
            final Double lo = asDouble( conf.get( 0 ) );
            final Double hi = asDouble( conf.get( 1 ) );
            if ( lo != null ) {
                min = BigDecimal.valueOf( lo );
            }
            if ( hi != null ) {
                max = BigDecimal.valueOf( hi );
            }
        }
        node.getNodeData().setDate( new Date( "", BigDecimal.valueOf( value ), min, max, "year" ) );
    }

    /** Store a discrete trait's posterior distribution as the {@code nextstrain:<trait>_set} + {@code _set_prob} pair
     *  (the same brace-list shape the ancestral-state pie renderer reads), so the pies come essentially for free. */
    private static void applyTraitConfidence( final PhylogenyNode node, final String trait,
                                              final Map<String, Object> conf ) {
        final StringBuilder states = new StringBuilder();
        final StringBuilder probs = new StringBuilder();
        int n = 0;
        for ( final Map.Entry<String, Object> e : conf.entrySet() ) {
            final Double p = asDouble( e.getValue() );
            if ( ( p == null ) || ForesterUtil.isEmpty( e.getKey() ) ) {
                continue;
            }
            if ( n++ > 0 ) {
                states.append( ',' );
                probs.append( ',' );
            }
            // quote the state name so a comma/space in it (e.g. "Korea, Republic of") does not corrupt the brace list
            // (parseBraceList unquotes); drop any embedded double-quote (rare) so the quoting stays well-formed
            states.append( '"' ).append( e.getKey().replace( "\"", "" ) ).append( '"' );
            probs.append( numberToString( p ) );
        }
        if ( n > 0 ) {
            addProperty( node, PREFIX + trait + "_set", "{" + states + "}" );
            addProperty( node, PREFIX + trait + "_set_prob", "{" + probs + "}" );
        }
    }

    /** Default branch lengths = successive {@code num_date} differences (a time tree); a node/parent without a date gets
     *  the default (unset) length. Root length is 0. */
    private static void setTimeBranchLengths( final PhylogenyNode node, final Double parent_date ) {
        final Double d = nodeDate( node );
        if ( node.isRoot() ) {
            node.setDistanceToParent( 0.0 );
        }
        else if ( ( parent_date != null ) && ( d != null ) ) {
            node.setDistanceToParent( Math.max( 0.0, d - parent_date ) ); // clamp a (spurious) negative branch to 0
        }
        for ( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            setTimeBranchLengths( node.getChildNode( i ), d );
        }
    }

    private static Double nodeDate( final PhylogenyNode node ) {
        if ( node.getNodeData().isHasDate() && ( node.getNodeData().getDate().getValue() != null ) ) {
            return node.getNodeData().getDate().getValue().doubleValue();
        }
        return null;
    }

    /** True if any node in the subtree carries a {@code num_date} value (a time-resolved build); false for a
     *  divergence-only build (no node has a date). */
    private static boolean hasAnyDate( final PhylogenyNode node ) {
        if ( nodeDate( node ) != null ) {
            return true;
        }
        for ( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            if ( hasAnyDate( node.getChildNode( i ) ) ) {
                return true;
            }
        }
        return false;
    }

    /** Fallback branch lengths from the cumulative {@code nextstrain:div} property (successive differences, a
     *  divergence tree); a node/parent without a div gets the default (unset) length. Root length is 0. */
    private static void setDivBranchLengths( final PhylogenyNode node, final Double parent_div ) {
        final Double d = nodeDiv( node );
        if ( node.isRoot() ) {
            node.setDistanceToParent( 0.0 );
        }
        else if ( ( parent_div != null ) && ( d != null ) ) {
            node.setDistanceToParent( Math.max( 0.0, d - parent_div ) ); // clamp a (spurious) negative branch to 0
        }
        for ( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            setDivBranchLengths( node.getChildNode( i ), d );
        }
    }

    private static Double nodeDiv( final PhylogenyNode node ) {
        final PropertiesList pl = node.getNodeData().getProperties();
        if ( pl != null ) {
            for ( final Property p : pl.getProperties() ) {
                if ( ( PREFIX + "div" ).equals( p.getRef() ) ) {
                    try {
                        return Double.valueOf( p.getValue() );
                    }
                    catch ( final NumberFormatException e ) {
                        return null;
                    }
                }
            }
        }
        return null;
    }

    private static void addProperty( final PhylogenyNode node, final String ref, final String value ) {
        if ( ForesterUtil.isEmpty( value ) ) {
            return;
        }
        PropertiesList pl = node.getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            node.getNodeData().setProperties( pl );
        }
        final String datatype = isNumeric( value ) ? "xsd:decimal" : "xsd:string";
        pl.addProperty( new Property( ref, value, "", datatype, AppliesTo.NODE ) );
    }

    private static boolean isNumeric( final String s ) {
        if ( ForesterUtil.isEmpty( s ) ) {
            return false;
        }
        try {
            Double.parseDouble( s );
            return true;
        }
        catch ( final NumberFormatException e ) {
            return false;
        }
    }

    private static String readSource( final Object source ) throws IOException {
        final StringBuilder sb = new StringBuilder();
        try ( final BufferedReader r = ParserUtils.createReader( source, "UTF-8" ) ) {
            final char[] buf = new char[ 8192 ];
            int n;
            while ( ( n = r.read( buf ) ) != -1 ) {
                sb.append( buf, 0, n );
            }
        }
        return sb.toString();
    }

    // ---- small typed accessors over the JsonParser value model (Map / List / String / Double / Boolean / null) ----

    @SuppressWarnings( "unchecked" )
    private static Map<String, Object> asObject( final Object o ) {
        return ( o instanceof Map ) ? (Map<String, Object>) o : null;
    }

    @SuppressWarnings( "unchecked" )
    private static List<Object> asArray( final Object o ) {
        return ( o instanceof List ) ? (List<Object>) o : null;
    }

    private static String asString( final Object o ) {
        return ( o instanceof String ) ? (String) o : null;
    }

    private static Double asDouble( final Object o ) {
        return ( o instanceof Double ) ? (Double) o : null;
    }

    private static boolean isScalar( final Object o ) {
        return ( o instanceof String ) || ( o instanceof Double ) || ( o instanceof Boolean );
    }

    private static String scalarToString( final Object o ) {
        if ( o instanceof Double ) {
            return numberToString( (Double) o );
        }
        return String.valueOf( o );
    }

    /** A compact string for a JSON number: a whole value drops the trailing ".0" (a clean categorical/integer property),
     *  otherwise the plain decimal. */
    private static String numberToString( final double d ) {
        if ( ( d == Math.rint( d ) ) && !Double.isInfinite( d ) && ( Math.abs( d ) < 1e15 ) ) {
            return Long.toString( (long) d );
        }
        return Double.toString( d );
    }
}
