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

import java.io.IOException;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;

/** Headless tests for {@link AuspiceJsonParser}: a small synthetic Auspice v2 dataset -> topology, dates + confidence,
 *  divergence + trait properties, trait-confidence as the pie-shaped {@code _set}/{@code _set_prob} pair, time branch
 *  lengths, and rejection of non-v2 / malformed input. */
public final class AuspiceJsonParserTest {

    private static final String DATASET = "{"
            + "\"version\":\"v2\","
            + "\"meta\":{\"title\":\"Test flu\",\"updated\":\"2024-01-01\"},"
            + "\"tree\":{"
            + "  \"name\":\"ROOT\","
            + "  \"node_attrs\":{\"div\":0,\"num_date\":{\"value\":2019.0,\"confidence\":[2018.8,2019.1]}},"
            + "  \"branch_attrs\":{\"labels\":{\"clade\":\"20A\"}},"
            + "  \"children\":["
            + "    {\"name\":\"A\",\"node_attrs\":{\"div\":0.002,"
            + "        \"num_date\":{\"value\":2020.0,\"confidence\":[2019.9,2020.2]},"
            + "        \"country\":{\"value\":\"USA\",\"confidence\":{\"USA\":0.9,\"Canada\":0.1}}}},"
            + "    {\"name\":\"B\",\"node_attrs\":{\"div\":0.003,\"num_date\":{\"value\":2020.5},"
            + "        \"country\":{\"value\":\"Canada\"}}}"
            + "  ]"
            + "}}";

    public static void main( final String[] args ) {
        System.out.println( "AuspiceJsonParser: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        try {
            final AuspiceJsonParser parser = new AuspiceJsonParser();
            parser.setSource( new StringBuffer( DATASET ) ); // a StringBuffer source is read as content
            final Phylogeny[] phys = parser.parse();
            if ( ( phys == null ) || ( phys.length != 1 ) ) {
                return fail( "expected exactly one phylogeny" );
            }
            final Phylogeny phy = phys[ 0 ];
            if ( !"Test flu".equals( phy.getName() ) ) {
                return fail( "phylogeny name from meta.title, got " + phy.getName() );
            }
            if ( phy.getNumberOfExternalNodes() != 2 ) {
                return fail( "expected 2 tips, got " + phy.getNumberOfExternalNodes() );
            }
            final PhylogenyNode root = phy.getRoot();
            if ( !"ROOT".equals( root.getName() ) ) {
                return fail( "root name" );
            }
            // root date value + 95% confidence -> Date value/min/max
            if ( !root.getNodeData().isHasDate() || ( dv( root ) == null ) || ( Math.abs( dv( root ) - 2019.0 ) > 1e-9 )
                    || ( Math.abs( root.getNodeData().getDate().getMin().doubleValue() - 2018.8 ) > 1e-9 )
                    || ( Math.abs( root.getNodeData().getDate().getMax().doubleValue() - 2019.1 ) > 1e-9 ) ) {
                return fail( "root num_date value/confidence -> Date value/min/max" );
            }
            // clade label -> a property on the root
            if ( !"20A".equals( prop( root, AuspiceJsonParser.PREFIX + "clade_label" ) ) ) {
                return fail( "branch_attrs.labels.clade -> nextstrain:clade_label" );
            }
            final PhylogenyNode a = named( phy, "A" );
            final PhylogenyNode b = named( phy, "B" );
            if ( ( a == null ) || ( b == null ) ) {
                return fail( "tips A and B" );
            }
            // trait value + divergence -> properties
            if ( !"USA".equals( prop( a, AuspiceJsonParser.PREFIX + "country" ) ) ) {
                return fail( "country trait value -> nextstrain:country" );
            }
            if ( !"0.002".equals( prop( a, AuspiceJsonParser.PREFIX + "div" ) ) ) {
                return fail( "div -> nextstrain:div, got " + prop( a, AuspiceJsonParser.PREFIX + "div" ) );
            }
            // trait confidence -> the pie-shaped set + set_prob pair (state names quoted so a comma in one survives)
            if ( !"{\"USA\",\"Canada\"}".equals( prop( a, AuspiceJsonParser.PREFIX + "country_set" ) )
                    || !"{0.9,0.1}".equals( prop( a, AuspiceJsonParser.PREFIX + "country_set_prob" ) ) ) {
                return fail( "trait confidence -> nextstrain:country_set / _set_prob (got "
                        + prop( a, AuspiceJsonParser.PREFIX + "country_set" ) + " / "
                        + prop( a, AuspiceJsonParser.PREFIX + "country_set_prob" ) + ")" );
            }
            // B has a value but NO confidence -> a plain property, no _set pair
            if ( !"Canada".equals( prop( b, AuspiceJsonParser.PREFIX + "country" ) )
                    || ( prop( b, AuspiceJsonParser.PREFIX + "country_set" ) != null ) ) {
                return fail( "a value-only trait must not fabricate a _set distribution" );
            }
            // default branch lengths = successive num_date differences
            if ( ( Math.abs( a.getDistanceToParent() - 1.0 ) > 1e-9 )
                    || ( Math.abs( b.getDistanceToParent() - 1.5 ) > 1e-9 ) ) {
                return fail( "time branch lengths: A=" + a.getDistanceToParent() + " B=" + b.getDistanceToParent() );
            }
            if ( Math.abs( root.getDistanceToParent() - 0.0 ) > 1e-9 ) {
                return fail( "root branch length must be 0" );
            }
            // a TIP keeps its point date but NOT the interval (that would misfire the fossil range bars); an INTERNAL
            // node keeps its interval (root's min/max were asserted above)
            if ( ( dv( a ) == null ) || ( a.getNodeData().getDate().getMin() != null )
                    || ( a.getNodeData().getDate().getMax() != null ) ) {
                return fail( "a tip's date interval must be stripped, keeping only the point date" );
            }
            // non-v2 / malformed input must throw
            if ( !throwsOn( "{\"version\":\"v1\",\"tree\":{}}" ) || !throwsOn( "{\"tree\":{}}" )
                    || !throwsOn( "not json" ) || !throwsOn( "[1,2,3]" ) ) {
                return fail( "a non-Auspice-v2 / malformed input must throw" );
            }
            if ( !commaStateNameQuoted() ) {
                return false;
            }
            if ( !divergenceOnlyBranchLengths() ) {
                return false;
            }
            return true;
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    /** A trait state name containing a comma (e.g. "Korea, Republic of") must be quoted in the {@code _set} brace list
     *  so it survives parseBraceList as ONE state (else it splits into two and the ancestral pie breaks). */
    private static boolean commaStateNameQuoted() throws Exception {
        final String ds = "{\"version\":\"v2\",\"tree\":{\"name\":\"R\",\"node_attrs\":{"
                + "\"num_date\":{\"value\":2020.0},"
                + "\"country\":{\"value\":\"Korea, Republic of\","
                + "  \"confidence\":{\"Korea, Republic of\":0.8,\"Japan\":0.2}}},"
                + "\"children\":[{\"name\":\"t\",\"node_attrs\":{\"num_date\":{\"value\":2020.5}}}]}}";
        final AuspiceJsonParser p = new AuspiceJsonParser();
        p.setSource( new StringBuffer( ds ) );
        final PhylogenyNode root = p.parse()[ 0 ].getRoot();
        final String set = prop( root, AuspiceJsonParser.PREFIX + "country_set" );
        if ( !"{\"Korea, Republic of\",\"Japan\"}".equals( set ) ) {
            return fail( "a comma-containing state name must be quoted, got " + set );
        }
        return true;
    }

    /** A divergence-only build (no num_date on any node) must fall back to branch lengths derived from the cumulative
     *  {@code div} property (successive differences), with the distance unit set to substitutions/site. */
    private static boolean divergenceOnlyBranchLengths() throws Exception {
        final String ds = "{\"version\":\"v2\",\"tree\":{\"name\":\"R\",\"node_attrs\":{\"div\":0.0},\"children\":["
                + "{\"name\":\"A\",\"node_attrs\":{\"div\":0.01}},"
                + "{\"name\":\"B\",\"node_attrs\":{\"div\":0.03}}]}}";
        final AuspiceJsonParser p = new AuspiceJsonParser();
        p.setSource( new StringBuffer( ds ) );
        final Phylogeny phy = p.parse()[ 0 ];
        final PhylogenyNode root = phy.getRoot();
        final PhylogenyNode a = named( phy, "A" );
        final PhylogenyNode b = named( phy, "B" );
        if ( ( a == null ) || ( b == null ) ) {
            return fail( "divergence-only: tips A and B" );
        }
        if ( ( Math.abs( root.getDistanceToParent() - 0.0 ) > 1e-12 )
                || ( Math.abs( a.getDistanceToParent() - 0.01 ) > 1e-12 )
                || ( Math.abs( b.getDistanceToParent() - 0.03 ) > 1e-12 ) ) {
            return fail( "divergence-only branch lengths from div: A=" + a.getDistanceToParent() + " B="
                    + b.getDistanceToParent() );
        }
        if ( !"subs/site".equals( phy.getDistanceUnit() ) ) {
            return fail( "divergence-only distance unit, got " + phy.getDistanceUnit() );
        }
        // no node has a date in a divergence-only build
        if ( dv( root ) != null ) {
            return fail( "a divergence-only build must carry no dates" );
        }
        return true;
    }

    private static boolean throwsOn( final String json ) {
        try {
            final AuspiceJsonParser p = new AuspiceJsonParser();
            p.setSource( new StringBuffer( json ) );
            p.parse();
            return false;
        }
        catch ( final IOException e ) {
            return true;
        }
    }

    private static Double dv( final PhylogenyNode n ) {
        return ( n.getNodeData().isHasDate() && ( n.getNodeData().getDate().getValue() != null ) )
                ? n.getNodeData().getDate().getValue().doubleValue() : null;
    }

    private static PhylogenyNode named( final Phylogeny phy, final String name ) {
        for ( final java.util.Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( name.equals( n.getName() ) ) {
                return n;
            }
        }
        return null;
    }

    private static String prop( final PhylogenyNode n, final String ref ) {
        if ( n.getNodeData().getProperties() == null ) {
            return null;
        }
        for ( final Property p : n.getNodeData().getProperties().getProperties() ) {
            if ( ref.equals( p.getRef() ) ) {
                return p.getValue();
            }
        }
        return null;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [AuspiceJsonParserTest] " + msg );
        return false;
    }

    private AuspiceJsonParserTest() {
    }
}
