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

package org.forester.io.parsers.nhx;

import java.io.File;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Headless unit tests for {@link BeastAnnotationParser}: the field mapping (posterior&rarr;confidence,
 * height/HPD&rarr;date, rate/traits&rarr;properties), the top-level tokenizer + interval/number helpers, and an
 * end-to-end parse through {@link NHXParser} (the wiring seam).
 */
public final class BeastAnnotationParserTest {

    public static boolean test() {
        return testFieldMapping() && testHelpers() && testMalformedTolerance() && testEndToEnd() && testMoreFields()
                && testLengthLedBlob() && testOptionOffKeepsComment() && testRoundTrip();
    }

    /** A blob whose FIRST field is 'length' (TreeAnnotator emits field order varies) must still route to the
     *  structured parser -- guards the removed '[&length' routing exclusion in NHXParser. */
    private static boolean testLengthLedBlob() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(A[&length=1.0,posterior=0.9,height=1.2,height_95%_HPD={0.8,1.5},rate=0.01]:1.0,B:1.0);" );
            final PhylogenyNode a = p.parse()[ 0 ].getNode( "A" );
            if ( !a.getNodeData().isHasDate() || ( a.getNodeData().getDate().getMin() == null )
                    || !a.getBranchData().isHasConfidences() || ( prop( a, "beast:rate" ) == null ) ) {
                return fail( "a length-led BEAST blob must still be parsed (date interval + posterior + rate)" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "length-led parse threw: " + e );
        }
    }

    /** refKey sanitization ('%' -> '_'), the height_range interval fallback (and 95%-HPD precedence), and the
     *  'date' desc + plain 'height' value branches. */
    private static boolean testMoreFields() {
        final PhylogenyNode n = new PhylogenyNode();
        BeastAnnotationParser.apply( "&rate_95%_HPD={0.001,0.004}", n );
        final Property hpd = prop( n, "beast:rate_95_HPD" );
        if ( ( hpd == null ) || !"{0.001,0.004}".equals( hpd.getValue() ) || !"xsd:string".equals( hpd.getDataType() ) ) {
            return fail( "a 'rate_95%_HPD' key must sanitize to a valid beast:rate_95_HPD string property" );
        }
        final PhylogenyNode n2 = new PhylogenyNode();
        BeastAnnotationParser.apply( "&height=1.5,height_range={1.0,2.0}", n2 );
        final Date d2 = n2.getNodeData().getDate();
        if ( ( d2 == null ) || ( d2.getMin() == null ) || ( Math.abs( d2.getMin().doubleValue() - 1.0 ) > 1e-9 )
                || ( d2.getMax() == null ) || ( Math.abs( d2.getMax().doubleValue() - 2.0 ) > 1e-9 ) ) {
            return fail( "height_range must fill the date interval when no 95% HPD is present" );
        }
        final PhylogenyNode n3 = new PhylogenyNode();
        BeastAnnotationParser.apply( "&height_95%_HPD={0.8,1.2},height_range={0.5,1.5}", n3 );
        if ( Math.abs( n3.getNodeData().getDate().getMin().doubleValue() - 0.8 ) > 1e-9 ) {
            return fail( "height_95%_HPD must take precedence over height_range" );
        }
        final PhylogenyNode n4 = new PhylogenyNode();
        BeastAnnotationParser.apply( "&date=2014-06-10,height=3.0", n4 );
        final Date d4 = n4.getNodeData().getDate();
        if ( ( d4 == null ) || !"2014-06-10".equals( d4.getDesc() ) || ( d4.getValue() == null )
                || ( Math.abs( d4.getValue().doubleValue() - 3.0 ) > 1e-9 ) ) {
            return fail( "date -> desc and plain height -> value; got desc='" + ( d4 == null ? "?" : d4.getDesc() )
                    + "'" );
        }
        return true;
    }

    /** With the option OFF the raw [&...] blob is kept as an nh:comment (backward compatible), not structured. */
    private static boolean testOptionOffKeepsComment() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( false );
            p.setSource( "(A[&posterior=0.9,height=1.2]:1.0,B:1.0);" );
            final PhylogenyNode a = p.parse()[ 0 ].getNode( "A" );
            if ( a.getNodeData().isHasDate() || a.getBranchData().isHasConfidences() ) {
                return fail( "with the option OFF, BEAST fields must NOT be parsed into structured data" );
            }
            if ( prop( a, "nh:comment" ) == null ) {
                return fail( "with the option OFF, the raw [&...] blob must be kept as an nh:comment" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "option-off parse threw: " + e );
        }
    }

    /** A parsed BEAST tree, saved as phyloXML and reloaded, reproduces the date intervals, posteriors and rate
     *  properties (the <date>/<confidence>/<property> writers + parsers round-trip). */
    private static boolean testRoundTrip() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "((A[&height=0.0,rate=0.003]:1.0,B[&height=0.0,rate=0.002]:1.0)"
                    + "[&posterior=0.95,height=1.0,height_95%_HPD={0.8,1.3},rate=0.0025]:0.5,"
                    + "C[&height=0.0,rate=0.004]:1.5)[&posterior=1.0,height=1.5,height_95%_HPD={1.2,1.9},rate=0.003];" );
            final Phylogeny phy = p.parse()[ 0 ];
            final File tmp = File.createTempFile( "beast_roundtrip", ".xml" );
            try {
                new PhylogenyWriter().toPhyloXML( phy, 0, tmp );
                final Phylogeny back = ParserBasedPhylogenyFactory.getInstance()
                        .create( tmp, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
                int intervals = 0;
                int posteriors = 0;
                int rates = 0;
                for( final PhylogenyNodeIterator it = back.iteratorPreorder(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if ( n.getNodeData().isHasDate() && ( n.getNodeData().getDate().getMin() != null )
                            && ( n.getNodeData().getDate().getMax() != null ) ) {
                        intervals++;
                    }
                    if ( n.getBranchData().isHasConfidences() ) {
                        posteriors++;
                    }
                    if ( prop( n, "beast:rate" ) != null ) {
                        rates++;
                    }
                }
                if ( ( intervals != 2 ) || ( posteriors != 2 ) || ( rates != 5 ) ) {
                    return fail( "phyloXML round-trip lost data: intervals=" + intervals + " posteriors=" + posteriors
                            + " rates=" + rates );
                }
                return true;
            }
            finally {
                tmp.delete();
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "round-trip threw: " + e );
        }
    }

    private static boolean testFieldMapping() {
        final PhylogenyNode n = new PhylogenyNode();
        BeastAnnotationParser.apply(
                "&posterior=0.99,height_median=1.44,height_mean=1.40,height_95%_HPD={1.435,1.465},"
                        + "rate=0.0031,location=\"Africa\"",
                n );
        if ( !n.getBranchData().isHasConfidences()
                || ( Math.abs( n.getBranchData().getConfidence( 0 ).getValue() - 0.99 ) > 1e-9 )
                || !"posterior".equals( n.getBranchData().getConfidence( 0 ).getType() ) ) {
            return fail( "posterior must become a Confidence of type 'posterior'" );
        }
        final Date d = n.getNodeData().getDate();
        if ( d == null ) {
            return fail( "a height field must produce a <date>" );
        }
        if ( ( d.getValue() == null ) || ( Math.abs( d.getValue().doubleValue() - 1.44 ) > 1e-9 ) ) {
            return fail( "date value must be height_median (1.44), got " + d.getValue() );
        }
        if ( ( d.getMin() == null ) || ( Math.abs( d.getMin().doubleValue() - 1.435 ) > 1e-9 )
                || ( d.getMax() == null ) || ( Math.abs( d.getMax().doubleValue() - 1.465 ) > 1e-9 ) ) {
            return fail( "date min/max must be the 95% HPD bounds {1.435,1.465}, got " + d.getMin() + "/" + d.getMax() );
        }
        final Property rate = prop( n, "beast:rate" );
        if ( ( rate == null ) || !"xsd:decimal".equals( rate.getDataType() ) || !"0.0031".equals( rate.getValue() ) ) {
            return fail( "rate must become a numeric (xsd:decimal) beast:rate property" );
        }
        final Property loc = prop( n, "beast:location" );
        if ( ( loc == null ) || !"xsd:string".equals( loc.getDataType() ) || !"Africa".equals( loc.getValue() ) ) {
            return fail( "location must become a categorical property with quotes stripped, got "
                    + ( loc == null ? "null" : loc.getValue() ) );
        }
        // height_mean is used only when there is no height_median
        final PhylogenyNode n2 = new PhylogenyNode();
        BeastAnnotationParser.apply( "&height_mean=2.5", n2 );
        if ( ( n2.getNodeData().getDate() == null ) || ( n2.getNodeData().getDate().getValue() == null )
                || ( Math.abs( n2.getNodeData().getDate().getValue().doubleValue() - 2.5 ) > 1e-9 ) ) {
            return fail( "height_mean must be the date value when no median is present" );
        }
        return true;
    }

    private static boolean testHelpers() {
        final List<String> toks = BeastAnnotationParser.splitTopLevel( "a=1,b={2,3},c=\"x,y\",d=4" );
        if ( ( toks.size() != 4 ) || !toks.get( 1 ).equals( "b={2,3}" ) || !toks.get( 2 ).equals( "c=\"x,y\"" ) ) {
            return fail( "splitTopLevel must not split commas inside {} or \"\": " + toks );
        }
        final String[] iv = BeastAnnotationParser.parseInterval( "{1.2,3.4}" );
        if ( ( iv == null ) || !"1.2".equals( iv[ 0 ] ) || !"3.4".equals( iv[ 1 ] ) ) {
            return fail( "parseInterval must return the two raw numbers of {1.2,3.4}" );
        }
        if ( ( BeastAnnotationParser.parseInterval( "{1.2}" ) != null )
                || ( BeastAnnotationParser.parseInterval( "5" ) != null )
                || ( BeastAnnotationParser.parseInterval( "{a,b}" ) != null ) ) {
            return fail( "parseInterval must reject non-two-number sets" );
        }
        if ( ( BeastAnnotationParser.parseNumber( "0.003" ) == null )
                || ( BeastAnnotationParser.parseNumber( "abc" ) != null )
                || ( BeastAnnotationParser.parseNumber( "NaN" ) != null )
                || ( BeastAnnotationParser.parseNumber( "{1,2}" ) != null ) ) {
            return fail( "parseNumber must accept finite numbers only" );
        }
        return true;
    }

    private static boolean testMalformedTolerance() {
        final PhylogenyNode n = new PhylogenyNode();
        // a bare token, an empty value, and a good field mixed together -- the good fields must still parse
        BeastAnnotationParser.apply( "&garbage,empty=,posterior=0.9,rate=0.01", n );
        if ( !n.getBranchData().isHasConfidences()
                || ( Math.abs( n.getBranchData().getConfidence( 0 ).getValue() - 0.9 ) > 1e-9 ) ) {
            return fail( "malformed fields must be skipped while good fields (posterior) still parse" );
        }
        if ( prop( n, "beast:rate" ) == null ) {
            return fail( "rate must still parse alongside malformed fields" );
        }
        // null / empty input must be a safe no-op
        BeastAnnotationParser.apply( null, n );
        BeastAnnotationParser.apply( "&", new PhylogenyNode() );
        return true;
    }

    private static boolean testEndToEnd() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "((A[&height=0.0,rate=0.003]:1.0,B[&height=0.0,rate=0.002]:1.0)"
                    + "[&posterior=0.95,height=1.0,height_95%_HPD={0.8,1.3},rate=0.0025]:0.5,"
                    + "C[&height=0.0,rate=0.004]:1.5)[&posterior=1.0,height=1.5,height_95%_HPD={1.2,1.9},rate=0.003];" );
            final Phylogeny[] phys = p.parse();
            if ( ( phys.length != 1 ) || ( phys[ 0 ].getNumberOfExternalNodes() != 3 ) ) {
                return fail( "should parse one 3-tip tree, got " + phys.length + " tree(s)" );
            }
            int internal_intervals = 0;
            int internal_posteriors = 0;
            int rate_props = 0;
            for( final PhylogenyNodeIterator it = phys[ 0 ].iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( !n.isExternal() && n.getNodeData().isHasDate() && ( n.getNodeData().getDate().getMin() != null )
                        && ( n.getNodeData().getDate().getMax() != null ) ) {
                    internal_intervals++;
                }
                if ( !n.isExternal() && n.getBranchData().isHasConfidences() ) {
                    internal_posteriors++;
                }
                if ( prop( n, "beast:rate" ) != null ) {
                    rate_props++;
                }
            }
            if ( internal_intervals != 2 ) {
                return fail( "both internal nodes must carry an HPD interval, got " + internal_intervals );
            }
            if ( internal_posteriors != 2 ) {
                return fail( "both internal nodes must carry a posterior confidence, got " + internal_posteriors );
            }
            if ( rate_props != 5 ) {
                return fail( "all 5 nodes carry a beast:rate property, got " + rate_props );
            }
            // and NO opaque nh:comment blob is left behind
            for( final PhylogenyNodeIterator it = phys[ 0 ].iteratorPreorder(); it.hasNext(); ) {
                if ( prop( it.next(), "nh:comment" ) != null ) {
                    return fail( "the structured parse must leave no opaque nh:comment blob" );
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "end-to-end parse threw: " + e );
        }
    }

    private static Property prop( final PhylogenyNode node, final String ref ) {
        if ( node.getNodeData().getProperties() == null ) {
            return null;
        }
        final List<Property> ps = node.getNodeData().getProperties().getProperties( ref );
        return ps.isEmpty() ? null : ps.get( 0 );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "BeastAnnotationParser test failed: " + msg );
        return false;
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }
}
