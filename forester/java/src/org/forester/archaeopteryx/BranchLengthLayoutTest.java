package org.forester.archaeopteryx;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

import java.math.BigDecimal;

/**
 * {@link BranchLengthLayout}: laying a dated tree out by time or by divergence, for the two shapes that carry both.
 * The point of the class is that ONE implementation serves an Auspice tree (which RECORDS divergence) and a BEAST
 * tree (where it is DERIVED from a clock rate) -- so both are exercised here, as is the refusal to offer the toggle
 * to a tree that only has one of the two.
 */
public class BranchLengthLayoutTest {

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }

    private static void date( final PhylogenyNode n, final double v ) {
        final Date d = new Date();
        d.setValue( new BigDecimal( String.valueOf( v ) ) );
        n.getNodeData().setDate( d );
    }

    private static void prop( final PhylogenyNode n, final String ref, final double v ) {
        if ( n.getNodeData().getProperties() == null ) {
            n.getNodeData().setProperties( new PropertiesList() );
        }
        n.getNodeData().getProperties()
                .addProperty( new Property( ref, String.valueOf( v ), "", "xsd:decimal", AppliesTo.NODE ) );
    }

    /** root(date 10) -> a(date 6), b(date 2); ages before present, as BEAST writes them. */
    private static Phylogeny twoTips() {
        final PhylogenyNode root = new PhylogenyNode();
        date( root, 10 );
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "a" );
        date( a, 6 );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "b" );
        date( b, 2 );
        root.addAsChild( a );
        root.addAsChild( b );
        final Phylogeny p = new Phylogeny();
        p.setRoot( root );
        p.externalNodesHaveChanged();
        return p;
    }

    private static double lengthOf( final Phylogeny p, final String name ) {
        for( final PhylogenyNodeIterator it = p.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( name.equals( n.getName() ) ) {
                return n.getDistanceToParent();
            }
        }
        return Double.NaN;
    }

    private static boolean fail( final String m ) {
        System.out.println( "  [BranchLengthLayoutTest] " + m );
        return false;
    }

    private static boolean eq( final double a, final double b ) {
        return Math.abs( a - b ) < 1e-9;
    }

    public static boolean test() {
        try {
            // --- a tree with dates but NO divergence source is not offered the toggle ---
            final Phylogeny dates_only = twoTips();
            if ( BranchLengthLayout.divergenceSource( dates_only ) != BranchLengthLayout.DIVERGENCE_SOURCE.NONE ) {
                return fail( "dates alone are not a divergence source" );
            }
            if ( BranchLengthLayout.isApplicable( dates_only ) ) {
                return fail( "a tree with only a time signal must NOT be offered the toggle" );
            }
            // --- BEAST shape: a per-branch clock rate, so divergence is DERIVED as rate x time ---
            final Phylogeny beast = twoTips();
            for( final PhylogenyNodeIterator it = beast.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( !n.isRoot() ) {
                    prop( n, BranchLengthLayout.RATE_PROPERTY_REF, 0.005 );
                }
            }
            if ( BranchLengthLayout.divergenceSource( beast ) != BranchLengthLayout.DIVERGENCE_SOURCE.CLOCK_RATE ) {
                return fail( "a per-branch rate must be recognised as a CLOCK_RATE divergence source" );
            }
            if ( !BranchLengthLayout.isApplicable( beast ) ) {
                return fail( "dates + a clock rate must offer the toggle (this is the BEAST case)" );
            }
            BranchLengthLayout.applyTime( beast );
            if ( !eq( lengthOf( beast, "a" ), 4 ) || !eq( lengthOf( beast, "b" ), 8 ) ) {
                return fail( "time lengths come from the date differences; got a=" + lengthOf( beast, "a" ) + " b="
                        + lengthOf( beast, "b" ) );
            }
            BranchLengthLayout.applyDivergence( beast );
            if ( !eq( lengthOf( beast, "a" ), 0.02 ) || !eq( lengthOf( beast, "b" ), 0.04 ) ) {
                return fail( "derived divergence must be rate x time; got a=" + lengthOf( beast, "a" ) + " b="
                        + lengthOf( beast, "b" ) );
            }
            BranchLengthLayout.applyTime( beast ); // and back again, exactly
            if ( !eq( lengthOf( beast, "a" ), 4 ) || !eq( lengthOf( beast, "b" ), 8 ) ) {
                return fail( "the toggle must be exactly reversible" );
            }
            // --- Auspice shape: a RECORDED cumulative divergence wins over any rate, and is differenced ---
            final Phylogeny auspice = twoTips();
            int i = 0;
            for( final PhylogenyNodeIterator it = auspice.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                prop( n, BranchLengthLayout.DIV_PROPERTY_REF, n.isRoot() ? 0.0 : ( 0.001 * ++i ) );
                if ( !n.isRoot() ) {
                    prop( n, BranchLengthLayout.RATE_PROPERTY_REF, 99 ); // must be IGNORED: recorded beats derived
                }
            }
            if ( BranchLengthLayout.divergenceSource( auspice ) != BranchLengthLayout.DIVERGENCE_SOURCE.STORED ) {
                return fail( "a recorded cumulative divergence must win over a derivable one" );
            }
            BranchLengthLayout.applyDivergence( auspice );
            if ( !eq( lengthOf( auspice, "a" ), 0.001 ) || !eq( lengthOf( auspice, "b" ), 0.002 ) ) {
                return fail( "recorded divergence is the successive difference; got a=" + lengthOf( auspice, "a" )
                        + " b=" + lengthOf( auspice, "b" ) );
            }
            // --- the label distinguishes a derived number from a recorded one ---
            if ( !"Divergence".equals( BranchLengthLayout.label( BranchLengthLayout.MODE.DIVERGENCE,
                                                                 BranchLengthLayout.DIVERGENCE_SOURCE.STORED ) ) ) {
                return fail( "a recorded divergence is just 'Divergence'" );
            }
            if ( !"Divergence (from clock rate)".equals( BranchLengthLayout
                    .label( BranchLengthLayout.MODE.DIVERGENCE, BranchLengthLayout.DIVERGENCE_SOURCE.CLOCK_RATE ) ) ) {
                return fail( "a DERIVED divergence must say so, so a reader does not quote an inference as a measurement" );
            }
            // --- a mostly-undated tree must not be offered a time layout ---
            final Phylogeny sparse = twoTips();
            for( final PhylogenyNodeIterator it = sparse.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( !n.isRoot() ) {
                    n.getNodeData().setDate( null );
                    prop( n, BranchLengthLayout.RATE_PROPERTY_REF, 0.005 );
                }
            }
            if ( BranchLengthLayout.isTimeDerivable( sparse ) ) {
                return fail( "a tree whose branches have no dated endpoints has no time layout" );
            }
            if ( BranchLengthLayout.isApplicable( sparse ) ) {
                return fail( "a rate without dates must NOT offer the toggle -- there is nothing to multiply" );
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
}
