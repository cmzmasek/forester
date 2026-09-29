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

    private static void propText( final PhylogenyNode n, final String ref, final String v ) {
        if ( n.getNodeData().getProperties() == null ) {
            n.getNodeData().setProperties( new PropertiesList() );
        }
        n.getNodeData().getProperties().addProperty( new Property( ref, v, "", "xsd:string", AppliesTo.NODE ) );
    }

    private static PhylogenyNode node( final Phylogeny p, final String name ) {
        for( final PhylogenyNodeIterator it = p.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( name.equals( n.getName() ) ) {
                return n;
            }
        }
        throw new IllegalStateException( "no node named " + name );
    }

    /** {@link #twoTips()} with "a" rated 0.005 and "b" stating {@code b_rate} as written (null: no rate at all). */
    private static Phylogeny ratedPair( final String b_rate ) {
        final Phylogeny p = twoTips();
        prop( node( p, "a" ), BranchLengthLayout.RATE_PROPERTY_REF, 0.005 );
        if ( b_rate != null ) {
            propText( node( p, "b" ), BranchLengthLayout.RATE_PROPERTY_REF, b_rate );
        }
        return p;
    }

    private static PhylogenyNode named( final String name, final double date ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        date( n, date );
        return n;
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
            // --- EVERY branch must state a rate (joint with Archaeopteryx.js). The pair differs by ONE rate. ---
            final Phylogeny both_rated = ratedPair( "0.0025" );
            if ( ( BranchLengthLayout.divergenceSource( both_rated ) != BranchLengthLayout.DIVERGENCE_SOURCE.CLOCK_RATE )
                    || !BranchLengthLayout.isApplicable( both_rated ) ) {
                return fail( "the rated half of the pair must be offered the toggle, or the refusal below proves nothing" );
            }
            final Phylogeny one_unrated = ratedPair( null );
            if ( !BranchLengthLayout.isTimeDerivable( one_unrated ) ) {
                return fail( "the unrated half of the pair must still have a time layout: the refusal is for the RATE" );
            }
            if ( BranchLengthLayout.divergenceSource( one_unrated ) != BranchLengthLayout.DIVERGENCE_SOURCE.NONE ) {
                return fail( "one branch without a rate: a clock rate is no divergence source, got "
                        + BranchLengthLayout.divergenceSource( one_unrated ) );
            }
            if ( BranchLengthLayout.isApplicable( one_unrated ) ) {
                return fail( "one branch without a rate must NOT be offered the toggle" );
            }
            // ...and a tree with no divergence source is left exactly as it was, never laid out at 0
            BranchLengthLayout.applyTime( one_unrated );
            BranchLengthLayout.applyDivergence( one_unrated );
            if ( !eq( lengthOf( one_unrated, "a" ), 4 ) || !eq( lengthOf( one_unrated, "b" ), 8 ) ) {
                return fail( "a tree with no divergence source must keep its lengths; got a="
                        + lengthOf( one_unrated, "a" ) + " b=" + lengthOf( one_unrated, "b" ) );
            }
            // --- what counts as a rate: finite and not negative, as a number is written ---
            for( final String good : new String[] { "0", "0.0", "1e-3", " 0.005 ", "12" } ) {
                if ( BranchLengthLayout.divergenceSource( ratedPair( good ) ) != BranchLengthLayout.DIVERGENCE_SOURCE.CLOCK_RATE ) {
                    return fail( "\"" + good + "\" is a rate" );
                }
            }
            for( final String bad : new String[] { "-0.001", "-1e-9", "fast", "NaN", "Infinity", "-Infinity", "", " " } ) {
                if ( BranchLengthLayout.divergenceSource( ratedPair( bad ) ) != BranchLengthLayout.DIVERGENCE_SOURCE.NONE ) {
                    return fail( "\"" + bad + "\" is no rate: the tree must have no divergence source" );
                }
            }
            // --- time BELOW a clade only: the clade's own branch and everything outside it are left alone ---
            // root(10) -> x(6) -> [ a(2), u(undated) ] ; root -> b(3); every length starts at 99
            final Phylogeny nested = new Phylogeny();
            final PhylogenyNode nr = new PhylogenyNode();
            date( nr, 10 );
            final PhylogenyNode nx = named( "x", 6 );
            final PhylogenyNode na = named( "a", 2 );
            final PhylogenyNode nu = named( "u", 0 );
            nu.getNodeData().setDate( null );
            final PhylogenyNode nb = named( "b", 3 );
            nr.addAsChild( nx );
            nr.addAsChild( nb );
            nx.addAsChild( na );
            nx.addAsChild( nu );
            nested.setRoot( nr );
            nested.externalNodesHaveChanged();
            for( final PhylogenyNode n : new PhylogenyNode[] { nx, na, nu, nb } ) {
                n.setDistanceToParent( 99 );
            }
            BranchLengthLayout.applyTimeBelow( nx );
            if ( !eq( na.getDistanceToParent(), 4 ) ) {
                return fail( "below the clade a branch spans its date gap (6 - 2); got " + na.getDistanceToParent() );
            }
            if ( !eq( nu.getDistanceToParent(), 0 ) ) {
                return fail( "below the clade an undated branch goes to 0; got " + nu.getDistanceToParent() );
            }
            if ( !eq( nx.getDistanceToParent(), 99 ) || !eq( nb.getDistanceToParent(), 99 ) ) {
                return fail( "the clade's own branch and the branches outside it must not be touched; got x="
                        + nx.getDistanceToParent() + " b=" + nb.getDistanceToParent() );
            }
            final PhylogenyNode deep = named( "deep", 1 );
            na.addAsChild( deep );
            deep.setDistanceToParent( 99 );
            BranchLengthLayout.applyTimeBelow( nx );
            if ( !eq( deep.getDistanceToParent(), 1 ) ) {
                return fail( "below the clade means all the way down (2 - 1); got " + deep.getDistanceToParent() );
            }
            BranchLengthLayout.applyTimeBelow( null ); // nothing to do, and nothing thrown
            // --- a rate written -0.0 is not negative, so it is a rate; its length must still be a plain 0 ---
            final Phylogeny minus_zero = ratedPair( "-0.0" );
            if ( BranchLengthLayout.divergenceSource( minus_zero ) != BranchLengthLayout.DIVERGENCE_SOURCE.CLOCK_RATE ) {
                return fail( "\"-0.0\" is not below zero: it is a rate" );
            }
            BranchLengthLayout.applyDivergence( minus_zero );
            if ( Double.doubleToRawLongBits( node( minus_zero, "b" ).getDistanceToParent() ) != 0L ) {
                return fail( "a rate of -0.0 must give a length of 0, not of -0.0 (it would be written \"-0.0\")" );
            }
            // --- the root has no branch: what it states as a rate is never asked for ---
            final Phylogeny bad_root = ratedPair( "0.0025" );
            propText( bad_root.getRoot(), BranchLengthLayout.RATE_PROPERTY_REF, "fast" );
            if ( BranchLengthLayout.divergenceSource( bad_root ) != BranchLengthLayout.DIVERGENCE_SOURCE.CLOCK_RATE ) {
                return fail( "the root's rate must not be asked for: it has no branch" );
            }
            // --- a tree of one node has no branch at all: "every branch" of none is not a clock model ---
            final Phylogeny lone = new Phylogeny();
            final PhylogenyNode lone_root = new PhylogenyNode();
            date( lone_root, 10 );
            prop( lone_root, BranchLengthLayout.RATE_PROPERTY_REF, 0.005 );
            lone.setRoot( lone_root );
            lone.externalNodesHaveChanged();
            if ( BranchLengthLayout.divergenceSource( lone ) != BranchLengthLayout.DIVERGENCE_SOURCE.NONE ) {
                return fail( "a tree without a branch has no divergence source" );
            }
            // --- a recorded divergence still wins, whatever the rates say (here: one branch has none) ---
            final Phylogeny stored_unrated = ratedPair( null );
            for( final PhylogenyNodeIterator it = stored_unrated.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                prop( n, BranchLengthLayout.DIV_PROPERTY_REF, n.isRoot() ? 0.0 : 0.003 );
            }
            if ( ( BranchLengthLayout.divergenceSource( stored_unrated ) != BranchLengthLayout.DIVERGENCE_SOURCE.STORED )
                    || !BranchLengthLayout.isApplicable( stored_unrated ) ) {
                return fail( "a recorded divergence needs no rates at all" );
            }
            // --- a date of exactly 0 IS a date: every contemporaneous BEAST tip sits at height 0 ---
            final Phylogeny zero_tips = twoTips();
            date( node( zero_tips, "a" ), 0 );
            date( node( zero_tips, "b" ), 0 );
            prop( node( zero_tips, "a" ), BranchLengthLayout.RATE_PROPERTY_REF, 0.005 );
            prop( node( zero_tips, "b" ), BranchLengthLayout.RATE_PROPERTY_REF, 0.002 );
            if ( node( zero_tips, "a" ).getNodeData().getDate().getValue().signum() != 0 ) {
                return fail( "fixture: the tip date must BE zero" );
            }
            if ( !BranchLengthLayout.isTimeDerivable( zero_tips ) ) {
                return fail( "tips dated 0 are dated: the tree has a time layout" );
            }
            if ( !BranchLengthLayout.isApplicable( zero_tips ) ) {
                return fail( "tips dated 0 and a rate on every branch must be offered the toggle" );
            }
            BranchLengthLayout.applyTime( zero_tips );
            if ( !eq( lengthOf( zero_tips, "a" ), 10 ) || !eq( lengthOf( zero_tips, "b" ), 10 ) ) {
                return fail( "a tip at 0 under a root at 10 spans 10; got a=" + lengthOf( zero_tips, "a" ) + " b="
                        + lengthOf( zero_tips, "b" ) );
            }
            BranchLengthLayout.applyDivergence( zero_tips );
            if ( !eq( lengthOf( zero_tips, "a" ), 0.05 ) || !eq( lengthOf( zero_tips, "b" ), 0.02 ) ) {
                return fail( "rate x time with a tip at 0; got a=" + lengthOf( zero_tips, "a" ) + " b="
                        + lengthOf( zero_tips, "b" ) );
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
}
