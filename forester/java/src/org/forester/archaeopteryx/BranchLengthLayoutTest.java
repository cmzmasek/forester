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

    /** root(10) -> x(6) -> [ t1(2), t2(1) ]; root -> y(5) -> t3(0)... as ages: FIVE branches (x, t1, t2, y, t3),
     *  every node dated, every branch rated 0.002, every length 99 until a layout is applied. */
    private static Phylogeny fiveBranches() {
        final PhylogenyNode root = new PhylogenyNode();
        date( root, 10 );
        final PhylogenyNode x = named( "x", 6 );
        final PhylogenyNode y = named( "y", 5 );
        final PhylogenyNode t1 = named( "t1", 2 );
        final PhylogenyNode t2 = named( "t2", 1 );
        final PhylogenyNode t3 = named( "t3", 0 );
        root.addAsChild( x );
        root.addAsChild( y );
        x.addAsChild( t1 );
        x.addAsChild( t2 );
        y.addAsChild( t3 );
        for( final PhylogenyNode n : new PhylogenyNode[] { x, y, t1, t2, t3 } ) {
            n.setDistanceToParent( 99 );
            prop( n, BranchLengthLayout.RATE_PROPERTY_REF, 0.002 );
        }
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
            // --- a VIEW is laid out the way the tree it is a view of runs. A cherry, calendar dates, one of its two
            //     spans running backwards: by itself it is a tie and reads as ages; told the whole tree's way, its
            //     spans carry the signs the whole tree gives them
            final Phylogeny cherry = twoTips();
            date( cherry.getRoot(), 2004 );
            date( node( cherry, "a" ), 2003 );
            date( node( cherry, "b" ), 2009 );
            for( final String tip : new String[] { "a", "b" } ) {
                prop( node( cherry, tip ), BranchLengthLayout.RATE_PROPERTY_REF, 0.002 );
            }
            if ( BranchLengthLayout.datesIncreaseTowardTips( cherry ) ) {
                return fail( "fixture: one span up, one down -- the cherry cannot say which way its dates run" );
            }
            BranchLengthLayout.applyTime( cherry, true );
            if ( !eq( lengthOf( cherry, "a" ), -1 ) || !eq( lengthOf( cherry, "b" ), 5 ) ) {
                return fail( "told the dates increase toward the tips: a = 2003 - 2004, b = 2009 - 2004; got a="
                        + lengthOf( cherry, "a" ) + " b=" + lengthOf( cherry, "b" ) );
            }
            BranchLengthLayout.applyDivergence( cherry, true );
            if ( ( Double.doubleToRawLongBits( lengthOf( cherry, "a" ) ) != 0L ) || !eq( lengthOf( cherry, "b" ), 0.01 ) ) {
                return fail( "...and in divergence a is 0 and b is 5 x 0.002; got a=" + lengthOf( cherry, "a" ) + " b="
                        + lengthOf( cherry, "b" ) );
            }
            BranchLengthLayout.applyTime( cherry, false );
            if ( !eq( lengthOf( cherry, "a" ), 1 ) || !eq( lengthOf( cherry, "b" ), -5 ) ) {
                return fail( "told they are ages, every span has the other sign; got a=" + lengthOf( cherry, "a" )
                        + " b=" + lengthOf( cherry, "b" ) );
            }
            // left to itself, a tree is laid out the way ITS dates run
            BranchLengthLayout.applyTime( cherry );
            if ( !eq( lengthOf( cherry, "a" ), 1 ) || !eq( lengthOf( cherry, "b" ), -5 ) ) {
                return fail( "a tie reads as ages; got a=" + lengthOf( cherry, "a" ) + " b=" + lengthOf( cherry, "b" ) );
            }
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
            // ===== THE JOINT RULE: both layouts can state EVERY branch =====
            // --- time: every node states a date, the root included. One tip undated of a tree of 5 branches was a
            //     majority and was offered; it is not any more.
            final Phylogeny five = fiveBranches();
            if ( !BranchLengthLayout.isTimeDerivable( five ) || !BranchLengthLayout.isApplicable( five ) ) {
                return fail( "fixture: five branches, every node dated and every branch rated, must be offered" );
            }
            for( final String undated : new String[] { "t1", "x", "root" } ) {
                final Phylogeny p = fiveBranches();
                ( "root".equals( undated ) ? p.getRoot() : node( p, undated ) ).getNodeData().setDate( null );
                if ( BranchLengthLayout.isTimeDerivable( p ) || BranchLengthLayout.isApplicable( p ) ) {
                    return fail( "\"" + undated + "\" states no date: the time layout cannot state every branch" );
                }
                // ...and the tree stays in the layout it arrived in
                final double before = lengthOf( p, "t3" );
                BranchLengthLayout.applyDivergence( p );
                if ( lengthOf( p, "t3" ) != before ) {
                    return fail( "a tree that is not offered the switch must keep its lengths (\"" + undated
                            + "\" undated); t3 " + before + " -> " + lengthOf( p, "t3" ) );
                }
            }
            // --- a date that is STATED is stated, whether zero or negative
            final Phylogeny zero_and_negative = fiveBranches();
            date( node( zero_and_negative, "t1" ), 0 );
            date( node( zero_and_negative, "t2" ), -3 );
            if ( !BranchLengthLayout.isApplicable( zero_and_negative ) ) {
                return fail( "a date of 0 and a date of -3 are stated: the tree is offered the switch" );
            }
            // --- recorded divergence: every node states it, the root included; stated in part, there is no
            //     divergence layout, and no falling back on the rates
            final Phylogeny recorded = fiveBranches();
            for( final PhylogenyNodeIterator it = recorded.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                prop( n, BranchLengthLayout.DIV_PROPERTY_REF, n.isRoot() ? 0.0 : 0.004 );
            }
            if ( ( BranchLengthLayout.divergenceSource( recorded ) != BranchLengthLayout.DIVERGENCE_SOURCE.STORED )
                    || !BranchLengthLayout.isApplicable( recorded ) ) {
                return fail( "fixture: a divergence recorded on every node must be offered" );
            }
            for( final String without : new String[] { "t2", "x", "root" } ) {
                final Phylogeny p = fiveBranches();
                for( final PhylogenyNodeIterator it = p.iteratorPreorder(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if ( !( n.isRoot() ? "root" : n.getName() ).equals( without ) ) {
                        prop( n, BranchLengthLayout.DIV_PROPERTY_REF, n.isRoot() ? 0.0 : 0.004 );
                    }
                }
                if ( BranchLengthLayout.divergenceSource( p ) != BranchLengthLayout.DIVERGENCE_SOURCE.NONE ) {
                    return fail( "\"" + without + "\" records no divergence: no divergence layout, and the rates on"
                            + " every branch are not fallen back on; got " + BranchLengthLayout.divergenceSource( p ) );
                }
                if ( BranchLengthLayout.isApplicable( p ) ) {
                    return fail( "\"" + without + "\" records no divergence: not offered" );
                }
            }
            // ===== A SPAN MAY BE NEGATIVE =====
            // --- ages (largest at the root): t1 is dated OLDER than its parent x, 7 against 6
            final Phylogeny older = fiveBranches();
            date( node( older, "t1" ), 7 );
            if ( BranchLengthLayout.datesIncreaseTowardTips( older ) ) {
                return fail( "four of five spans run down: the dates are ages" );
            }
            BranchLengthLayout.applyTime( older );
            if ( !eq( lengthOf( older, "t1" ), -1 ) || !eq( lengthOf( older, "t2" ), 5 ) || !eq( lengthOf( older, "x" ), 4 ) ) {
                return fail( "time keeps the sign: t1 is 1 BEFORE its parent; got t1=" + lengthOf( older, "t1" )
                        + " t2=" + lengthOf( older, "t2" ) + " x=" + lengthOf( older, "x" ) );
            }
            // summed from the root, signed spans land every node on its own date: x at 10-6, t1 at 10-7
            if ( !eq( lengthOf( older, "x" ) + lengthOf( older, "t1" ), 10 - 7 ) ) {
                return fail( "root to t1 must span root date - t1 date" );
            }
            BranchLengthLayout.applyDivergence( older );
            if ( Double.doubleToRawLongBits( lengthOf( older, "t1" ) ) != 0L ) {
                return fail( "divergence draws a span that runs backwards at 0; got " + lengthOf( older, "t1" ) );
            }
            if ( !eq( lengthOf( older, "t2" ), 5 * 0.002 ) ) {
                return fail( "...and every other branch at rate x span; got t2=" + lengthOf( older, "t2" ) );
            }
            BranchLengthLayout.applyTime( older );
            if ( !eq( lengthOf( older, "t1" ), -1 ) ) {
                return fail( "and back in time the span is negative again; got " + lengthOf( older, "t1" ) );
            }
            // --- calendar dates (increasing toward the tips): the same tree, mirrored about year 2000
            final Phylogeny calendar = fiveBranches();
            for( final PhylogenyNodeIterator it = calendar.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                date( n, 2010 - n.getNodeData().getDate().getValue().doubleValue() );
            }
            date( node( calendar, "t1" ), 2003 ); // its parent x is dated 2004
            if ( !BranchLengthLayout.datesIncreaseTowardTips( calendar ) ) {
                return fail( "four of five spans run up: the dates are calendar dates" );
            }
            BranchLengthLayout.applyTime( calendar );
            if ( !eq( lengthOf( calendar, "t1" ), -1 ) || !eq( lengthOf( calendar, "t2" ), 5 )
                    || !eq( lengthOf( calendar, "x" ), 4 ) ) {
                return fail( "calendar time keeps the sign too; got t1=" + lengthOf( calendar, "t1" ) + " t2="
                        + lengthOf( calendar, "t2" ) + " x=" + lengthOf( calendar, "x" ) );
            }
            // --- a recorded divergence that FALLS along a branch is drawn at 0 as well
            final Phylogeny falling = fiveBranches();
            for( final PhylogenyNodeIterator it = falling.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                prop( n, BranchLengthLayout.DIV_PROPERTY_REF,
                      n.isRoot() ? 0.0 : ( "t1".equals( n.getName() ) ? 0.003 : ( "x".equals( n.getName() ) ? 0.004 : 0.009 ) ) );
            }
            BranchLengthLayout.applyDivergence( falling );
            if ( ( Double.doubleToRawLongBits( lengthOf( falling, "t1" ) ) != 0L ) || !eq( lengthOf( falling, "t2" ), 0.005 ) ) {
                return fail( "a recorded divergence falling from 0.004 to 0.003 is drawn at 0; got t1="
                        + lengthOf( falling, "t1" ) + " t2=" + lengthOf( falling, "t2" ) );
            }
            // --- which way the dates run: the majority of the pairs that differ; none, or a tie, reads as ages
            final Phylogeny level = twoTips();
            date( node( level, "a" ), 10 );
            date( node( level, "b" ), 10 );
            if ( BranchLengthLayout.datesIncreaseTowardTips( level ) ) {
                return fail( "no pair differs: ages" );
            }
            final Phylogeny tie = twoTips();
            date( node( tie, "a" ), 12 );
            date( node( tie, "b" ), 8 );
            if ( BranchLengthLayout.datesIncreaseTowardTips( tie ) ) {
                return fail( "one pair up, one down: ages" );
            }
            BranchLengthLayout.applyTime( level );
            if ( Double.doubleToRawLongBits( lengthOf( level, "a" ) ) != 0L ) {
                return fail( "a span of nothing is a plain 0, never -0.0" );
            }
            // a pair with EQUAL dates does not vote. Real builds are full of them (1405 of the 9205 branches of one
            // Nextstrain tree); counted as running down, four of them would outvote the one pair that runs up here,
            // and a calendar tree would be read as ages, every span with the wrong sign
            final Phylogeny mostly_level = fiveBranches();
            for( final PhylogenyNodeIterator it = mostly_level.iteratorPreorder(); it.hasNext(); ) {
                date( it.next(), 2000 );
            }
            date( node( mostly_level, "t2" ), 2001 );
            if ( !BranchLengthLayout.datesIncreaseTowardTips( mostly_level ) ) {
                return fail( "one pair runs up, four are level: the dates increase toward the tips" );
            }
            BranchLengthLayout.applyTime( mostly_level );
            if ( !eq( lengthOf( mostly_level, "t2" ), 1 ) || !eq( lengthOf( mostly_level, "t1" ), 0 ) ) {
                return fail( "...and the one span is +1; got t2=" + lengthOf( mostly_level, "t2" ) + " t1="
                        + lengthOf( mostly_level, "t1" ) );
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
            if ( BranchLengthLayout.isTimeDerivable( lone ) ) {
                return fail( "a tree without a branch has no time layout either" );
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
