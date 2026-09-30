package org.forester.archaeopteryx;

import org.forester.archaeopteryx.BranchLengthLayout.TimeLengths;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

import java.math.BigDecimal;

/**
 * {@link BranchLengthLayout}: laying a dated tree out by time or by divergence, for the two shapes that carry both.
 * ONE implementation serves an Auspice tree (which RECORDS divergence) and a BEAST tree (where it is DERIVED from a
 * clock rate) -- so both are exercised here, as are the refusals: the joint rule that both layouts must be able to
 * state every branch and have some depth, what counts as a rate, and that TIME is the lengths the tree had, given
 * back to the digit.
 */
public class BranchLengthLayoutTest {

    private static final BranchLengthLayout.DIVERGENCE_SOURCE STORED = BranchLengthLayout.DIVERGENCE_SOURCE.STORED;
    private static final BranchLengthLayout.DIVERGENCE_SOURCE CLOCK  = BranchLengthLayout.DIVERGENCE_SOURCE.CLOCK_RATE;
    private static final BranchLengthLayout.DIVERGENCE_SOURCE NONE   = BranchLengthLayout.DIVERGENCE_SOURCE.NONE;
    private static final String                               RATE   = BranchLengthLayout.RATE_PROPERTY_REF;
    private static final String                               DIV    = BranchLengthLayout.DIV_PROPERTY_REF;

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }

    private static void date( final PhylogenyNode n, final double v ) {
        final Date d = new Date();
        d.setValue( new BigDecimal( String.valueOf( v ) ) );
        n.getNodeData().setDate( d );
    }

    private static void prop( final PhylogenyNode n, final String ref, final double v ) {
        propText( n, ref, String.valueOf( v ) );
    }

    private static void propText( final PhylogenyNode n, final String ref, final String v ) {
        if ( n.getNodeData().getProperties() == null ) {
            n.getNodeData().setProperties( new PropertiesList() );
        }
        n.getNodeData().getProperties().addProperty( new Property( ref, v, "", "xsd:string", AppliesTo.NODE ) );
    }

    private static PhylogenyNode named( final String name, final double date, final double length ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        date( n, date );
        n.setDistanceToParent( length );
        return n;
    }

    private static Phylogeny tree( final PhylogenyNode root ) {
        final Phylogeny p = new Phylogeny();
        p.setRoot( root );
        p.externalNodesHaveChanged();
        return p;
    }

    /** root(date 10) -> a(date 6), b(date 2); ages before present, as BEAST writes them. The file states a length of
     *  4.5 for a -- NOT the gap between its dates, 4 -- and of 8 for b. No rate, no divergence. */
    private static Phylogeny twoTips() {
        final PhylogenyNode root = new PhylogenyNode();
        date( root, 10 );
        root.addAsChild( named( "a", 6, 4.5 ) );
        root.addAsChild( named( "b", 2, 8 ) );
        return tree( root );
    }

    /** {@link #twoTips()} with "a" rated 0.005 and "b" stating {@code b_rate} as written (null: no rate at all). */
    private static Phylogeny ratedPair( final String b_rate ) {
        final Phylogeny p = twoTips();
        prop( node( p, "a" ), RATE, 0.005 );
        if ( b_rate != null ) {
            propText( node( p, "b" ), RATE, b_rate );
        }
        return p;
    }

    /** root(10) -> x(6) -> [ t1(2), t2(1) ]; root -> y(5) -> t3(0), as ages: FIVE branches (x, t1, t2, y, t3), every
     *  node dated, every branch rated 0.002 and stating as its length the gap between its dates. */
    private static Phylogeny fiveBranches() {
        final PhylogenyNode root = new PhylogenyNode();
        date( root, 10 );
        final PhylogenyNode x = named( "x", 6, 4 );
        final PhylogenyNode y = named( "y", 5, 5 );
        final PhylogenyNode t1 = named( "t1", 2, 4 );
        final PhylogenyNode t2 = named( "t2", 1, 5 );
        final PhylogenyNode t3 = named( "t3", 0, 5 );
        root.addAsChild( x );
        root.addAsChild( y );
        x.addAsChild( t1 );
        x.addAsChild( t2 );
        y.addAsChild( t3 );
        for( final PhylogenyNode n : new PhylogenyNode[] { x, y, t1, t2, t3 } ) {
            prop( n, RATE, 0.002 );
        }
        return tree( root );
    }

    /** {@link #fiveBranches()} recording a divergence on every node but {@code without} ("root" for the root),
     *  0.004 deeper along every branch, and stating no rate. */
    private static Phylogeny recording( final String without ) {
        final Phylogeny p = fiveBranches();
        for( final PhylogenyNodeIterator it = p.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            n.getNodeData().setProperties( null );
            if ( !( n.isRoot() ? "root" : n.getName() ).equals( without ) ) {
                prop( n, DIV, n.isRoot() ? 0.0 : ( n.getParent().isRoot() ? 0.004 : 0.008 ) );
            }
        }
        return p;
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

    private static double lengthOf( final Phylogeny p, final String name ) {
        return node( p, name ).getDistanceToParent();
    }

    private static boolean plainZero( final double d ) {
        return Double.doubleToRawLongBits( d ) == 0L;
    }

    private static boolean fail( final String m ) {
        System.out.println( "  [BranchLengthLayoutTest] " + m );
        return false;
    }

    private static boolean eq( final double a, final double b ) {
        return Math.abs( a - b ) < 1e-9;
    }

    private static boolean eq( final Double a, final double b ) {
        return ( a != null ) && eq( a.doubleValue(), b );
    }

    public static boolean test() {
        try {
            return twoShapesOk() && timeIsWhatTheTreeHadOk() && everyBranchOk() && whatIsARateOk() && depthOk()
                    && negativeOk() && noLengthKeptOk() && directionOk() && arrivalOk() && whatIsARecordedDivergenceOk()
                    && keptAcrossADeleteOk() && aPartOk();
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
    }

    // ---- a length kept is kept WITH ITS PARENT: the tree may lose nodes while divergence is on screen ----------

    private static boolean keptAcrossADeleteOk() {
        final Phylogeny p = fiveBranches();
        node( p, "x" ).setDistanceToParent( 4.25 ); // neither is its date gap (4), so a sum is not a gap
        node( p, "t1" ).setDistanceToParent( 4.5 );
        final TimeLengths kept = TimeLengths.onScreen( p );
        PhylogenyMethods.removeNode( node( p, "x" ), p ); // as Delete Node does: t1 and t2 now hang from the root
        if ( ( node( p, "t1" ).getParent() != p.getRoot() ) || ( node( p, "t2" ).getParent() != p.getRoot() ) ) {
            return fail( "kept across a delete: fixture -- t1 and t2 must hang from the root" );
        }
        if ( !eq( kept.of( node( p, "t1" ) ), 8.75 ) || !eq( kept.of( node( p, "t2" ) ), 9.25 ) ) {
            return fail( "a node removed since: its children span its kept length and their own, t1 8.75 t2 9.25; got"
                    + " t1=" + kept.of( node( p, "t1" ) ) + " t2=" + kept.of( node( p, "t2" ) ) );
        }
        // divergence PIECE BY PIECE: t1's own piece at its rate, x's piece at the rate x stated (x rated 0.004)
        final Phylogeny pr = fiveBranches();
        node( pr, "x" ).setDistanceToParent( 4.25 );
        node( pr, "t1" ).setDistanceToParent( 4.5 );
        node( pr, "x" ).getNodeData().setProperties( null );
        propText( node( pr, "x" ), RATE, "0.004" );
        final TimeLengths kept_pr = TimeLengths.onScreen( pr );
        PhylogenyMethods.removeNode( node( pr, "x" ), pr );
        if ( !eq( kept_pr.divergence( node( pr, "t1" ), 0.002 ), ( 0.002 * 4.5 ) + ( 0.004 * 4.25 ) ) ) {
            return fail( "a merged branch's divergence is its pieces', each at its own rate: "
                    + ( ( 0.002 * 4.5 ) + ( 0.004 * 4.25 ) ) + "; got " + kept_pr.divergence( node( pr, "t1" ), 0.002 ) );
        }
        BranchLengthLayout.applyDivergence( pr, kept_pr );
        if ( !eq( lengthOf( pr, "t1" ), ( 0.002 * 4.5 ) + ( 0.004 * 4.25 ) ) ) {
            return fail( "Div lays a merged branch out piece by piece; got " + lengthOf( pr, "t1" ) );
        }
        // back in time, the pieces are remembered while the branch on screen is still their sum...
        BranchLengthLayout.applyTime( pr, kept_pr );
        final TimeLengths refreshed = kept_pr.refreshedBy( pr );
        if ( !eq( refreshed.divergence( node( pr, "t1" ), 0.002 ), ( 0.002 * 4.5 ) + ( 0.004 * 4.25 ) ) ) {
            return fail( "refreshed while the branch on screen is the sum of its pieces: the pieces stay" );
        }
        // ...and given up when the branch is edited to something else: then it is one piece, at its own rate
        node( pr, "t1" ).setDistanceToParent( 9 );
        if ( !eq( kept_pr.refreshedBy( pr ).divergence( node( pr, "t1" ), 0.002 ), 0.002 * 9 )
                || !eq( kept_pr.refreshedBy( pr ).of( node( pr, "t1" ) ), 9 ) ) {
            return fail( "a merged branch edited to another length is one piece of that length" );
        }
        // a NEGATIVE piece: clamped by itself, max(0, span) x rate -- not the total
        final Phylogeny neg = fiveBranches();
        node( neg, "x" ).setDistanceToParent( -0.25 );
        node( neg, "t1" ).setDistanceToParent( 4.5 );
        final TimeLengths kept_neg = TimeLengths.onScreen( neg );
        PhylogenyMethods.removeNode( node( neg, "x" ), neg );
        if ( !eq( kept_neg.divergence( node( neg, "t1" ), 0.002 ), 0.002 * 4.5 ) ) {
            return fail( "a negative piece adds nothing, the others count in full: " + ( 0.002 * 4.5 ) + "; got "
                    + kept_neg.divergence( node( neg, "t1" ), 0.002 ) + " (clamping the total gives " + ( 0.002 * 4.25 ) + ")" );
        }
        // the tree code's merge dropped the negative piece: t1 is 4.5 on screen; settled, it is the signed 4.25
        if ( !eq( lengthOf( neg, "t1" ), 4.5 ) ) {
            return fail( "settle: fixture -- the tree code's merge drops a negative length, t1 4.5; got " + lengthOf( neg, "t1" ) );
        }
        if ( ( BranchLengthLayout.settleMergedBranches( neg, kept_neg ) != 2 ) || !eq( lengthOf( neg, "t1" ), 4.25 )
                || !eq( lengthOf( neg, "t2" ), 4.75 ) ) {
            return fail( "both branches that took in x are given the signed sum of their pieces: t1 4.25, t2 4.75; got t1="
                    + lengthOf( neg, "t1" ) + " t2=" + lengthOf( neg, "t2" ) );
        }
        if ( BranchLengthLayout.settleMergedBranches( neg, kept_neg ) != 0 ) {
            return fail( "settling twice changes nothing" );
        }
        // a merged branch edited to another length is not touched
        node( neg, "t1" ).setDistanceToParent( 7 );
        if ( ( BranchLengthLayout.settleMergedBranches( neg, kept_neg ) != 0 ) || !eq( lengthOf( neg, "t1" ), 7 ) ) {
            return fail( "a branch that is no longer the tree code's merge of its pieces is left as it is" );
        }
        // a branch that states no length on screen has none, whatever was remembered for it
        node( pr, "t2" ).setDistanceToParent( PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT );
        if ( kept_pr.refreshedBy( pr ).of( node( pr, "t2" ) ) != null ) {
            return fail( "a branch stating no length on screen has none after a refresh" );
        }
        // a removed node that stated no rate: its piece at the rate of the branch it joined
        final Phylogeny unrated_x = fiveBranches();
        node( unrated_x, "x" ).getNodeData().setProperties( null );
        final TimeLengths kept_ux = TimeLengths.onScreen( unrated_x );
        PhylogenyMethods.removeNode( node( unrated_x, "x" ), unrated_x );
        if ( !eq( kept_ux.divergence( node( unrated_x, "t1" ), 0.002 ), 0.002 * 8 ) ) {
            return fail( "a removed piece with no rate goes at the rate of its branch" );
        }
        if ( !eq( kept.of( node( p, "t3" ) ), 5 ) ) {
            return fail( "a branch whose parent is the one it was kept to is its own kept length" );
        }
        // a node added since (a fresh id) has nothing kept, and is completed from its dates; the kept set is not
        // touched by completing it
        final PhylogenyNode z = named( "z", 3, 99 );
        p.getRoot().addAsChild( z );
        if ( kept.of( z ) != null ) {
            return fail( "a node added since has no kept length" );
        }
        if ( !eq( kept.completedFromDates( false, p ).of( z ), 7 ) || ( kept.of( z ) != null ) ) {
            return fail( "a node added since is completed from its dates (7), in a NEW set" );
        }
        // a branch MOVED elsewhere: its way up to its new parent is not all kept, so its dates decide
        final Phylogeny q = fiveBranches();
        node( q, "t3" ).setDistanceToParent( 5.5 );
        final TimeLengths kept_q = TimeLengths.onScreen( q );
        final PhylogenyNode t3 = node( q, "t3" );
        node( q, "y" ).removeChildNode( t3 );
        node( q, "x" ).addAsChild( t3 );
        if ( kept_q.of( t3 ) != null ) {
            return fail( "a branch moved under a node that is not above its old place has no kept length" );
        }
        // over(): the other set's lengths for the given tree's branches, these for the rest
        final Phylogeny r = fiveBranches();
        final TimeLengths mine = TimeLengths.onScreen( r );
        node( r, "t1" ).setDistanceToParent( 7 );
        node( r, "t3" ).setDistanceToParent( 8 );
        final TimeLengths theirs = TimeLengths.onScreen( r );
        final PhylogenyNode y_clade = node( r, "y" );
        final Phylogeny just_y = r.copy( y_clade );
        final TimeLengths merged = mine.over( theirs, just_y );
        if ( !eq( merged.of( node( r, "t3" ) ), 8 ) || !eq( merged.of( node( r, "t1" ) ), 4 )
                || !eq( mine.of( node( r, "t3" ) ), 5 ) ) {
            return fail( "over(): t3 (in the tree given) theirs, 8; t1 (outside it) mine, 4; mine untouched; got t3="
                    + merged.of( node( r, "t3" ) ) + " t1=" + merged.of( node( r, "t1" ) ) );
        }
        return true;
    }

    /** A CLADE of a tree in divergence is laid out without a depth of its own: identical sequences have none. */
    private static boolean aPartOk() {
        final Phylogeny flat = fiveBranches();
        for( final PhylogenyNodeIterator it = flat.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isRoot() ) {
                n.getNodeData().setProperties( null );
                propText( n, RATE, "0" );
            }
        }
        if ( BranchLengthLayout.isApplicable( flat ) ) {
            return fail( "a part: fixture -- rates of 0 everywhere, no depth, not offered by itself" );
        }
        final TimeLengths time = TimeLengths.onScreen( flat );
        BranchLengthLayout.applyDivergence( flat, time );
        if ( lengthOf( flat, "t2" ) != 5 ) {
            return fail( "a tree by itself without depth is left as it is" );
        }
        BranchLengthLayout.applyDivergenceToPart( flat, time );
        for( final String n : new String[] { "x", "y", "t1", "t2", "t3" } ) {
            if ( !plainZero( lengthOf( flat, n ) ) ) {
                return fail( "a part of a tree in divergence is laid out in divergence, depth or none: " + n + "="
                        + lengthOf( flat, n ) );
            }
        }
        final Phylogeny unrated = fiveBranches();
        node( unrated, "t2" ).getNodeData().setProperties( null );
        BranchLengthLayout.applyDivergenceToPart( unrated, TimeLengths.onScreen( unrated ) );
        if ( lengthOf( unrated, "t1" ) != 4 ) {
            return fail( "a part that cannot state every branch is left as it is" );
        }
        return true;
    }

    // ---- the two shapes, and the tree that has only one signal ------------------------------------------------

    private static boolean twoShapesOk() {
        final Phylogeny dates_only = twoTips();
        if ( BranchLengthLayout.divergenceSource( dates_only ) != NONE ) {
            return fail( "dates alone are not a divergence source" );
        }
        if ( BranchLengthLayout.isApplicable( dates_only ) ) {
            return fail( "a tree with only a time signal must NOT be offered the toggle" );
        }
        // BEAST shape: a rate on every branch, so divergence is DERIVED as rate x the length in time
        final Phylogeny beast = ratedPair( "0.005" );
        if ( BranchLengthLayout.divergenceSource( beast ) != CLOCK ) {
            return fail( "a per-branch rate must be recognised as a CLOCK_RATE divergence source" );
        }
        if ( !BranchLengthLayout.isApplicable( beast ) ) {
            return fail( "dates + a clock rate must offer the toggle (this is the BEAST case)" );
        }
        // Auspice shape: a RECORDED cumulative divergence wins over any rate, and is differenced
        final Phylogeny auspice = ratedPair( "99" ); // the rates must be IGNORED: recorded beats derived
        prop( auspice.getRoot(), DIV, 0.0 );
        prop( node( auspice, "a" ), DIV, 0.001 );
        prop( node( auspice, "b" ), DIV, 0.002 );
        if ( BranchLengthLayout.divergenceSource( auspice ) != STORED ) {
            return fail( "a recorded cumulative divergence must win over a derivable one" );
        }
        BranchLengthLayout.applyDivergence( auspice, TimeLengths.onScreen( auspice ) );
        if ( !eq( lengthOf( auspice, "a" ), 0.001 ) || !eq( lengthOf( auspice, "b" ), 0.002 ) ) {
            return fail( "recorded divergence is the successive difference; got a=" + lengthOf( auspice, "a" )
                    + " b=" + lengthOf( auspice, "b" ) );
        }
        // the label distinguishes a derived number from a recorded one
        if ( !"Divergence".equals( BranchLengthLayout.label( BranchLengthLayout.MODE.DIVERGENCE, STORED ) ) ) {
            return fail( "a recorded divergence is just 'Divergence'" );
        }
        if ( !"Divergence (from clock rate)".equals( BranchLengthLayout.label( BranchLengthLayout.MODE.DIVERGENCE, CLOCK ) ) ) {
            return fail( "a DERIVED divergence must say so, so a reader does not quote an inference as a measurement" );
        }
        return true;
    }

    // ---- TIME is the lengths the tree had: kept on leaving time, given back to the digit -----------------------

    private static boolean timeIsWhatTheTreeHadOk() {
        final Phylogeny beast = ratedPair( "0.005" );
        final TimeLengths kept = TimeLengths.onScreen( beast );
        BranchLengthLayout.applyDivergence( beast, kept );
        // rate x the length the FILE states for a, 4.5 -- not x the gap between its dates, 4
        if ( !eq( lengthOf( beast, "a" ), 0.0225 ) || !eq( lengthOf( beast, "b" ), 0.04 ) ) {
            return fail( "derived divergence must be rate x the length in time; got a=" + lengthOf( beast, "a" )
                    + " b=" + lengthOf( beast, "b" ) );
        }
        // laid out a second time from the same kept lengths (two trees of a tab share their nodes): the same
        BranchLengthLayout.applyDivergence( beast, kept );
        if ( !eq( lengthOf( beast, "a" ), 0.0225 ) || !eq( lengthOf( beast, "b" ), 0.04 ) ) {
            return fail( "laid out twice, divergence must not be taken for time; got a=" + lengthOf( beast, "a" )
                    + " b=" + lengthOf( beast, "b" ) );
        }
        BranchLengthLayout.applyTime( beast, kept );
        if ( ( lengthOf( beast, "a" ) != 4.5 ) || ( lengthOf( beast, "b" ) != 8 ) ) {
            return fail( "time gives back the lengths the tree had, to the digit; got a=" + lengthOf( beast, "a" )
                    + " b=" + lengthOf( beast, "b" ) );
        }
        // what is kept is kept by node id, which a COPY of the tree shares with it (an undo snapshot is one)
        final Phylogeny snapshot = beast.copy();
        BranchLengthLayout.applyDivergence( snapshot, kept );
        if ( !eq( lengthOf( snapshot, "a" ), 0.0225 ) || ( lengthOf( beast, "a" ) != 4.5 ) ) {
            return fail( "a copy is laid out from the lengths kept of its original, and alone; got copy a="
                    + lengthOf( snapshot, "a" ) + ", original a=" + lengthOf( beast, "a" ) );
        }
        BranchLengthLayout.applyTime( snapshot, kept );
        if ( lengthOf( snapshot, "a" ) != 4.5 ) {
            return fail( "...and back; got " + lengthOf( snapshot, "a" ) );
        }
        // kept from two trees, the one kept LAST counts
        final Phylogeny edited = beast.copy();
        node( edited, "a" ).setDistanceToParent( 4.75 );
        if ( TimeLengths.onScreen( beast, edited ).of( node( beast, "a" ) ).doubleValue() != 4.75 ) {
            return fail( "of two trees holding a node of the same id, the one kept last counts" );
        }
        return true;
    }

    // ---- THE JOINT RULE: both layouts can state EVERY branch ---------------------------------------------------

    private static boolean everyBranchOk() {
        final Phylogeny five = fiveBranches();
        if ( !BranchLengthLayout.isTimeDerivable( five ) || !BranchLengthLayout.isApplicable( five ) ) {
            return fail( "fixture: five branches, every node dated, every branch rated and long, must be offered" );
        }
        // time: a date on every node, the root included. One undated of five was a majority, and was offered once.
        for( final String undated : new String[] { "t1", "x", "root" } ) {
            final Phylogeny p = fiveBranches();
            ( "root".equals( undated ) ? p.getRoot() : node( p, undated ) ).getNodeData().setDate( null );
            if ( BranchLengthLayout.isTimeDerivable( p ) || BranchLengthLayout.isApplicable( p ) ) {
                return fail( "\"" + undated + "\" states no date: the time layout cannot state every branch" );
            }
            // ...and the tree stays in the layout it arrived in
            BranchLengthLayout.applyDivergence( p, TimeLengths.onScreen( p ) );
            if ( lengthOf( p, "t3" ) != 5 ) {
                return fail( "a tree that is not offered the switch must keep its lengths (\"" + undated
                        + "\" undated); t3 5 -> " + lengthOf( p, "t3" ) );
            }
        }
        // time: a LENGTH on every branch -- it is what the time layout shows, and what divergence is a rate times
        final Phylogeny no_length = fiveBranches();
        node( no_length, "t2" ).setDistanceToParent( org.forester.phylogeny.data.PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT );
        if ( BranchLengthLayout.isTimeDerivable( no_length ) || BranchLengthLayout.isApplicable( no_length ) ) {
            return fail( "t2 states no length: the time layout cannot state every branch" );
        }
        if ( TimeLengths.onScreen( no_length ).of( node( no_length, "t2" ) ) != null ) {
            return fail( "a length that is not stated is not kept" );
        }
        // a date or a length that is STATED is stated, whether zero or negative
        final Phylogeny zero_and_negative = fiveBranches();
        date( node( zero_and_negative, "t1" ), 0 );
        date( node( zero_and_negative, "t2" ), -3 );
        node( zero_and_negative, "t1" ).setDistanceToParent( 0 );
        node( zero_and_negative, "t2" ).setDistanceToParent( -2 );
        if ( !BranchLengthLayout.isApplicable( zero_and_negative ) ) {
            return fail( "dates of 0 and -3, lengths of 0 and -2, are stated: the tree is offered the switch" );
        }
        // a clock rate on EVERY branch: the pair differs by ONE rate
        final Phylogeny one_unrated = ratedPair( null );
        if ( !BranchLengthLayout.isTimeDerivable( one_unrated ) ) {
            return fail( "the unrated half of the pair must still have a time layout: the refusal is for the RATE" );
        }
        if ( ( BranchLengthLayout.divergenceSource( one_unrated ) != NONE ) || BranchLengthLayout.isApplicable( one_unrated ) ) {
            return fail( "one branch without a rate: no divergence source, not offered; got "
                    + BranchLengthLayout.divergenceSource( one_unrated ) );
        }
        BranchLengthLayout.applyDivergence( one_unrated, TimeLengths.onScreen( one_unrated ) );
        if ( ( lengthOf( one_unrated, "a" ) != 4.5 ) || ( lengthOf( one_unrated, "b" ) != 8 ) ) {
            return fail( "a tree with no divergence source must keep its lengths; got a=" + lengthOf( one_unrated, "a" )
                    + " b=" + lengthOf( one_unrated, "b" ) );
        }
        // the root has no branch: what it states as a rate is never asked for
        final Phylogeny bad_root = ratedPair( "0.0025" );
        propText( bad_root.getRoot(), RATE, "fast" );
        if ( BranchLengthLayout.divergenceSource( bad_root ) != CLOCK ) {
            return fail( "the root's rate must not be asked for: it has no branch" );
        }
        // a tree of one node has no branch at all: "every branch" of none is no layout
        final PhylogenyNode lone_root = new PhylogenyNode();
        date( lone_root, 10 );
        prop( lone_root, RATE, 0.005 );
        final Phylogeny lone = tree( lone_root );
        if ( ( BranchLengthLayout.divergenceSource( lone ) != NONE ) || BranchLengthLayout.isTimeDerivable( lone ) ) {
            return fail( "a tree without a branch has no divergence source and no time layout" );
        }
        // recorded divergence: every node states it, the root included; stated in part there is no divergence
        // layout, and no falling back on the rates
        final Phylogeny recorded = recording( "nobody" );
        if ( ( BranchLengthLayout.divergenceSource( recorded ) != STORED ) || !BranchLengthLayout.isApplicable( recorded ) ) {
            return fail( "fixture: a divergence recorded on every node must be offered" );
        }
        for( final String without : new String[] { "t2", "x", "root" } ) {
            final Phylogeny p = recording( without );
            for( final PhylogenyNodeIterator it = p.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( !n.isRoot() ) {
                    prop( n, RATE, 0.002 );
                }
            }
            if ( ( BranchLengthLayout.divergenceSource( p ) != NONE ) || BranchLengthLayout.isApplicable( p ) ) {
                return fail( "\"" + without + "\" records no divergence: no divergence layout, and the rates on every"
                        + " branch are not fallen back on; got " + BranchLengthLayout.divergenceSource( p ) );
            }
        }
        // ...a recorded divergence needs no rates at all
        if ( !BranchLengthLayout.isApplicable( recording( "nobody" ) ) ) {
            return fail( "a recorded divergence needs no rates" );
        }
        return true;
    }

    // ---- what counts as a rate: a plain decimal, finite, not negative (JOINT) -----------------------------------

    private static boolean whatIsARateOk() {
        for( final String good : new String[] { "0.0025", "1e-3", "1E-3", "2.5e-3", " 0.005 ", "12", "+0.5", ".5", "5.",
                "-0.0", "0", "0.0" } ) {
            if ( BranchLengthLayout.divergenceSource( ratedPair( good ) ) != CLOCK ) {
                return fail( "\"" + good + "\" is a rate" );
            }
        }
        for( final String bad : new String[] { "-0.001", "-1e-9", "fast", "NaN", "Infinity", "-Infinity", "", " ", "1e400",
                "0.005d", "0.005f", "1d", "0x1p-8", "0x10", "0.01abc", "1,5", "1_000", "0.5 per year", "1e", ".", "+" } ) {
            if ( BranchLengthLayout.divergenceSource( ratedPair( bad ) ) != NONE ) {
                return fail( "\"" + bad + "\" is no rate: the tree must have no divergence source" );
            }
        }
        // a rate written -0.0 is not negative, so it is a rate; its length must still be a plain 0
        final Phylogeny minus_zero = ratedPair( "-0.0" );
        BranchLengthLayout.applyDivergence( minus_zero, TimeLengths.onScreen( minus_zero ) );
        if ( !eq( lengthOf( minus_zero, "a" ), 0.0225 ) || !plainZero( lengthOf( minus_zero, "b" ) ) ) {
            return fail( "a rate of -0.0 must give a length of 0, not of -0.0 (it would be written \"-0.0\"); got a="
                    + lengthOf( minus_zero, "a" ) + " b=" + lengthOf( minus_zero, "b" ) );
        }
        return true;
    }

    // ---- what counts as a recorded divergence: a plain decimal, finite -- negative too (JOINT) -------------------

    /** {@link #recording} with t2's divergence written {@code text}. */
    private static Phylogeny recordingT2( final String text ) {
        final Phylogeny p = recording( "t2" );
        propText( node( p, "t2" ), DIV, text );
        return p;
    }

    private static boolean whatIsARecordedDivergenceOk() {
        for( final String good : new String[] { "0.008", "8e-3", "8E-3", " 0.008 ", "+0.008", ".008", "0", "0.0", "-0.0",
                "5.", "-0.001", "-1e-9" } ) {
            if ( BranchLengthLayout.divergenceSource( recordingT2( good ) ) != STORED ) {
                return fail( "\"" + good + "\" is a recorded divergence" );
            }
        }
        for( final String bad : new String[] { "0.008d", "0.008f", "1d", "0x1p-8", "0x10", "0.01abc", "1,5", "1_000",
                "NaN", "Infinity", "-Infinity", "", " ", "1e400", "1e", ".", "+", "abc" } ) {
            if ( BranchLengthLayout.divergenceSource( recordingT2( bad ) ) != NONE ) {
                return fail( "\"" + bad + "\" is no recorded divergence: the tree must have no divergence source" );
            }
        }
        // a NEGATIVE divergence is a value: the branch along which it falls is drawn at 0 (Christian, joint)
        final Phylogeny falls = recordingT2( "-0.001" );
        if ( !BranchLengthLayout.isApplicable( falls ) ) {
            return fail( "a tree recording a negative divergence on one node is offered the switch" );
        }
        BranchLengthLayout.applyDivergence( falls, TimeLengths.onScreen( falls ) );
        if ( !plainZero( lengthOf( falls, "t2" ) ) || !eq( lengthOf( falls, "t1" ), 0.004 ) ) {
            return fail( "a branch along which the recorded divergence falls is drawn at 0; got t2=" + lengthOf( falls, "t2" )
                    + " t1=" + lengthOf( falls, "t1" ) );
        }
        return true;
    }

    // ---- each picture has some depth, measured at the tips (JOINT) ---------------------------------------------

    private static boolean depthOk() {
        // every rate 0: divergence is 0 along every branch, a picture of nothing
        final Phylogeny all_zero = fiveBranches();
        for( final PhylogenyNodeIterator it = all_zero.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isRoot() ) {
                n.getNodeData().setProperties( null );
                propText( n, RATE, "0" );
            }
        }
        if ( BranchLengthLayout.divergenceSource( all_zero ) != CLOCK ) {
            return fail( "fixture: a rate of 0 on every branch is a rate on every branch" );
        }
        if ( BranchLengthLayout.isApplicable( all_zero ) ) {
            return fail( "divergence of 0 along every branch: not offered" );
        }
        BranchLengthLayout.applyDivergence( all_zero, TimeLengths.onScreen( all_zero ) );
        if ( lengthOf( all_zero, "t3" ) != 5 ) {
            return fail( "...and the tree keeps its lengths; t3 5 -> " + lengthOf( all_zero, "t3" ) );
        }
        // ONE branch with a rate above 0 is depth enough
        node( all_zero, "t2" ).getNodeData().setProperties( null );
        propText( node( all_zero, "t2" ), RATE, "0.002" );
        if ( !BranchLengthLayout.isApplicable( all_zero ) ) {
            return fail( "one branch with divergence above 0: offered" );
        }
        // rates above 0, but no length above 0 to multiply: no depth either
        final Phylogeny no_length_above_zero = fiveBranches();
        for( final PhylogenyNodeIterator it = no_length_above_zero.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isRoot() ) {
                n.setDistanceToParent( n.isExternal() ? -1 : 0 );
            }
        }
        if ( BranchLengthLayout.isApplicable( no_length_above_zero ) ) {
            return fail( "every length 0 or negative: divergence is 0 along every branch, not offered" );
        }
        // a recorded divergence that is the same on every node has no depth
        final Phylogeny flat = recording( "nobody" );
        for( final PhylogenyNodeIterator it = flat.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            n.getNodeData().setProperties( null );
            prop( n, DIV, 0.004 );
        }
        if ( ( BranchLengthLayout.divergenceSource( flat ) != STORED ) || BranchLengthLayout.isApplicable( flat ) ) {
            return fail( "a divergence recorded the same on every node has no depth: not offered" );
        }
        // TIME has depth when some TIP is dated differently from the root -- not when some branch spans time
        final Phylogeny tips_on_the_root = fiveBranches();
        for( final String tip : new String[] { "t1", "t2", "t3" } ) {
            date( node( tips_on_the_root, tip ), 10 );
        }
        if ( BranchLengthLayout.isApplicable( tips_on_the_root ) ) {
            return fail( "every tip on the root's date, the inner nodes not: the time picture has no depth at the tips" );
        }
        date( node( tips_on_the_root, "t3" ), 9.5 );
        if ( !BranchLengthLayout.isApplicable( tips_on_the_root ) ) {
            return fail( "one tip off the root's date: offered" );
        }
        // a CONSTANT rate draws the same picture at another scale: a necessary condition is not failed by it
        if ( !BranchLengthLayout.isApplicable( fiveBranches() ) ) {
            return fail( "a constant rate has depth in both pictures" );
        }
        return true;
    }

    // ---- a length may be negative: given back as stated in time, 0 in divergence -------------------------------

    private static boolean negativeOk() {
        final Phylogeny backwards = fiveBranches();
        node( backwards, "t1" ).setDistanceToParent( -1 ); // as a file states a node before its parent
        final TimeLengths kept = TimeLengths.onScreen( backwards );
        if ( !BranchLengthLayout.isApplicable( backwards, kept ) ) {
            return fail( "a negative length is stated: the tree is offered the switch" );
        }
        BranchLengthLayout.applyDivergence( backwards, kept );
        if ( !plainZero( lengthOf( backwards, "t1" ) ) ) {
            return fail( "divergence states 0 for a length that is negative; got " + lengthOf( backwards, "t1" ) );
        }
        if ( !eq( lengthOf( backwards, "t2" ), 5 * 0.002 ) ) {
            return fail( "...and rate x length for every other branch; got t2=" + lengthOf( backwards, "t2" ) );
        }
        BranchLengthLayout.applyTime( backwards, kept );
        if ( lengthOf( backwards, "t1" ) != -1 ) {
            return fail( "and back in time the length is -1 again; got " + lengthOf( backwards, "t1" ) );
        }
        // a recorded divergence that FALLS along a branch is 0 as well
        final Phylogeny falling = recording( "nobody" );
        node( falling, "t1" ).getNodeData().setProperties( null );
        prop( node( falling, "t1" ), DIV, 0.003 ); // its parent x records 0.004
        BranchLengthLayout.applyDivergence( falling, TimeLengths.onScreen( falling ) );
        if ( !plainZero( lengthOf( falling, "t1" ) ) || !eq( lengthOf( falling, "t2" ), 0.004 ) ) {
            return fail( "a recorded divergence falling from 0.004 to 0.003 is 0; got t1=" + lengthOf( falling, "t1" )
                    + " t2=" + lengthOf( falling, "t2" ) );
        }
        return true;
    }

    // ---- a branch NO length was kept for: the gap between its dates, signed; no dates either, 0 ----------------

    private static boolean noLengthKeptOk() {
        final TimeLengths nothing = new TimeLengths();
        // ages (largest at the root): t1 is dated OLDER than its parent x, 7 against 6
        final Phylogeny older = fiveBranches();
        date( node( older, "t1" ), 7 );
        for( final String n : new String[] { "x", "y", "t1", "t2", "t3" } ) {
            node( older, n ).setDistanceToParent( 99 ); // what divergence left behind
        }
        if ( nothing.of( node( older, "t1" ) ) != null ) {
            return fail( "fixture: nothing was kept" );
        }
        BranchLengthLayout.applyTime( older, nothing.completedFromDates( false, older ) );
        if ( !eq( lengthOf( older, "t1" ), -1 ) || !eq( lengthOf( older, "t2" ), 5 ) || !eq( lengthOf( older, "x" ), 4 ) ) {
            return fail( "no length kept: the gap between the dates, signed; got t1=" + lengthOf( older, "t1" ) + " t2="
                    + lengthOf( older, "t2" ) + " x=" + lengthOf( older, "x" ) );
        }
        // told the dates run the other way, every gap has the other sign
        BranchLengthLayout.applyTime( older, nothing.completedFromDates( true, older ) );
        if ( !eq( lengthOf( older, "t1" ), 1 ) || !eq( lengthOf( older, "t2" ), -5 ) ) {
            return fail( "the gap is signed the way it is TOLD the dates run; got t1=" + lengthOf( older, "t1" ) + " t2="
                    + lengthOf( older, "t2" ) );
        }
        // completing leaves what it completes as it was: nothing is kept that was not
        if ( nothing.of( node( older, "t2" ) ) != null ) {
            return fail( "completing from the dates must not change the lengths it was asked to complete" );
        }
        // calendar dates (increasing toward the tips): the same tree, mirrored about year 2000
        final Phylogeny calendar = fiveBranches();
        for( final PhylogenyNodeIterator it = calendar.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            date( n, 2010 - n.getNodeData().getDate().getValue().doubleValue() );
            n.setDistanceToParent( 99 );
        }
        date( node( calendar, "t1" ), 2003 ); // its parent x is dated 2004
        BranchLengthLayout.applyTime( calendar, nothing.completedFromDates( true, calendar ) );
        if ( !eq( lengthOf( calendar, "t1" ), -1 ) || !eq( lengthOf( calendar, "t2" ), 5 ) || !eq( lengthOf( calendar, "x" ), 4 ) ) {
            return fail( "no length kept, calendar dates: signed the same; got t1=" + lengthOf( calendar, "t1" )
                    + " t2=" + lengthOf( calendar, "t2" ) + " x=" + lengthOf( calendar, "x" ) );
        }
        // a length that WAS kept wins over the dates; the others still come from them
        final Phylogeny one_kept = fiveBranches();
        node( one_kept, "t2" ).setDistanceToParent( 5.25 );
        for( final String n : new String[] { "x", "y", "t1", "t3" } ) {
            node( one_kept, n ).setDistanceToParent( PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ); // stating none
        }
        final TimeLengths only_t2 = TimeLengths.onScreen( one_kept );
        for( final String n : new String[] { "x", "y", "t1", "t2", "t3" } ) {
            node( one_kept, n ).setDistanceToParent( 99 );
        }
        BranchLengthLayout.applyTime( one_kept, only_t2.completedFromDates( false, one_kept ) );
        if ( ( lengthOf( one_kept, "t2" ) != 5.25 ) || !eq( lengthOf( one_kept, "t1" ), 4 ) ) {
            return fail( "the kept length wins, the gap serves the rest; got t2=" + lengthOf( one_kept, "t2" ) + " t1="
                    + lengthOf( one_kept, "t1" ) );
        }
        // neither a kept length nor dates: no length in time, and laid out at 0 -- never the divergence it holds
        final Phylogeny undated = fiveBranches();
        node( undated, "t3" ).getNodeData().setDate( null );
        node( undated, "t3" ).setDistanceToParent( 0.01 );
        final TimeLengths completed = nothing.completedFromDates( false, undated );
        if ( ( completed.of( node( undated, "t3" ) ) != null ) || ( completed.of( node( undated, "t1" ) ) == null ) ) {
            return fail( "a branch with no dates cannot be completed from them; one with dates can" );
        }
        BranchLengthLayout.applyTime( undated, completed );
        if ( !plainZero( lengthOf( undated, "t3" ) ) ) {
            return fail( "no length in time: 0; got " + lengthOf( undated, "t3" ) );
        }
        // and the offer can be asked of a tree that shows divergence: its time lengths are the completed ones
        final Phylogeny showing_div = fiveBranches();
        final TimeLengths kept = TimeLengths.onScreen( showing_div );
        BranchLengthLayout.applyDivergence( showing_div, kept );
        if ( !BranchLengthLayout.isApplicable( showing_div, nothing.completedFromDates( false, showing_div ) ) ) {
            return fail( "a tree showing divergence, nothing kept: offered on the gaps between its dates" );
        }
        return true;
    }

    /** A two-node tree whose one branch is a copy of {@code n}'s: same id, same length. */
    // ---- which way the dates run: the majority of the pairs that differ; none, or a tie, reads as ages ----------

    private static boolean directionOk() {
        if ( BranchLengthLayout.datesIncreaseTowardTips( fiveBranches() ) ) {
            return fail( "five spans run down: the dates are ages" );
        }
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
        node( level, "a" ).setDistanceToParent( 99 );
        final TimeLengths level_gaps = new TimeLengths().completedFromDates( false, level );
        if ( ( level_gaps.of( node( level, "a" ) ) == null ) || !plainZero( level_gaps.of( node( level, "a" ) ).doubleValue() ) ) {
            return fail( "a gap of nothing is a plain 0, never -0.0 and never missing; got " + level_gaps.of( node( level, "a" ) ) );
        }
        // a pair with EQUAL dates does not vote. Real builds are full of them (1405 of the 9205 branches of one
        // Nextstrain tree); counted as running down, four of them would outvote the one pair that runs up here,
        // and a calendar tree would be read as ages, every gap with the wrong sign
        final Phylogeny mostly_level = fiveBranches();
        for( final PhylogenyNodeIterator it = mostly_level.iteratorPreorder(); it.hasNext(); ) {
            date( it.next(), 2000 );
        }
        date( node( mostly_level, "t2" ), 2001 );
        if ( !BranchLengthLayout.datesIncreaseTowardTips( mostly_level ) ) {
            return fail( "one pair runs up, four are level: the dates increase toward the tips" );
        }
        BranchLengthLayout.applyTime( mostly_level, new TimeLengths()
                .completedFromDates( BranchLengthLayout.datesIncreaseTowardTips( mostly_level ), mostly_level ) );
        if ( !eq( lengthOf( mostly_level, "t2" ), 1 ) || !plainZero( lengthOf( mostly_level, "t1" ) ) ) {
            return fail( "...and the one gap is +1; got t2=" + lengthOf( mostly_level, "t2" ) + " t1="
                    + lengthOf( mostly_level, "t1" ) );
        }
        return true;
    }

    // ---- the layout a tree ARRIVES in: asked once, of a tree both layouts can state ----------------------------

    private static boolean arrivalOk() {
        // a recording tree whose lengths are the gaps between its dates is showing time
        final Phylogeny recording = recording( "nobody" );
        if ( BranchLengthLayout.arrivesShowingDivergence( recording ) ) {
            return fail( "lengths that are the date gaps: the tree arrives showing time" );
        }
        // ...laid out by divergence and handed over like that (a tree saved from the Div view): showing divergence
        final Phylogeny saved_in_div = recording( "nobody" );
        BranchLengthLayout.applyDivergence( saved_in_div, TimeLengths.onScreen( saved_in_div ) );
        if ( !eq( lengthOf( saved_in_div, "t1" ), 0.004 ) ) {
            return fail( "fixture: the tree must BE in divergence lengths" );
        }
        if ( !BranchLengthLayout.arrivesShowingDivergence( saved_in_div ) ) {
            return fail( "lengths that are its divergence: the tree arrives showing divergence" );
        }
        // its time is what its dates say, its divergence what it records -- and back
        final TimeLengths from_dates = new TimeLengths().completedFromDates( false, saved_in_div );
        BranchLengthLayout.applyTime( saved_in_div, from_dates );
        if ( !eq( lengthOf( saved_in_div, "t1" ), 4 ) || !eq( lengthOf( saved_in_div, "y" ), 5 ) ) {
            return fail( "arrived in divergence: time is laid out from the dates; got t1=" + lengthOf( saved_in_div, "t1" )
                    + " y=" + lengthOf( saved_in_div, "y" ) );
        }
        // 19 branches in 20 (JOINT): of five, ONE off its date gap is one too many; to a millionth it is not off
        final Phylogeny one_off = recording( "nobody" );
        node( one_off, "t3" ).setDistanceToParent( 5.01 );
        if ( !BranchLengthLayout.arrivesShowingDivergence( one_off ) ) {
            return fail( "four lengths of five on their date gap is under 19 in 20: not showing time" );
        }
        node( one_off, "t3" ).setDistanceToParent( 5.000001 );
        if ( BranchLengthLayout.arrivesShowingDivergence( one_off ) ) {
            return fail( "a length within a millionth of its date gap IS its date gap" );
        }
        // the gap is SIGNED: a node stated 1 before its parent, dated 1 before it, is on its gap
        final Phylogeny backwards = recording( "nobody" );
        date( node( backwards, "t1" ), 7 ); // older than its parent x at 6
        node( backwards, "t1" ).setDistanceToParent( -1 );
        if ( BranchLengthLayout.arrivesShowingDivergence( backwards ) ) {
            return fail( "a negative length on a gap that runs backwards is its date gap" );
        }
        node( backwards, "t1" ).setDistanceToParent( 1 );
        if ( !BranchLengthLayout.arrivesShowingDivergence( backwards ) ) {
            return fail( "...and +1 on a gap of -1 is not" );
        }
        // a CLOCK tree's lengths are not its date gaps even in time, so that is not asked of it: it shows time...
        final Phylogeny clock = ratedPair( "0.005" ); // a is 4.5 long, its dates 4 apart
        if ( BranchLengthLayout.arrivesShowingDivergence( clock ) ) {
            return fail( "a clock tree whose lengths are not its date gaps is showing time all the same" );
        }
        // ...unless its file says its lengths are divergence, by the unit the divergence layout stamps
        clock.setDistanceUnit( BranchLengthLayout.distanceUnit( BranchLengthLayout.MODE.DIVERGENCE, "year" ) );
        if ( !BranchLengthLayout.arrivesShowingDivergence( clock ) ) {
            return fail( "a clock tree whose unit says subs/site arrives showing divergence" );
        }
        clock.setDistanceUnit( "year" );
        if ( BranchLengthLayout.arrivesShowingDivergence( clock ) ) {
            return fail( "a clock tree in years arrives showing time" );
        }
        // a tree that cannot be laid out both ways is not asked: it shows what its file states
        final Phylogeny one_unrated = ratedPair( null );
        one_unrated.setDistanceUnit( "subs/site" );
        final Phylogeny undated = recording( "nobody" );
        node( undated, "t2" ).getNodeData().setDate( null );
        node( undated, "t2" ).setDistanceToParent( 0.004 );
        if ( BranchLengthLayout.arrivesShowingDivergence( one_unrated ) || BranchLengthLayout.arrivesShowingDivergence( undated )
                || BranchLengthLayout.arrivesShowingDivergence( twoTips() ) ) {
            return fail( "a tree that has not both layouts is not asked which it arrives in" );
        }
        return true;
    }
}
