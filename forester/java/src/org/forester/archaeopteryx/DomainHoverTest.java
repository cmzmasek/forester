// The rollover for a protein domain: what it says, and where it decides a domain is.
//
// See TreePanel.domainAt, RenderableDomainArchitecture.domainAtX and NodeHoverText.domainRows.

package org.forester.archaeopteryx;

import java.awt.GraphicsEnvironment;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.phylogeny.data.RenderableDomainArchitecture;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;

public class DomainHoverTest {

    private static String rowsOf( final List<NodeHoverText.Row> rows ) {
        final StringBuilder sb = new StringBuilder();
        for( final NodeHoverText.Row r : rows ) {
            sb.append( r ).append( " | " );
        }
        return sb.toString();
    }

    /**
     * FIVE tips, architectures on two of them: a strong SH3 with a Pfam accession, a strong Pkinase, and a
     * WEAK domain above the e-value threshold that must never be drawn or hit.
     * <p>
     * Five rather than two on purpose. With two tips the circular layout puts them at near-symmetric angles,
     * and a hit-test that rotates the point the WRONG WAY still lands on a box -- that mutation survived
     * until the fixture had spokes whose angles a sign flip actually moves.
     */
    private static Phylogeny fixture() throws Exception {
        final Phylogeny p = Phylogeny
                .createInstanceFromNhxString( "(prot_a:0.2,prot_b:0.3,prot_c:0.25,prot_d:0.35,prot_e:0.15)" );
        for( final int which : new int[] { 0, 2 } ) {
            final List<PhylogenyData> ds = new ArrayList<PhylogenyData>();
            ds.add( new ProteinDomain( "SH3", 10, 60, "PF00018", 1e-6 ) );
            ds.add( new ProteinDomain( "Pkinase", 185, 445, "PF00069", 1e-30 ) );
            ds.add( new ProteinDomain( "Weak", 460, 490, 1.0 ) ); // above the 1e-3 threshold: not drawn
            final Sequence seq = new Sequence();
            seq.setName( p.getExternalNodes().get( which ).getName() );
            seq.setDomainArchitecture( new DomainArchitecture( ds, 500 ) );
            p.getExternalNodes().get( which ).getNodeData().addSequence( seq );
        }
        return p;
    }

    /** Lay out and paint {@code tp} once at this size, offscreen. */
    private static void paintOnce( final TreePanel tp, final int w, final int h ) {
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        tp.printAll( g );
        g.dispose();
    }

    /** {@code n} tips, every one carrying the same two-domain architecture -- a tree dense enough that
     *  "Auto-hide Labels" starts dropping tips. Every tip carries one, so any tip the rollover refuses was
     *  refused because it was not DRAWN, never because it had nothing to draw. */
    private static Phylogeny denseFixture( final int n ) throws Exception {
        final StringBuilder sb = new StringBuilder( "(" );
        for( int i = 0; i < n; ++i ) {
            sb.append( i > 0 ? "," : "" ).append( "prot_" ).append( i ).append( ":0.2" );
        }
        final Phylogeny p = Phylogeny.createInstanceFromNhxString( sb.append( ")" ).toString() );
        for( final PhylogenyNode tip : p.getExternalNodes() ) {
            final List<PhylogenyData> ds = new ArrayList<PhylogenyData>();
            ds.add( new ProteinDomain( "SH3", 10, 60, "PF00018", 1e-6 ) );
            ds.add( new ProteinDomain( "Pkinase", 185, 445, "PF00069", 1e-30 ) );
            final Sequence seq = new Sequence();
            seq.setName( tip.getName() );
            seq.setDomainArchitecture( new DomainArchitecture( ds, 500 ) );
            tip.getNodeData().addSequence( seq );
        }
        return p;
    }

    /**
     * Repaints, then requires that the pixels the hit-test calls a domain are pixels something was DRAWN on.
     *
     * Two things this deliberately does NOT do, both learned the hard way. It does not match a domain's
     * nominal colour: a box is shaded, and {@code colorFor} hands out palette entries in first-seen order, so
     * the colour is not a stable label. And it does not take a corner pixel as the background: a standalone
     * run picks up the developer's saved preferences, which may be a DARK theme, and then "not the corner
     * colour" means nothing. The background is the image's most common colour, which is the panel's whatever
     * the theme -- and what remains is the question that matters: if the track's origin is wrong, the
     * hit-test answers over empty canvas.
     */
    private static boolean alignedAtOrigin( final TreePanel tp, final String what, final int w, final int h ) {
        return alignedAtOrigin( tp, what, w, h, 1 );
    }

    /** As above, sampling every {@code step}-th pixel -- a dense tree makes a full sweep needlessly slow. */
    private static boolean alignedAtOrigin( final TreePanel tp, final String what, final int w, final int h,
                                            final int step ) {
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        tp.printAll( g );
        g.dispose();
        final java.util.Map<Integer, Integer> counts = new java.util.HashMap<Integer, Integer>();
        for( int y = 0; y < h; y += 2 ) {
            for( int x = 0; x < w; x += 2 ) {
                final Integer k = Integer.valueOf( img.getRGB( x, y ) & 0xFFFFFF );
                counts.merge( k, Integer.valueOf( 1 ), Integer::sum );
            }
        }
        int background = 0;
        int most = -1;
        for( final java.util.Map.Entry<Integer, Integer> e : counts.entrySet() ) {
            if ( e.getValue().intValue() > most ) {
                most = e.getValue().intValue();
                background = e.getKey().intValue();
            }
        }
        int claimed = 0;
        int on_ink = 0;
        String first_bad = null;
        for( int y = 0; y < h; y += step ) {
            for( int x = 0; x < w; x += step ) {
                final TreePanel.DomainHit hit = tp.domainAt( x, y );
                if ( hit == null ) {
                    continue;
                }
                ++claimed;
                if ( ( img.getRGB( x, y ) & 0xFFFFFF ) != background ) {
                    ++on_ink;
                }
                else if ( first_bad == null ) {
                    first_bad = "(" + x + "," + y + ") is empty canvas but hit-tests as "
                            + hit.domain().getName();
                }
            }
        }
        if ( claimed < 50 ) {
            System.out.println( what + ": the hit-test claims almost nothing (" + claimed
                    + " px), so alignment cannot be checked" );
            return false;
        }
        if ( on_ink < ( ( claimed * 9 ) / 10 ) ) {
            System.out.println( what + ": the hit-test answers over empty canvas for " + ( claimed - on_ink )
                    + " of " + claimed + " pixels it claims. " + first_bad );
            return false;
        }
        return true;
    }

    public static boolean test() {
        try {
            // (1) what the card says
            final PhylogenyNode tip = fixture().getExternalNodes().get( 0 );
            final ProteinDomain sh3 = new ProteinDomain( "SH3", 10, 60, "PF00018", 1e-6 );
            final List<NodeHoverText.Row> rows = NodeHoverText.domainRows( tip, sh3, 500 );
            final String s = rowsOf( rows );
            for( final String expected : new String[] { "Domain: SH3", "E-value: 1E-6", "Residues: 10–60",
                                                        "(51 aa)", "Protein length: 500 aa", "Tip: prot_a",
                                                        "Accession: PF00018" } ) {
                if ( !s.contains( expected ) ) {
                    System.out.println( "the domain card is missing [" + expected + "]: " + s );
                    return false;
                }
            }
            // the domain's NAME is data, so it must be a value row -- a heading is drawn small-caps
            // uppercase, which would misrepresent a case-sensitive name like "Bcl-2" or "wnt"
            for( final NodeHoverText.Row r : rows ) {
                if ( r.isHeading() ) {
                    System.out.println( "the domain card should have no heading rows: " + r );
                    return false;
                }
            }
            // (2) deliberate non-behaviour: no accession, no Pfam row; no e-value, no E-value row
            final String bare = rowsOf( NodeHoverText
                    .domainRows( tip, new ProteinDomain( "Bcl-2", 90, 188 ), 0 ) );
            if ( bare.contains( "Accession" ) || bare.contains( "Pfam entry" ) ) {
                System.out.println( "a domain with no accession must not claim one: " + bare );
                return false;
            }
            // an accession from ANOTHER source is shown as what it is, and does NOT earn the Pfam entry
            // link -- calling a SMART accession "Pfam" would be wrong the same way linking a name is
            final String smart = rowsOf( NodeHoverText
                    .domainRows( tip, new ProteinDomain( "CARD", 6, 90, "SM00114", 1e-9 ), 400 ) );
            if ( !smart.contains( "Accession: SM00114" ) ) {
                System.out.println( "a non-Pfam accession should still be shown: " + smart );
                return false;
            }
            if ( smart.contains( "the Pfam entry" ) ) {
                System.out.println( "a non-Pfam accession must not promise a Pfam entry: " + smart );
                return false;
            }
            if ( TreePanel.domainUrl( new ProteinDomain( "CARD", 6, 90, "SM00114", 1e-9 ) )
                    .contains( "/entry/pfam/" ) ) {
                System.out.println( "a non-Pfam accession must not be linked as one" );
                return false;
            }
            if ( bare.contains( "E-value" ) ) {
                System.out.println( "a domain with no e-value must not print one: " + bare );
                return false;
            }
            if ( bare.contains( "Protein length" ) ) {
                System.out.println( "an unknown protein length must not print as 0: " + bare );
                return false;
            }
            if ( !bare.contains( "Domain: Bcl-2" ) || !bare.contains( "Tip: prot_a" ) ) {
                System.out.println( "the name and the tip are always there: " + bare );
                return false;
            }
            // (3) the accession comes from the phyloXML id, and only when it looks like one
            if ( !"PF00452".equals( NodeHoverText
                    .pfamAccession( new ProteinDomain( "Bcl-2", 1, 2, "PF00452", 1e-9 ) ) ) ) {
                System.out.println( "a Pfam accession in the id was not found" );
                return false;
            }
            if ( !"PF00452".equals( NodeHoverText
                    .pfamAccession( new ProteinDomain( "x", 1, 2, "PF00452.25", 1e-9 ) ) ) ) {
                System.out.println( "a versioned accession (PF00452.25) should still yield PF00452" );
                return false;
            }
            // a six-digit accession must not be truncated to five
            if ( !"PF123456".equals( NodeHoverText
                    .pfamAccession( new ProteinDomain( "x", 1, 2, "PF123456", 1e-9 ) ) ) ) {
                System.out.println( "a longer accession was truncated: "
                        + NodeHoverText.pfamAccession( new ProteinDomain( "x", 1, 2, "PF123456", 1e-9 ) ) );
                return false;
            }
            if ( NodeHoverText.pfamAccession( new ProteinDomain( "x", 1, 2, "SM00109", 1e-9 ) ) != null ) {
                System.out.println( "a non-Pfam identifier must not be offered as a Pfam accession" );
                return false;
            }
            if ( NodeHoverText.pfamAccession( new ProteinDomain( "x", 1, 2 ) ) != null ) {
                System.out.println( "a domain with no id has no accession" );
                return false;
            }
            // (3b) the two link routes, chosen by the data. An ACCESSION addresses an InterPro entry; an
            // IDENTIFIER does not -- /entry/pfam/NB-ARC/ is a 404 while /entry/pfam/PF00931/ is a 200 -- so a
            // named-only domain must get a SEARCH, and must say so rather than promising an entry.
            final ProteinDomain with_acc = new ProteinDomain( "NB-ARC", 1, 2, "PF00931", 1e-9 );
            final ProteinDomain name_only = new ProteinDomain( "NB-ARC", 1, 2 );
            if ( !TreePanel.domainUrl( with_acc ).contains( "/entry/pfam/PF00931/" ) ) {
                System.out.println( "an accession should link to its entry: " + TreePanel.domainUrl( with_acc ) );
                return false;
            }
            if ( TreePanel.domainUrl( name_only ).contains( "/entry/pfam/" ) ) {
                System.out.println( "an identifier must NOT be linked as an accession (it 404s): "
                        + TreePanel.domainUrl( name_only ) );
                return false;
            }
            if ( !TreePanel.domainUrl( name_only ).contains( "/search/text/NB-ARC" ) ) {
                System.out.println( "a named-only domain should fall back to a search: "
                        + TreePanel.domainUrl( name_only ) );
                return false;
            }
            if ( TreePanel.domainUrl( new ProteinDomain( "", 1, 2 ) ) != null ) {
                System.out.println( "a domain with neither a name nor an accession links nowhere" );
                return false;
            }
            // A name needing encoding must be escaped FOR A PATH SEGMENT. Asserting only "no raw space"
            // was the classic check of a value against its own encoded self: URLEncoder's "Bcl-2+like"
            // contains no space and passes it, yet '+' in a path is a literal plus, so InterPro searched
            // for the string "Bcl-2+like" and found nothing.
            final String spaced = TreePanel.domainUrl( new ProteinDomain( "Bcl-2 like", 1, 2 ) );
            if ( !spaced.contains( "Bcl-2%20like" ) ) {
                System.out.println( "a space must be escaped as %20 for a path segment: " + spaced );
                return false;
            }
            if ( spaced.contains( "+" ) || spaced.contains( " " ) ) {
                System.out.println( "the search URL is form-encoded, not path-encoded: " + spaced );
                return false;
            }
            // and the card names the route, because the two are different promises
            final String acc_rows = rowsOf( NodeHoverText.domainRows( tip, with_acc, 100 ) );
            final String nam_rows = rowsOf( NodeHoverText.domainRows( tip, name_only, 100 ) );
            if ( !acc_rows.contains( "Click: the Pfam entry" ) ) {
                System.out.println( "an accession should promise the entry: " + acc_rows );
                return false;
            }
            if ( nam_rows.contains( "the Pfam entry" ) || !nam_rows.contains( "InterPro" ) ) {
                System.out.println( "a named-only domain should promise a lookup, not an entry: " + nam_rows );
                return false;
            }
            // (3c) the real-world shape, from the Archaeopteryx.js session's shared fixture: ONE protein
            // whose domains do not all carry an accession, so the two routes are offered side by side on
            // the same architecture rather than on two different files.
            final ProteinDomain card = new ProteinDomain( "CARD", 6, 90, "PF00619", 7.0E-26 );
            final ProteinDomain nbarc = new ProteinDomain( "NB-ARC", 109, 414, "PF00931", 7.2E-117 );
            final ProteinDomain wd40 = new ProteinDomain( "WD40", 1168, 1204, 0.3 ); // no id
            if ( !TreePanel.domainUrl( card ).contains( "/entry/pfam/PF00619/" )
                    || !TreePanel.domainUrl( nbarc ).contains( "/entry/pfam/PF00931/" ) ) {
                System.out.println( "an accessioned domain of a mixed architecture lost its entry link" );
                return false;
            }
            if ( !TreePanel.domainUrl( wd40 ).contains( "/search/text/WD40" ) ) {
                System.out.println( "the un-accessioned domain beside them should still be looked up: "
                        + TreePanel.domainUrl( wd40 ) );
                return false;
            }
            final String wd_rows = rowsOf( NodeHoverText.domainRows( tip, wd40, 1248 ) );
            if ( wd_rows.contains( "Accession" ) || wd_rows.contains( "the Pfam entry" ) ) {
                System.out.println( "the un-accessioned domain must not borrow its neighbours' promise: "
                        + wd_rows );
                return false;
            }
            if ( !wd_rows.contains( "Protein length: 1248 aa" ) ) {
                System.out.println( "the protein length is a property of the protein, not the domain: "
                        + wd_rows );
                return false;
            }
            // (4) the link is the InterPro location; pfam.xfam.org is retired
            if ( !"https://www.ebi.ac.uk/interpro/entry/pfam/PF00018/".equals( TreePanel.pfamUrl( "PF00018" ) ) ) {
                System.out.println( "wrong Pfam URL: " + TreePanel.pfamUrl( "PF00018" ) );
                return false;
            }
            // (5) an e-value reads as a power of ten, a plain number does not
            if ( !"1E-6".equals( NodeHoverText.formatEValue( 1e-6 ) ) ) {
                System.out.println( "e-value formatting: " + NodeHoverText.formatEValue( 1e-6 ) );
                return false;
            }
            if ( !"0.5".equals( NodeHoverText.formatEValue( 0.5 ) ) ) {
                System.out.println( "a plain value should not be exponential: "
                        + NodeHoverText.formatEValue( 0.5 ) );
                return false;
            }
            // (5b) the click gate. Opening a browser is the one irreversible thing this feature does, so
            // what does NOT open one is worth pinning: a right-click (which belongs to the context menu), any
            // modified click (ctrl/shift/alt/meta all mean something else here), and the SECOND event of a
            // double click -- two MOUSE_CLICKED events would otherwise open two tabs.
            final javax.swing.JPanel src = new javax.swing.JPanel();
            final int CLICKED = java.awt.event.MouseEvent.MOUSE_CLICKED;
            final java.awt.event.MouseEvent plain = new java.awt.event.MouseEvent( src, CLICKED, 0L, 0, 5, 5, 1,
                    false, java.awt.event.MouseEvent.BUTTON1 );
            if ( !TreePanel.isPlainSingleLeftClickForTest( plain ) ) {
                System.out.println( "a plain left click must open the domain's page" );
                return false;
            }
            final Object[][] refused = {
                    { "a right click", Integer.valueOf( 0 ), Integer.valueOf( 1 ),
                            Integer.valueOf( java.awt.event.MouseEvent.BUTTON3 ) },
                    { "a middle click", Integer.valueOf( 0 ), Integer.valueOf( 1 ),
                            Integer.valueOf( java.awt.event.MouseEvent.BUTTON2 ) },
                    { "the second event of a double click", Integer.valueOf( 0 ), Integer.valueOf( 2 ),
                            Integer.valueOf( java.awt.event.MouseEvent.BUTTON1 ) },
                    { "a ctrl-click", Integer.valueOf( java.awt.event.InputEvent.CTRL_DOWN_MASK ),
                            Integer.valueOf( 1 ), Integer.valueOf( java.awt.event.MouseEvent.BUTTON1 ) },
                    { "a shift-click", Integer.valueOf( java.awt.event.InputEvent.SHIFT_DOWN_MASK ),
                            Integer.valueOf( 1 ), Integer.valueOf( java.awt.event.MouseEvent.BUTTON1 ) },
                    { "an alt-click", Integer.valueOf( java.awt.event.InputEvent.ALT_DOWN_MASK ),
                            Integer.valueOf( 1 ), Integer.valueOf( java.awt.event.MouseEvent.BUTTON1 ) },
                    { "a meta-click", Integer.valueOf( java.awt.event.InputEvent.META_DOWN_MASK ),
                            Integer.valueOf( 1 ), Integer.valueOf( java.awt.event.MouseEvent.BUTTON1 ) }, };
            for( final Object[] c : refused ) {
                final java.awt.event.MouseEvent e = new java.awt.event.MouseEvent( src, CLICKED, 0L,
                        ( (Integer) c[ 1 ] ).intValue(), 5, 5, ( (Integer) c[ 2 ] ).intValue(), false,
                        ( (Integer) c[ 3 ] ).intValue() );
                if ( TreePanel.isPlainSingleLeftClickForTest( e ) ) {
                    System.out.println( c[ 0 ] + " must not open a browser tab" );
                    return false;
                }
            }
            // (6) the hit-test agrees with the drawing, including which domains are drawn at all
            if ( GraphicsEnvironment.isHeadless() ) {
                return true;
            }
            final Phylogeny phy = fixture();
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, conf, "domhover" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    mp.getControlPanel().setCheckbox( DisplayOption.SHOW_DOMAIN_ARCHITECTURES, true );
                    tp.setSize( 1400, 800 );
                    tp.calcParametersForPainting( 1400, 800 );
                    final BufferedImage img = new BufferedImage( 1400, 800, BufferedImage.TYPE_INT_ARGB );
                    final Graphics2D g = img.createGraphics();
                    tp.printAll( g );
                    g.dispose();
                    int hits = 0;
                    boolean saw_sh3 = false;
                    boolean saw_weak = false;
                    for( int y = 0; y < 800; y += 2 ) {
                        for( int x = 0; x < 1400; x += 2 ) {
                            final TreePanel.DomainHit h = tp.domainAt( x, y );
                            if ( h != null ) {
                                ++hits;
                                if ( "SH3".equals( h.domain().getName() ) ) {
                                    saw_sh3 = true;
                                }
                                if ( "Weak".equals( h.domain().getName() ) ) {
                                    saw_weak = true;
                                }
                            }
                        }
                    }
                    if ( hits < 1 ) {
                        System.out.println( "the domain hit-test never fires anywhere on a drawn architecture" );
                        ok[ 0 ] = false;
                        return;
                    }
                    if ( !saw_sh3 ) {
                        System.out.println( "a drawn domain was never hit" );
                        ok[ 0 ] = false;
                        return;
                    }
                    // A domain above the e-value threshold is NOT drawn, so the rollover must not claim it:
                    // that is the whole reason the hit-test asks the renderable rather than the architecture.
                    if ( saw_weak ) {
                        System.out.println( "the hit-test found a domain that the threshold keeps off screen" );
                        ok[ 0 ] = false;
                        return;
                    }
                    // (7) the hit-test must agree with the PIXELS, not merely fire somewhere. Compared
                    // against what was painted rather than against a second copy of the geometry: an offset
                    // applied to both would cancel out, and this is the one failure a shared formula is
                    // supposed to make impossible -- so it is worth proving rather than assuming.
                    // Run for EVERY origin the painter can choose. The three are computed in different
                    // places -- alignedPhylogramDomainColumnX, the tips-plus-longest-label sum, and
                    // verticalDomainColumnStart -- so one of them being wrong is invisible to a test that
                    // exercises only another. That is exactly how the vertical origin shipped wrong.
                    if ( !alignedAtOrigin( tp, "cladogram", 1400, 800 ) ) {
                        ok[ 0 ] = false;
                        return;
                    }
                    mp.getControlPanel()
                            .setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM );
                    if ( !alignedAtOrigin( tp, "phylogram", 1400, 800 ) ) {
                        ok[ 0 ] = false;
                        return;
                    }
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    if ( !alignedAtOrigin( tp, "root-top", 1400, 800 ) ) {
                        ok[ 0 ] = false;
                        return;
                    }
                    // CIRCULAR draws a concentric domain ring, and CLAUDE.md calls it the most important
                    // layout, so it is not an approved exception -- it gets the same rollover. Domain boxes
                    // are only drawn there with RADIAL labels (domainBoxesDrawnInCurrentLayout), which is
                    // the condition the hit-test has to share rather than guess at.
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    // TRIANGULAR is drawn by the SAME rectangular paint loop (paintPhylogeny takes that
                    // branch for everything that is not CIRCULAR or UNROOTED), so the boxes are there and
                    // the rollover has to answer on them. It did not: the hit-test listed the three types it
                    // expected instead of asking what the painter does, and the feature was silently dead in
                    // one of the five display types.
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR );
                    if ( !alignedAtOrigin( tp, "triangular", 1400, 800 ) ) {
                        ok[ 0 ] = false;
                        return;
                    }
                    mp.getOptions().setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.RADIAL );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    if ( !alignedAtOrigin( tp, "circular", 2000, 1600 ) ) {
                        ok[ 0 ] = false;
                        return;
                    }
                    // UNROOTED draws the same radial track, pivoting on each tip rather than on a common
                    // centre, so it gets the rollover too -- all five display types at parity, no exception.
                    mp.getOptions().setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.RADIAL );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                    if ( !alignedAtOrigin( tp, "unrooted", 2000, 1600 ) ) {
                        ok[ 0 ] = false;
                        return;
                    }
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    // and the gate itself: with labels not radial, nothing is drawn, so nothing may be hit
                    mp.getOptions().setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.HORIZONTAL );
                    tp.setSize( 2000, 1600 );
                    tp.calcParametersForPainting( 2000, 1600 );
                    int claimed_without_boxes = 0;
                    for( int yy = 0; yy < 1600; yy += 3 ) {
                        for( int xx = 0; xx < 2000; xx += 3 ) {
                            if ( tp.domainAt( xx, yy ) != null ) {
                                ++claimed_without_boxes;
                            }
                        }
                    }
                    if ( claimed_without_boxes > 0 ) {
                        System.out.println( "circular with non-radial labels draws no domain boxes, yet the "
                                + "hit-test claims " + claimed_without_boxes + " points" );
                        ok[ 0 ] = false;
                    }
                }
                catch ( final Throwable t ) {
                    t.printStackTrace( System.out );
                    ok[ 0 ] = false;
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
            if ( !ok[ 0 ] ) {
                return false;
            }
            // (8) "Auto-hide Labels": on a dense tree the painter keeps only every n-th tip and draws nothing for
            // the rest, so the rollover must not report a domain on a tip whose track was dropped. Reported from
            // use -- the geometry alone still finds every track, drawn or not. Same instrument as (7): what the
            // hit-test claims must be pixels something was drawn on. Under hiding the drawn boxes are one row
            // tall but n rows apart, so a hidden tip's band is empty canvas and a wrong claim has nowhere to hide.
            final Phylogeny dense = denseFixture( 60 );
            final MainFrame[] mf2 = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf2[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { dense }, conf, "domhover-dense" ) );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final MainPanel mp = mf2[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    mp.getControlPanel().setCheckbox( DisplayOption.SHOW_DOMAIN_ARCHITECTURES, true );
                    mp.getControlPanel().setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
                    // Find a canvas where hiding actually bites AND the rows are still far enough apart that a
                    // hidden tip's band is not covered by its drawn neighbours' boxes (the box height follows the
                    // row spacing, clamped at DOMAIN_STRUCTURE_HEIGHT_MIN). The font decides where that is, so it
                    // is measured rather than assumed: a fixture that quietly fails its own precondition pins
                    // nothing, and this one would pass for the wrong reason if nothing were hidden at all.
                    // Hiding starts when the labels no longer fit the row spacing, so a BIG label font puts
                    // that threshold at a comfortable row spacing rather than a cramped one. It has to be
                    // comfortable: a domain box is as tall as one row but no shorter than
                    // DOMAIN_STRUCTURE_HEIGHT_MIN (6 px), so at a row spacing under 6 the DRAWN boxes spill over
                    // their neighbours' rows and cover exactly the empty canvas this case has to see. That is
                    // what the first version of this fixture did, and every mutation survived it.
                    mp.getTreeFontSet().setBaseFont( new java.awt.Font( "SansSerif", java.awt.Font.PLAIN, 28 ) );
                    final int w = 1200;
                    int h = 0;
                    final StringBuilder tried = new StringBuilder();
                    for( final int candidate : new int[] { 900, 1100, 1300, 1500, 1700, 1900, 2100 } ) {
                        tp.setSize( w, candidate );
                        tp.calcParametersForPainting( w, candidate );
                        tried.append( " " ).append( candidate ).append( "px:y-dist=" ).append( tp.getYdistance() )
                                .append( ",font=" ).append( tp.getLargeFontHeight() ).append( ",hiding=" )
                                .append( tp.labelsDynamicallyHidden() );
                        if ( tp.labelsDynamicallyHidden() && ( tp.getYdistance() >= 7.0f ) ) {
                            h = candidate;
                            break;
                        }
                    }
                    if ( h == 0 ) {
                        System.out.println( "auto-hide: no canvas both hid labels and kept the rows more than a "
                                + "domain box apart, so the fixture pins nothing. Tried:" + tried );
                        ok[ 0 ] = false;
                        return;
                    }
                    // Paint TALL first, where nothing is hidden and every tip's track is drawn, and only then
                    // dense. What the rollover may report has to be what the LAST paint drew, so a record that
                    // is added to but never cleared -- the shape of staleness this codebase has shipped before --
                    // shows up as the zoomed-out tips still answering after the zoom back in.
                    paintOnce( tp, w, 4000 );
                    if ( !alignedAtOrigin( tp, "auto-hide", w, h, 2 ) ) {
                        ok[ 0 ] = false;
                    }
                }
                catch ( final Throwable t ) {
                    t.printStackTrace( System.out );
                    ok[ 0 ] = false;
                }
                finally {
                    ( (JFrame) mf2[ 0 ] ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
    }

    public static void main( final String[] args ) {
        if ( test() ) {
            System.out.println( "DomainHoverTest: OK." );
        }
        else {
            System.out.println( "DomainHoverTest: FAILED." );
        }
    }
}
