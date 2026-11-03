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

import java.awt.GraphicsEnvironment;
import java.io.File;
import java.util.EnumMap;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * A saved figure wins over the load-time guesses ({@link AptxUtil#lookAtSomeTreePropertiesForAptxControlSettings}):
 * every setting a tree's figure names opens as it was saved, whatever the data would switch on, and a setting the
 * figure does not name is guessed as for any tree. The one exception: a tree without branch lengths opens as a
 * cladogram whatever its figure says.
 * <p>
 * Every tree here is opened through File &gt; Open's own path ({@link AptxUtil#addPhylogeniesToTabs}), the step the
 * older figure tests stepped around -- they added the tab directly, or asked only about settings no guess touches,
 * which is how a tree saved with its confidence values off went on reopening with them on. Each check reads the tab
 * the moment that path returns, before any later fit could carry a widget onto it.
 * <p>
 * Fixture: {@code forester/demo/saved-figure-none.xml}, a tree that reaches every guess (support values, events,
 * branch colour and width on its first tip, domains, an alignment, lengths on every branch), with figures written
 * here; and the shipped {@code saved-figure.xml} as it is.
 */
public final class SavedFigureWinsTest {

    /** The on/off settings a load-time guess decides (the phylogram choice is the seventh, a three-way one). */
    static final DisplayOption[] GUESSED_SWITCHES = { DisplayOption.WRITE_CONFIDENCE_VALUES, DisplayOption.WRITE_EVENTS,
            DisplayOption.USE_STYLE, DisplayOption.WIDTH_BRANCHES, DisplayOption.SHOW_DOMAIN_ARCHITECTURES,
            DisplayOption.SHOW_MSA };

    /** Every on/off Display option: what the manual promises a saved figure restores ("every Display-panel
     *  checkbox"). The phylogram choice is a three-way one, carried by {@code displaytype}. */
    static final DisplayOption[] ALL_SWITCHES = java.util.Arrays.stream( DisplayOption.values() )
            .filter( o -> o != DisplayOption.DISPLAY_AS_PHYLOGRAM ).toArray( DisplayOption[]::new );

    /** A scrollbar's width, at most: the fit reads the viewport before fitWidth, which may add or drop one. */
    private static final int SCROLLBAR_ALLOWANCE = 24;

    private static final String FIXTURE = "forester/demo/saved-figure-none.xml";
    private static final String DEMO    = "forester/demo/saved-figure.xml";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SavedFigureWins: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        boolean ok = names();
        if ( GraphicsEnvironment.isHeadless() ) {
            return ok;
        }
        try {
            final Options.PHYLOGENY_DISPLAY_TYPE guess = guessedType();
            final Options.PHYLOGENY_DISPLAY_TYPE other = ( guess == Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM )
                    ? Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM : Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM;
            ok &= everyNamedSettingWins( other );
            ok &= eachUnnamedSettingIsGuessed( guess, other );
            ok &= theOtherDirection();
            ok &= noLengthsOpensACladogram( other );
            ok &= domainsAndLabelDirection();
            ok &= saveAndReopen( guess, other );
            ok &= everyCheckboxComesBack();
            ok &= domainTreeReopens();
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected " + e );
        }
        return ok;
    }

    // ---- what a figure names, by the test applyTo uses --------------------------------------------------------------

    private static boolean names() {
        boolean ok = true;
        final FigureSpec empty_value = FigureSpec.parse( "v1;show.WRITE_EVENTS=" );
        ok &= ck( empty_value.names( DisplayOption.WRITE_EVENTS ), "an empty show. value is named (applyTo sets it off)" );
        ok &= ck( Boolean.FALSE.equals( empty_value.shows( DisplayOption.WRITE_EVENTS ) ), "...and reads as off" );
        ok &= ck( !empty_value.names( DisplayOption.WRITE_CONFIDENCE_VALUES )
                && ( empty_value.shows( DisplayOption.WRITE_CONFIDENCE_VALUES ) == null ), "an absent key is not named" );
        ok &= ck( Boolean.TRUE.equals( FigureSpec.parse( "v1;show.SHOW_MSA=true" ).shows( DisplayOption.SHOW_MSA ) ),
                  "show.SHOW_MSA=true reads as on" );
        final FigureSpec aligned = FigureSpec.parse( "v1;displaytype=ALIGNED_PHYLOGRAM" );
        ok &= ck( aligned.names( DisplayOption.DISPLAY_AS_PHYLOGRAM ), "a known displaytype names the phylogram choice" );
        ok &= ck( aligned.shows( DisplayOption.DISPLAY_AS_PHYLOGRAM ) == null,
                  "the phylogram choice is three-way: shows() has no on/off for it" );
        ok &= ck( !FigureSpec.parse( "v1;displaytype=SIDEWAYS" ).names( DisplayOption.DISPLAY_AS_PHYLOGRAM ),
                  "an unknown displaytype names nothing (applyTo skips it)" );
        ok &= ck( !FigureSpec.parse( "v1;displaytype=" ).names( DisplayOption.DISPLAY_AS_PHYLOGRAM ),
                  "an empty displaytype names nothing" );
        ok &= ck( !FigureSpec.parse( "v1;show.DISPLAY_AS_PHYLOGRAM=true" ).names( DisplayOption.DISPLAY_AS_PHYLOGRAM ),
                  "the phylogram choice is named by displaytype only" );
        ok &= ck( FigureSpec.parse( "v1;labeldirection=HORIZONTAL" ).namesLabelDirection(), "labeldirection=HORIZONTAL" );
        ok &= ck( !FigureSpec.parse( "v1;labeldirection=DIAGONAL" ).namesLabelDirection(),
                  "an unknown label direction names nothing" );
        ok &= ck( !FigureSpec.parse( "v1;layout=CIRCULAR" ).namesLabelDirection(), "no labeldirection key" );
        return ok;
    }

    // ---- the figure wins -------------------------------------------------------------------------------------------

    /** Every guessed setting named the other way opens as named -- the fixture, and the shipped demo as it is. */
    private static boolean everyNamedSettingWins( final Options.PHYLOGENY_DISPLAY_TYPE other ) throws Exception {
        boolean ok = true;
        // + a show. key for the phylogram choice, which no writer writes: displaytype carries it, so it is not applied
        final Opened o = open( fixture( otherWay( other, null ) + ";show.DISPLAY_AS_PHYLOGRAM=false" ), null );
        ok &= ck( o._type == other, "a figure's phylogram choice must open as saved: " + other + ", got " + o._type );
        ok &= ck( !o._phylogram_as_boolean, "show.DISPLAY_AS_PHYLOGRAM must not become a per-tab boolean (what "
                + "applyTo sets and what names() reports are one test)" );
        for( final DisplayOption opt : GUESSED_SWITCHES ) {
            ok &= ck( !o._shows.get( opt ), "a figure saved with " + opt + " off must open with it off" );
        }
        ok &= ck( o.widgetsAgree(), "the checkboxes must show what the tab draws: " + o._shows + " vs " + o._widgets );
        // ticking the box only re-fits the tree, so the track's width is set as the tree opens, on or off
        ok &= ck( o.domainTrackFitted(), "a domain track the figure keeps off still gets the load's width, got "
                + o._domain_width + " for " + o._viewport_width );
        final Opened demo = open( FigureRenderer.readTrees( new File( System.getProperty( "user.dir" ), DEMO ) )[ 0 ],
                                  null );
        ok &= ck( demo._type == Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM,
                  "saved-figure.xml opens as an aligned phylogram, got " + demo._type );
        for( final DisplayOption opt : GUESSED_SWITCHES ) {
            ok &= ck( !demo._shows.get( opt ), "saved-figure.xml opens with " + opt + " off" );
        }
        return ok;
    }

    /**
     * Each setting the figure does NOT name is guessed, the rest stay as named. Also the precondition of every check
     * above: the fixture reaches each guess, and the guess differs from what the figure says. And each guess reaches
     * the TAB by itself -- with the confidence guess skipped, nothing after the colour or width guess pushed the
     * checkboxes onto the tab.
     */
    private static boolean eachUnnamedSettingIsGuessed( final Options.PHYLOGENY_DISPLAY_TYPE guess,
                                                        final Options.PHYLOGENY_DISPLAY_TYPE other ) throws Exception {
        boolean ok = true;
        final Opened phylogram = open( fixture( otherWay( other, DisplayOption.DISPLAY_AS_PHYLOGRAM ) ), null );
        ok &= ck( phylogram._type == guess, "unnamed, the phylogram choice is guessed: " + guess + ", got "
                + phylogram._type );
        for( final DisplayOption opt : GUESSED_SWITCHES ) {
            ok &= ck( !phylogram._shows.get( opt ), "...while " + opt + " stays off as named" );
        }
        for( final DisplayOption unnamed : GUESSED_SWITCHES ) {
            // every switch OFF in the window first: a new tab adopts the checkboxes for what its figure does not name,
            // and Visual Styles is ON by default -- "guessed on" must come from the guess, not from the default
            final Opened o = open( fixture( otherWay( other, unnamed ) ), SavedFigureWinsTest::allSwitchesOff );
            ok &= ck( o._shows.get( unnamed ), "unnamed, " + unnamed + " is guessed on for this tree" );
            ok &= ck( o._type == other, "...while the phylogram choice stays " + other + " (" + unnamed + "), got "
                    + o._type );
            for( final DisplayOption opt : GUESSED_SWITCHES ) {
                if ( opt != unnamed ) {
                    ok &= ck( !o._shows.get( opt ), "...while " + opt + " stays off as named (" + unnamed + ")" );
                }
            }
            ok &= ck( o.widgetsAgree(), "the checkboxes must show what the tab draws (" + unnamed + "): " + o._shows
                    + " vs " + o._widgets );
            if ( unnamed == DisplayOption.SHOW_DOMAIN_ARCHITECTURES ) {
                ok &= ck( o.domainTrackFitted(), "a guessed domain track is fitted to a quarter of the window, got "
                        + o._domain_width + " for " + o._viewport_width );
            }
        }
        return ok;
    }

    /** The two-way guesses, the other way: a figure naming confidence values and events ON keeps them on a tree
     *  that has none; a figure's phylogram on a tree measured too thinly for the guess to draw one is kept. */
    private static boolean theOtherDirection() throws Exception {
        boolean ok = true;
        final String on = "v1;layout=RECTANGULAR;show.WRITE_CONFIDENCE_VALUES=true;show.WRITE_EVENTS=true";
        final Opened named = open( withoutSupportAndEvents( fixture( on ) ), null );
        final Opened guessed = open( withoutSupportAndEvents( fixture( "v1;layout=RECTANGULAR" ) ), null );
        ok &= ck( !guessed._shows.get( DisplayOption.WRITE_CONFIDENCE_VALUES )
                && !guessed._shows.get( DisplayOption.WRITE_EVENTS ),
                  "precondition: a tree without support values or events has both guessed off" );
        ok &= ck( named._shows.get( DisplayOption.WRITE_CONFIDENCE_VALUES ),
                  "a figure saved with confidence values on keeps them on" );
        ok &= ck( named._shows.get( DisplayOption.WRITE_EVENTS ), "a figure saved with events on keeps them on" );
        final Opened thin_named = open( thinlyMeasured( fixture( "v1;layout=RECTANGULAR;displaytype=UNALIGNED_PHYLOGRAM" ) ),
                                        null );
        final Opened thin_guessed = open( thinlyMeasured( fixture( "v1;layout=RECTANGULAR" ) ), null );
        ok &= ck( thin_guessed._type == Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM,
                  "precondition: a tree with lengths on its tips only is guessed a cladogram, got " + thin_guessed._type );
        ok &= ck( thin_named._type == Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM,
                  "...and a figure saved as a phylogram keeps it, got " + thin_named._type );
        ok &= ck( thin_named._phylogram_enabled, "...with the phylogram buttons live (the tree has lengths)" );
        return ok;
    }

    /** The exception: with no branch lengths there is no phylogram to draw, whatever the figure says. */
    private static boolean noLengthsOpensACladogram( final Options.PHYLOGENY_DISPLAY_TYPE other ) throws Exception {
        final Phylogeny phy = fixture( "v1;layout=RECTANGULAR;displaytype=" + other.name() );
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            it.next().setDistanceToParent( PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT );
        }
        final Opened o = open( phy, null );
        boolean ok = ck( o._type == Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM,
                         "a tree without branch lengths opens as a cladogram whatever its figure says, got " + o._type );
        ok &= ck( !o._phylogram_enabled, "...with the phylogram buttons greyed" );
        return ok;
    }

    /**
     * Domains in a radial figure: a figure that names a label direction keeps it (the domain guess switches to
     * radial labels only when no direction is named), and a figure that names the track ON still gets its fit.
     */
    private static boolean domainsAndLabelDirection() throws Exception {
        boolean ok = true;
        final String circular_domains = "v1;layout=CIRCULAR;show.SHOW_DOMAIN_ARCHITECTURES=true";
        final Opened guessed = open( fixture( circular_domains ),
                                     mf -> mf.setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.HORIZONTAL ) );
        ok &= ck( guessed._direction == Options.NODE_LABEL_DIRECTION.RADIAL,
                  "precondition: with no direction named, a radial tree with domains switches to radial labels, got "
                          + guessed._direction );
        final Opened named = open( fixture( circular_domains + ";labeldirection=HORIZONTAL" ),
                                   mf -> mf.setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.RADIAL ) );
        ok &= ck( named._direction == Options.NODE_LABEL_DIRECTION.HORIZONTAL,
                  "a figure's label direction must stand, got " + named._direction );
        ok &= ck( named._direction_menu_radial == ( named._direction == Options.NODE_LABEL_DIRECTION.RADIAL ),
                  "the Radial Labels menu item must show the direction drawn" );
        ok &= ck( named._shows.get( DisplayOption.SHOW_DOMAIN_ARCHITECTURES ), "the figure's domain track is on" );
        ok &= ck( named.domainTrackFitted(), "a figure that names the track on still gets its initial fit, got "
                + named._domain_width + " for " + named._viewport_width );
        return ok;
    }

    /** The real thing: what File > Save writes (MainFrame: syncFigureToTree, then the phyloXML writer), read back
     *  and opened through File > Open's path. */
    private static boolean saveAndReopen( final Options.PHYLOGENY_DISPLAY_TYPE guess,
                                          final Options.PHYLOGENY_DISPLAY_TYPE other ) throws Exception {
        boolean ok = true;
        final boolean[] drive = { false, false };
        final Phylogeny saved = saved( fixture( null ), mf -> {
            // opened without a figure: every guess fires
            final Opened before = new Opened( mf );
            drive[ 0 ] = ( before._type == guess );
            for( final DisplayOption opt : GUESSED_SWITCHES ) {
                drive[ 0 ] &= before._shows.get( opt );
            }
            // what the user does: the other phylogram, every guessed setting off
            final ControlPanel cp = mf.getMainPanel().getControlPanel();
            cp.setTreeDisplayType( other );
            for( final DisplayOption opt : GUESSED_SWITCHES ) {
                cp.setCheckbox( opt, false );
            }
            final Opened changed = new Opened( mf );
            drive[ 1 ] = ( changed._type == other );
            for( final DisplayOption opt : GUESSED_SWITCHES ) {
                drive[ 1 ] &= !changed._shows.get( opt );
            }
        } );
        ok &= ck( drive[ 0 ], "precondition: opened without a figure, every guess fires" );
        ok &= ck( drive[ 1 ], "precondition: the user's changes took" );
        final Opened back = open( saved, null );
        ok &= ck( back._type == other, "saved as " + other + ", reopened as " + back._type );
        for( final DisplayOption opt : GUESSED_SWITCHES ) {
            ok &= ck( !back._shows.get( opt ), "saved with " + opt + " off, reopened with it on" );
        }
        return ok;
    }

    /**
     * The manual's promise, "every Display-panel checkbox" comes back: every on/off option set one way, saved as
     * File &gt; Save does and reopened in a window whose checkboxes all say the other way -- all off, then all on. A
     * figure key that never reached the tab fails here (Show Internal / External Data went into a store the paint
     * does not read, and came back on), and so does a later load-time guess that forgets to ask the figure.
     */
    private static boolean everyCheckboxComesBack() throws Exception {
        boolean ok = true;
        for( final boolean v : new boolean[] { false, true } ) {
            final boolean[] took = { false };
            final Phylogeny saved = saved( fixture( null ), mf -> {
                setEverySwitch( mf, v );
                took[ 0 ] = true;
                for( final DisplayOption opt : ALL_SWITCHES ) {
                    took[ 0 ] &= ( mf.getMainPanel().getCurrentTreePanel().shows( opt ) == v );
                }
            } );
            ok &= ck( took[ 0 ], "precondition: every switch " + ( v ? "on" : "off" ) + " in the tab that is saved" );
            final Opened back = open( saved, mf -> setEverySwitch( mf, !v ) );
            for( final DisplayOption opt : ALL_SWITCHES ) {
                ok &= ck( back._shows.get( opt ) == v, "saved with " + opt + ( v ? " on" : " off" )
                        + ", reopened with it " + ( back._shows.get( opt ) ? "on" : "off" ) );
            }
            ok &= ck( back.widgetsAgree(), "the checkboxes must show what the tab draws (all " + v + "): "
                    + back._shows + " vs " + back._widgets );
        }
        return ok;
    }

    /**
     * A tree saved with its domain track showing -- as every domain tree opens -- reopens. Its figure turns the
     * track on while it is applied, before the tab's first paint, and a label measure that assumed the paint had
     * come first threw a ClassCastException (0.11.198 and earlier): the file could not be opened again.
     */
    private static boolean domainTreeReopens() throws Exception {
        final boolean[] on = { false };
        final Phylogeny saved = saved( fixture( null ),
                                       mf -> on[ 0 ] = new Opened( mf )._shows.get( DisplayOption.SHOW_DOMAIN_ARCHITECTURES ) );
        boolean ok = ck( on[ 0 ] && Boolean.TRUE.equals( FigureSpec.readFrom( saved ).shows(
                DisplayOption.SHOW_DOMAIN_ARCHITECTURES ) ), "precondition: saved with the domain track on" );
        try {
            final Opened back = open( saved, null );
            ok &= ck( back._shows.get( DisplayOption.SHOW_DOMAIN_ARCHITECTURES ), "reopened with the domain track on" );
        }
        catch ( final java.lang.reflect.InvocationTargetException e ) {
            return fail( "a tree saved with its domains showing must reopen, got " + e.getCause() );
        }
        return ok;
    }

    /** Opens {@code phy} as File > Open does, lets {@code change} act on the window, saves the tab as File > Save
     *  does (its figure stamped, the phyloXML writer) and reads the file back. */
    private static Phylogeny saved( final Phylogeny phy, final Before change ) throws Exception {
        final File out = File.createTempFile( "saved-figure-wins", ".xml" );
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { placeholder() }, new Configuration(), "saved figure wins" ) );
            try {
                SwingUtilities.invokeAndWait( () -> {
                    AptxUtil.addPhylogeniesToTabs( new Phylogeny[] { phy }, "", "", mf[ 0 ].getConfiguration(),
                                                   mf[ 0 ].getMainPanel() );
                    change.run( mf[ 0 ] );
                    final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                    tp.syncFigureToTree();
                    try {
                        PhylogenyWriter.createPhylogenyWriter().toPhyloXML( out, tp.getPhylogeny(), 0 );
                    }
                    catch ( final Exception e ) {
                        throw new RuntimeException( e );
                    }
                } );
            }
            finally {
                SwingUtilities.invokeAndWait( () -> ( ( JFrame ) mf[ 0 ] ).dispose() );
            }
            return FigureRenderer.readTrees( out )[ 0 ];
        }
        finally {
            out.delete();
        }
    }

    // ---- helpers ---------------------------------------------------------------------------------------------------

    /** What one tab shows the moment File > Open's path has added it. */
    private static final class Opened {

        final Options.PHYLOGENY_DISPLAY_TYPE _type;
        final Map<DisplayOption, Boolean>    _shows   = new EnumMap<>( DisplayOption.class );
        final Map<DisplayOption, Boolean>    _widgets = new EnumMap<>( DisplayOption.class );
        final boolean                        _phylogram_enabled;
        final Options.NODE_LABEL_DIRECTION   _direction;
        final boolean                        _direction_menu_radial;
        final double                         _domain_width;
        final int                            _viewport_width;
        final boolean                        _phylogram_as_boolean;

        Opened( final MainFrame mf ) {
            final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
            final ControlPanel cp = mf.getMainPanel().getControlPanel();
            _type = tp.getTreeDisplayTypeForThisTab();
            for( final DisplayOption opt : ALL_SWITCHES ) {
                _shows.put( opt, tp.shows( opt ) );
                _widgets.put( opt, cp.currentValueOf( opt ) );
            }
            _phylogram_enabled = cp.getDisplayAsCladogramRb().isEnabled();
            _direction = mf.getOptions().getNodeLabelDirection();
            _direction_menu_radial = mf._label_direction_cbmi.isSelected();
            _domain_width = tp.domainStructureWidthForTest();
            _viewport_width = mf.getMainPanel().getSizeOfViewport().width;
            _phylogram_as_boolean = tp.hasOwnDisplayState( DisplayOption.DISPLAY_AS_PHYLOGRAM );
        }

        boolean widgetsAgree() {
            return _shows.equals( _widgets );
        }

        /** The load's fit: a quarter of the window (ControlPanel.INITIAL_DOMAIN_WIDTH_SCREEN_FRACTION), never the
         *  default -- asserted on a window where the two differ by more than a scrollbar can move the fit. */
        boolean domainTrackFitted() {
            final double slack = 0.25 * SCROLLBAR_ALLOWANCE;
            return ( _viewport_width > 0 )
                    && ( Math.abs( ( 0.25 * _viewport_width ) - AptxConstants.DOMAIN_STRUCTURE_DEFAULT_WIDTH ) > slack )
                    && ( Math.abs( _domain_width - ( 0.25 * _viewport_width ) ) <= slack );
        }
    }

    private interface Before {

        void run( MainFrame mf );
    }

    /** Opens {@code phy} as File > Open does, in a fresh window (no earlier tab's checkboxes carry over); {@code
     *  before} acts on the window first. */
    private static Opened open( final Phylogeny phy, final Before before ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { placeholder() }, new Configuration(), "saved figure wins" ) );
        try {
            final Opened[] o = new Opened[ 1 ];
            SwingUtilities.invokeAndWait( () -> {
                if ( before != null ) {
                    before.run( mf[ 0 ] );
                }
                AptxUtil.addPhylogeniesToTabs( new Phylogeny[] { phy }, "", "", mf[ 0 ].getConfiguration(),
                                               mf[ 0 ].getMainPanel() );
                if ( mf[ 0 ].getMainPanel().getCurrentTreePanel().getPhylogeny() != phy ) {
                    throw new IllegalStateException( "the opened tree is not the current tab" );
                }
                o[ 0 ] = new Opened( mf[ 0 ] );
            } );
            return o[ 0 ];
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( ( JFrame ) mf[ 0 ] ).dispose() );
        }
    }

    /** Every guessed switch off in the window's current tab (and so in the checkboxes a new tab adopts). */
    private static void allSwitchesOff( final MainFrame mf ) {
        for( final DisplayOption opt : GUESSED_SWITCHES ) {
            mf.getMainPanel().getControlPanel().setCheckbox( opt, false );
        }
        for( final DisplayOption opt : GUESSED_SWITCHES ) {
            if ( mf.getMainPanel().getControlPanel().currentValueOf( opt ) ) {
                throw new IllegalStateException( "precondition: could not switch " + opt + " off" );
            }
        }
    }

    /** Every on/off Display option set to {@code on} in the window's current tab (and so in the checkboxes a new
     *  tab adopts for what its figure does not name). */
    private static void setEverySwitch( final MainFrame mf, final boolean on ) {
        final ControlPanel cp = mf.getMainPanel().getControlPanel();
        for( final DisplayOption opt : ALL_SWITCHES ) {
            cp.setCheckbox( opt, on );
        }
        for( final DisplayOption opt : ALL_SWITCHES ) {
            if ( cp.currentValueOf( opt ) != on ) {
                throw new IllegalStateException( "precondition: could not switch " + opt + ( on ? " on" : " off" ) );
            }
        }
    }

    /** The phylogram choice the guess makes for the fixture in a window opened with this machine's preferences. */
    private static Options.PHYLOGENY_DISPLAY_TYPE guessedType() throws Exception {
        final Opened o = open( fixture( "v1;layout=RECTANGULAR" ), null );
        return o._type;
    }

    /** A figure naming the phylogram choice as {@code type} and every guessed switch OFF -- each the other way from
     *  what the fixture's guesses give -- except {@code except} (null: none). */
    private static String otherWay( final Options.PHYLOGENY_DISPLAY_TYPE type, final DisplayOption except ) {
        final StringBuilder sb = new StringBuilder( "v1;layout=RECTANGULAR" );
        if ( except != DisplayOption.DISPLAY_AS_PHYLOGRAM ) {
            sb.append( ";displaytype=" ).append( type.name() );
        }
        for( final DisplayOption opt : GUESSED_SWITCHES ) {
            if ( opt != except ) {
                sb.append( ";show." ).append( opt.name() ).append( "=false" );
            }
        }
        return sb.toString();
    }

    /** A fresh copy of the fixture tree carrying {@code figure} (null: none). */
    private static Phylogeny fixture( final String figure ) throws Exception {
        final Phylogeny phy = FigureRenderer.readTrees( new File( System.getProperty( "user.dir" ), FIXTURE ) )[ 0 ];
        if ( FigureSpec.readFrom( phy ) != null ) {
            throw new IllegalStateException( FIXTURE + " must carry no figure of its own" );
        }
        if ( figure != null ) {
            FigureSpec.writeToTree( phy, FigureSpec.parse( figure ) );
        }
        return phy;
    }

    private static Phylogeny withoutSupportAndEvents( final Phylogeny phy ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            n.getBranchData().getConfidences().clear();
            n.getNodeData().setEvent( null );
        }
        if ( AptxUtil.isHasAtLeastOneBranchWithSupportValues( phy ) || AptxUtil.isHasAtLeastNodeWithEvent( phy ) ) {
            throw new IllegalStateException( "the stripped fixture still carries support values or events" );
        }
        return phy;
    }

    /** Lengths on the tips only: a tree that has lengths, but too few for the guess to draw a phylogram. */
    private static Phylogeny thinlyMeasured( final Phylogeny phy ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() ) {
                n.setDistanceToParent( PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT );
            }
        }
        if ( !AptxUtil.isHasAtLeastOneBranchLengthLargerThanZero( phy ) || AptxUtil.isMostlyMeasuredForPhylogram( phy ) ) {
            throw new IllegalStateException( "the thinly measured fixture is not" );
        }
        return phy;
    }

    /** The window's first tab: two tips with lengths and no data, so it reaches no guess of its own. */
    private static Phylogeny placeholder() {
        final PhylogenyNode root = new PhylogenyNode();
        for( final String name : new String[] { "a", "b" } ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( name );
            tip.setDistanceToParent( 0.1 );
            root.addAsChild( tip );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean ck( final boolean b, final String msg ) {
        return b || fail( msg );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [SavedFigureWinsTest] " + msg );
        return false;
    }
}
