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
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.LinkedHashMap;
import java.util.Collections;
import java.util.Map;
import java.util.Set;

import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;
import org.forester.ws.seqdb.TaxonomicLineageService;

/**
 * The shared "which rank, from which source" controls of "Annotate Clades by Rank" and "Colorize Subtrees via
 * Taxonomic Rank": a rank chooser (plus optional extra-level choosers) whose entries state, per rank, how many nodes
 * are annotated at it, what share of the tips the CHOSEN {@link TaxonSource} places right now, and -- for an online
 * source -- how many names it would still look up; and the three source radio buttons. Changing the source rebuilds
 * every chooser's entries, keeping the selected rank, so the number the user reads is always the number the chosen
 * source will deliver. No network: the figures come from the tree and the service's cache.
 */
final class RankSourcePanel {

    /** The "no rank" entry of an extra-level chooser. */
    static final String            NONE = "(none)";
    // The major Linnaean ranks always offered (a database can place a tip at any of them), even when no node in the
    // tree is annotated there. NCBI calls the domain level "superkingdom" in older and "domain" in newer dumps.
    private static final Set<String> OFFERABLE_RANKS = new HashSet<String>( Arrays.asList( "domain", "superkingdom",
            "kingdom", "phylum", "class", "order", "family", "genus", "species" ) );

    private final Phylogeny                _phy;
    private final TaxonomicLineageService  _svc;
    private final String                   _first_entry;
    private final JComboBox<String>        _rank_box;
    private final List<JComboBox<String>>  _extra_boxes = new ArrayList<JComboBox<String>>();
    private final JRadioButton[]           _radios      = new JRadioButton[ TaxonSource.values().length ];
    private final JPanel                   _source_panel;
    private final JCheckBox                _write_lineages;
    private final List<Runnable>           _listeners   = new ArrayList<Runnable>();

    /**
     * @param phy          the tree the figures are computed on
     * @param svc          the lineage cache consulted (never fetched); may be null
     * @param initial      the source the radios open on
     * @param first_entry  an entry placed before the ranks in the main chooser (e.g. "stop drawing"), or null
     * @param extra_levels how many extra-level choosers to build (each opens on {@link #NONE})
     */
    RankSourcePanel( final Phylogeny phy, final TaxonomicLineageService svc, final TaxonSource initial,
                     final String first_entry, final int extra_levels ) {
        _phy = phy;
        _svc = svc;
        _first_entry = first_entry;
        final ButtonGroup group = new ButtonGroup();
        _source_panel = new JPanel( new java.awt.GridLayout( 0, 1, 0, 0 ) );
        for( final TaxonSource s : TaxonSource.values() ) {
            final JRadioButton rb = new JRadioButton( s.getLabel(), s == initial );
            rb.setToolTipText( tooltipOf( s ) );
            rb.addActionListener( e -> rebuildChoices() );
            group.add( rb );
            _source_panel.add( rb );
            _radios[ s.ordinal() ] = rb;
        }
        _write_lineages = new JCheckBox( "Also write each looked-up lineage into the tree as metadata fields"
                + " (ncbi:phylum \u2026 ncbi:species; undoable)", true );
        _write_lineages.setToolTipText( "After an online source has looked tips up, keep what it found IN the tree: one"
                + " field per rank, so the next annotation, Color by, search and table work offline and survive a save."
                + " Replaces a tip's earlier ncbi: fields; tips the database could not place are left untouched." );
        _write_lineages.setEnabled( initial.isOnline() );
        final List<Coverage> cov = coverage( phy, svc, initial );
        final List<String> entries = new ArrayList<String>();
        for( final Coverage c : cov ) {
            entries.add( entry( c.label, c.annotated, c.placed, c.tips, c.to_look_up, initial.isOnline() ) );
        }
        final String[] choices = entries.toArray( new String[ 0 ] );
        _rank_box = new JComboBox<String>( withFirst( choices, first_entry ) );
        for( int i = 0; i < extra_levels; ++i ) {
            final JComboBox<String> box = new JComboBox<String>( withFirst( choices, NONE ) );
            box.setSelectedIndex( 0 );
            _extra_boxes.add( box );
        }
        // open on the finest rank naming 2..20 taxa (the joint default); a caller with a remembered rank selects
        // over it afterwards
        final String open_on = defaultRank( cov, 20 );
        if ( open_on != null ) {
            selectRank( open_on );
        }
    }

    private static String tooltipOf( final TaxonSource s ) {
        switch ( s ) {
            case TREE:
                return "Only what the tree already carries: a tip's own taxonomy at the rank, a field named for the rank in any"
                        + " namespace (vipr:Genus, ncbi:family, gtdb:family), else the nearest ancestor naming one. A tip with"
                        + " none does not vote: it is swept into the clade around it. Nothing is looked up.";
            case DATABASE:
                return "Every tip is looked up (its NCBI id, else its name, else its label) and its own identity wins;"
                        + " an ancestor's annotation only fills a gap the database cannot. Needs an internet connection.";
            default:
                return "The tree places what it can; only the tips neither their own nor an ancestor's annotation"
                        + " places are looked up. Needs an internet connection for those.";
        }
    }

    /** {@code choices} preceded by {@code first} when it is not null; {@code choices} itself otherwise, so a chooser
     *  with nothing to remove shows only ranks. */
    static String[] withFirst( final String[] choices, final String first ) {
        if ( first == null ) {
            return choices;
        }
        final String[] out = new String[ choices.length + 1 ];
        out[ 0 ] = first;
        System.arraycopy( choices, 0, out, 1, choices.length );
        return out;
    }

    /** A rank's coverage on the tree under a source: see {@link #coverage}. */
    static final class Coverage {

        final String key;        // the rank's key (RankVocabulary)
        final String label;      // what to show: the canonical spelling, else the file's own
        final int    annotated;  // nodes whose <taxonomy> is at this rank
        final int    placed;     // tips the source places now
        final int    tips;
        final int    taxa;       // distinct taxa among the placed tips
        final int    to_look_up; // names an online source would still look up

        Coverage( final String key, final String label, final int annotated, final int placed, final int tips,
                  final int taxa, final int to_look_up ) {
            this.key = key;
            this.label = label;
            this.annotated = annotated;
            this.placed = placed;
            this.tips = tips;
            this.taxa = taxa;
            this.to_look_up = to_look_up;
        }

        /** The rank as the placement code takes it: the canonical spelling of a listed rank, else the key. */
        String rank() {
            return ( RankVocabulary.label( key ) != null ) ? RankVocabulary.label( key ) : key;
        }
    }

    /**
     * The ranks the tree names anywhere -- in a {@code <taxonomy>}'s rank, or as a property's local name in any
     * namespace ({@link RankVocabulary#propertyRank}) -- plus, for an online source, the major Linnaean ranks
     * ({@link #OFFERABLE_RANKS}: a database can place a tip at a rank no node literally carries); in canonical
     * order, each with what {@code source} places now. The joint offer rule (Archaeopteryx.js {@code rankCoverage}):
     * a rank is worth offering when it places a tip and its taxa repeat -- at most three taxa per five placed tips,
     * the Color-by's own near-unique line -- so a strain-like field drops out by itself (its clades would be single
     * tips, nothing to colour or to mark); an online source's rank with names still to look up is offered too.
     */
    static List<Coverage> coverage( final Phylogeny phy, final TaxonomicLineageService svc, final TaxonSource source ) {
        final List<Coverage> out = new ArrayList<Coverage>();
        if ( ( phy == null ) || phy.isEmpty() ) {
            return out;
        }
        final Map<String, String> named = new LinkedHashMap<String, String>(); // key -> label, the file's first spelling
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy() ) {
                final String raw = n.getNodeData().getTaxonomy().getRank();
                final String key = RankVocabulary.key( raw );
                if ( key != null ) {
                    named.putIfAbsent( key, ( RankVocabulary.label( key ) != null ) ? RankVocabulary.label( key ) : raw.trim() );
                }
            }
            if ( n.getNodeData().getProperties() != null ) {
                for( final Property p : n.getNodeData().getProperties().getProperties() ) {
                    final String key = RankVocabulary.propertyRank( p.getRef() );
                    // a scheme of the file's own counts only when a TIP carries it (an internal-only branch label,
                    // the Auspice reader's nextstrain:clade_label, is no rank); a listed rank counts wherever it is
                    if ( ( key != null ) && ( RankVocabulary.isListed( key ) || n.isExternal() ) ) {
                        named.putIfAbsent( key, ( RankVocabulary.label( key ) != null ) ? RankVocabulary.label( key )
                                : RankVocabulary.localName( p.getRef() ) );
                    }
                }
            }
        }
        if ( source.isOnline() ) {
            for( final String rank : OFFERABLE_RANKS ) {
                named.putIfAbsent( RankVocabulary.key( rank ), rank );
            }
        }
        final List<String> keys = new ArrayList<String>( named.keySet() );
        Collections.sort( keys, RankVocabulary::compare );
        final Map<String, Integer> annotated = AptxUtil.getRankCounts( phy );
        final int tips = phy.getNumberOfExternalNodes();
        for( final String key : keys ) {
            final String label = named.get( key );
            final String rank = ( RankVocabulary.label( key ) != null ) ? RankVocabulary.label( key ) : key;
            final Map<PhylogenyNode, TaxonPlacement> placements = TreePanelUtil.placeNodes( phy, rank, svc, source );
            final int placed = TreePanelUtil.placedTipCount( placements );
            final int taxa = TreePanelUtil.placedTaxonCount( placements );
            final int to_look_up = source.isOnline() ? TreePanelUtil.unresolvedTipTaxa( phy, rank, svc, source ).size() : 0;
            final Integer ann = annotated.get( rank );
            final boolean repeats = ( placed > 0 ) && ( ( taxa * 5 ) <= ( placed * 3 ) );
            if ( repeats || ( source.isOnline() && ( to_look_up > 0 ) ) ) {
                out.add( new Coverage( key, label, ( ann == null ) ? 0 : ann, placed, tips, taxa, to_look_up ) );
            }
        }
        return out;
    }

    /** The rank entries for {@code source} ({@link #coverage}), each reading {@code "family (19 annotated; 29%
     *  placed)"}, and for an online source {@code "family (19 annotated; 29% placed now; 136 to look up)"}. */
    static String[] choices( final Phylogeny phy, final TaxonomicLineageService svc, final TaxonSource source ) {
        final List<String> out = new ArrayList<String>();
        for( final Coverage c : coverage( phy, svc, source ) ) {
            out.add( entry( c.label, c.annotated, c.placed, c.tips, c.to_look_up, source.isOnline() ) );
        }
        return out.toArray( new String[ 0 ] );
    }

    /** The rank to open a chooser on (the joint rule, Archaeopteryx.js {@code defaultRank}): the FINEST listed rank
     *  naming 2 to {@code max_taxa} taxa (a legend's worth: a genus level on a family tree); else, among the file's
     *  own clade / lineage schemes, the one naming the most taxa within that cap; else the rank naming the most
     *  taxa; else null -- always one of the offered ranks. Returns the entry's label. */
    static String defaultRank( final List<Coverage> coverage, final int max_taxa ) {
        if ( ( coverage == null ) || coverage.isEmpty() ) {
            return null;
        }
        final int cap = ( max_taxa > 1 ) ? max_taxa : 20;
        for( int i = coverage.size() - 1; i >= 0; --i ) {
            final Coverage c = coverage.get( i );
            if ( RankVocabulary.isListed( c.key ) && ( c.taxa >= 2 ) && ( c.taxa <= cap ) ) {
                return c.label;
            }
        }
        Coverage scheme = null;
        for( final Coverage c : coverage ) {
            if ( !RankVocabulary.isListed( c.key ) && ( c.taxa >= 2 ) && ( c.taxa <= cap ) && ( ( scheme == null ) || ( c.taxa > scheme.taxa ) ) ) {
                scheme = c;
            }
        }
        if ( scheme != null ) {
            return scheme.label;
        }
        Coverage best = coverage.get( 0 );
        for( final Coverage c : coverage ) {
            if ( c.taxa > best.taxa ) {
                best = c;
            }
        }
        return best.label;
    }

    /** One chooser entry; see {@link #choices}. */
    static String entry( final String rank, final int annotated, final int placed, final int tips,
                         final int to_look_up, final boolean online ) {
        final StringBuilder sb = new StringBuilder( rank ).append( " (" )
                .append( ( annotated > 0 ) ? ( annotated + " annotated; " ) : "not annotated; " )
                .append( AptxUtil.coveragePercentLabel( placed, tips ) ).append( " placed" );
        if ( online ) {
            sb.append( " now; " ).append( to_look_up ).append( " to look up" );
        }
        return sb.append( ')' ).toString();
    }

    /** The bare rank of a chooser entry -- the text before its "(" -- or null for {@link #NONE} and for the first
     *  entry a caller supplied. */
    static String bareRank( final String entry ) {
        if ( ForesterUtil.isEmpty( entry ) || NONE.equals( entry ) ) {
            return null;
        }
        final int paren = entry.indexOf( '(' );
        final String bare = ( paren > 0 ) ? entry.substring( 0, paren ).trim() : entry.trim();
        return bare.isEmpty() ? null : bare;
    }

    /** Rebuilds every chooser's entries for the selected source, keeping each chooser's selected rank. */
    void rebuildChoices() {
        final TaxonSource source = selectedSource();
        _write_lineages.setEnabled( source.isOnline() ); // only an online source has something looked up to keep
        final String[] choices = choices( _phy, _svc, source );
        reselect( _rank_box, withFirst( choices, _first_entry ) );
        for( final JComboBox<String> box : _extra_boxes ) {
            reselect( box, withFirst( choices, NONE ) );
        }
        for( final Runnable r : _listeners ) {
            r.run();
        }
    }

    private void reselect( final JComboBox<String> box, final String[] entries ) {
        final String was = (String) box.getSelectedItem();
        final String bare = ( ( _first_entry != null ) && _first_entry.equals( was ) ) ? null : bareRank( was );
        box.setModel( new DefaultComboBoxModel<String>( entries ) );
        // a selected first entry / NONE has no bare rank and sits at index 0 of the new model too
        final int idx = AptxUtil.indexOfRank( entries, bare );
        box.setSelectedIndex( ( idx >= 0 ) ? idx : 0 );
    }

    JComboBox<String> rankBox() {
        return _rank_box;
    }

    List<JComboBox<String>> extraRankBoxes() {
        return _extra_boxes;
    }

    JPanel sourcePanel() {
        return _source_panel;
    }

    /** The "write each looked-up lineage into the tree" checkbox: on by default, enabled for an online source. */
    JCheckBox writeLineagesBox() {
        return _write_lineages;
    }

    /** True when the lineages the online source looked up are to be written into the tree as metadata fields --
     *  never for the tree source, which looks nothing up. */
    boolean isWriteLineages() {
        return selectedSource().isOnline() && _write_lineages.isSelected();
    }

    TaxonSource selectedSource() {
        for( final TaxonSource s : TaxonSource.values() ) {
            if ( _radios[ s.ordinal() ].isSelected() ) {
                return s;
            }
        }
        return TaxonSource.TREE;
    }

    void setSelectedSource( final TaxonSource s ) {
        _radios[ s.ordinal() ].setSelected( true );
        rebuildChoices();
    }

    /** The main chooser's bare rank, or null when its first entry is selected. */
    String selectedRank() {
        final String item = (String) _rank_box.getSelectedItem();
        if ( ( _first_entry != null ) && _first_entry.equals( item ) ) {
            return null;
        }
        return bareRank( item );
    }

    /** True when the caller-supplied first entry (e.g. "stop drawing") is selected. */
    boolean isFirstEntrySelected() {
        return ( _first_entry != null ) && _first_entry.equals( _rank_box.getSelectedItem() );
    }

    /** The bare rank of extra chooser {@code i}, or null for {@link #NONE}. */
    String selectedExtraRank( final int i ) {
        return bareRank( (String) _extra_boxes.get( i ).getSelectedItem() );
    }

    /** Selects {@code bare_rank} in the main chooser; a rank it does not offer leaves the selection alone. */
    void selectRank( final String bare_rank ) {
        final String[] entries = new String[ _rank_box.getItemCount() ];
        for( int i = 0; i < entries.length; ++i ) {
            entries[ i ] = _rank_box.getItemAt( i );
        }
        final int idx = AptxUtil.indexOfRank( entries, bare_rank );
        if ( idx >= 0 ) {
            _rank_box.setSelectedIndex( idx );
        }
    }

    void setSelectedIndex( final int idx ) {
        _rank_box.setSelectedIndex( idx );
    }

    /** Runs {@code r} after every rebuild (a source change). */
    void addSourceListener( final Runnable r ) {
        _listeners.add( r );
    }

    JRadioButton radioForTest( final TaxonSource s ) {
        return _radios[ s.ordinal() ];
    }
}
