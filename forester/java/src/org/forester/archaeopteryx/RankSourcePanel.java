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
import java.util.Map;
import java.util.Set;

import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterUtil;
import org.forester.util.TaxonomyUtil;
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
        _rank_box = new JComboBox<String>( withFirst( choices( phy, svc, initial ), first_entry ) );
        for( int i = 0; i < extra_levels; ++i ) {
            final JComboBox<String> box = new JComboBox<String>( withFirst( choices( phy, svc, initial ), NONE ) );
            box.setSelectedIndex( 0 );
            _extra_boxes.add( box );
        }
    }

    private static String tooltipOf( final TaxonSource s ) {
        switch ( s ) {
            case TREE:
                return "Only what the tree already carries: a tip's own taxonomy at the rank, a per-rank metadata field"
                        + " (ncbi:family, gtdb:family), else the nearest ancestor annotated at the rank. Nothing is looked up.";
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

    /**
     * The rank entries for {@code source}, in canonical rank order: every rank some node is annotated at, plus the
     * major Linnaean ranks ({@link #OFFERABLE_RANKS}) even when none is -- a database can place a tip at a rank no
     * node literally carries. Each reads {@code "family (19 annotated; 29% placed)"}, and for an online source
     * {@code "family (19 annotated; 29% placed now; 136 to look up)"}.
     */
    static String[] choices( final Phylogeny phy, final TaxonomicLineageService svc, final TaxonSource source ) {
        final List<String> out = new ArrayList<String>();
        if ( ( phy == null ) || phy.isEmpty() ) {
            return new String[ 0 ];
        }
        final Map<String, Integer> annotated = AptxUtil.getRankCounts( phy );
        final int tips = phy.getNumberOfExternalNodes();
        for( final String rank : TaxonomyUtil.TAXONOMY_RANKS_LIST ) {
            if ( rank.equals( TaxonomyUtil.UNKNOWN ) || rank.equals( TaxonomyUtil.OTHER ) ) {
                continue;
            }
            final Integer n = annotated.get( rank );
            final boolean present = ( n != null ) && ( n > 0 );
            if ( !present && !OFFERABLE_RANKS.contains( rank ) ) {
                continue;
            }
            final int placed = TreePanelUtil.placedTipCount( TreePanelUtil.placeNodes( phy, rank, svc, source ) );
            final int to_look_up = source.isOnline() ? TreePanelUtil.unresolvedTipTaxa( phy, rank, svc, source ).size() : 0;
            out.add( entry( rank, present ? n : 0, placed, tips, to_look_up, source.isOnline() ) );
        }
        return out.toArray( new String[ 0 ] );
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
