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

/**
 * The "Display Data" options offered as checkboxes in the left control panel. Replaces the old
 * positional {@code String[][] Configuration.display_options} array (title + a stringly-typed
 * "yes"/"no"/"?" default) and its parallel {@code final static int} index constants: each option is
 * now a type-safe enum constant carrying its display title and its {@link DefaultState}. Because the
 * identity is an enum, adding or removing an option no longer shifts any index (the fragility the old
 * array had, guarded by a now-removed alignment test).
 */
public enum DisplayOption {

    DISPLAY_AS_PHYLOGRAM( "Phylogram", DefaultState.GUESS ),
    SHOW_NODE_NAMES( "Node Name", DefaultState.ON ),
    SHOW_TAX_CODE( "Taxonomy Code", DefaultState.ON ),
    WRITE_CONFIDENCE_VALUES( "Confidence Values", DefaultState.GUESS ),
    WRITE_EVENTS( "Node Events", DefaultState.GUESS ),
    USE_STYLE( "Visual Styles/Branch Colors", DefaultState.ON ),
    WIDTH_BRANCHES( "Branch Widths", DefaultState.OFF ),
    SHOW_DOMAIN_ARCHITECTURES( "Domain Architectures", DefaultState.OFF ),
    SHOW_SEQ_NAMES( "Seq Name", DefaultState.OFF ),
    SHOW_SEQUENCE_ACC( "Seq Accession", DefaultState.OFF ),
    DISPLAY_INTERNAL_DATA( "Show Internal Data", DefaultState.ON ),
    // "Auto-hide Labels": renamed from the historical "Dyna Hide", which named the mechanism (dynamic hiding)
    // rather than what the user sees (labels disappearing on their own when the tree is too dense)
    DYNAMICALLY_HIDE_DATA( "Auto-hide Labels", DefaultState.ON ),
    SHOW_TAXONOMY_SCIENTIFIC_NAMES( "Taxonomy Scientific", DefaultState.ON ),
    SHOW_TAXONOMY_COMMON_NAMES( "Taxonomy Common", DefaultState.OFF ),
    SHOW_SEQ_SYMBOLS( "Seq Symbol", DefaultState.OFF ),
    NODE_DATA_POPUP( "Rollover", DefaultState.ON ),
    SHOW_PROPERTIES( "Properties", DefaultState.OFF ),
    SHOW_GENE_NAMES( "Gene Name", DefaultState.OFF ),
    WRITE_BRANCH_LENGTH_VALUES( "Branch Length Values", DefaultState.OFF ),
    SHOW_TAX_RANK( "Taxonomy Rank", DefaultState.OFF ),
    DISPLAY_EXTERNAL_DATA( "Show External Data", DefaultState.ON ),
    SHORTEN_LABELS( "Shorten Labels", DefaultState.ON );

    /** The default state of a display checkbox: ON / OFF, or GUESS (derived from the tree at load
     *  time, e.g. phylogram-vs-cladogram or whether confidence/event data is present). Replaces the
     *  old "yes" / "no" / "?" third column. */
    public enum DefaultState {
        ON, OFF, GUESS
    }

    private final String       _title;
    private final DefaultState _default_state;

    DisplayOption( final String title, final DefaultState default_state ) {
        _title = title;
        _default_state = default_state;
    }

    /** The checkbox label. */
    public String title() {
        return _title;
    }

    public DefaultState defaultState() {
        return _default_state;
    }

    /** True if the checkbox is on by default (the old "yes"). */
    public boolean isCheckedByDefault() {
        return _default_state == DefaultState.ON;
    }

    /** True if the default is derived from the tree at load time (the old "?"). */
    public boolean isGuessedByDefault() {
        return _default_state == DefaultState.GUESS;
    }
}
