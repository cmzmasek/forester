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

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.DefaultTableModel;

import org.forester.archaeopteryx.tools.TipDateExtractor;
import org.forester.archaeopteryx.tools.TipDateExtractor.DayMonthOrder;
import org.forester.archaeopteryx.tools.TipDateExtractor.Summary;
import org.forester.archaeopteryx.tools.TipDateExtractor.TipDate;
import org.forester.phylogeny.Phylogeny;

/**
 * The "Extract Dates from Labels…" preview dialog: shows what date {@link TipDateExtractor} finds in each tip label
 * (a per-tip table + a roll-up summary) BEFORE anything is written, so the user can confirm the auto-detection. When
 * an ambiguous numeric date (e.g. {@code 05/03/2021}) is present it exposes a <b>day-first / month-first</b> toggle
 * (re-previewed live); a "skip tips that already have a date" checkbox is always shown. Apply returns the chosen
 * options to {@link MainFrame} (which writes the dates, undoably); Cancel writes nothing.
 */
final class TipDateExtractionDialog extends JDialog {

    private final Phylogeny      _phylogeny;
    private DayMonthOrder        _order          = DayMonthOrder.DAY_FIRST;
    private boolean              _skip_existing  = true;
    private boolean              _applied        = false;
    private final JLabel         _summary_label  = new JLabel();
    private JPanel               _order_row;     // the day/month toggle row -- shown only when a date is ambiguous
    private final DefaultTableModel _table_model = new DefaultTableModel( new Object[] { "Tip label", "Date found",
            "Decimal year" }, 0 ) {
        @Override
        public boolean isCellEditable( final int r, final int c ) {
            return false;
        }
    };

    TipDateExtractionDialog( final MainFrame parent, final Phylogeny phylogeny ) {
        super( parent, "Extract Dates from Labels", true );
        _phylogeny = phylogeny;
        buildUi();
        refreshPreview();
        pack();
        setSize( new Dimension( Math.max( 560, getWidth() ), Math.min( 560, Math.max( 380, getHeight() ) ) ) );
        setLocationRelativeTo( parent );
    }

    private void buildUi() {
        final JPanel content = new JPanel( new BorderLayout( 8, 8 ) );
        content.setBorder( BorderFactory.createEmptyBorder( 10, 12, 10, 12 ) );

        _summary_label.setBorder( BorderFactory.createEmptyBorder( 0, 0, 4, 0 ) );
        content.add( _summary_label, BorderLayout.NORTH );

        final JTable table = new JTable( _table_model );
        table.setRowSelectionAllowed( false );
        table.getTableHeader().setReorderingAllowed( false );
        final JScrollPane sp = new JScrollPane( table );
        sp.setPreferredSize( new Dimension( 540, 300 ) );
        content.add( sp, BorderLayout.CENTER );

        // controls: the adaptive day/month toggle + the skip-existing checkbox
        final JPanel controls = new JPanel();
        controls.setLayout( new BoxLayout( controls, BoxLayout.PAGE_AXIS ) );

        _order_row = new JPanel( new FlowLayout( FlowLayout.LEFT, 6, 0 ) );
        _order_row.add( new JLabel( "Day/month order:" ) );
        final JComboBox<String> order_combo = new JComboBox<>( new String[] { "Day first (e.g. 05/03 = 5 March)",
                "Month first (e.g. 05/03 = 3 May)" } );
        order_combo.addActionListener( e -> {
            _order = ( order_combo.getSelectedIndex() == 1 ) ? DayMonthOrder.MONTH_FIRST : DayMonthOrder.DAY_FIRST;
            refreshPreview();
        } );
        _order_row.add( order_combo );
        // visibility (only when a genuinely ambiguous numeric date is present) is set from the preview in refreshPreview
        controls.add( _order_row );

        final JCheckBox skip = new JCheckBox( "Skip tips that already have a date", true );
        skip.addActionListener( e -> {
            _skip_existing = skip.isSelected();
            refreshPreview(); // the effective "will write" count changes -> update the summary + skip marks
        } );
        final JPanel skip_row = new JPanel( new FlowLayout( FlowLayout.LEFT, 6, 0 ) );
        skip_row.add( skip );
        controls.add( skip_row );

        final JPanel south = new JPanel( new BorderLayout() );
        south.add( controls, BorderLayout.CENTER );

        final JPanel buttons = new JPanel( new FlowLayout( FlowLayout.RIGHT ) );
        final JButton cancel = new JButton( "Cancel" );
        cancel.addActionListener( e -> {
            _applied = false;
            dispose();
        } );
        final JButton apply = new JButton( "Apply" );
        apply.addActionListener( e -> {
            _applied = true;
            dispose();
        } );
        buttons.add( cancel );
        buttons.add( apply );
        south.add( buttons, BorderLayout.SOUTH );
        content.add( south, BorderLayout.SOUTH );

        setContentPane( content );
        getRootPane().setDefaultButton( apply );
    }

    /** Re-run the preview for the current day/month order and rebuild the table + summary. */
    private void refreshPreview() {
        final List<TipDate> rows = TipDateExtractor.preview( _phylogeny, _order );
        _table_model.setRowCount( 0 );
        boolean any_ambiguous = false;
        int will_write = 0;
        for ( final TipDate t : rows ) {
            if ( t.match() == null ) {
                _table_model.addRow( new Object[] { t.label(), "—", "—" } );
                continue;
            }
            any_ambiguous |= t.match().ambiguous();
            final boolean skipped = _skip_existing && t.alreadyDated();
            if ( !skipped ) {
                will_write++;
            }
            _table_model.addRow( new Object[] { t.label(), t.match().matchedText(),
                    TipDateExtractor.decimalYearString( t.match().decimalYear() )
                            + ( skipped ? "  (already dated → skip)" : "" ) } );
        }
        // the day/month toggle is relevant only when a genuinely ambiguous numeric date is present (derived from the
        // preview we just built, so no separate scan)
        _order_row.setVisible( any_ambiguous );
        _summary_label.setText( summaryText( TipDateExtractor.summarize( rows ), will_write ) );
    }

    /** The header HTML: detected format, matched/unmatched counts, how many will actually be written (some matched tips
     *  may already have a date and be skipped), and the matched date range. */
    static String summaryText( final Summary s, final int will_write ) {
        final StringBuilder sb = new StringBuilder( "<html>" );
        if ( s.matched() == 0 ) {
            return sb.append( "No dates were recognised in the tip labels.</html>" ).toString();
        }
        sb.append( "Detected format: <b>" ).append( s.dominantFormat() ).append( "</b> &nbsp;&middot;&nbsp; matched <b>" )
                .append( s.matched() ).append( "</b> / " ).append( s.matched() + s.unmatched() ).append( " tips" );
        if ( s.ambiguous() > 0 ) {
            sb.append( " &nbsp;&middot;&nbsp; " ).append( s.ambiguous() )
                    .append( " ambiguous (using the day/month setting below)" );
        }
        if ( will_write < s.matched() ) {
            sb.append( " &nbsp;&middot;&nbsp; will write <b>" ).append( will_write ).append( "</b> (" )
                    .append( s.matched() - will_write ).append( " already dated, skipped)" );
        }
        sb.append( "<br>Range: <b>" ).append( TipDateExtractor.decimalYearString( s.minYear() ) ).append( "</b> &ndash; <b>" )
                .append( TipDateExtractor.decimalYearString( s.maxYear() ) ).append( "</b> (decimal years)" );
        return sb.append( "</html>" ).toString();
    }

    boolean isApplied() {
        return _applied;
    }

    DayMonthOrder getOrder() {
        return _order;
    }

    boolean isSkipExisting() {
        return _skip_existing;
    }

    /** Test hook: the current summary header text (HTML). */
    String summaryTextForTest() {
        return _summary_label.getText();
    }
}
