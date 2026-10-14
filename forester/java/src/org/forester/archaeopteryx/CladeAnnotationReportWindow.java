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
import java.awt.Frame;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.function.Consumer;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.table.DefaultTableModel;

/**
 * The modeless window that shows a {@link CladeAnnotationReport}: the per-rank summary on top, the rows as a
 * sortable table, and an "Export TSV…" button. Opened by "Annotate Clades by Rank" (report-only, or from the
 * result dialog's "Show Report…") and by "Colorize Subtrees via Taxonomic Rank". Reading only; it never touches
 * the tree.
 */
final class CladeAnnotationReportWindow extends JFrame {

    private static final long      serialVersionUID = 1L;
    private final CladeAnnotationReport _report;
    private final JTable           _table;
    private final JTextArea        _summary;
    private final File             _initial_dir;
    private final Consumer<File>   _remember_dir;

    /**
     * @param owner        the main window (for placement); may be null
     * @param title        the window title
     * @param report       what to show
     * @param initial_dir  where the export file chooser opens; may be null
     * @param remember_dir told the directory of a written export; may be null
     */
    CladeAnnotationReportWindow( final Frame owner, final String title, final CladeAnnotationReport report,
                                 final File initial_dir, final Consumer<File> remember_dir ) {
        super( title );
        _report = report;
        _initial_dir = initial_dir;
        _remember_dir = remember_dir;
        setDefaultCloseOperation( DISPOSE_ON_CLOSE );
        _summary = new JTextArea( report.getSummary() );
        _summary.setEditable( false );
        _summary.setLineWrap( true );
        _summary.setWrapStyleWord( true );
        _summary.setFont( FormWidgets.monoFont( _summary.getFont() ) );
        _summary.setRows( Math.max( 2, Math.min( 6, report.getSummaryLines().size() + 1 ) ) );
        final JScrollPane summary_scroll = new JScrollPane( _summary );
        final DefaultTableModel model = new DefaultTableModel( CladeAnnotationReport.COLUMNS, 0 ) {

            private static final long serialVersionUID = 1L;

            @Override
            public boolean isCellEditable( final int row, final int column ) {
                return false;
            }
        };
        for( final CladeAnnotationReport.Row r : report.getRows() ) {
            model.addRow( new Object[] { r.getRank(), r.getNode(), r.getTaxon(), r.getSource(), r.getNote() } );
        }
        _table = new JTable( model );
        _table.setAutoCreateRowSorter( true );
        _table.setFillsViewportHeight( true );
        _table.getColumnModel().getColumn( 0 ).setPreferredWidth( 70 );
        _table.getColumnModel().getColumn( 1 ).setPreferredWidth( 220 );
        _table.getColumnModel().getColumn( 2 ).setPreferredWidth( 160 );
        _table.getColumnModel().getColumn( 3 ).setPreferredWidth( 200 );
        _table.getColumnModel().getColumn( 4 ).setPreferredWidth( 360 );
        final JScrollPane table_scroll = new JScrollPane( _table );
        table_scroll.setPreferredSize( new Dimension( 1040, 420 ) );
        final JButton export = new JButton( "Export TSV…" );
        export.setToolTipText( "Save the table as a tab-separated file" );
        export.addActionListener( e -> exportInteractively() );
        final JButton close = new JButton( "Close" );
        close.addActionListener( e -> dispose() );
        final JPanel buttons = new JPanel( new FlowLayout( FlowLayout.RIGHT ) );
        buttons.add( export );
        buttons.add( close );
        getContentPane().setLayout( new BorderLayout( 6, 6 ) );
        getContentPane().add( summary_scroll, BorderLayout.NORTH );
        getContentPane().add( table_scroll, BorderLayout.CENTER );
        getContentPane().add( buttons, BorderLayout.SOUTH );
        pack();
        setLocationRelativeTo( owner );
    }

    private void exportInteractively() {
        final JFileChooser fc = new JFileChooser();
        if ( _initial_dir != null ) {
            fc.setCurrentDirectory( _initial_dir );
        }
        fc.setSelectedFile( new File( _initial_dir, "clade-annotation-report.tsv" ) );
        if ( fc.showSaveDialog( this ) != JFileChooser.APPROVE_OPTION ) {
            return;
        }
        final File f = fc.getSelectedFile();
        if ( f.exists() && ( JOptionPane.showConfirmDialog( this, "Overwrite " + f.getName() + "?", "Export TSV",
                JOptionPane.OK_CANCEL_OPTION ) != JOptionPane.OK_OPTION ) ) {
            return;
        }
        try {
            exportTo( f );
            if ( ( _remember_dir != null ) && ( f.getParentFile() != null ) ) {
                _remember_dir.accept( f.getParentFile() );
            }
        }
        catch ( final IOException e ) {
            JOptionPane.showMessageDialog( this, "Could not write " + f + ":\n" + e.getMessage(), "Export TSV",
                                           JOptionPane.ERROR_MESSAGE );
        }
    }

    /** Writes the report as TSV to {@code f}. */
    void exportTo( final File f ) throws IOException {
        Files.write( f.toPath(), _report.toTsv().getBytes( StandardCharsets.UTF_8 ) );
    }

    CladeAnnotationReport getReport() {
        return _report;
    }

    JTable tableForTest() {
        return _table;
    }

    String summaryForTest() {
        return _summary.getText();
    }
}
