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
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.io.File;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.swing.AbstractAction;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JToolBar;
import javax.swing.JViewport;
import javax.swing.KeyStroke;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.forester.phylogeny.Phylogeny;

/**
 * A stand-alone window hosting a {@link TanglegramPanel}. Deliberately a separate top-level frame (not a tab): the
 * main workspace's tab machinery is index-coupled to {@code TreePanel}, so a non-{@code TreePanel} tab would desync
 * it. A separate window also fits a comparison view -- it carries only a minimal toolbar (undo/redo of rotations,
 * zoom/fit) plus a live summary of the link (field, connectors, crossings, unmatched), and none of the tree-editing
 * control panel. Clicking a clade's vertical bar flips it to untangle; the crossing count updates live.
 */
final class TanglegramFrame extends JFrame {

    private static final long     serialVersionUID = 1L;
    private final TanglegramPanel _panel;
    private final String          _left_name;
    private final String          _right_name;
    private final JLabel          _summary         = new JLabel();
    private final JButton         _undo_button     = new JButton( "Undo" );
    private final JButton         _redo_button     = new JButton( "Redo" );
    private final JCheckBox       _phylogram_cb    = new JCheckBox( "Branch lengths" );
    private JComboBox<String>     _color_selector;

    TanglegramFrame( final Phylogeny left, final Phylogeny right, final TanglegramLinker.LinkField field,
                     final String left_name, final String right_name ) {
        this( left, right, field, field, null, left_name, right_name );
    }

    /** Link the two trees on possibly-DIFFERENT fields holding the same value (see {@link TanglegramLinker#link}). */
    TanglegramFrame( final Phylogeny left, final Phylogeny right, final TanglegramLinker.LinkField left_field,
                     final TanglegramLinker.LinkField right_field, final String left_name, final String right_name ) {
        this( left, right, left_field, right_field, null, left_name, right_name );
    }

    /** Link the two trees through an external association table (different names per tree; see
     *  {@link TanglegramLinker#linkByAssociation}); a null {@code associations} falls back to the value join. */
    TanglegramFrame( final Phylogeny left, final Phylogeny right, final TanglegramLinker.LinkField left_field,
                     final TanglegramLinker.LinkField right_field, final Map<String, List<String>> associations,
                     final String left_name, final String right_name ) {
        super( "Tanglegram: " + left_name + "  ↔  " + right_name );
        _left_name = left_name;
        _right_name = right_name;
        _panel = new TanglegramPanel( left, right, left_field, right_field, associations );
        _panel.setChangeListener( this::refresh );
        final JScrollPane scroll = new JScrollPane( _panel );
        scroll.getViewport().setScrollMode( JViewport.SIMPLE_SCROLL_MODE );
        setLayout( new BorderLayout() );
        add( buildToolbar(), BorderLayout.NORTH );
        add( scroll, BorderLayout.CENTER );
        installUndoRedoKeys();
        refresh();
        setDefaultCloseOperation( DISPOSE_ON_CLOSE );
        setSize( 1000, 700 );
        setLocationRelativeTo( null );
    }

    TanglegramPanel getTanglegramPanel() {
        return _panel;
    }

    // ---- test hooks --------------------------------------------------------------------------------------------

    boolean isUndoButtonEnabledForTest() {
        return _undo_button.isEnabled();
    }

    boolean isRedoButtonEnabledForTest() {
        return _redo_button.isEnabled();
    }

    String summaryTextForTest() {
        return _summary.getText();
    }

    void clickUndoForTest() {
        _undo_button.doClick();
    }

    void clickRedoForTest() {
        _redo_button.doClick();
    }

    int colorItemCountForTest() {
        return _color_selector.getItemCount();
    }

    void selectColorForTest( final int index ) {
        _color_selector.setSelectedIndex( index );
    }

    boolean isPhylogramCheckboxEnabledForTest() {
        return _phylogram_cb.isEnabled();
    }

    void clickPhylogramForTest() {
        _phylogram_cb.doClick();
    }

    private JToolBar buildToolbar() {
        final JToolBar bar = new JToolBar();
        bar.setFloatable( false );
        _undo_button.setToolTipText( "Undo the last clade flip" );
        _redo_button.setToolTipText( "Redo" );
        _undo_button.addActionListener( e -> _panel.undo() );
        _redo_button.addActionListener( e -> _panel.redo() );
        final JButton untangle = new JButton( "Auto-untangle" );
        untangle.setToolTipText( "Flip clades in both trees to reduce connector crossings" );
        untangle.addActionListener( e -> _panel.autoUntangle() );
        final JButton zoom_in = new JButton( "Zoom In" );
        final JButton zoom_out = new JButton( "Zoom Out" );
        final JButton fit = new JButton( "Fit" );
        zoom_in.addActionListener( e -> _panel.zoomIn() );
        zoom_out.addActionListener( e -> _panel.zoomOut() );
        fit.addActionListener( e -> _panel.fit() );
        // aligned-phylogram toggle: only meaningful when a tree carries branch lengths, else greyed out
        _phylogram_cb.setToolTipText( _panel.hasBranchLengths()
                ? "Draw the branches to scale (aligned phylogram) instead of a lined-up cladogram"
                : "No branch lengths in these trees -- a cladogram is the only option" );
        _phylogram_cb.setEnabled( _panel.hasBranchLengths() );
        _phylogram_cb.addActionListener( e -> _panel.setPhylogram( _phylogram_cb.isSelected() ) );
        final JButton export = new JButton( "Export..." );
        export.setToolTipText( "Export the tanglegram figure as PDF, SVG, EPS, or PNG" );
        export.addActionListener( e -> exportFigure() );
        bar.add( untangle );
        bar.addSeparator();
        bar.add( _undo_button );
        bar.add( _redo_button );
        bar.addSeparator();
        bar.add( zoom_in );
        bar.add( zoom_out );
        bar.add( fit );
        bar.add( _phylogram_cb );
        bar.add( export );
        bar.addSeparator();
        bar.add( new JLabel( "Colour: " ) );
        bar.add( buildColorSelector() );
        bar.addSeparator();
        bar.add( _summary );
        bar.add( Box.createHorizontalGlue() );
        return bar;
    }

    /** The "Colour connectors by:" selector: Uniform, Crossings, then one entry per available tip attribute. */
    private JComboBox<String> buildColorSelector() {
        final List<TanglegramColoring.Field> fields = _panel.availableColorFields();
        _color_selector = new JComboBox<>();
        _color_selector.addItem( "Uniform" );
        _color_selector.addItem( "Crossings" );
        for( final TanglegramColoring.Field field : fields ) {
            _color_selector.addItem( field.label() );
        }
        _color_selector
                .setToolTipText( "Colour the connectors uniformly, highlight the crossings, or colour by a tip attribute" );
        _color_selector.setMaximumSize( _color_selector.getPreferredSize() ); // don't let the toolbar stretch it
        _color_selector.addActionListener( e -> {
            final int index = _color_selector.getSelectedIndex();
            if ( index <= 0 ) {
                _panel.setUniformColoring();
            }
            else if ( index == 1 ) {
                _panel.setCrossingColoring();
            }
            else {
                _panel.setFieldColoring( fields.get( index - 2 ) );
            }
        } );
        return _color_selector;
    }

    /** Save the tanglegram as PDF / SVG / EPS / PNG via a file chooser (the format follows the chosen filter). */
    private void exportFigure() {
        final JFileChooser chooser = new JFileChooser();
        chooser.setDialogTitle( "Export Tanglegram" );
        chooser.setAcceptAllFileFilterUsed( false );
        final Map<FileFilter, TanglegramExporter.Format> by_filter = new LinkedHashMap<>();
        for( final TanglegramExporter.Format format : TanglegramExporter.Format.values() ) {
            final FileNameExtensionFilter filter = new FileNameExtensionFilter( format.label(), format.extension() );
            chooser.addChoosableFileFilter( filter );
            by_filter.put( filter, format );
        }
        if ( chooser.showSaveDialog( this ) != JFileChooser.APPROVE_OPTION ) {
            return;
        }
        final TanglegramExporter.Format format = by_filter.getOrDefault( chooser.getFileFilter(),
                                                                         TanglegramExporter.Format.PDF );
        File file = chooser.getSelectedFile();
        if ( !file.getName().toLowerCase().endsWith( "." + format.extension() ) ) {
            file = new File( file.getParentFile(), file.getName() + "." + format.extension() );
        }
        if ( file.exists() && ( JOptionPane.showConfirmDialog( this, file.getName() + " already exists. Overwrite?",
                                                              "Export Tanglegram",
                                                              JOptionPane.YES_NO_OPTION ) != JOptionPane.YES_OPTION ) ) {
            return;
        }
        final Dimension size = _panel.exportSize();
        try {
            final String msg = TanglegramExporter.write( file, format, _panel, size.width, size.height );
            JOptionPane.showMessageDialog( this, "Wrote " + msg, "Export Tanglegram",
                                           JOptionPane.INFORMATION_MESSAGE );
        }
        catch ( final Throwable ex ) {
            // Throwable, not Exception: a large export can throw OutOfMemoryError (an Error), which must not escape
            // onto the EDT -- report it instead.
            JOptionPane.showMessageDialog( this, "Could not export the tanglegram: " + ex.getMessage(),
                                           "Export Tanglegram", JOptionPane.ERROR_MESSAGE );
        }
    }

    private void installUndoRedoKeys() {
        final int mod = Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx();
        final int wifw = JComponent.WHEN_IN_FOCUSED_WINDOW;
        getRootPane().getInputMap( wifw ).put( KeyStroke.getKeyStroke( KeyEvent.VK_Z, mod ), "tangle-undo" );
        getRootPane().getInputMap( wifw )
                .put( KeyStroke.getKeyStroke( KeyEvent.VK_Z, mod | InputEvent.SHIFT_DOWN_MASK ), "tangle-redo" );
        getRootPane().getActionMap().put( "tangle-undo", new AbstractAction() {

            @Override
            public void actionPerformed( final ActionEvent e ) {
                _panel.undo();
            }
        } );
        getRootPane().getActionMap().put( "tangle-redo", new AbstractAction() {

            @Override
            public void actionPerformed( final ActionEvent e ) {
                _panel.redo();
            }
        } );
    }

    /** Refreshes the summary text (crossings change on each flip) and the undo/redo button state. */
    private void refresh() {
        _summary.setText( summary() );
        _undo_button.setEnabled( _panel.canUndo() );
        _redo_button.setEnabled( _panel.canRedo() );
    }

    private String summary() {
        return "Link: " + linkDescription() + "   •   " + _panel.getResult().getLinks().size()
                + " connectors, " + _panel.getCrossingCount() + " crossings (entanglement "
                + String.format( java.util.Locale.ROOT, "%.2f", _panel.getEntanglement() ) + ")   •   " + _left_name
                + " (" + _panel.getLeftTipCount() + " tips) ↔ " + _right_name + " (" + _panel.getRightTipCount()
                + " tips)   •   " + _panel.getUnmatchedCount() + " unmatched";
    }

    /** One field label when both trees link on the same field, else "leftLabel ↔ rightLabel"; notes an association
     *  file when the trees were linked through a mapping table rather than by equal value. */
    private String linkDescription() {
        final TanglegramLinker.LinkField left = _panel.getLeftField();
        final TanglegramLinker.LinkField right = _panel.getRightField();
        final String fields = ( left == right ) ? left.label() : ( left.label() + " ↔ " + right.label() );
        return _panel.isAssociationLinked() ? ( fields + " via association file" ) : fields;
    }
}
