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
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JColorChooser;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.NodeVisualData.FontType;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;

/**
 * A small modal editor for the per-node visual style ({@link NodeVisualData}): font style/size/colour and node
 * mark shape/fill/size/colour. Each attribute has an "Apply" checkbox -- only ticked attributes are written -- so
 * a bulk edit over the selected/found set changes just the chosen attributes and leaves the rest of each node's
 * style untouched. Opened for ONE node from the "Node Style" click-to action (pre-filled from that node), or for
 * the whole selected/found set from Tools -> "Node Style for Selected Nodes" (all attributes start un-ticked). The
 * spec it builds is applied via {@link TreePanel#applyNodeStyleEdit(List, NodeStyleEditor.Spec)} (undo + provenance).
 */
final class NodeStyleDialog extends JDialog {

    private static final long serialVersionUID = 1L;

    private static final String[] FONT_STYLES = { "Plain", "Bold", "Italic", "Bold Italic" };
    private static final String[] SHAPES      = { "Default", "Circle", "Rectangle" };
    private static final String[] FILLS       = { "Default", "Solid", "None", "Gradient" };

    private final TreePanel           _tp;
    private final List<PhylogenyNode> _targets;

    private final JCheckBox        _apply_font_style = new JCheckBox();
    private final JCheckBox        _apply_font_size  = new JCheckBox();
    private final JCheckBox        _apply_font_color = new JCheckBox();
    private final JCheckBox        _apply_shape      = new JCheckBox();
    private final JCheckBox        _apply_fill       = new JCheckBox();
    private final JCheckBox        _apply_node_size  = new JCheckBox();
    private final JCheckBox        _apply_node_color = new JCheckBox();

    private final JComboBox<String> _font_style_combo = new JComboBox<>( FONT_STYLES );
    private final JSpinner          _font_size_spinner;
    private final JButton           _font_color_button = new JButton( "Choose..." );
    private final JComboBox<String> _shape_combo       = new JComboBox<>( SHAPES );
    private final JComboBox<String> _fill_combo        = new JComboBox<>( FILLS );
    private final JSpinner          _node_size_spinner;
    private final JButton           _node_color_button = new JButton( "Choose..." );

    private Color _font_color = Color.BLACK;
    private Color _node_color = Color.BLACK;

    NodeStyleDialog( final TreePanel tp, final List<PhylogenyNode> targets, final boolean single_node ) {
        super( SwingUtilities.getWindowAncestor( tp ), "Node Style" );
        setModalityType( ModalityType.APPLICATION_MODAL );
        _tp = tp;
        _targets = targets;
        final int default_font_size = 12; // the tree default; only the initial spinner value when "Apply size" is ticked
        final int default_node_size = ( tp.getOptions() != null ) ? tp.getOptions().getDefaultNodeShapeSize() : 5;
        _font_size_spinner = new JSpinner( new SpinnerNumberModel( clamp( default_font_size, 2, 64 ), 2, 64, 1 ) );
        _node_size_spinner = new JSpinner( new SpinnerNumberModel( clamp( default_node_size, 0, 100 ), 0, 100, 1 ) );

        final JPanel grid = new JPanel( new GridBagLayout() );
        grid.setBorder( BorderFactory.createEmptyBorder( 10, 12, 6, 12 ) );
        int row = 0;
        addSectionHeader( grid, row++, "Font" );
        addRow( grid, row++, _apply_font_style, "Style:", _font_style_combo );
        addRow( grid, row++, _apply_font_size, "Size:", _font_size_spinner );
        _font_color_button.addActionListener( e -> chooseColor( true ) );
        addRow( grid, row++, _apply_font_color, "Color:", _font_color_button );
        addSectionHeader( grid, row++, "Node mark" );
        addRow( grid, row++, _apply_shape, "Shape:", _shape_combo );
        addRow( grid, row++, _apply_fill, "Fill:", _fill_combo );
        addRow( grid, row++, _apply_node_size, "Size:", _node_size_spinner );
        _node_color_button.addActionListener( e -> chooseColor( false ) );
        addRow( grid, row++, _apply_node_color, "Color:", _node_color_button );

        // each Apply checkbox enables/disables its own control, so an un-ticked (won't-be-written) attribute reads
        // as inactive
        wireEnable( _apply_font_style, _font_style_combo );
        wireEnable( _apply_font_size, _font_size_spinner );
        wireEnable( _apply_font_color, _font_color_button );
        wireEnable( _apply_shape, _shape_combo );
        wireEnable( _apply_fill, _fill_combo );
        wireEnable( _apply_node_size, _node_size_spinner );
        wireEnable( _apply_node_color, _node_color_button );

        if ( single_node && ( targets != null ) && !targets.isEmpty() ) {
            prefillFrom( targets.get( 0 ) );
        }
        syncAllEnabled();

        final int count = ( targets == null ) ? 0 : targets.size();
        final JLabel header = new JLabel( single_node ? "Editing the clicked node." : ( "Editing " + count
                + ( count == 1 ? " selected / found node." : " selected / found nodes." ) ) );
        header.setBorder( BorderFactory.createEmptyBorder( 10, 12, 0, 12 ) );
        final JLabel hint = new JLabel(
                "<html><body style='width: 300px'><i>Tick an attribute to change it; unticked attributes are left"
                        + " as they are.</i></body></html>" );
        hint.setBorder( BorderFactory.createEmptyBorder( 2, 12, 6, 12 ) );
        final JPanel north = new JPanel( new BorderLayout() );
        north.add( header, BorderLayout.NORTH );
        north.add( hint, BorderLayout.SOUTH );

        final JButton ok = new JButton( "OK" );
        ok.addActionListener( e -> {
            applyAndClose();
        } );
        final JButton cancel = new JButton( "Cancel" );
        cancel.addActionListener( e -> setVisible( false ) );
        final JPanel south = new JPanel( new FlowLayout( FlowLayout.RIGHT ) );
        south.add( cancel );
        south.add( ok );

        setLayout( new BorderLayout() );
        add( north, BorderLayout.NORTH );
        add( grid, BorderLayout.CENTER );
        add( south, BorderLayout.SOUTH );
        getRootPane().setDefaultButton( ok );
        pack();
        setLocationRelativeTo( SwingUtilities.getWindowAncestor( tp ) );
    }

    private void applyAndClose() {
        final NodeStyleEditor.Spec spec = buildSpec();
        if ( !spec.isEmpty() ) {
            _tp.applyNodeStyleEdit( _targets, spec );
        }
        setVisible( false );
    }

    /** Builds the {@link NodeStyleEditor.Spec} from the current controls: each attribute is included only when its
     *  Apply checkbox is ticked (else null = leave unchanged). */
    NodeStyleEditor.Spec buildSpec() {
        final Integer font_style = _apply_font_style.isSelected()
                ? Integer.valueOf( fontStyleFromLabel( (String) _font_style_combo.getSelectedItem() ) ) : null;
        final Integer font_size = _apply_font_size.isSelected()
                ? Integer.valueOf( ( (Number) _font_size_spinner.getValue() ).intValue() ) : null;
        final Color font_color = _apply_font_color.isSelected() ? _font_color : null;
        final NodeShape shape = _apply_shape.isSelected() ? shapeFromLabel( (String) _shape_combo.getSelectedItem() )
                : null;
        final NodeFill fill = _apply_fill.isSelected() ? fillFromLabel( (String) _fill_combo.getSelectedItem() ) : null;
        final Float node_size = _apply_node_size.isSelected()
                ? Float.valueOf( ( (Number) _node_size_spinner.getValue() ).floatValue() ) : null;
        final Color node_color = _apply_node_color.isSelected() ? _node_color : null;
        return new NodeStyleEditor.Spec( font_style, font_size, font_color, shape, fill, node_size, node_color );
    }

    private void prefillFrom( final PhylogenyNode node ) {
        final NodeVisualData vis = node.getNodeData().getNodeVisualData();
        if ( vis == null ) {
            return;
        }
        if ( vis.getFontStyle() != FontType.PLAIN ) {
            _font_style_combo.setSelectedItem( labelForFontStyle( vis.getFontStyle() ) );
            _apply_font_style.setSelected( true );
        }
        if ( vis.getFontSize() > 0 ) {
            // expand the spinner range to fit the stored value rather than clamping it -- a clamp here would pre-tick
            // "Apply size" AND emit a SMALLER value, silently shrinking a font size the user never touched
            setSpinnerExpanding( _font_size_spinner, vis.getFontSize() );
            _apply_font_size.setSelected( true );
        }
        if ( vis.getFontColor() != null ) {
            setFontColor( vis.getFontColor() );
            _apply_font_color.setSelected( true );
        }
        if ( vis.getShape() != NodeShape.DEFAULT ) {
            _shape_combo.setSelectedItem( ( vis.getShape() == NodeShape.RECTANGLE ) ? "Rectangle" : "Circle" );
            _apply_shape.setSelected( true );
        }
        if ( vis.getFillType() != NodeFill.DEFAULT ) {
            _fill_combo.setSelectedItem( labelForFill( vis.getFillType() ) );
            _apply_fill.setSelected( true );
        }
        if ( vis.getSize() >= 0 ) {
            setSpinnerExpanding( _node_size_spinner, Math.round( vis.getSize() ) );
            _apply_node_size.setSelected( true );
        }
        if ( vis.getNodeColor() != null ) {
            setNodeColor( vis.getNodeColor() );
            _apply_node_color.setSelected( true );
        }
    }

    private void chooseColor( final boolean font ) {
        final Color chosen = JColorChooser.showDialog( this, font ? "Font Color" : "Node Color",
                                                       font ? _font_color : _node_color );
        if ( chosen != null ) {
            if ( font ) {
                setFontColor( chosen );
            }
            else {
                setNodeColor( chosen );
            }
        }
    }

    private void setFontColor( final Color c ) {
        _font_color = c;
        _font_color_button.setBackground( c );
    }

    private void setNodeColor( final Color c ) {
        _node_color = c;
        _node_color_button.setBackground( c );
    }

    // ---- layout helpers ----------------------------------------------------------------------------------------

    private static void addSectionHeader( final JPanel grid, final int row, final String text ) {
        final GridBagConstraints gc = new GridBagConstraints();
        gc.gridx = 0;
        gc.gridy = row;
        gc.gridwidth = 3;
        gc.anchor = GridBagConstraints.WEST;
        gc.insets = new Insets( row == 0 ? 0 : 8, 0, 2, 0 );
        final JLabel l = new JLabel( text );
        l.setFont( l.getFont().deriveFont( Font.BOLD ) );
        grid.add( l, gc );
    }

    private static void addRow( final JPanel grid, final int row, final JCheckBox apply, final String label,
                                final Component control ) {
        final GridBagConstraints gc = new GridBagConstraints();
        gc.gridy = row;
        gc.insets = new Insets( 2, 0, 2, 6 );
        gc.gridx = 0;
        gc.anchor = GridBagConstraints.WEST;
        apply.setToolTipText( "Apply this attribute" );
        grid.add( apply, gc );
        gc.gridx = 1;
        grid.add( new JLabel( label ), gc );
        gc.gridx = 2;
        gc.anchor = GridBagConstraints.WEST;
        grid.add( control, gc );
    }

    private static void wireEnable( final JCheckBox apply, final Component control ) {
        apply.addItemListener( e -> control.setEnabled( apply.isSelected() ) );
    }

    private void syncAllEnabled() {
        _font_style_combo.setEnabled( _apply_font_style.isSelected() );
        _font_size_spinner.setEnabled( _apply_font_size.isSelected() );
        _font_color_button.setEnabled( _apply_font_color.isSelected() );
        _shape_combo.setEnabled( _apply_shape.isSelected() );
        _fill_combo.setEnabled( _apply_fill.isSelected() );
        _node_size_spinner.setEnabled( _apply_node_size.isSelected() );
        _node_color_button.setEnabled( _apply_node_color.isSelected() );
    }

    private static int clamp( final int v, final int lo, final int hi ) {
        return Math.max( lo, Math.min( hi, v ) );
    }

    /** Sets a spinner to {@code value}, first widening its min/max to include it, so a valid stored style value that
     *  lies outside the spinner's usual range is shown exactly (never clamped to a different value). */
    private static void setSpinnerExpanding( final JSpinner spinner, final int value ) {
        final SpinnerNumberModel m = (SpinnerNumberModel) spinner.getModel();
        if ( ( (Number) m.getMaximum() ).intValue() < value ) {
            m.setMaximum( Integer.valueOf( value ) );
        }
        if ( ( (Number) m.getMinimum() ).intValue() > value ) {
            m.setMinimum( Integer.valueOf( value ) );
        }
        spinner.setValue( Integer.valueOf( value ) );
    }

    // ---- label <-> enum mapping --------------------------------------------------------------------------------

    private static int fontStyleFromLabel( final String s ) {
        if ( "Bold".equals( s ) ) {
            return Font.BOLD;
        }
        if ( "Italic".equals( s ) ) {
            return Font.ITALIC;
        }
        if ( "Bold Italic".equals( s ) ) {
            return Font.BOLD + Font.ITALIC;
        }
        return Font.PLAIN;
    }

    private static String labelForFontStyle( final FontType t ) {
        switch ( t ) {
            case BOLD:
                return "Bold";
            case ITALIC:
                return "Italic";
            case BOLD_ITALIC:
                return "Bold Italic";
            default:
                return "Plain";
        }
    }

    private static NodeShape shapeFromLabel( final String s ) {
        if ( "Circle".equals( s ) ) {
            return NodeShape.CIRCLE;
        }
        if ( "Rectangle".equals( s ) ) {
            return NodeShape.RECTANGLE;
        }
        return NodeShape.DEFAULT;
    }

    private static NodeFill fillFromLabel( final String s ) {
        if ( "Solid".equals( s ) ) {
            return NodeFill.SOLID;
        }
        if ( "None".equals( s ) ) {
            return NodeFill.NONE;
        }
        if ( "Gradient".equals( s ) ) {
            return NodeFill.GRADIENT;
        }
        return NodeFill.DEFAULT;
    }

    private static String labelForFill( final NodeFill f ) {
        switch ( f ) {
            case SOLID:
                return "Solid";
            case NONE:
                return "None";
            case GRADIENT:
                return "Gradient";
            default:
                return "Default";
        }
    }

    // ---- test hooks --------------------------------------------------------------------------------------------

    /** Ticks "Apply font colour" and sets the colour (as if the user picked it in the chooser). */
    void setApplyFontColorForTest( final Color c ) {
        setFontColor( c );
        _apply_font_color.setSelected( true );
    }

    /** Ticks "Apply node shape" and selects the shape. */
    void setApplyShapeForTest( final NodeShape shape ) {
        _shape_combo.setSelectedItem( ( shape == NodeShape.RECTANGLE ) ? "Rectangle"
                : ( shape == NodeShape.CIRCLE ) ? "Circle" : "Default" );
        _apply_shape.setSelected( true );
    }

    boolean isApplyFontStyleTickedForTest() {
        return _apply_font_style.isSelected();
    }

    void applyForTest() {
        applyAndClose();
    }
}
