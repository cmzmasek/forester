// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: www.phylosoft.org/

package org.forester.archaeopteryx;

import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.math.BigDecimal;
import java.net.URL;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JEditorPane;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTree;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.text.Position;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import org.forester.archaeopteryx.tools.ImageLoader;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Distribution;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.MultipleUris;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Point;
import org.forester.phylogeny.data.Reference;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Uri;
import org.forester.util.FailedConditionCheckException;
import org.forester.util.ForesterUtil;

class NodeEditPanel extends JPanel {

    private static final long                            serialVersionUID = 5120159904388100771L;
    private final JTree                                  _tree;
    private final JEditorPane                            _pane;
    private final PhylogenyNode                          _my_node;
    private final TreePanel                              _tree_panel;
    private final Map<DefaultMutableTreeNode, TagNumber> _map;

    public NodeEditPanel( final PhylogenyNode phylogeny_node, final TreePanel tree_panel ) {
        _map = new HashMap<DefaultMutableTreeNode, TagNumber>();
        _my_node = phylogeny_node;
        _tree_panel = tree_panel;
        String node_name = "";
        if ( !ForesterUtil.isEmpty( phylogeny_node.getName() ) ) {
            node_name = phylogeny_node.getName() + " ";
        }
        final DefaultMutableTreeNode top = new DefaultMutableTreeNode( "Node " + node_name );
        createNodes( top, phylogeny_node );
        _tree = new JTree( top );
        getJTree().setEditable( true );
        getJTree().setFocusable( true );
        getJTree().setToggleClickCount( 1 );
        getJTree().setInvokesStopCellEditing( true );
        final JScrollPane tree_view = new JScrollPane( getJTree() );
        _pane = new JEditorPane();
        _pane.setEditable( true );
        final JScrollPane data_view = new JScrollPane( _pane );
        final JSplitPane split_pane = new JSplitPane( JSplitPane.VERTICAL_SPLIT );
        split_pane.setTopComponent( tree_view );
        // split_pane.setBottomComponent( data_view );
        data_view.setMinimumSize( Constants.NODE_PANEL_SPLIT_MINIMUM_SIZE );
        tree_view.setMinimumSize( Constants.NODE_PANEL_SPLIT_MINIMUM_SIZE );
        // split_pane.setDividerLocation( 400 );
        split_pane.setPreferredSize( Constants.NODE_PANEL_SIZE );
        add( split_pane );
        getJTree().getSelectionModel().setSelectionMode( TreeSelectionModel.SINGLE_TREE_SELECTION );
        getJTree().addKeyListener( new KeyListener() {

            @Override
            public void keyPressed( final KeyEvent e ) {
                keyEvent( e );
            }

            @Override
            public void keyReleased( final KeyEvent e ) {
                keyEvent( e );
            }

            @Override
            public void keyTyped( final KeyEvent e ) {
                keyEvent( e );
            }
        } );
        for( int i = 0; i < getJTree().getRowCount(); i++ ) {
            getJTree().expandRow( i );
        }
        collapsePath( NodePanel.BASIC );
        collapsePath( NodePanel.TAXONOMY );
        collapsePath( NodePanel.SEQUENCE );
        collapsePath( NodePanel.EVENTS );
        collapsePath( NodePanel.DATE );
        collapsePath( NodePanel.DISTRIBUTION );
        collapsePath( NodePanel.LIT_REFERENCE );
        getJTree().addTreeSelectionListener( new TreeSelectionListener() {

            @Override
            public void valueChanged( final TreeSelectionEvent e ) {
                final TreePath new_path = e.getNewLeadSelectionPath();
                final TreePath old_path = e.getOldLeadSelectionPath();
                if ( new_path != null ) {
                    writeBack( ( DefaultMutableTreeNode ) new_path.getLastPathComponent() );
                }
                if ( old_path != null ) {
                    writeBack( ( DefaultMutableTreeNode ) old_path.getLastPathComponent() );
                }
            }
        } );
    }

    private void addBasics( final DefaultMutableTreeNode top, final PhylogenyNode phylogeny_node, final String name ) {
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelementEditable( category, NodePanel.NODE_NAME, phylogeny_node.getName(), PHYLOXML_TAG.NODE_NAME );
        String bl = "";
        if ( phylogeny_node.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
            bl = ForesterUtil.FORMATTER_6.format( phylogeny_node.getDistanceToParent() );
        }
        addSubelementEditable( category, NodePanel.NODE_BRANCH_LENGTH, bl, PHYLOXML_TAG.NODE_BRANCH_LENGTH );
        int counter = 0;
        if ( phylogeny_node.getBranchData().isHasConfidences() ) {
            for( int i = phylogeny_node.getBranchData().getConfidences().size() - 1; i >= 0; i-- ) {
                if ( phylogeny_node.getBranchData().getConfidences().get( i ).getValue() == Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                    phylogeny_node.getBranchData().getConfidences().remove( i );
                }
            }
            for( final PhylogenyData conf : phylogeny_node.getBranchData().getConfidences() ) {
                final Confidence my_conf = ( Confidence ) ( conf );
                addSubelementEditable( category,
                                       NodePanel.CONFIDENCE + " [" + counter + "]",
                                       ForesterUtil.FORMATTER_6.format( my_conf.getValue() ),
                                       PHYLOXML_TAG.CONFIDENCE_VALUE,
                                       NodePanel.CONFIDENCE_TYPE,
                                       my_conf.getType(),
                                       PHYLOXML_TAG.CONFIDENCE_TYPE,
                                       counter++ );
            }
        }
        addSubelementEditable( category,
                               NodePanel.CONFIDENCE + " [" + counter + "]",
                               "",
                               PHYLOXML_TAG.CONFIDENCE_VALUE,
                               NodePanel.CONFIDENCE_TYPE,
                               "",
                               PHYLOXML_TAG.CONFIDENCE_TYPE,
                               counter );
        String bw = "1";
        if ( ( phylogeny_node.getBranchData().getBranchWidth() != null )
                && ( phylogeny_node.getBranchData().getBranchWidth().getValue() != BranchWidth.BRANCH_WIDTH_DEFAULT_VALUE ) ) {
            bw = ForesterUtil.FORMATTER_3.format( phylogeny_node.getBranchData().getBranchWidth().getValue() );
        }
        addSubelementEditable( category, NodePanel.NODE_BRANCH_WIDTH, bw, PHYLOXML_TAG.NODE_BRANCH_WIDTH );
    }

    //    private void addAnnotation( final DefaultMutableTreeNode top, final Annotation ann, final String name ) {
    //        DefaultMutableTreeNode category;
    //        category = new DefaultMutableTreeNode( name );
    //        top.add( category );
    //        addSubelementEditable( category, "Reference", ann.getRef() , PHYLOXML_TAG.);
    //        addSubelementEditable( category, "Description", ann.getDesc() , PHYLOXML_TAG.);
    //        addSubelementEditable( category, "Source", ann.getSource(), PHYLOXML_TAG. );
    //        addSubelementEditable( category, "Type", ann.getType(), PHYLOXML_TAG. );
    //        addSubelementEditable( category, "Evidence", ann.getEvidence() , PHYLOXML_TAG.);
    //        if ( ann.getConfidence() != null ) {
    //            addSubelementEditable( category, "Confidence", ann.getConfidence().asText().toString() , PHYLOXML_TAG.);
    //        }
    //        if ( ann.getProperties() != null ) {
    //            addProperties( category, ann.getProperties(), "Properties", PHYLOXML_TAG. );
    //        }
    //    }
    //    private void addAnnotations( final DefaultMutableTreeNode top,
    //                                 final List<PhylogenyData> annotations,
    //                                 final DefaultMutableTreeNode category ) {
    //        if ( ( annotations != null ) && ( annotations.size() > 0 ) ) {
    //            category.add( new DefaultMutableTreeNode( "Annotations" ) );
    //            final DefaultMutableTreeNode last = top.getLastLeaf();
    //            int i = 0;
    //            for( final PhylogenyData ann : annotations ) {
    //                addAnnotation( last, ( Annotation ) ann, "Annotation " + ( i++ ) );
    //            }
    //        }
    //    }
    private void addDate( final DefaultMutableTreeNode top, Date date, final String name ) {
        if ( date == null ) {
            date = new Date();
        }
        DefaultMutableTreeNode category;
        category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelementEditable( category, NodePanel.DATE_DESCRIPTION, date.getDesc(), PHYLOXML_TAG.DATE_DESCRIPTION );
        addSubelementEditable( category,
                               NodePanel.DATE_VALUE,
                               String.valueOf( date.getValue() != null ? date.getValue() : "" ),
                               PHYLOXML_TAG.DATE_VALUE );
        addSubelementEditable( category,
                               NodePanel.DATE_MIN,
                               String.valueOf( date.getMin() != null ? date.getMin() : "" ),
                               PHYLOXML_TAG.DATE_MIN );
        addSubelementEditable( category,
                               NodePanel.DATE_MAX,
                               String.valueOf( date.getMax() != null ? date.getMax() : "" ),
                               PHYLOXML_TAG.DATE_MAX );
        addSubelementEditable( category, NodePanel.DATE_UNIT, date.getUnit(), PHYLOXML_TAG.DATE_UNIT );
    }

    private void addDistribution( final DefaultMutableTreeNode top, Distribution dist, final String name ) {
        if ( dist == null ) {
            dist = new Distribution( "" );
        }
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        Point p0 = null;
        if ( ( dist.getPoints() != null ) && ( dist.getPoints().size() > 0 ) ) {
            p0 = dist.getPoints().get( 0 );
        }
        else {
            p0 = new Point();
        }
        addSubelementEditable( category, NodePanel.DIST_DESCRIPTION, dist.getDesc(), PHYLOXML_TAG.DIST_DESC );
        addSubelementEditable( category,
                               NodePanel.DIST_GEODETIC_DATUM,
                               p0.getGeodeticDatum(),
                               PHYLOXML_TAG.DIST_GEODETIC );
        addSubelementEditable( category,
                               NodePanel.DIST_LATITUDE,
                               String.valueOf( p0.getLatitude() != null ? p0.getLatitude() : "" ),
                               PHYLOXML_TAG.DIST_LAT );
        addSubelementEditable( category,
                               NodePanel.DIST_LONGITUDE,
                               String.valueOf( p0.getLongitude() != null ? p0.getLongitude() : "" ),
                               PHYLOXML_TAG.DIST_LONG );
        addSubelementEditable( category,
                               NodePanel.DIST_ALTITUDE,
                               String.valueOf( p0.getAltitude() != null ? p0.getAltitude() : "" ),
                               PHYLOXML_TAG.DIST_ALT );
        addSubelementEditable( category,
                               NodePanel.DIST_ALT_UNIT,
                               String.valueOf( p0.getAltiudeUnit() != null ? p0.getAltiudeUnit() : "" ),
                               PHYLOXML_TAG.DIST_ALT_UNIT );
    }

    private void addEvents( final DefaultMutableTreeNode top, Event events, final String name ) {
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        if ( events == null ) {
            events = new Event();
        }
        top.add( category );
        addSubelementEditable( category,
                               NodePanel.EVENTS_DUPLICATIONS,
                               String.valueOf( events.getNumberOfDuplications() >= 0 ? events.getNumberOfDuplications()
                                       : 0 ),
                               PHYLOXML_TAG.EVENTS_DUPLICATIONS );
        addSubelementEditable( category,
                               NodePanel.EVENTS_SPECIATIONS,
                               String.valueOf( events.getNumberOfSpeciations() >= 0 ? events.getNumberOfSpeciations()
                                       : 0 ),
                               PHYLOXML_TAG.EVENTS_SPECIATIONS );
        addSubelementEditable( category,
                               NodePanel.EVENTS_GENE_LOSSES,
                               String.valueOf( events.getNumberOfGeneLosses() >= 0 ? events.getNumberOfGeneLosses() : 0 ),
                               PHYLOXML_TAG.EVENTS_GENE_LOSSES );
    }

    private void addMapping( final DefaultMutableTreeNode mtn, final TagNumber tag ) {
        if ( getMap().containsKey( mtn ) ) {
            throw new IllegalArgumentException( "key " + mtn + " already present" );
        }
        if ( getMap().containsValue( tag ) ) {
            throw new IllegalArgumentException( "value " + tag + " already present" );
        }
        getMap().put( mtn, tag );
    }

    private void addReference( final DefaultMutableTreeNode top, Reference ref, final String name ) {
        if ( ref == null ) {
            ref = new Reference( "" );
        }
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelementEditable( category,
                               NodePanel.LIT_REFERENCE_DESC,
                               ref.getDescription(),
                               PHYLOXML_TAG.LIT_REFERENCE_DESC );
        addSubelementEditable( category, NodePanel.LIT_REFERENCE_DOI, ref.getDoi(), PHYLOXML_TAG.LIT_REFERENCE_DOI );
    }

    private void addSequence( final DefaultMutableTreeNode top, Sequence seq, final String name ) {
        if ( seq == null ) {
            seq = new Sequence();
        }
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        Accession acc = seq.getAccession();
        if ( acc == null ) {
            acc = new Accession( "", "" );
        }
        addSubelementEditable( category, NodePanel.SEQ_NAME, seq.getName(), PHYLOXML_TAG.SEQ_NAME );
        addSubelementEditable( category, NodePanel.SEQ_SYMBOL, seq.getSymbol(), PHYLOXML_TAG.SEQ_SYMBOL );
        addSubelementEditable( category,
                               NodePanel.SEQ_ACCESSION,
                               acc.getValue(),
                               PHYLOXML_TAG.SEQ_ACC_VALUE,
                               "Source",
                               acc.getSource(),
                               PHYLOXML_TAG.SEQ_ACC_SOURCE );
        addSubelementEditable( category, NodePanel.SEQ_LOCATION, seq.getLocation(), PHYLOXML_TAG.SEQ_LOCATION );
        addSubelementEditable( category, NodePanel.SEQ_TYPE, seq.getType(), PHYLOXML_TAG.SEQ_TYPE );
        addSubelementEditable( category, NodePanel.SEQ_MOL_SEQ, seq.getMolecularSequence(), PHYLOXML_TAG.SEQ_MOL_SEQ );
        int uri_counter = 0;
        if ( seq.getUris() != null ) {
            for( final Uri uri : seq.getUris() ) {
                if ( uri != null ) {
                    addSubelementEditable( category, NodePanel.SEQ_URI + " [" + uri_counter + "]", uri.getValue()
                            .toString(), PHYLOXML_TAG.SEQ_URI, uri_counter++ );
                }
            }
        }
        addSubelementEditable( category,
                               NodePanel.SEQ_URI + " [" + uri_counter + "]",
                               "",
                               PHYLOXML_TAG.SEQ_URI,
                               uri_counter );
        //  addAnnotations( top, seq.getAnnotations(), category );
    }

    private void addSubelementEditable( final DefaultMutableTreeNode node,
                                        final String name,
                                        final String value,
                                        final PHYLOXML_TAG phyloxml_tag ) {
        addSubelementEditable( node, name, value, phyloxml_tag, 0 );
    }

    private void addSubelementEditable( final DefaultMutableTreeNode node,
                                        final String name,
                                        final String value,
                                        final PHYLOXML_TAG phyloxml_tag,
                                        final int number ) {
        String my_value = value;
        if ( ForesterUtil.isEmpty( my_value ) ) {
            my_value = "";
        }
        final DefaultMutableTreeNode name_node = new DefaultMutableTreeNode( name );
        final DefaultMutableTreeNode value_node = new DefaultMutableTreeNode( my_value );
        name_node.add( value_node );
        node.add( name_node );
        addMapping( name_node, new TagNumber( phyloxml_tag, number ) );
    }

    private void addSubelementEditable( final DefaultMutableTreeNode node,
                                        final String name,
                                        final String value,
                                        final PHYLOXML_TAG phyloxml_value_tag,
                                        final String source_name,
                                        final String source_value,
                                        final PHYLOXML_TAG phyloxml_source_tag ) {
        addSubelementEditable( node, name, value, phyloxml_value_tag, source_name, source_value, phyloxml_source_tag, 0 );
    }

    private void addSubelementEditable( final DefaultMutableTreeNode node,
                                        final String name,
                                        final String value,
                                        final PHYLOXML_TAG phyloxml_value_tag,
                                        final String source_name,
                                        final String source_value,
                                        final PHYLOXML_TAG phyloxml_source_tag,
                                        final int number ) {
        String my_value = value;
        if ( ForesterUtil.isEmpty( my_value ) ) {
            my_value = "";
        }
        String my_source_value = source_value;
        if ( ForesterUtil.isEmpty( my_source_value ) ) {
            my_source_value = "";
        }
        final DefaultMutableTreeNode name_node = new DefaultMutableTreeNode( name );
        final DefaultMutableTreeNode source_name_node = new DefaultMutableTreeNode( source_name );
        final DefaultMutableTreeNode source_value_node = new DefaultMutableTreeNode( my_source_value );
        final DefaultMutableTreeNode value_node = new DefaultMutableTreeNode( my_value );
        name_node.add( source_name_node );
        source_name_node.add( source_value_node );
        name_node.add( value_node );
        node.add( name_node );
        addMapping( name_node, new TagNumber( phyloxml_value_tag, number ) );
        addMapping( source_name_node, new TagNumber( phyloxml_source_tag, number ) );
    }

    private void addTaxonomy( final DefaultMutableTreeNode top, Taxonomy tax, final String name ) {
        if ( tax == null ) {
            tax = new Taxonomy();
        }
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        Identifier id = tax.getIdentifier();
        if ( id == null ) {
            id = new Identifier();
        }
        addSubelementEditable( category,
                               NodePanel.TAXONOMY_IDENTIFIER,
                               id.getValue(),
                               PHYLOXML_TAG.TAXONOMY_ID_VALUE,
                               "Provider",
                               id.getProvider(),
                               PHYLOXML_TAG.TAXONOMY_ID_PROVIDER );
        addSubelementEditable( category, NodePanel.TAXONOMY_CODE, tax.getTaxonomyCode(), PHYLOXML_TAG.TAXONOMY_CODE );
        addSubelementEditable( category,
                               NodePanel.TAXONOMY_SCIENTIFIC_NAME,
                               tax.getScientificName(),
                               PHYLOXML_TAG.TAXONOMY_SCIENTIFIC_NAME );
        addSubelementEditable( category,
                               NodePanel.TAXONOMY_AUTHORITY,
                               tax.getAuthority(),
                               PHYLOXML_TAG.TAXONOMY_AUTHORITY );
        addSubelementEditable( category,
                               NodePanel.TAXONOMY_COMMON_NAME,
                               tax.getCommonName(),
                               PHYLOXML_TAG.TAXONOMY_COMMON_NAME );
        for( int i = tax.getSynonyms().size() - 1; i >= 0; i-- ) {
            if ( ForesterUtil.isEmpty( tax.getSynonyms().get( i ) ) ) {
                tax.getSynonyms().remove( i );
            }
        }
        int syn_counter = 0;
        for( final String syn : tax.getSynonyms() ) {
            addSubelementEditable( category,
                                   NodePanel.TAXONOMY_SYNONYM + " [" + syn_counter + "]",
                                   syn,
                                   PHYLOXML_TAG.TAXONOMY_SYNONYM,
                                   syn_counter++ );
        }
        addSubelementEditable( category,
                               NodePanel.TAXONOMY_SYNONYM + " [" + syn_counter + "]",
                               "",
                               PHYLOXML_TAG.TAXONOMY_SYNONYM,
                               syn_counter );
        addSubelementEditable( category, NodePanel.TAXONOMY_RANK, tax.getRank(), PHYLOXML_TAG.TAXONOMY_RANK );
        int uri_counter = 0;
        if ( tax.getUris() != null ) {
            for( final Uri uri : tax.getUris() ) {
                if ( uri != null ) {
                    addSubelementEditable( category, NodePanel.TAXONOMY_URI + " [" + uri_counter + "]", uri.getValue()
                            .toString(), PHYLOXML_TAG.TAXONOMY_URI, uri_counter++ );
                }
            }
        }
        addSubelementEditable( category,
                               NodePanel.TAXONOMY_URI + " [" + uri_counter + "]",
                               "",
                               PHYLOXML_TAG.TAXONOMY_URI,
                               uri_counter );
    }

    private void collapsePath( final String name ) {
        final TreePath tp = getJTree().getNextMatch( name, 0, Position.Bias.Forward );
        if ( tp != null ) {
            getJTree().collapsePath( tp );
        }
    }

    private void createNodes( final DefaultMutableTreeNode top, final PhylogenyNode phylogeny_node ) {
        if ( !phylogeny_node.getNodeData().isHasTaxonomy() ) {
            phylogeny_node.getNodeData().addTaxonomy( new Taxonomy() );
        }
        if ( !phylogeny_node.getNodeData().isHasSequence() ) {
            phylogeny_node.getNodeData().addSequence( new Sequence() );
        }
        if ( !phylogeny_node.getNodeData().isHasDistribution() ) {
            phylogeny_node.getNodeData().addDistribution( new Distribution( "" ) );
        }
        if ( !phylogeny_node.getNodeData().isHasReference() ) {
            phylogeny_node.getNodeData().addReference( new Reference( "" ) );
        }
        addBasics( top, phylogeny_node, NodePanel.BASIC );
        addTaxonomy( top, phylogeny_node.getNodeData().getTaxonomy(), NodePanel.TAXONOMY );
        addSequence( top, phylogeny_node.getNodeData().getSequence(), NodePanel.SEQUENCE );
        if ( !phylogeny_node.isExternal() ) {
            addEvents( top, phylogeny_node.getNodeData().getEvent(), NodePanel.EVENTS );
        }
        addDate( top, phylogeny_node.getNodeData().getDate(), NodePanel.DATE );
        addDistribution( top, phylogeny_node.getNodeData().getDistribution(), NodePanel.DISTRIBUTION );
        addReference( top, phylogeny_node.getNodeData().getReference(), NodePanel.LIT_REFERENCE );
        //  addProperties( top, phylogeny_node.getNodeData().getProperties(), "Properties" );
    }

    private void formatError( final DefaultMutableTreeNode mtn, final PhyloXmlDataFormatException e ) {
        JOptionPane.showMessageDialog( this, e.getMessage(), "Format error", JOptionPane.ERROR_MESSAGE );
        mtn.setUserObject( "" );
        getJTree().repaint();
    }

    private JTree getJTree() {
        return _tree;
    }

    private Map<DefaultMutableTreeNode, TagNumber> getMap() {
        return _map;
    }

    private TagNumber getMapping( final DefaultMutableTreeNode mtn ) {
        return getMap().get( mtn );
    }

    PhylogenyNode getMyNode() {
        return _my_node;
    }

    private DefaultMutableTreeNode getSelectedTreeNode() {
        final TreePath selectionPath = getJTree().getSelectionPath();
        if ( selectionPath != null ) {
            final Object[] path = selectionPath.getPath();
            if ( path.length > 0 ) {
                return ( DefaultMutableTreeNode ) path[ path.length - 1 ]; // Last node
            }
        }
        return null;
    }

    private TreePanel getTreePanel() {
        return _tree_panel;
    }

    private void keyEvent( final KeyEvent e ) {
        if ( e.getKeyCode() == KeyEvent.VK_ENTER ) {
            writeBack( getSelectedTreeNode() );
        }
    }

    private List<Point> obtainPoints() {
        AptxUtil.ensurePresenceOfDistribution( getMyNode() );
        Distribution d = getMyNode().getNodeData().getDistribution();
        if ( d.getPoints() == null ) {
            d = new Distribution( d.getDesc(), new ArrayList<Point>(), d.getPolygons() );
            getMyNode().getNodeData().setDistribution( d );
        }
        final List<Point> ps = d.getPoints();
        if ( ps.isEmpty() ) {
            ps.add( new Point() );
        }
        else if ( ps.get( 0 ) == null ) {
            ps.set( 0, new Point() );
        }
        return ps;
    }

    private BigDecimal parseBigDecimal( final DefaultMutableTreeNode mtn, final String value ) {
        if ( ForesterUtil.isEmpty( value ) ) {
            return new BigDecimal( 0 );
        }
        BigDecimal i = null;
        try {
            i = new BigDecimal( value );
        }
        catch ( final NumberFormatException e ) {
            JOptionPane.showMessageDialog( this, "illegal value: " + value, "Error", JOptionPane.ERROR_MESSAGE );
            mtn.setUserObject( "" );
        }
        return i;
    }

    private int parsePositiveInt( final DefaultMutableTreeNode mtn, final String value ) {
        if ( ForesterUtil.isEmpty( value ) ) {
            return 0;
        }
        int i = -1;
        try {
            i = ForesterUtil.parseInt( value );
        }
        catch ( final ParseException e ) {
            JOptionPane.showMessageDialog( this, "illegal value: " + value, "Error", JOptionPane.ERROR_MESSAGE );
            mtn.setUserObject( "" );
        }
        if ( i < 0 ) {
            JOptionPane.showMessageDialog( this, "illegal value: " + value, "Error", JOptionPane.ERROR_MESSAGE );
            mtn.setUserObject( "" );
        }
        return i;
    }

    void writeAll() {
        for( int i = 0; i < getJTree().getRowCount(); i++ ) {
            final TreePath p = getJTree().getPathForRow( i );
            writeBack( ( DefaultMutableTreeNode ) p.getLastPathComponent() );
        }
    }

    private void writeBack( final DefaultMutableTreeNode mtn ) {
        if ( !getMap().containsKey( mtn ) ) {
            final DefaultMutableTreeNode parent = ( DefaultMutableTreeNode ) mtn.getParent();
            if ( getMap().containsKey( parent ) ) {
                writeBack( mtn, getMapping( parent ) );
            }
        }
    }

    private void writeBack( final DefaultMutableTreeNode mtn, final TagNumber tag_number ) {
        if ( tag_number == null ) {
            return;
        }
        String value = mtn.toString();
        if ( value == null ) {
            value = "";
        }
        value = value.replaceAll( "\\s+", " " );
        value = value.trim();
        mtn.setUserObject( value );
        getJTree().repaint();
        final PHYLOXML_TAG tag = tag_number.getTag();
        final int number = tag_number.getNumber();
        switch ( tag ) {
            case NODE_NAME:
                getMyNode().setName( value );
                break;
            case NODE_BRANCH_LENGTH:
                if ( ForesterUtil.isEmpty( value ) ) {
                    getMyNode().setDistanceToParent( PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT );
                }
                else {
                    try {
                        getMyNode().setDistanceToParent( ForesterUtil.parseDouble( value ) );
                    }
                    catch ( final ParseException e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "failed to parse branch length from: " + value,
                                                       "Error",
                                                       JOptionPane.ERROR_MESSAGE );
                        mtn.setUserObject( "" );
                    }
                }
                break;
            case NODE_BRANCH_WIDTH:
                if ( ForesterUtil.isEmpty( value ) || value.equals( "1" ) ) {
                    if ( getMyNode().getBranchData().getBranchWidth() != null ) {
                        getMyNode().getBranchData().setBranchWidth( new BranchWidth() );
                    }
                }
                else {
                    try {
                        final double bw = ForesterUtil.parseDouble( value );
                        if ( bw >= 0 ) {
                            getMyNode().getBranchData().setBranchWidth( new BranchWidth( bw ) );
                        }
                    }
                    catch ( final ParseException e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "failed to parse branch width from: " + value,
                                                       "Error",
                                                       JOptionPane.ERROR_MESSAGE );
                        mtn.setUserObject( "" );
                    }
                }
                break;
            case CONFIDENCE_VALUE:
                double confidence = Confidence.CONFIDENCE_DEFAULT_VALUE;
                if ( !ForesterUtil.isEmpty( value ) ) {
                    try {
                        confidence = ForesterUtil.parseDouble( value );
                    }
                    catch ( final ParseException e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "failed to parse confidence value from: " + value,
                                                       "Error",
                                                       JOptionPane.ERROR_MESSAGE );
                        mtn.setUserObject( "" );
                        break;
                    }
                }
                if ( getMyNode().getBranchData().getConfidences().size() < number ) {
                    throw new FailedConditionCheckException();
                }
                else if ( getMyNode().getBranchData().getConfidences().size() == number ) {
                    if ( confidence >= 0 ) {
                        getMyNode().getBranchData().getConfidences().add( new Confidence( confidence, "unknown" ) );
                    }
                }
                else {
                    final String type = getMyNode().getBranchData().getConfidences().get( number ).getType();
                    final double sd = getMyNode().getBranchData().getConfidences().get( number ).getStandardDeviation();
                    getMyNode().getBranchData().getConfidences().set( number, new Confidence( confidence, type, sd ) );
                }
                break;
            case CONFIDENCE_TYPE:
                if ( getMyNode().getBranchData().getConfidences().size() < number ) {
                    throw new FailedConditionCheckException();
                }
                else if ( getMyNode().getBranchData().getConfidences().size() == number ) {
                    if ( !ForesterUtil.isEmpty( value ) ) {
                        getMyNode().getBranchData().getConfidences().add( new Confidence( 0, value ) );
                    }
                }
                else {
                    final double v = getMyNode().getBranchData().getConfidences().get( number ).getValue();
                    final double sd = getMyNode().getBranchData().getConfidences().get( number ).getStandardDeviation();
                    getMyNode().getBranchData().getConfidences().set( number, new Confidence( v, value, sd ) );
                }
                break;
            case TAXONOMY_CODE:
                AptxUtil.ensurePresenceOfTaxonomy( getMyNode() );
                try {
                    getMyNode().getNodeData().getTaxonomy().setTaxonomyCode( value );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    formatError( mtn, e );
                    break;
                }
                break;
            case TAXONOMY_SCIENTIFIC_NAME:
                AptxUtil.ensurePresenceOfTaxonomy( getMyNode() );
                getMyNode().getNodeData().getTaxonomy().setScientificName( value );
                break;
            case TAXONOMY_COMMON_NAME:
                AptxUtil.ensurePresenceOfTaxonomy( getMyNode() );
                getMyNode().getNodeData().getTaxonomy().setCommonName( value );
                break;
            case TAXONOMY_RANK:
                AptxUtil.ensurePresenceOfTaxonomy( getMyNode() );
                try {
                    getMyNode().getNodeData().getTaxonomy().setRank( value.toLowerCase() );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    formatError( mtn, e );
                    break;
                }
                break;
            case TAXONOMY_AUTHORITY:
                AptxUtil.ensurePresenceOfTaxonomy( getMyNode() );
                getMyNode().getNodeData().getTaxonomy().setAuthority( value );
                break;
            case TAXONOMY_URI: {
                Uri uri = null;
                if ( !ForesterUtil.isEmpty( value ) ) {
                    try {
                        uri = new Uri( new URL( value ).toURI() );
                    }
                    catch ( final Exception e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "failed to parse URL from: " + value,
                                                       "Error",
                                                       JOptionPane.ERROR_MESSAGE );
                        mtn.setUserObject( "" );
                    }
                }
                if ( uri != null ) {
                    AptxUtil.ensurePresenceOfTaxonomy( getMyNode() );
                }
                addUri( mtn, uri, number, getMyNode().getNodeData().getTaxonomy() );
                break;
            }
            case TAXONOMY_SYNONYM:
                if ( getMyNode().getNodeData().getTaxonomy().getSynonyms().size() < number ) {
                    throw new FailedConditionCheckException();
                }
                else if ( getMyNode().getNodeData().getTaxonomy().getSynonyms().size() == number ) {
                    if ( !ForesterUtil.isEmpty( value ) ) {
                        AptxUtil.ensurePresenceOfTaxonomy( getMyNode() );
                        getMyNode().getNodeData().getTaxonomy().getSynonyms().add( value );
                    }
                }
                else {
                    getMyNode().getNodeData().getTaxonomy().getSynonyms().set( number, value );
                }
                break;
            case TAXONOMY_ID_VALUE:
                AptxUtil.ensurePresenceOfTaxonomy( getMyNode() );
                if ( getMyNode().getNodeData().getTaxonomy().getIdentifier() == null ) {
                    getMyNode().getNodeData().getTaxonomy().setIdentifier( new Identifier( value ) );
                }
                else {
                    final String provider = getMyNode().getNodeData().getTaxonomy().getIdentifier().getProvider();
                    getMyNode().getNodeData().getTaxonomy().setIdentifier( new Identifier( value, provider ) );
                }
                break;
            case TAXONOMY_ID_PROVIDER:
                AptxUtil.ensurePresenceOfTaxonomy( getMyNode() );
                if ( getMyNode().getNodeData().getTaxonomy().getIdentifier() == null ) {
                    getMyNode().getNodeData().getTaxonomy().setIdentifier( new Identifier( "", value ) );
                }
                else {
                    final String v = getMyNode().getNodeData().getTaxonomy().getIdentifier().getValue();
                    getMyNode().getNodeData().getTaxonomy().setIdentifier( new Identifier( v, value ) );
                }
                break;
            case SEQ_LOCATION:
                AptxUtil.ensurePresenceOfSequence( getMyNode() );
                getMyNode().getNodeData().getSequence().setLocation( value );
                break;
            case SEQ_MOL_SEQ:
                AptxUtil.ensurePresenceOfSequence( getMyNode() );
                getMyNode().getNodeData().getSequence().setMolecularSequence( value.replaceAll( "[^a-zA-Z-]", "" ) );
                break;
            case SEQ_NAME:
                AptxUtil.ensurePresenceOfSequence( getMyNode() );
                getMyNode().getNodeData().getSequence().setName( value );
                break;
            case SEQ_SYMBOL:
                AptxUtil.ensurePresenceOfSequence( getMyNode() );
                try {
                    getMyNode().getNodeData().getSequence().setSymbol( value );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    formatError( mtn, e );
                    break;
                }
                break;
            case SEQ_TYPE:
                AptxUtil.ensurePresenceOfSequence( getMyNode() );
                try {
                    getMyNode().getNodeData().getSequence().setType( value.toLowerCase() );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    formatError( mtn, e );
                    break;
                }
                break;
            case SEQ_ACC_SOURCE:
                AptxUtil.ensurePresenceOfSequence( getMyNode() );
                if ( getMyNode().getNodeData().getSequence().getAccession() == null ) {
                    getMyNode().getNodeData().getSequence().setAccession( new Accession( "", value ) );
                }
                else {
                    final String v = getMyNode().getNodeData().getSequence().getAccession().getValue();
                    getMyNode().getNodeData().getSequence().setAccession( new Accession( v, value ) );
                }
                break;
            case SEQ_ACC_VALUE:
                AptxUtil.ensurePresenceOfSequence( getMyNode() );
                if ( getMyNode().getNodeData().getSequence().getAccession() == null ) {
                    getMyNode().getNodeData().getSequence().setAccession( new Accession( value, "" ) );
                }
                else {
                    final String source = getMyNode().getNodeData().getSequence().getAccession().getSource();
                    getMyNode().getNodeData().getSequence().setAccession( new Accession( value, source ) );
                }
                break;
            case SEQ_URI: {
                Uri uri = null;
                if ( !ForesterUtil.isEmpty( value ) ) {
                    try {
                        uri = new Uri( new URL( value ).toURI() );
                    }
                    catch ( final Exception e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "failed to parse URL from: " + value,
                                                       "Error",
                                                       JOptionPane.ERROR_MESSAGE );
                        mtn.setUserObject( "" );
                    }
                }
                if ( uri != null ) {
                    AptxUtil.ensurePresenceOfSequence( getMyNode() );
                }
                addUri( mtn, uri, number, getMyNode().getNodeData().getSequence() );
                break;
            }
            case LIT_REFERENCE_DESC:
                if ( !getMyNode().getNodeData().isHasReference() ) {
                    getMyNode().getNodeData().setReference( new Reference( "" ) );
                }
                getMyNode().getNodeData().getReference().setValue( value );
                break;
            case LIT_REFERENCE_DOI:
                if ( !getMyNode().getNodeData().isHasReference() ) {
                    getMyNode().getNodeData().setReference( new Reference( "" ) );
                }
                try {
                    getMyNode().getNodeData().getReference().setDoi( value );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    formatError( mtn, e );
                    break;
                }
                break;
            case EVENTS_DUPLICATIONS:
                if ( !getMyNode().getNodeData().isHasEvent() ) {
                    getMyNode().getNodeData().setEvent( new Event() );
                }
                getMyNode().getNodeData().getEvent().setDuplications( parsePositiveInt( mtn, value ) );
                break;
            case EVENTS_SPECIATIONS:
                if ( !getMyNode().getNodeData().isHasEvent() ) {
                    getMyNode().getNodeData().setEvent( new Event() );
                }
                getMyNode().getNodeData().getEvent().setSpeciations( parsePositiveInt( mtn, value ) );
                break;
            case EVENTS_GENE_LOSSES:
                if ( !getMyNode().getNodeData().isHasEvent() ) {
                    getMyNode().getNodeData().setEvent( new Event() );
                }
                getMyNode().getNodeData().getEvent().setGeneLosses( parsePositiveInt( mtn, value ) );
                break;
            case DATE_DESCRIPTION:
                AptxUtil.ensurePresenceOfDate( getMyNode() );
                getMyNode().getNodeData().getDate().setDesc( value );
                break;
            case DATE_MAX:
                AptxUtil.ensurePresenceOfDate( getMyNode() );
                getMyNode().getNodeData().getDate().setMax( parseBigDecimal( mtn, value ) );
                break;
            case DATE_MIN:
                AptxUtil.ensurePresenceOfDate( getMyNode() );
                getMyNode().getNodeData().getDate().setMin( parseBigDecimal( mtn, value ) );
                break;
            case DATE_UNIT:
                AptxUtil.ensurePresenceOfDate( getMyNode() );
                getMyNode().getNodeData().getDate().setUnit( value );
                break;
            case DATE_VALUE:
                AptxUtil.ensurePresenceOfDate( getMyNode() );
                getMyNode().getNodeData().getDate().setValue( parseBigDecimal( mtn, value ) );
                break;
            case DIST_ALT: {
                final BigDecimal new_value = parseBigDecimal( mtn, value );
                if ( new_value != null ) {
                    final List<Point> ps = obtainPoints();
                    final Point p = ps.get( 0 );
                    final Point p_new = new Point( p.getGeodeticDatum(),
                                                   p.getLatitude(),
                                                   p.getLongitude(),
                                                   new_value,
                                                   ForesterUtil.isEmpty( p.getAltiudeUnit() ) ? "?"
                                                           : p.getAltiudeUnit() );
                    ps.set( 0, p_new );
                }
                break;
            }
            case DIST_DESC: {
                AptxUtil.ensurePresenceOfDistribution( getMyNode() );
                final Distribution d = getMyNode().getNodeData().getDistribution();
                getMyNode().getNodeData().setDistribution( new Distribution( value, d.getPoints(), d.getPolygons() ) );
                break;
            }
            case DIST_GEODETIC: {
                if ( !ForesterUtil.isEmpty( value ) ) {
                    final List<Point> ps = obtainPoints();
                    final Point p = ps.get( 0 );
                    final Point p_new = new Point( value,
                                                   p.getLatitude(),
                                                   p.getLongitude(),
                                                   p.getAltitude(),
                                                   p.getAltiudeUnit() );
                    ps.set( 0, p_new );
                }
                break;
            }
            case DIST_ALT_UNIT: {
                if ( !ForesterUtil.isEmpty( value ) ) {
                    final List<Point> ps = obtainPoints();
                    final Point p = ps.get( 0 );
                    final Point p_new = new Point( p.getGeodeticDatum(),
                                                   p.getLatitude(),
                                                   p.getLongitude(),
                                                   p.getAltitude(),
                                                   value );
                    ps.set( 0, p_new );
                }
                break;
            }
            case DIST_LAT: {
                final BigDecimal new_value = parseBigDecimal( mtn, value );
                if ( new_value != null ) {
                    final List<Point> ps = obtainPoints();
                    final Point p = ps.get( 0 );
                    final Point p_new = new Point( p.getGeodeticDatum(),
                                                   new_value,
                                                   p.getLongitude(),
                                                   p.getAltitude(),
                                                   p.getAltiudeUnit() );
                    ps.set( 0, p_new );
                }
                break;
            }
            case DIST_LONG: {
                final BigDecimal new_value = parseBigDecimal( mtn, value );
                if ( new_value != null ) {
                    final List<Point> ps = obtainPoints();
                    final Point p = ps.get( 0 );
                    final Point p_new = new Point( p.getGeodeticDatum(),
                                                   p.getLatitude(),
                                                   new_value,
                                                   p.getAltitude(),
                                                   p.getAltiudeUnit() );
                    ps.set( 0, p_new );
                }
                break;
            }
            default:
                throw new IllegalArgumentException( "unknown: " + tag );
        }
        getJTree().repaint();
        getTreePanel().setEdited( true );
        getTreePanel().repaint();
    }

    private void addUri( final DefaultMutableTreeNode mtn, final Uri uri, final int number, final MultipleUris mu ) {
        if ( uri != null ) {
            if ( mu.getUris() == null ) {
                mu.setUris( new ArrayList<Uri>() );
            }
        }
        if ( ( uri != null ) && ( mu.getUris() == null ) ) {
            mu.setUris( new ArrayList<Uri>() );
        }
        if ( ( uri != null ) && ( mu.getUris().size() == number ) ) {
            mu.getUris().add( uri );
        }
        if ( ( mu.getUris() != null ) && ( mu.getUris().size() != number ) ) {
            mu.getUris().set( number, uri );
        }
        final ImageLoader il = new ImageLoader( getTreePanel() );
        new Thread( il ).start();
    }

    private enum PHYLOXML_TAG {
        NODE_NAME,
        NODE_BRANCH_LENGTH,
        NODE_BRANCH_WIDTH,
        TAXONOMY_CODE,
        TAXONOMY_SCIENTIFIC_NAME,
        TAXONOMY_AUTHORITY,
        TAXONOMY_COMMON_NAME,
        TAXONOMY_SYNONYM,
        TAXONOMY_RANK,
        TAXONOMY_URI,
        SEQ_SYMBOL,
        SEQ_NAME,
        SEQ_LOCATION,
        SEQ_TYPE,
        SEQ_MOL_SEQ,
        SEQ_URI,
        DATE_DESCRIPTION,
        DATE_VALUE,
        DATE_MIN,
        DATE_MAX,
        DATE_UNIT,
        TAXONOMY_ID_VALUE,
        TAXONOMY_ID_PROVIDER,
        SEQ_ACC_VALUE,
        SEQ_ACC_SOURCE,
        CONFIDENCE_VALUE,
        CONFIDENCE_TYPE,
        LIT_REFERENCE_DESC,
        LIT_REFERENCE_DOI,
        EVENTS_DUPLICATIONS,
        EVENTS_SPECIATIONS,
        EVENTS_GENE_LOSSES,
        DIST_DESC,
        DIST_GEODETIC,
        DIST_LAT,
        DIST_LONG,
        DIST_ALT,
        DIST_ALT_UNIT
    }

    private class TagNumber {

        final private PHYLOXML_TAG _tag;
        final private int          _number;

        TagNumber( final PHYLOXML_TAG tag, final int number ) {
            _tag = tag;
            _number = number;
        }

        int getNumber() {
            return _number;
        }

        PHYLOXML_TAG getTag() {
            return _tag;
        }

        @Override
        public String toString() {
            return getTag() + "_" + getNumber();
        }
    }
}
