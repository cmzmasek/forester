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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.archaeopteryx;

import java.awt.Color;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;

import javax.swing.JEditorPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTree;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.text.Position;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;

import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.BinaryCharacters;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Distribution;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Point;
import org.forester.phylogeny.data.PropertiesMap;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Reference;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Uri;
import org.forester.util.ForesterUtil;

class NodePanel extends JPanel implements TreeSelectionListener {

    static final String         BASIC                    = "Basic";
    static final String         BINARY_CHARACTERS        = "Binary characters";
    static final String         CONFIDENCE               = "Confidence";
    static final String         CONFIDENCE_TYPE          = "type";
    static final String         DATE                     = "Date";
    static final String         DATE_DESCRIPTION         = "Description";
    static final String         DATE_MAX                 = "Max";
    static final String         DATE_MIN                 = "Min";
    static final String         DATE_UNIT                = "Unit";
    static final String         DATE_VALUE               = "Value";
    static final String         DIST_ALT_UNIT            = "Altitude unit";
    static final String         DIST_ALTITUDE            = "Altitude";
    static final String         DIST_DESCRIPTION         = "Description";
    static final String         DIST_GEODETIC_DATUM      = "Geodetic datum";
    static final String         DIST_LATITUDE            = "Latitude";
    static final String         DIST_LONGITUDE           = "Longitude";
    static final String         DISTRIBUTION             = "Distribution";
    static final String         EVENTS                   = "Events";
    static final String         EVENTS_DUPLICATIONS      = "Duplications";
    static final String         EVENTS_GENE_LOSSES       = "Gene losses";
    static final String         EVENTS_SPECIATIONS       = "Speciations";
    static final String         LIT_REFERENCE            = "Reference";
    static final String         LIT_REFERENCE_DESC       = "Description";
    static final String         LIT_REFERENCE_DOI        = "DOI";
    static final String         NODE_BRANCH_COLOR        = "Branch color";
    static final String         NODE_BRANCH_LENGTH       = "Branch length";
    static final String         NODE_BRANCH_WIDTH        = "Branch width";
    static final String         NODE_NAME                = "Name";
    static final String         PROP                     = "Properties";
    static final String         REFERENCE                = "Reference";
    static final String         SEQ_ACCESSION            = "Accession";
    static final String         SEQ_LOCATION             = "Location";
    static final String         SEQ_MOL_SEQ              = "Mol seq";
    static final String         SEQ_NAME                 = "Name";
    static final String         SEQ_SYMBOL               = "Symbol";
    static final String         SEQ_TYPE                 = "Type";
    static final String         SEQ_URI                  = "URI";
    static final String         SEQUENCE                 = "Sequence";
    static final String         TAXONOMY                 = "Taxonomy";
    static final String         TAXONOMY_AUTHORITY       = "Authority";
    static final String         TAXONOMY_CODE            = "Code";
    static final String         TAXONOMY_COMMON_NAME     = "Common name";
    static final String         TAXONOMY_IDENTIFIER      = "Identifier";
    static final String         TAXONOMY_RANK            = "Rank";
    static final String         TAXONOMY_SCIENTIFIC_NAME = "Scientific name";
    static final String         TAXONOMY_SYNONYM         = "Synonym";
    static final String         TAXONOMY_URI             = "URI";
    private static final String SEQ_GENE_NAME            = "Gene name";
    private static final long   serialVersionUID         = 5120159904388100771L;
    private final JEditorPane   _pane;
    private final JTree         _tree;

    public NodePanel( final PhylogenyNode phylogeny_node ) {
        String node_name = "";
        if ( !ForesterUtil.isEmpty( phylogeny_node.getName() ) ) {
            node_name = phylogeny_node.getName() + " ";
        }
        final DefaultMutableTreeNode top = new DefaultMutableTreeNode( "Node " + node_name );
        createNodes( top, phylogeny_node );
        _tree = new JTree( top );
        _tree.setEditable( false );
        getJTree().setToggleClickCount( 1 );
        expandPath( BASIC );
        expandPath( TAXONOMY );
        expandPath( SEQUENCE );
        expandPath( EVENTS );
        final JScrollPane tree_view = new JScrollPane( getJTree() );
        _pane = new JEditorPane();
        _pane.setEditable( false );
        final JScrollPane data_view = new JScrollPane( _pane );
        final JSplitPane split_pane = new JSplitPane( JSplitPane.VERTICAL_SPLIT );
        split_pane.setTopComponent( tree_view );
        split_pane.setBottomComponent( data_view );
        data_view.setMinimumSize( Constants.NODE_PANEL_SPLIT_MINIMUM_SIZE );
        tree_view.setMinimumSize( Constants.NODE_PANEL_SPLIT_MINIMUM_SIZE );
        split_pane.setDividerLocation( 400 );
        split_pane.setPreferredSize( Constants.NODE_PANEL_SIZE );
        add( split_pane );
    }

    @Override
    public void valueChanged( final TreeSelectionEvent e ) {
        // Do nothing.
    }

    private void expandPath( final String name ) {
        final TreePath tp = getJTree().getNextMatch( name, 0, Position.Bias.Forward );
        if ( tp != null ) {
            getJTree().expandPath( tp );
        }
    }

    private JTree getJTree() {
        return _tree;
    }

    private static void addAnnotation( final DefaultMutableTreeNode top, final Annotation ann, final String name ) {
        DefaultMutableTreeNode category;
        category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelement( category, "Source", ann.getSource() );
        addSubelement( category, "Type", ann.getType() );
        addSubelement( category, "Evidence", ann.getEvidence() );
        if ( ann.getConfidence() != null ) {
            addSubelement( category, CONFIDENCE, ann.getConfidence().asText().toString() );
        }
        if ( ann.getProperties() != null ) {
            addProperties( category, ann.getProperties(), PROP );
        }
    }

    private static void addAnnotations( final DefaultMutableTreeNode top,
                                        final SortedSet<Annotation> annotations,
                                        final DefaultMutableTreeNode category ) {
        if ( ( annotations != null ) && ( annotations.size() > 0 ) ) {
            category.add( new DefaultMutableTreeNode( "Annotations" ) );
            final DefaultMutableTreeNode last = top.getLastLeaf();
            for( final Annotation ann : annotations ) {
                addAnnotation( last, ann, ann.asText().toString() );
            }
        }
    }

    private static void addBasics( final DefaultMutableTreeNode top,
                                   final PhylogenyNode phylogeny_node,
                                   final String name ) {
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelement( category, NODE_NAME, phylogeny_node.getName() );
        if ( phylogeny_node.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
            addSubelement( category,
                           NODE_BRANCH_LENGTH,
                           ForesterUtil.FORMATTER_6.format( phylogeny_node.getDistanceToParent() ) );
        }
        if ( phylogeny_node.getBranchData().isHasConfidences() ) {
            for( final PhylogenyData conf : phylogeny_node.getBranchData().getConfidences() ) {
                addSubelement( category, CONFIDENCE, conf.asText().toString() );
            }
        }
        if ( !phylogeny_node.isExternal() ) {
            addSubelement( category, "Children", String.valueOf( phylogeny_node.getNumberOfDescendants() ) );
            addSubelement( category,
                           "External children",
                           String.valueOf( phylogeny_node.getAllExternalDescendants().size() ) );
            final SortedMap<Taxonomy, Integer> distinct_tax = PhylogenyMethods
                    .obtainDistinctTaxonomyCounts( phylogeny_node );
            if ( distinct_tax != null ) {
                final int no_tax = PhylogenyMethods.calculateNumberOfExternalNodesWithoutTaxonomy( phylogeny_node );
                final int tax_count = distinct_tax.size();
                addSubelement( category, "Distinct external taxonomies", String.valueOf( tax_count ) );
                if ( no_tax > 0 ) {
                    addSubelement( category, "External nodes without taxonomy", String.valueOf( no_tax ) );
                }
            }
        }
        if ( !phylogeny_node.isRoot() ) {
            addSubelement( category, "Depth", String.valueOf( phylogeny_node.calculateDepth() ) );
            final double d = phylogeny_node.calculateDistanceToRoot();
            if ( d > 0 ) {
                addSubelement( category, "Distance to root", String.valueOf( ForesterUtil.FORMATTER_6.format( d ) ) );
            }
        }
        if ( ( phylogeny_node.getBranchData().getBranchWidth() != null )
                && ( phylogeny_node.getBranchData().getBranchWidth().getValue() != BranchWidth.BRANCH_WIDTH_DEFAULT_VALUE ) ) {
            addSubelement( category,
                           NODE_BRANCH_WIDTH,
                           ForesterUtil.FORMATTER_3.format( phylogeny_node.getBranchData().getBranchWidth().getValue() ) );
        }
        if ( ( phylogeny_node.getBranchData().getBranchColor() != null ) ) {
            final Color c = phylogeny_node.getBranchData().getBranchColor().getValue();
            addSubelement( category, NODE_BRANCH_COLOR, c.getRed() + ", " + c.getGreen() + ", " + c.getBlue() );
        }
    }

    private static void addBinaryCharacters( final DefaultMutableTreeNode top,
                                             final BinaryCharacters bc,
                                             final String name ) {
        DefaultMutableTreeNode category;
        category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelement( category, "Gained", String.valueOf( bc.getGainedCount() ) );
        addSubelement( category, "Lost", String.valueOf( bc.getLostCount() ) );
        addSubelement( category, "Present", String.valueOf( bc.getPresentCount() ) );
        final DefaultMutableTreeNode chars = new DefaultMutableTreeNode( "Lists" );
        category.add( chars );
        addSubelement( chars, "Gained", bc.getGainedCharactersAsStringBuffer().toString() );
        addSubelement( chars, "Lost", bc.getLostCharactersAsStringBuffer().toString() );
        addSubelement( chars, "Present", bc.getPresentCharactersAsStringBuffer().toString() );
    }

    private static void addCrossReference( final DefaultMutableTreeNode top, final Accession x, final String name ) {
        DefaultMutableTreeNode category;
        category = new DefaultMutableTreeNode( name );
        top.add( category );
    }

    private static void addCrossReferences( final DefaultMutableTreeNode top,
                                            final SortedSet<Accession> xs,
                                            final DefaultMutableTreeNode category ) {
        if ( ( xs != null ) && ( xs.size() > 0 ) ) {
            category.add( new DefaultMutableTreeNode( "Cross references" ) );
            final DefaultMutableTreeNode last = top.getLastLeaf();
            for( final Accession x : xs ) {
                addCrossReference( last, x, x.asText().toString() );
            }
        }
    }

    private static void addDate( final DefaultMutableTreeNode top, final Date date, final String name ) {
        DefaultMutableTreeNode category;
        category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelement( category, DATE_DESCRIPTION, date.getDesc() );
        addSubelement( category, DATE_VALUE, String.valueOf( date.getValue() ) );
        addSubelement( category, DATE_MIN, String.valueOf( date.getMin() ) );
        addSubelement( category, DATE_MAX, String.valueOf( date.getMax() ) );
        addSubelement( category, DATE_UNIT, date.getUnit() );
    }

    private static void addDistribution( final DefaultMutableTreeNode top, final Distribution dist, final String name ) {
        DefaultMutableTreeNode category;
        category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelement( category, DIST_DESCRIPTION, dist.getDesc() );
        if ( ( dist.getPoints() != null ) && ( dist.getPoints().size() > 0 ) ) {
            final Point p0 = dist.getPoints().get( 0 );
            if ( ( p0 != null ) && !Point.isSeemsEmpty( p0 ) ) {
                addSubelement( category, DIST_GEODETIC_DATUM, p0.getGeodeticDatum() );
                addSubelement( category, DIST_LATITUDE, String.valueOf( p0.getLatitude() ) );
                addSubelement( category, DIST_LONGITUDE, String.valueOf( p0.getLongitude() ) );
                String alt_unit = p0.getAltiudeUnit();
                if ( ForesterUtil.isEmpty( alt_unit ) ) {
                    alt_unit = "?";
                }
                addSubelement( category, DIST_ALTITUDE, String.valueOf( p0.getAltitude() ) + alt_unit );
            }
        }
    }

    private static void addEvents( final DefaultMutableTreeNode top, final Event events, final String name ) {
        DefaultMutableTreeNode category;
        category = new DefaultMutableTreeNode( name );
        top.add( category );
        if ( events.getNumberOfDuplications() > 0 ) {
            addSubelement( category, EVENTS_DUPLICATIONS, String.valueOf( events.getNumberOfDuplications() ) );
        }
        if ( events.getNumberOfSpeciations() > 0 ) {
            addSubelement( category, EVENTS_SPECIATIONS, String.valueOf( events.getNumberOfSpeciations() ) );
        }
        if ( events.getNumberOfGeneLosses() > 0 ) {
            addSubelement( category, EVENTS_GENE_LOSSES, String.valueOf( events.getNumberOfGeneLosses() ) );
        }
        addSubelement( category, "Type", events.getEventType().toString() );
        if ( events.getConfidence() != null ) {
            addSubelement( category, CONFIDENCE, events.getConfidence().asText().toString() );
        }
    }

    private static void addLineage( final DefaultMutableTreeNode top,
                                    final List<String> lineage,
                                    final DefaultMutableTreeNode category ) {
        if ( ( lineage != null ) && ( lineage.size() > 0 ) ) {
            final StringBuilder sb = new StringBuilder();
            for( final String lin : lineage ) {
                if ( !ForesterUtil.isEmpty( lin ) ) {
                    sb.append( lin );
                    sb.append( " > " );
                }
            }
            String str = null;
            if ( sb.length() > 1 ) {
                str = sb.substring( 0, sb.length() - 3 );
            }
            if ( !ForesterUtil.isEmpty( str ) ) {
                addSubelement( category, "Lineage", str );
            }
        }
    }

    private static void addProperties( final DefaultMutableTreeNode top,
                                       final PropertiesMap properties,
                                       final String string ) {
        final SortedMap<String, Property> properties_map = properties.getProperties();
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( "Properties " );
        top.add( category );
        for( final String key : properties_map.keySet() ) {
            final Property prop = properties_map.get( key );
            category.add( new DefaultMutableTreeNode( prop.getRef() + "=" + prop.getValue() + " " + prop.getUnit()
                    + " [" + prop.getAppliesTo().toString() + "]" ) );
        }
    }

    private static void addReference( final DefaultMutableTreeNode top, final Reference ref, final String name ) {
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelement( category, LIT_REFERENCE_DOI, ref.getDoi() );
        addSubelement( category, LIT_REFERENCE_DESC, ref.getDescription() );
    }

    private static void addSequence( final DefaultMutableTreeNode top, final Sequence seq, final String name ) {
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelement( category, SEQ_NAME, seq.getName() );
        addSubelement( category, SEQ_SYMBOL, seq.getSymbol() );
        addSubelement( category, SEQ_GENE_NAME, seq.getGeneName() );
        if ( seq.getAccession() != null ) {
            addSubelement( category, SEQ_ACCESSION, seq.getAccession().asText().toString() );
        }
        addSubelement( category, SEQ_LOCATION, seq.getLocation() );
        addSubelement( category, SEQ_TYPE, seq.getType() );
        addSubelement( category, SEQ_MOL_SEQ, seq.getMolecularSequence() );
        if ( ( seq.getAnnotations() != null ) && !seq.getAnnotations().isEmpty() ) {
            addAnnotations( top, seq.getAnnotations(), category );
        }
        if ( ( seq.getCrossReferences() != null ) && !seq.getCrossReferences().isEmpty() ) {
            addCrossReferences( top, seq.getCrossReferences(), category );
        }
        if ( ( seq.getUris() != null ) && !seq.getUris().isEmpty() ) {
            addUris( top, seq.getUris(), category );
        }
    }

    private static void addSubelement( final DefaultMutableTreeNode node, final String name, final String value ) {
        if ( !ForesterUtil.isEmpty( value ) ) {
            node.add( new DefaultMutableTreeNode( name + ": " + value ) );
        }
    }

    private static void addTaxonomy( final DefaultMutableTreeNode top, final Taxonomy tax, final String name ) {
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        if ( tax.getIdentifier() != null ) {
            addSubelement( category, TAXONOMY_IDENTIFIER, tax.getIdentifier().asText().toString() );
        }
        addSubelement( category, TAXONOMY_CODE, tax.getTaxonomyCode() );
        addSubelement( category, TAXONOMY_SCIENTIFIC_NAME, tax.getScientificName() );
        addSubelement( category, TAXONOMY_AUTHORITY, tax.getAuthority() );
        addSubelement( category, TAXONOMY_COMMON_NAME, tax.getCommonName() );
        for( final String syn : tax.getSynonyms() ) {
            addSubelement( category, TAXONOMY_SYNONYM, syn );
        }
        addSubelement( category, TAXONOMY_RANK, tax.getRank() );
        if ( ( tax.getUris() != null ) && !tax.getUris().isEmpty() ) {
            addUris( top, tax.getUris(), category );
        }
        if ( ( tax.getLineage() != null ) && !tax.getLineage().isEmpty() ) {
            addLineage( top, tax.getLineage(), category );
        }
    }

    private static void addUri( final DefaultMutableTreeNode top, final Uri uri, final String name ) {
        DefaultMutableTreeNode category;
        category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelement( category, "Description", uri.getDescription() );
        addSubelement( category, "Type", uri.getType() );
        addSubelement( category, "URI", uri.getValue().toString() );
    }

    private static void addUris( final DefaultMutableTreeNode top,
                                 final List<Uri> uris,
                                 final DefaultMutableTreeNode category ) {
        if ( ( uris != null ) && ( uris.size() > 0 ) ) {
            category.add( new DefaultMutableTreeNode( "URIs" ) );
            final DefaultMutableTreeNode last = top.getLastLeaf();
            int i = 0;
            for( final Uri uri : uris ) {
                if ( uri != null ) {
                    addUri( last, uri, "URI " + ( i++ ) );
                }
            }
        }
    }

    private static void createNodes( final DefaultMutableTreeNode top, final PhylogenyNode phylogeny_node ) {
        addBasics( top, phylogeny_node, BASIC );
        // Taxonomy
        if ( phylogeny_node.getNodeData().isHasTaxonomy() ) {
            addTaxonomy( top, phylogeny_node.getNodeData().getTaxonomy(), TAXONOMY );
        }
        // Sequence
        if ( phylogeny_node.getNodeData().isHasSequence() ) {
            addSequence( top, phylogeny_node.getNodeData().getSequence(), SEQUENCE );
        }
        // Events
        if ( phylogeny_node.getNodeData().isHasEvent() ) {
            addEvents( top, phylogeny_node.getNodeData().getEvent(), EVENTS );
        }
        // Date
        if ( phylogeny_node.getNodeData().isHasDate() ) {
            addDate( top, phylogeny_node.getNodeData().getDate(), DATE );
        }
        // Distribution
        if ( phylogeny_node.getNodeData().isHasDistribution() ) {
            addDistribution( top, phylogeny_node.getNodeData().getDistribution(), DISTRIBUTION );
        }
        // Reference
        if ( phylogeny_node.getNodeData().isHasReference() ) {
            addReference( top, phylogeny_node.getNodeData().getReference(), LIT_REFERENCE );
        }
        // BinaryCharacters
        if ( phylogeny_node.getNodeData().isHasBinaryCharacters() ) {
            addBinaryCharacters( top, phylogeny_node.getNodeData().getBinaryCharacters(), BINARY_CHARACTERS );
        }
        // Properties
        if ( phylogeny_node.getNodeData().isHasProperties() ) {
            addProperties( top, phylogeny_node.getNodeData().getProperties(), PROP );
        }
    }
}
