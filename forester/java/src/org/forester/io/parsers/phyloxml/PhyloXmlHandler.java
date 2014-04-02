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

package org.forester.io.parsers.phyloxml;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.forester.io.parsers.phyloxml.data.BinaryCharactersParser;
import org.forester.io.parsers.phyloxml.data.BranchWidthParser;
import org.forester.io.parsers.phyloxml.data.ColorParser;
import org.forester.io.parsers.phyloxml.data.ConfidenceParser;
import org.forester.io.parsers.phyloxml.data.DateParser;
import org.forester.io.parsers.phyloxml.data.DistributionParser;
import org.forester.io.parsers.phyloxml.data.EventParser;
import org.forester.io.parsers.phyloxml.data.IdentifierParser;
import org.forester.io.parsers.phyloxml.data.PropertyParser;
import org.forester.io.parsers.phyloxml.data.ReferenceParser;
import org.forester.io.parsers.phyloxml.data.SequenceParser;
import org.forester.io.parsers.phyloxml.data.SequenceRelationParser;
import org.forester.io.parsers.phyloxml.data.TaxonomyParser;
import org.forester.io.parsers.util.PhylogenyParserException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.BinaryCharacters;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Distribution;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.PropertiesMap;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Reference;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.SequenceRelation;
import org.forester.phylogeny.data.SequenceRelation.SEQUENCE_RELATION_TYPE;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.FailedConditionCheckException;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

public final class PhyloXmlHandler extends DefaultHandler {

    private static final String                              PHYLOXML               = "phyloxml";
    private String                                           _current_element_name;
    private Phylogeny                                        _current_phylogeny;
    private List<Phylogeny>                                  _phylogenies;
    private XmlElement                                       _current_xml_element;
    private PhylogenyNode                                    _current_node;
    private static Map<Phylogeny, HashMap<String, Sequence>> phylogenySequencesById = new HashMap<Phylogeny, HashMap<String, Sequence>>();

    PhyloXmlHandler() {
        // Constructor.
    }

    private void addNode() {
        final PhylogenyNode new_node = new PhylogenyNode();
        getCurrentNode().addAsChild( new_node );
        setCurrentNode( new_node );
    }

    @Override
    public void characters( final char[] chars, final int start_index, final int end_index ) {
        if ( ( ( getCurrentXmlElement() != null ) && ( getCurrentElementName() != null ) )
                && !getCurrentElementName().equals( PhyloXmlMapping.CLADE )
                && !getCurrentElementName().equals( PhyloXmlMapping.PHYLOGENY ) ) {
            if ( !ForesterUtil.isEmpty( getCurrentXmlElement().getValueAsString() ) ) {
                getCurrentXmlElement().appendValue( new String( chars, start_index, end_index ) );
            }
            else {
                getCurrentXmlElement().setValue( new String( chars, start_index, end_index ) );
            }
        }
    }

    @Override
    public void endElement( final String namespace_uri, final String local_name, final String qualified_name )
            throws SAXException {
        if ( ForesterUtil.isEmpty( namespace_uri ) || namespace_uri.startsWith( ForesterConstants.PHYLO_XML_LOCATION ) ) {
            if ( local_name.equals( PhyloXmlMapping.CLADE ) ) {
                try {
                    mapElementToPhylogenyNode( getCurrentXmlElement(), getCurrentNode() );
                    if ( !getCurrentNode().isRoot() ) {
                        setCurrentNode( getCurrentNode().getParent() );
                    }
                    getCurrentXmlElement().setValue( null );
                    setCurrentXmlElement( getCurrentXmlElement().getParent() );
                }
                catch ( final PhylogenyParserException ex ) {
                    throw new SAXException( ex.getMessage() );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    throw new SAXException( e.getMessage() );
                }
            }
            else if ( local_name.equals( PhyloXmlMapping.SEQUENCE_RELATION ) ) {
                try {
                    if ( getCurrentPhylogeny() != null ) {
                        final SequenceRelation seqRelation = ( SequenceRelation ) SequenceRelationParser
                                .getInstance( getCurrentPhylogeny() ).parse( getCurrentXmlElement() );
                        final Map<String, Sequence> sequencesById = getSequenceMapByIdForPhylogeny( getCurrentPhylogeny() );
                        final Sequence ref0 = sequencesById.get( seqRelation.getRef0().getSourceId() ), ref1 = sequencesById
                                .get( seqRelation.getRef1().getSourceId() );
                        if ( ref0 != null ) {
                            // check for reverse relation
                            boolean fFoundReverse = false;
                            for( final SequenceRelation sr : ref0.getSequenceRelations() ) {
                                if ( sr.getType().equals( seqRelation.getType() )
                                        && ( ( sr.getRef0().isEqual( ref1 ) && sr.getRef1().isEqual( ref0 ) ) || ( sr
                                                .getRef0().isEqual( ref0 ) && sr.getRef1().isEqual( ref1 ) ) ) ) {
                                    // in this case we don't need to re-add it, but we make sure we don't loose the confidence value
                                    fFoundReverse = true;
                                    if ( ( sr.getConfidence() == null ) && ( seqRelation.getConfidence() != null ) ) {
                                        sr.setConfidence( seqRelation.getConfidence() );
                                    }
                                }
                            }
                            if ( !fFoundReverse ) {
                                ref0.addSequenceRelation( seqRelation );
                            }
                        }
                        if ( ref1 != null ) {
                            // check for reverse relation
                            boolean fFoundReverse = false;
                            for( final SequenceRelation sr : ref1.getSequenceRelations() ) {
                                if ( sr.getType().equals( seqRelation.getType() )
                                        && ( ( sr.getRef0().isEqual( ref1 ) && sr.getRef1().isEqual( ref0 ) ) || ( sr
                                                .getRef0().isEqual( ref0 ) && sr.getRef1().isEqual( ref1 ) ) ) ) {
                                    // in this case we don't need to re-add it, but we make sure we don't loose the confidence value
                                    fFoundReverse = true;
                                    if ( ( sr.getConfidence() == null ) && ( seqRelation.getConfidence() != null ) ) {
                                        sr.setConfidence( seqRelation.getConfidence() );
                                    }
                                }
                            }
                            if ( !fFoundReverse ) {
                                ref1.addSequenceRelation( seqRelation );
                            }
                        }
                        // we add the type to the current phylogeny so we can know it needs to be displayed in the combo
                        final Collection<SEQUENCE_RELATION_TYPE> relationTypesForCurrentPhylogeny = getCurrentPhylogeny()
                                .getRelevantSequenceRelationTypes();
                        if ( !relationTypesForCurrentPhylogeny.contains( seqRelation.getType() ) ) {
                            relationTypesForCurrentPhylogeny.add( seqRelation.getType() );
                        }
                    }
                }
                catch ( final PhyloXmlDataFormatException ex ) {
                    throw new SAXException( ex.getMessage() );
                }
            }
            else if ( local_name.equals( PhyloXmlMapping.PHYLOGENY ) ) {
                try {
                    PhyloXmlHandler.mapElementToPhylogeny( getCurrentXmlElement(), getCurrentPhylogeny() );
                }
                catch ( final PhylogenyParserException e ) {
                    throw new SAXException( e.getMessage() );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    throw new SAXException( e.getMessage() );
                }
                finishPhylogeny();
                reset();
            }
            else if ( local_name.equals( PHYLOXML ) ) {
                // Do nothing.
            }
            else if ( ( getCurrentPhylogeny() != null ) && ( getCurrentXmlElement().getParent() != null ) ) {
                setCurrentXmlElement( getCurrentXmlElement().getParent() );
            }
            setCurrentElementName( null );
        }
    }

    private void finishPhylogeny() throws SAXException {
        getCurrentPhylogeny().recalculateNumberOfExternalDescendants( false );
        getPhylogenies().add( getCurrentPhylogeny() );
        final HashMap<String, Sequence> phyloSequences = phylogenySequencesById.get( getCurrentPhylogeny() );
        if ( phyloSequences != null ) {
            getCurrentPhylogeny().setSequenceRelationQueries( phyloSequences.values() );
            phylogenySequencesById.remove( getCurrentPhylogeny() );
        }
    }

    private String getCurrentElementName() {
        return _current_element_name;
    }

    private PhylogenyNode getCurrentNode() {
        return _current_node;
    }

    private Phylogeny getCurrentPhylogeny() {
        return _current_phylogeny;
    }

    private XmlElement getCurrentXmlElement() {
        return _current_xml_element;
    }

    List<Phylogeny> getPhylogenies() {
        return _phylogenies;
    }

    private void init() {
        reset();
        setPhylogenies( new ArrayList<Phylogeny>() );
    }

    private void initCurrentNode() {
        if ( getCurrentNode() != null ) {
            throw new FailedConditionCheckException( "attempt to create new current node when current node already exists" );
        }
        if ( getCurrentPhylogeny() == null ) {
            throw new FailedConditionCheckException( "attempt to create new current node for non-existing phylogeny" );
        }
        final PhylogenyNode node = new PhylogenyNode();
        getCurrentPhylogeny().setRoot( node );
        setCurrentNode( getCurrentPhylogeny().getRoot() );
    }

    private void mapElementToPhylogenyNode( final XmlElement xml_element, final PhylogenyNode node )
            throws PhylogenyParserException, PhyloXmlDataFormatException {
        if ( xml_element.isHasAttribute( PhyloXmlMapping.BRANCH_LENGTH ) ) {
            double d = 0;
            try {
                d = Double.parseDouble( xml_element.getAttribute( PhyloXmlMapping.BRANCH_LENGTH ) );
            }
            catch ( final NumberFormatException e ) {
                throw new PhylogenyParserException( "ill formatted distance in clade attribute ["
                        + xml_element.getAttribute( PhyloXmlMapping.BRANCH_LENGTH ) + "]: " + e.getMessage() );
            }
            node.setDistanceToParent( d );
        }
        if ( xml_element.isHasAttribute( PhyloXmlMapping.NODE_COLLAPSE ) ) {
            final String collapse_str = xml_element.getAttribute( PhyloXmlMapping.NODE_COLLAPSE );
            if ( !ForesterUtil.isEmpty( collapse_str ) && collapse_str.trim().equalsIgnoreCase( "true" ) ) {
                node.setCollapse( true );
            }
        }
        for( int i = 0; i < xml_element.getNumberOfChildElements(); ++i ) {
            final XmlElement element = xml_element.getChildElement( i );
            final String qualified_name = element.getQualifiedName();
            if ( qualified_name.equals( PhyloXmlMapping.BRANCH_LENGTH ) ) {
                if ( node.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
                    throw new PhylogenyParserException( "ill advised attempt to set distance twice for the same clade (probably via element and via attribute)" );
                }
                node.setDistanceToParent( element.getValueAsDouble() );
            }
            if ( qualified_name.equals( PhyloXmlMapping.NODE_NAME ) ) {
                node.setName( element.getValueAsString() );
            }
            //  else if ( qualified_name.equals( PhyloXmlMapping.NODE_IDENTIFIER ) ) {
            //      node.getNodeData().setNodeIdentifier( ( Identifier ) IdentifierParser.getInstance().parse( element ) );
            //  }
            else if ( qualified_name.equals( PhyloXmlMapping.TAXONOMY ) ) {
                node.getNodeData().addTaxonomy( ( Taxonomy ) TaxonomyParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.SEQUENCE ) ) {
                final Sequence sequence = ( Sequence ) SequenceParser.getInstance().parse( element );
                node.getNodeData().addSequence( sequence );
                // we temporarily store all sequences that have a source ID so we can access them easily when we need to attach relations to them
                final String sourceId = sequence.getSourceId();
                if ( ( getCurrentPhylogeny() != null ) && !ForesterUtil.isEmpty( sourceId ) ) {
                    getSequenceMapByIdForPhylogeny( getCurrentPhylogeny() ).put( sourceId, sequence );
                }
            }
            else if ( qualified_name.equals( PhyloXmlMapping.DISTRIBUTION ) ) {
                node.getNodeData().addDistribution( ( Distribution ) DistributionParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.CLADE_DATE ) ) {
                node.getNodeData().setDate( ( Date ) DateParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.REFERENCE ) ) {
                node.getNodeData().addReference( ( Reference ) ReferenceParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.BINARY_CHARACTERS ) ) {
                node.getNodeData().setBinaryCharacters( ( BinaryCharacters ) BinaryCharactersParser.getInstance()
                        .parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.COLOR ) ) {
                node.getBranchData().setBranchColor( ( BranchColor ) ColorParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.CONFIDENCE ) ) {
                node.getBranchData().addConfidence( ( Confidence ) ConfidenceParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.WIDTH ) ) {
                node.getBranchData().setBranchWidth( ( BranchWidth ) BranchWidthParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.EVENTS ) ) {
                node.getNodeData().setEvent( ( Event ) EventParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.PROPERTY ) ) {
                final Property prop = ( Property ) PropertyParser.getInstance().parse( element );
                if ( prop.getRef().startsWith( NodeVisualData.APTX_VISUALIZATION_REF ) ) {
                    if ( prop.getAppliesTo() == AppliesTo.NODE ) {
                        if ( node.getNodeData().getNodeVisualData() == null ) {
                            node.getNodeData().setNodeVisualData( new NodeVisualData() );
                        }
                        final NodeVisualData vd = node.getNodeData().getNodeVisualData();
                        vd.parseProperty( prop );
                    }
                    else {
                        System.err.println( "Do not know how to handle " + NodeVisualData.APTX_VISUALIZATION_REF
                                + " property applied to " + prop.getAppliesTo() );
                    }
                }
                else {
                    if ( !node.getNodeData().isHasProperties() ) {
                        node.getNodeData().setProperties( new PropertiesMap() );
                    }
                    node.getNodeData().getProperties().addProperty( prop );
                }
            }
        }
    }

    private void newClade() {
        if ( getCurrentNode() == null ) {
            initCurrentNode();
        }
        else {
            addNode();
        }
    }

    private void newPhylogeny() {
        setCurrentPhylogeny( new Phylogeny() );
    }

    private void reset() {
        setCurrentPhylogeny( null );
        setCurrentNode( null );
        setCurrentElementName( null );
        setCurrentXmlElement( null );
    }

    private void setCurrentElementName( final String element_name ) {
        _current_element_name = element_name;
    }

    private void setCurrentNode( final PhylogenyNode current_node ) {
        _current_node = current_node;
    }

    private void setCurrentPhylogeny( final Phylogeny phylogeny ) {
        _current_phylogeny = phylogeny;
    }

    private void setCurrentXmlElement( final XmlElement element ) {
        _current_xml_element = element;
    }

    private void setPhylogenies( final List<Phylogeny> phylogenies ) {
        _phylogenies = phylogenies;
    }

    @Override
    public void startDocument() throws SAXException {
        init();
    }

    @Override
    public void startElement( final String namespace_uri,
                              final String local_name,
                              final String qualified_name,
                              final Attributes attributes ) throws SAXException {
        if ( ForesterUtil.isEmpty( namespace_uri ) || namespace_uri.startsWith( ForesterConstants.PHYLO_XML_LOCATION ) ) {
            setCurrentElementName( local_name );
            if ( local_name.equals( PhyloXmlMapping.CLADE ) ) {
                final XmlElement element = new XmlElement( namespace_uri, local_name, local_name, attributes );
                getCurrentXmlElement().addChildElement( element );
                setCurrentXmlElement( element );
                newClade();
            }
            else if ( local_name.equals( PhyloXmlMapping.PHYLOGENY ) ) {
                setCurrentXmlElement( new XmlElement( "", "", "", null ) );
                newPhylogeny();
                final XmlElement element = new XmlElement( namespace_uri, local_name, local_name, attributes );
                if ( element.isHasAttribute( PhyloXmlMapping.PHYLOGENY_IS_REROOTABLE_ATTR ) ) {
                    getCurrentPhylogeny().setRerootable( Boolean.parseBoolean( element
                            .getAttribute( PhyloXmlMapping.PHYLOGENY_IS_REROOTABLE_ATTR ) ) );
                }
                if ( element.isHasAttribute( PhyloXmlMapping.PHYLOGENY_BRANCHLENGTH_UNIT_ATTR ) ) {
                    getCurrentPhylogeny()
                            .setDistanceUnit( element.getAttribute( PhyloXmlMapping.PHYLOGENY_BRANCHLENGTH_UNIT_ATTR ) );
                }
                if ( element.isHasAttribute( PhyloXmlMapping.PHYLOGENY_IS_ROOTED_ATTR ) ) {
                    getCurrentPhylogeny().setRooted( Boolean.parseBoolean( element
                            .getAttribute( PhyloXmlMapping.PHYLOGENY_IS_ROOTED_ATTR ) ) );
                }
                if ( element.isHasAttribute( PhyloXmlMapping.PHYLOGENY_TYPE_ATTR ) ) {
                    getCurrentPhylogeny().setType( ( element.getAttribute( PhyloXmlMapping.PHYLOGENY_TYPE_ATTR ) ) );
                }
            }
            else if ( local_name.equals( PHYLOXML ) ) {
            }
            else if ( getCurrentPhylogeny() != null ) {
                final XmlElement element = new XmlElement( namespace_uri, local_name, local_name, attributes );
                getCurrentXmlElement().addChildElement( element );
                setCurrentXmlElement( element );
            }
        }
    }

    public static boolean attributeEqualsValue( final XmlElement element,
                                                final String attributeName,
                                                final String attributeValue ) {
        final String attr = element.getAttribute( attributeName );
        return ( ( attr != null ) && attr.equals( attributeValue ) );
    }

    public static String getAtttributeValue( final XmlElement element, final String attributeName ) {
        final String attr = element.getAttribute( attributeName );
        if ( attr != null ) {
            return attr;
        }
        else {
            return "";
        }
    }

    static public Map<String, Sequence> getSequenceMapByIdForPhylogeny( final Phylogeny ph ) {
        HashMap<String, Sequence> seqMap = phylogenySequencesById.get( ph );
        if ( seqMap == null ) {
            seqMap = new HashMap<String, Sequence>();
            phylogenySequencesById.put( ph, seqMap );
        }
        return seqMap;
    }

    private static void mapElementToPhylogeny( final XmlElement xml_element, final Phylogeny phylogeny )
            throws PhylogenyParserException, PhyloXmlDataFormatException {
        for( int i = 0; i < xml_element.getNumberOfChildElements(); ++i ) {
            final XmlElement element = xml_element.getChildElement( i );
            final String qualified_name = element.getQualifiedName();
            if ( qualified_name.equals( PhyloXmlMapping.PHYLOGENY_NAME ) ) {
                phylogeny.setName( element.getValueAsString() );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.PHYLOGENY_DESCRIPTION ) ) {
                phylogeny.setDescription( element.getValueAsString() );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.IDENTIFIER ) ) {
                phylogeny.setIdentifier( ( Identifier ) IdentifierParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.CONFIDENCE ) ) {
                phylogeny.setConfidence( ( Confidence ) ConfidenceParser.getInstance().parse( element ) );
            }
        }
    }
}
