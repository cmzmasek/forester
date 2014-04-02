// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.util.ForesterUtil;

public class NodeData implements PhylogenyData {

    public enum NODE_DATA {
        NODE_NAME,
        EVENT,
        SEQUENCE_NAME,
        GENE_NAME,
        SEQUENCE_SYMBOL,
        SEQUENCE_MOL_SEQ,
        SEQUENCE_MOL_SEQ_FASTA,
        SEQUENCE_ACC,
        TAXONOMY_SCIENTIFIC_NAME,
        TAXONOMY_COMM0N_NAME,
        TAXONOMY_CODE,
        UNKNOWN;
    }
    private String             _node_name;
    private Event              _event;
    private List<Sequence>     _sequences;
    private List<Taxonomy>     _taxonomies;
    private List<Distribution> _distributions;
    private Date               _date;
    private BinaryCharacters   _binary_characters;
    private PropertiesMap      _properties;
    private List<Reference>    _references;
    private List<Double>       _vector;
    private NodeVisualData     _node_visual_data;

    public NodeData() {
        init();
    }

    private void init() {
        _node_name = "";
        _event = null;
        _sequences = null;
        _taxonomies = null;
        _distributions = null;
        _date = null;
        _binary_characters = null;
        _properties = null;
        _references = null;
        _vector = null;
        _node_visual_data = null;
    }

    public void addDistribution( final Distribution distribution ) {
        if ( _distributions == null ) {
            _distributions = new ArrayList<Distribution>();
        }
        _distributions.add( distribution );
    }

    public void addReference( final Reference reference ) {
        if ( _references == null ) {
            _references = new ArrayList<Reference>();
        }
        _references.add( reference );
    }

    public void addSequence( final Sequence sequence ) {
        if ( _sequences == null ) {
            _sequences = new ArrayList<Sequence>();
        }
        _sequences.add( sequence );
    }

    public void addTaxonomy( final Taxonomy taxonomy ) {
        if ( _taxonomies == null ) {
            _taxonomies = new ArrayList<Taxonomy>();
        }
        _taxonomies.add( taxonomy );
    }

    @Override
    public StringBuffer asSimpleText() {
        throw new UnsupportedOperationException();
    }

    @Override
    public StringBuffer asText() {
        throw new UnsupportedOperationException();
    }

    @Override
    public PhylogenyData copy() {
        final NodeData new_data = new NodeData();
        new_data.setNodeName( getNodeName() );
        if ( ( getSequences() != null ) && ( getSequences().size() > 0 ) ) {
            new_data.setSequences( new ArrayList<Sequence>() );
            for( final Sequence s : getSequences() ) {
                if ( s != null ) {
                    new_data.addSequence( ( Sequence ) s.copy() );
                }
            }
        }
        if ( isHasEvent() ) {
            new_data.setEvent( ( Event ) getEvent().copy() );
        }
        if ( ( getTaxonomies() != null ) && ( getTaxonomies().size() > 0 ) ) {
            new_data.setTaxonomies( new ArrayList<Taxonomy>() );
            for( final Taxonomy t : getTaxonomies() ) {
                if ( t != null ) {
                    new_data.addTaxonomy( ( Taxonomy ) t.copy() );
                }
            }
        }
        if ( isHasBinaryCharacters() ) {
            new_data.setBinaryCharacters( ( BinaryCharacters ) getBinaryCharacters().copy() );
        }
        if ( ( getReferences() != null ) && ( getReferences().size() > 0 ) ) {
            new_data.setReferences( new ArrayList<Reference>() );
            for( final Reference r : getReferences() ) {
                if ( r != null ) {
                    new_data.addReference( ( Reference ) r.copy() );
                }
            }
        }
        if ( ( getDistributions() != null ) && ( getDistributions().size() > 0 ) ) {
            new_data.setDistributions( new ArrayList<Distribution>() );
            for( final Distribution d : getDistributions() ) {
                if ( d != null ) {
                    new_data.addDistribution( ( Distribution ) d.copy() );
                }
            }
        }
        if ( ( getNodeVisualData() != null ) && !getNodeVisualData().isEmpty() ) {
            new_data.setNodeVisualData( ( NodeVisualData ) getNodeVisualData().copy() );
        }
        if ( isHasDate() ) {
            new_data.setDate( ( Date ) getDate().copy() );
        }
        if ( isHasProperties() ) {
            new_data.setProperties( ( PropertiesMap ) getProperties().copy() );
        }
        return new_data;
    }

    public BinaryCharacters getBinaryCharacters() {
        return _binary_characters;
    }

    public Date getDate() {
        return _date;
    }

    /**
     * Convenience method -- always returns the first Distribution.
     *  
     * @return Distribution
     */
    public Distribution getDistribution() {
        return getDistribution( 0 );
    }

    public Distribution getDistribution( final int index ) {
        return _distributions.get( index );
    }

    public List<Distribution> getDistributions() {
        return _distributions;
    }

    public Event getEvent() {
        return _event;
    }

    public PropertiesMap getProperties() {
        return _properties;
    }

    /**
     * Convenience method -- always returns the first Reference.
     * 
     *  @return Reference
     *  
     */
    public Reference getReference() {
        return getReference( 0 );
    }

    public Reference getReference( final int index ) {
        return _references.get( index );
    }

    public List<Reference> getReferences() {
        return _references;
    }

    /**
     * Convenience method -- always returns the first Sequence.
     * 
     * @return Sequence
     */
    public Sequence getSequence() {
        return getSequence( 0 );
    }

    public Sequence getSequence( final int index ) {
        return _sequences.get( index );
    }

    public List<Sequence> getSequences() {
        return _sequences;
    }

    public List<Taxonomy> getTaxonomies() {
        return _taxonomies;
    }

    /**
     * Convenience method -- always returns the first Taxonomy.
     * 
     * @return  Taxonomy
     * 
     */
    public Taxonomy getTaxonomy() {
        return getTaxonomy( 0 );
    }

    public Taxonomy getTaxonomy( final int index ) {
        return _taxonomies.get( index );
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new NoSuchMethodError();
    }

    public boolean isHasBinaryCharacters() {
        return getBinaryCharacters() != null;
    }

    public boolean isEmpty() {
        return ( ForesterUtil.isEmpty( _node_name ) && !isHasSequence() && !isHasTaxonomy() && !isHasBinaryCharacters()
                && !isHasDate() && !isHasDistribution() && !isHasEvent() && !isHasProperties() && !isHasReference() && ( ( _vector == null ) || _vector
                .isEmpty() ) );
    }

    public boolean isHasDate() {
        return ( getDate() != null )
                && ( !ForesterUtil.isEmpty( getDate().getDesc() ) || !ForesterUtil.isNull( getDate().getMax() )
                        || !ForesterUtil.isNull( getDate().getMin() ) || !ForesterUtil.isNull( getDate().getValue() ) || !ForesterUtil
                        .isEmpty( getDate().getUnit() ) );
    }

    public boolean isHasDistribution() {
        return ( ( ( getDistributions() != null ) && ( getDistributions().size() > 0 ) ) && ( ( !ForesterUtil
                .isEmpty( getDistribution().getDesc() ) )
                || ( ( getDistribution().getPoints() != null ) && ( getDistribution().getPoints().size() > 0 ) ) || ( ( getDistribution()
                .getPolygons() != null ) && ( getDistribution().getPolygons().size() > 0 ) ) ) );
    }

    public boolean isHasEvent() {
        return getEvent() != null;
    }

    public boolean isHasProperties() {
        return ( getProperties() != null ) && ( getProperties().size() > 0 );
    }

    public boolean isHasReference() {
        return ( ( getReferences() != null ) && ( getReferences().size() > 0 ) )
                && ( !ForesterUtil.isEmpty( getReference().getDoi() ) || !ForesterUtil.isEmpty( getReference()
                        .getDescription() ) );
    }

    public boolean isHasSequence() {
        return ( getSequences() != null ) && ( getSequences().size() > 0 ) && ( getSequences().get( 0 ) != null );
    }

    public boolean isHasTaxonomy() {
        return ( getTaxonomies() != null ) && ( getTaxonomies().size() > 0 ) && ( getTaxonomies().get( 0 ) != null );
    }

    public void setBinaryCharacters( final BinaryCharacters binary_characters ) {
        _binary_characters = binary_characters;
    }

    public void setDate( final Date date ) {
        _date = date;
    }

    /**
     * Convenience method -- always sets the first Distribution.
     * 
     */
    public void setDistribution( final Distribution distribution ) {
        if ( _distributions == null ) {
            _distributions = new ArrayList<Distribution>();
        }
        if ( _distributions.size() == 0 ) {
            _distributions.add( distribution );
        }
        else {
            _distributions.set( 0, distribution );
        }
    }

    public void setDistribution( final int index, final Distribution distribution ) {
        if ( _distributions == null ) {
            _distributions = new ArrayList<Distribution>();
        }
        _distributions.set( index, distribution );
    }

    private void setDistributions( final List<Distribution> distributions ) {
        _distributions = distributions;
    }

    public void setEvent( final Event event ) {
        _event = event;
    }

    public void setProperties( final PropertiesMap custom_data ) {
        _properties = custom_data;
    }

    public void setReference( final int index, final Reference reference ) {
        if ( _references == null ) {
            _references = new ArrayList<Reference>();
        }
        _references.set( index, reference );
    }

    /**
     * Convenience method -- always sets the first Reference.
     * 
     */
    public void setReference( final Reference reference ) {
        if ( _references == null ) {
            _references = new ArrayList<Reference>();
        }
        if ( _references.size() == 0 ) {
            _references.add( reference );
        }
        else {
            _references.set( 0, reference );
        }
    }

    private void setReferences( final List<Reference> references ) {
        _references = references;
    }

    public void setSequence( final int index, final Sequence sequence ) {
        if ( _sequences == null ) {
            _sequences = new ArrayList<Sequence>();
        }
        _sequences.set( index, sequence );
    }

    /**
     * Convenience method -- always sets the first Sequence.
     * 
     */
    public void setSequence( final Sequence sequence ) {
        if ( _sequences == null ) {
            _sequences = new ArrayList<Sequence>();
        }
        if ( _sequences.size() == 0 ) {
            _sequences.add( sequence );
        }
        else {
            _sequences.set( 0, sequence );
        }
    }

    private void setSequences( final List<Sequence> sequences ) {
        _sequences = sequences;
    }

    private void setTaxonomies( final List<Taxonomy> taxonomies ) {
        _taxonomies = taxonomies;
    }

    public void setTaxonomy( final int index, final Taxonomy taxonomy ) {
        if ( _taxonomies == null ) {
            _taxonomies = new ArrayList<Taxonomy>();
        }
        _taxonomies.set( index, taxonomy );
    }

    /**
     * Convenience method -- always sets the first Taxonomy.
     * 
     */
    public void setTaxonomy( final Taxonomy taxonomy ) {
        if ( _taxonomies == null ) {
            _taxonomies = new ArrayList<Taxonomy>();
        }
        if ( _taxonomies.size() == 0 ) {
            _taxonomies.add( taxonomy );
        }
        else {
            _taxonomies.set( 0, taxonomy );
        }
    }

    @Override
    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        if ( isHasTaxonomy() ) {
            sb.append( getTaxonomy().toNHX() );
        }
        if ( isHasSequence() ) {
            sb.append( getSequence().toNHX() );
        }
        if ( isHasEvent() ) {
            sb.append( getEvent().toNHX() );
        }
        return sb;
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( isHasTaxonomy() ) {
            for( final Taxonomy t : getTaxonomies() ) {
                if ( !t.isEmpty() ) {
                    t.toPhyloXML( writer, level, indentation );
                }
            }
        }
        if ( isHasSequence() ) {
            for( final Sequence s : getSequences() ) {
                if ( !s.isEmpty() ) {
                    s.toPhyloXML( writer, level, indentation );
                }
            }
        }
        if ( isHasEvent() ) {
            getEvent().toPhyloXML( writer, level, indentation );
        }
        if ( isHasBinaryCharacters() ) {
            getBinaryCharacters().toPhyloXML( writer, level, indentation );
        }
        if ( isHasDistribution() ) {
            for( final Distribution d : getDistributions() ) {
                d.toPhyloXML( writer, level, indentation );
            }
        }
        if ( isHasDate() ) {
            getDate().toPhyloXML( writer, level, indentation );
        }
        if ( isHasReference() ) {
            for( final Reference r : getReferences() ) {
                r.toPhyloXML( writer, level, indentation );
            }
        }
        if ( isHasProperties() ) {
            getProperties().toPhyloXML( writer, level, indentation.substring( 0, indentation.length() - 2 ) );
        }
        if ( ( level == 0 ) && ( getNodeVisualData() != null ) && !getNodeVisualData().isEmpty() ) {
            getNodeVisualData().toPhyloXML( writer, level, indentation.substring( 0, indentation.length() - 2 ) );
        }
        if ( ( getVector() != null )
                && !getVector().isEmpty()
                && ( ( getProperties() == null ) || getProperties()
                        .getPropertiesWithGivenReferencePrefix( PhyloXmlUtil.VECTOR_PROPERTY_REF ).isEmpty() ) ) {
            final List<Property> ps = vectorToProperties( getVector() );
            final String my_indent = indentation.substring( 0, indentation.length() - 2 );
            for( final Property p : ps ) {
                p.toPhyloXML( writer, level, my_indent );
            }
        }
    }

    private List<Property> vectorToProperties( final List<Double> vector ) {
        final List<Property> properties = new ArrayList<Property>();
        for( int i = 0; i < vector.size(); ++i ) {
            properties.add( new Property( PhyloXmlUtil.VECTOR_PROPERTY_REF + i,
                                          String.valueOf( vector.get( i ) ),
                                          "",
                                          PhyloXmlUtil.VECTOR_PROPERTY_TYPE,
                                          AppliesTo.NODE ) );
        }
        return properties;
    }

    public void setVector( final List<Double> vector ) {
        _vector = vector;
    }

    public List<Double> getVector() {
        return _vector;
    }

    public String getNodeName() {
        return _node_name;
    }

    public void setNodeName( final String node_name ) {
        _node_name = node_name;
    }

    public void setNodeVisualData( final NodeVisualData node_visual_data ) {
        _node_visual_data = node_visual_data;
    }

    public NodeVisualData getNodeVisualData() {
        return _node_visual_data;
    }
}
