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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.util.ForesterUtil;

public class Sequence implements PhylogenyData, MultipleUris, Comparable<Sequence> {

    private Accession              _accession;
    private SortedSet<Annotation>  _annotations;
    private DomainArchitecture     _da;
    private String                 _gene_name;
    private String                 _location;
    private String                 _mol_sequence;
    private boolean                _mol_sequence_is_aligned;
    private String                 _name;
    private List<SequenceRelation> _seq_relations;
    private String                 _source_id;
    private String                 _symbol;
    private String                 _type;
    private List<Uri>              _uris;
    private SortedSet<Accession>   _xrefs;

    public Sequence() {
        init();
    }

    public void addAnnotation( final Annotation annotation ) {
        getAnnotations().add( annotation );
    }

    public void addCrossReference( final Accession cross_reference ) {
        if ( getCrossReferences() == null ) {
            setCrossReferences( new TreeSet<Accession>() );
        }
        getCrossReferences().add( cross_reference );
    }

    public void addSequenceRelation( final SequenceRelation sr ) {
        getSequenceRelations().add( sr );
    }

    @Override
    public void addUri( final Uri uri ) {
        if ( getUris() == null ) {
            setUris( new ArrayList<Uri>() );
        }
        getUris().add( uri );
    }

    @Override
    public StringBuffer asSimpleText() {
        final StringBuffer sb = new StringBuffer();
        if ( getAccession() != null ) {
            sb.append( "[" );
            sb.append( getAccession() );
            sb.append( "] " );
        }
        if ( !ForesterUtil.isEmpty( getName() ) ) {
            sb.append( getName() );
            sb.append( " " );
        }
        if ( !ForesterUtil.isEmpty( getLocation() ) ) {
            sb.append( getLocation() );
        }
        return sb;
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public int compareTo( final Sequence o ) {
        if ( ( !ForesterUtil.isEmpty( getName() ) ) && ( !ForesterUtil.isEmpty( o.getName() ) ) ) {
            return getName().compareTo( o.getName() );
        }
        if ( ( !ForesterUtil.isEmpty( getSymbol() ) ) && ( !ForesterUtil.isEmpty( o.getSymbol() ) ) ) {
            return getSymbol().compareTo( o.getSymbol() );
        }
        if ( ( !ForesterUtil.isEmpty( getGeneName() ) ) && ( !ForesterUtil.isEmpty( o.getGeneName() ) ) ) {
            return getGeneName().compareTo( o.getGeneName() );
        }
        if ( ( getAccession() != null ) && ( o.getAccession() != null )
                && !ForesterUtil.isEmpty( getAccession().getValue() )
                && !ForesterUtil.isEmpty( o.getAccession().getValue() ) ) {
            return getAccession().getValue().compareTo( o.getAccession().getValue() );
        }
        if ( ( !ForesterUtil.isEmpty( getMolecularSequence() ) )
                && ( !ForesterUtil.isEmpty( o.getMolecularSequence() ) ) ) {
            return getMolecularSequence().compareTo( o.getMolecularSequence() );
        }
        return 0;
    }

    /**
     * Not a deep copy.
     * 
     */
    @Override
    public PhylogenyData copy() {
        final Sequence seq = new Sequence();
        seq.setAnnotations( getAnnotations() );
        seq.setName( getName() );
        seq.setGeneName( getGeneName() );
        try {
            seq.setSymbol( getSymbol() );
        }
        catch ( final PhyloXmlDataFormatException e ) {
            e.printStackTrace();
        }
        seq.setMolecularSequence( getMolecularSequence() );
        seq.setMolecularSequenceAligned( isMolecularSequenceAligned() );
        seq.setLocation( getLocation() );
        if ( getAccession() != null ) {
            seq.setAccession( ( Accession ) getAccession().copy() );
        }
        else {
            seq.setAccession( null );
        }
        try {
            seq.setType( getType() );
        }
        catch ( final PhyloXmlDataFormatException e ) {
            e.printStackTrace();
        }
        if ( getUris() != null ) {
            seq.setUris( new ArrayList<Uri>() );
            for( final Uri uri : getUris() ) {
                if ( uri != null ) {
                    seq.getUris().add( uri );
                }
            }
        }
        if ( getDomainArchitecture() != null ) {
            seq.setDomainArchitecture( ( DomainArchitecture ) getDomainArchitecture().copy() );
        }
        else {
            seq.setDomainArchitecture( null );
        }
        if ( getCrossReferences() != null ) {
            seq.setCrossReferences( new TreeSet<Accession>() );
            for( final Accession x : getCrossReferences() ) {
                if ( x != null ) {
                    seq.getCrossReferences().add( x );
                }
            }
        }
        return seq;
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            return false;
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return isEqual( ( Sequence ) o );
        }
    }

    public Accession getAccession() {
        return _accession;
    }

    public Annotation getAnnotation( final int i ) {
        return ( Annotation ) getAnnotations().toArray()[ i ];
    }

    public SortedSet<Annotation> getAnnotations() {
        if ( _annotations == null ) {
            _annotations = new TreeSet<Annotation>();
        }
        return _annotations;
    }

    public SortedSet<Accession> getCrossReferences() {
        return _xrefs;
    }

    public DomainArchitecture getDomainArchitecture() {
        return _da;
    }

    public String getGeneName() {
        return _gene_name;
    }

    public String getLocation() {
        return _location;
    }

    public String getMolecularSequence() {
        return _mol_sequence;
    }

    public String getName() {
        return _name;
    }

    public List<SequenceRelation> getSequenceRelations() {
        if ( _seq_relations == null ) {
            _seq_relations = new ArrayList<SequenceRelation>();
        }
        return _seq_relations;
    }

    public String getSourceId() {
        return _source_id;
    }

    public String getSymbol() {
        return _symbol;
    }

    public String getType() {
        return _type;
    }

    @Override
    public Uri getUri( final int index ) {
        return getUris().get( index );
    }

    @Override
    public List<Uri> getUris() {
        return _uris;
    }

    @Override
    public int hashCode() {
        if ( getAccession() != null ) {
            return getAccession().hashCode();
        }
        int result = getName().hashCode();
        if ( getSymbol().length() > 0 ) {
            result ^= getName().hashCode();
        }
        if ( getGeneName().length() > 0 ) {
            result ^= getGeneName().hashCode();
        }
        if ( getMolecularSequence().length() > 0 ) {
            result ^= getMolecularSequence().hashCode();
        }
        return result;
    }

    public boolean hasSequenceRelations() {
        return _seq_relations.size() > 0;
    }

    public void init() {
        setName( "" );
        setGeneName( "" );
        setMolecularSequence( "" );
        setMolecularSequenceAligned( false );
        setLocation( "" );
        setAccession( null );
        try {
            setSymbol( "" );
        }
        catch ( final PhyloXmlDataFormatException e ) {
            e.printStackTrace();
        }
        try {
            setType( "" );
        }
        catch ( final PhyloXmlDataFormatException e ) {
            e.printStackTrace();
        }
        setDomainArchitecture( null );
        setUris( null );
        setSequenceRelations( null );
        setSourceId( null );
        setCrossReferences( null );
        setAnnotations( null );
    }

    public boolean isEmpty() {
        return ( getAccession() == null ) && ForesterUtil.isEmpty( getName() ) && ForesterUtil.isEmpty( getSymbol() )
                && ForesterUtil.isEmpty( getGeneName() ) && ForesterUtil.isEmpty( getType() )
                && ForesterUtil.isEmpty( getLocation() ) && ForesterUtil.isEmpty( getSourceId() )
                && ForesterUtil.isEmpty( getMolecularSequence() ) && ( getDomainArchitecture() == null )
                && ForesterUtil.isEmpty( _annotations ) && ForesterUtil.isEmpty( _uris )
                && ForesterUtil.isEmpty( _seq_relations )
                && ( ( getCrossReferences() == null ) || getCrossReferences().isEmpty() );
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        if ( this == data ) {
            return true;
        }
        final Sequence s = ( Sequence ) data;
        if ( ( getAccession() != null ) && ( s.getAccession() != null ) ) {
            return getAccession().isEqual( s.getAccession() );
        }
        return s.getMolecularSequence().equals( getMolecularSequence() ) && s.getName().equals( getName() )
                && s.getSymbol().equals( getSymbol() ) && s.getGeneName().equals( getGeneName() );
    }

    public boolean isMolecularSequenceAligned() {
        return _mol_sequence_is_aligned;
    }

    public void setAccession( final Accession accession ) {
        _accession = accession;
    }

    public void setDomainArchitecture( final DomainArchitecture ds ) {
        _da = ds;
    }

    public void setGeneName( final String gene_name ) {
        _gene_name = gene_name;
    }

    public void setLocation( final String description ) {
        _location = description;
    }

    public void setMolecularSequence( final String mol_sequence ) {
        _mol_sequence = mol_sequence;
    }

    public void setMolecularSequenceAligned( final boolean aligned ) {
        _mol_sequence_is_aligned = aligned;
    }

    public void setName( final String name ) {
        _name = name;
    }

    public void setSourceId( final String source_id ) {
        _source_id = source_id;
    }

    public void setSymbol( final String symbol ) throws PhyloXmlDataFormatException {
        if ( !ForesterUtil.isEmpty( symbol ) && !PhyloXmlUtil.SEQUENCE_SYMBOL_PATTERN.matcher( symbol ).matches() ) {
            throw new PhyloXmlDataFormatException( "illegal sequence symbol: [" + symbol + "]" );
        }
        _symbol = symbol;
    }

    public void setType( final String type ) throws PhyloXmlDataFormatException {
        if ( !ForesterUtil.isEmpty( type ) && !PhyloXmlUtil.SEQUENCE_TYPES.contains( type ) ) {
            throw new PhyloXmlDataFormatException( "illegal sequence type: [" + type + "]" );
        }
        _type = type;
    }

    @Override
    public void setUris( final List<Uri> uris ) {
        _uris = uris;
    }

    @Override
    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        if ( getName().length() > 0 ) {
            sb.append( ":" );
            sb.append( NHXtags.GENE_NAME );
            sb.append( ForesterUtil.replaceIllegalNhxCharacters( getName() ) );
        }
        if ( getAccession() != null ) {
            getAccession().toNHX();
        }
        return sb;
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( isEmpty() ) {
            return;
        }
        final String my_ind = indentation + PhylogenyWriter.PHYLO_XML_INTENDATION_BASE;
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.SEQUENCE, PhyloXmlMapping.SEQUENCE_TYPE, getType() );
        if ( !ForesterUtil.isEmpty( getSymbol() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.SEQUENCE_SYMBOL, getSymbol(), indentation );
        }
        if ( ( getAccession() != null ) && !ForesterUtil.isEmpty( getAccession().getValue() ) ) {
            getAccession().toPhyloXML( writer, level, indentation );
        }
        if ( !ForesterUtil.isEmpty( getName() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.SEQUENCE_NAME, getName(), indentation );
        }
        if ( !ForesterUtil.isEmpty( getGeneName() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.SEQUENCE_GENE_NAME, getGeneName(), indentation );
        }
        if ( !ForesterUtil.isEmpty( getLocation() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.SEQUENCE_LOCATION, getLocation(), indentation );
        }
        if ( !ForesterUtil.isEmpty( getMolecularSequence() ) ) {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.SEQUENCE_MOL_SEQ,
                                             getMolecularSequence(),
                                             PhyloXmlMapping.SEQUENCE_MOL_SEQ_ALIGNED_ATTR,
                                             String.valueOf( isMolecularSequenceAligned() ),
                                             indentation );
        }
        if ( ( getUris() != null ) && !getUris().isEmpty() ) {
            for( final Uri uri : getUris() ) {
                if ( uri != null ) {
                    uri.toPhyloXML( writer, level, indentation );
                }
            }
        }
        if ( ( getAnnotations() != null ) && !getAnnotations().isEmpty() ) {
            for( final PhylogenyData annotation : getAnnotations() ) {
                annotation.toPhyloXML( writer, level, my_ind );
            }
        }
        if ( ( getCrossReferences() != null ) && !getCrossReferences().isEmpty() ) {
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( my_ind );
            PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.SEQUENCE_X_REFS );
            for( final PhylogenyData x : getCrossReferences() ) {
                x.toPhyloXML( writer, level, my_ind );
            }
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( my_ind );
            PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.SEQUENCE_X_REFS );
        }
        if ( getDomainArchitecture() != null ) {
            getDomainArchitecture().toPhyloXML( writer, level, my_ind );
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.SEQUENCE );
    }

    @Override
    public String toString() {
        return asText().toString();
    }

    private void setAnnotations( final SortedSet<Annotation> annotations ) {
        _annotations = annotations;
    }

    private void setCrossReferences( final TreeSet<Accession> cross_references ) {
        _xrefs = cross_references;
    }

    private void setSequenceRelations( final List<SequenceRelation> seq_relations ) {
        _seq_relations = seq_relations;
    }
}
