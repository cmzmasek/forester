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
// WWW: www.phylosoft.org/forester

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.util.ForesterUtil;

public class Taxonomy implements PhylogenyData, MultipleUris, Comparable<Taxonomy> {

    private String       _scientific_name;
    private String       _common_name;
    private List<String> _synonyms;
    private String       _authority;
    private Identifier   _identifier;
    private String       _taxonomy_code;
    private String       _rank;
    private List<Uri>    _uris;
    private List<String> _lineage;

    public Taxonomy() {
        init();
    }

    @Override
    public StringBuffer asSimpleText() {
        return asText();
    }

    @Override
    public Uri getUri( final int index ) {
        return getUris().get( index );
    }

    @Override
    public void addUri( final Uri uri ) {
        if ( getUris() == null ) {
            setUris( new ArrayList<Uri>() );
        }
        getUris().add( uri );
    }

    @Override
    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        if ( getIdentifier() != null ) {
            sb.append( "[" );
            sb.append( getIdentifier().asSimpleText() );
            sb.append( "]" );
        }
        if ( !ForesterUtil.isEmpty( getTaxonomyCode() ) ) {
            if ( sb.length() > 0 ) {
                sb.append( " " );
            }
            sb.append( "[" );
            sb.append( getTaxonomyCode() );
            sb.append( "]" );
        }
        if ( !ForesterUtil.isEmpty( getScientificName() ) ) {
            if ( sb.length() > 0 ) {
                sb.append( " " );
            }
            sb.append( getScientificName() );
            if ( !ForesterUtil.isEmpty( getAuthority() ) ) {
                sb.append( " (" );
                sb.append( getAuthority() );
                sb.append( ")" );
            }
        }
        if ( !ForesterUtil.isEmpty( getCommonName() ) ) {
            if ( sb.length() > 0 ) {
                sb.append( " " );
            }
            sb.append( getCommonName() );
        }
        return sb;
    }

    @Override
    public PhylogenyData copy() {
        final Taxonomy t = new Taxonomy();
        try {
            t.setTaxonomyCode( getTaxonomyCode() );
        }
        catch ( final PhyloXmlDataFormatException e ) {
            e.printStackTrace();
        }
        t.setScientificName( getScientificName() );
        t.setCommonName( getCommonName() );
        t.setAuthority( getAuthority() );
        for( final String syn : getSynonyms() ) {
            t.getSynonyms().add( syn );
        }
        if ( getIdentifier() != null ) {
            t.setIdentifier( ( Identifier ) getIdentifier().copy() );
        }
        else {
            t.setIdentifier( null );
        }
        try {
            t.setRank( new String( getRank() ) );
        }
        catch ( final PhyloXmlDataFormatException e ) {
            e.printStackTrace();
        }
        if ( getUris() != null ) {
            t.setUris( new ArrayList<Uri>() );
            for( final Uri uri : getUris() ) {
                if ( uri != null ) {
                    t.getUris().add( uri );
                }
            }
        }
        if ( getLineage() != null ) {
            t.setLineage( new ArrayList<String>() );
            for( final String l : getLineage() ) {
                if ( l != null ) {
                    t.getLineage().add( l );
                }
            }
        }
        return t;
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
            return isEqual( ( Taxonomy ) o );
        }
    }

    public String getAuthority() {
        return _authority;
    }

    public String getCommonName() {
        return _common_name;
    }

    public Identifier getIdentifier() {
        return _identifier;
    }

    public String getRank() {
        return _rank;
    }

    public String getScientificName() {
        return _scientific_name;
    }

    public List<String> getSynonyms() {
        if ( _synonyms == null ) {
            _synonyms = new ArrayList<String>();
        }
        return _synonyms;
    }

    public String getTaxonomyCode() {
        return _taxonomy_code;
    }

    @Override
    public List<Uri> getUris() {
        return _uris;
    }

    @Override
    public int hashCode() {
        if ( ( getIdentifier() != null ) && !ForesterUtil.isEmpty( getIdentifier().getValue() ) ) {
            return getIdentifier().hashCode();
        }
        else if ( !ForesterUtil.isEmpty( getTaxonomyCode() ) ) {
            return getTaxonomyCode().hashCode();
        }
        else if ( !ForesterUtil.isEmpty( getScientificName() ) ) {
            if ( !ForesterUtil.isEmpty( getAuthority() ) ) {
                return ( getScientificName().toLowerCase() + getAuthority().toLowerCase() ).hashCode();
            }
            return getScientificName().toLowerCase().hashCode();
        }
        else {
            return getCommonName().toLowerCase().hashCode();
        }
    }

    public void init() {
        setScientificName( "" );
        setCommonName( "" );
        setIdentifier( null );
        try {
            setRank( "" );
        }
        catch ( final PhyloXmlDataFormatException e ) {
            e.printStackTrace();
        }
        try {
            setTaxonomyCode( "" );
        }
        catch ( final PhyloXmlDataFormatException e ) {
            e.printStackTrace();
        }
        setAuthority( "" );
        setSynonyms( null );
        setUris( null );
        setLineage( null );
    }

    public boolean isEmpty() {
        return ( ( getIdentifier() == null ) && ForesterUtil.isEmpty( getTaxonomyCode() )
                && ForesterUtil.isEmpty( getCommonName() ) && ForesterUtil.isEmpty( getScientificName() )
                && ForesterUtil.isEmpty( getRank() ) && ForesterUtil.isEmpty( _uris )
                && ForesterUtil.isEmpty( getAuthority() ) && ForesterUtil.isEmpty( _synonyms ) && ForesterUtil
                    .isEmpty( _lineage ) );
    }

    /**
     * 
     * If this and taxonomy 'data' has an identifier, comparison will be based on that.
     * Otherwise,  if this and taxonomy 'data' has a code, comparison will be based on that.
     * Otherwise,  if Taxonomy 'data' has a scientific name, comparison will be
     * based on that (case insensitive!).
     * Otherwise,  if Taxonomy 'data' has a common  name, comparison will be
     * based on that (case insensitive!).
     * (Note. This is important and should not be change without a very good reason.)
     * 
     */
    @Override
    public boolean isEqual( final PhylogenyData data ) {
        if ( this == data ) {
            return true;
        }
        final Taxonomy tax = ( Taxonomy ) data;
        if ( ( getIdentifier() != null ) && ( tax.getIdentifier() != null )
                && !ForesterUtil.isEmpty( getIdentifier().getValue() )
                && !ForesterUtil.isEmpty( tax.getIdentifier().getValue() ) ) {
            return getIdentifier().isEqual( tax.getIdentifier() );
        }
        else if ( !ForesterUtil.isEmpty( getTaxonomyCode() ) && !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
            return getTaxonomyCode().equals( tax.getTaxonomyCode() );
        }
        else if ( !ForesterUtil.isEmpty( getScientificName() ) && !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
            if ( !ForesterUtil.isEmpty( getAuthority() ) && !ForesterUtil.isEmpty( tax.getAuthority() ) ) {
                return ( getScientificName().equalsIgnoreCase( tax.getScientificName() ) )
                        && ( getAuthority().equalsIgnoreCase( tax.getAuthority() ) );
            }
            return getScientificName().equalsIgnoreCase( tax.getScientificName() );
        }
        else if ( !ForesterUtil.isEmpty( getCommonName() ) && !ForesterUtil.isEmpty( tax.getCommonName() ) ) {
            return getCommonName().equalsIgnoreCase( tax.getCommonName() );
        }
        else if ( !ForesterUtil.isEmpty( getScientificName() ) && !ForesterUtil.isEmpty( tax.getCommonName() ) ) {
            return getScientificName().equalsIgnoreCase( tax.getCommonName() );
        }
        else if ( !ForesterUtil.isEmpty( getCommonName() ) && !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
            return getCommonName().equalsIgnoreCase( tax.getScientificName() );
        }
        throw new RuntimeException( "comparison not possible with empty fields" );
    }

    public void setAuthority( final String authority ) {
        _authority = authority;
    }

    public void setCommonName( final String common_name ) {
        _common_name = common_name;
    }

    public void setIdentifier( final Identifier identifier ) {
        _identifier = identifier;
    }

    public void setRank( final String rank ) throws PhyloXmlDataFormatException {
        if ( !ForesterUtil.isEmpty( rank ) && !PhyloXmlUtil.TAXONOMY_RANKS_SET.contains( rank ) ) {
            throw new PhyloXmlDataFormatException( "illegal rank: [" + rank + "]" );
        }
        _rank = rank;
    }

    public void setScientificName( final String scientific_name ) {
        _scientific_name = scientific_name;
    }

    private void setSynonyms( final List<String> synonyms ) {
        _synonyms = synonyms;
    }

    public void setTaxonomyCode( final String taxonomy_code ) throws PhyloXmlDataFormatException {
        if ( !ForesterUtil.isEmpty( taxonomy_code )
                && !PhyloXmlUtil.TAXOMONY_CODE_PATTERN.matcher( taxonomy_code ).matches() ) {
            throw new PhyloXmlDataFormatException( "illegal taxonomy code: [" + taxonomy_code + "]" );
        }
        _taxonomy_code = taxonomy_code;
    }

    @Override
    public void setUris( final List<Uri> uris ) {
        _uris = uris;
    }

    @Override
    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        if ( getIdentifier() != null ) {
            sb.append( ':' + NHXtags.TAXONOMY_ID );
            sb.append( ForesterUtil.replaceIllegalNhxCharacters( getIdentifier().getValue() ) );
        }
        final StringBuffer species = new StringBuffer();
        if ( !ForesterUtil.isEmpty( getTaxonomyCode() ) ) {
            species.append( ForesterUtil.replaceIllegalNhxCharacters( getTaxonomyCode() ) );
        }
        if ( !ForesterUtil.isEmpty( getScientificName() ) ) {
            ForesterUtil.appendSeparatorIfNotEmpty( species, '|' );
            species.append( ForesterUtil.replaceIllegalNhxCharacters( getScientificName() ) );
        }
        if ( !ForesterUtil.isEmpty( getCommonName() ) ) {
            ForesterUtil.appendSeparatorIfNotEmpty( species, '|' );
            species.append( ForesterUtil.replaceIllegalNhxCharacters( getCommonName() ) );
        }
        if ( species.length() > 0 ) {
            sb.append( ':' + NHXtags.SPECIES_NAME );
            sb.append( species );
        }
        return sb;
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( isEmpty() ) {
            return;
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.TAXONOMY );
        if ( ( getIdentifier() != null ) && !ForesterUtil.isEmpty( getIdentifier().getValue() ) ) {
            getIdentifier().toPhyloXML( writer, level, indentation );
        }
        if ( !ForesterUtil.isEmpty( getTaxonomyCode() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.TAXONOMY_CODE, getTaxonomyCode(), indentation );
        }
        if ( !ForesterUtil.isEmpty( getScientificName() ) ) {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.TAXONOMY_SCIENTIFIC_NAME,
                                             getScientificName(),
                                             indentation );
        }
        if ( !ForesterUtil.isEmpty( getAuthority() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.TAXONOMY_AUTHORITY, getAuthority(), indentation );
        }
        if ( !ForesterUtil.isEmpty( getCommonName() ) ) {
            PhylogenyDataUtil
                    .appendElement( writer, PhyloXmlMapping.TAXONOMY_COMMON_NAME, getCommonName(), indentation );
        }
        if ( _synonyms != null ) {
            for( final String syn : getSynonyms() ) {
                if ( !ForesterUtil.isEmpty( syn ) ) {
                    PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.TAXONOMY_SYNONYM, syn, indentation );
                }
            }
        }
        if ( !ForesterUtil.isEmpty( getRank() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.TAXONOMY_RANK, getRank(), indentation );
        }
        if ( getUris() != null ) {
            for( final Uri uri : getUris() ) {
                if ( uri != null ) {
                    uri.toPhyloXML( writer, level, indentation );
                }
            }
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.TAXONOMY );
    }

    @Override
    public String toString() {
        return asText().toString();
    }

    @Override
    public int compareTo( final Taxonomy o ) {
        if ( equals( o ) ) {
            return 0;
        }
        else if ( !ForesterUtil.isEmpty( getScientificName() ) && !ForesterUtil.isEmpty( o.getScientificName() ) ) {
            return getScientificName().compareToIgnoreCase( o.getScientificName() );
        }
        else if ( !ForesterUtil.isEmpty( getCommonName() ) && !ForesterUtil.isEmpty( o.getCommonName() ) ) {
            return getCommonName().compareToIgnoreCase( o.getCommonName() );
        }
        else if ( !ForesterUtil.isEmpty( getTaxonomyCode() ) && !ForesterUtil.isEmpty( o.getTaxonomyCode() ) ) {
            return getTaxonomyCode().compareToIgnoreCase( o.getTaxonomyCode() );
        }
        return 0;
    }

    public void setLineage( final List<String> lineage ) {
        _lineage = lineage;
    }

    public List<String> getLineage() {
        return _lineage;
    }
}
