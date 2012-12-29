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

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.util.ForesterUtil;

public class Annotation implements PhylogenyData, MultipleUris, Comparable<Annotation> {

    private String        _desc;
    private String        _type;
    private String        _source;
    private final String  _ref;
    private String        _evidence;
    private Confidence    _confidence;
    private PropertiesMap _properties;
    private List<Uri>     _uris;

    public Annotation( final String ref ) {
        if ( ForesterUtil.isEmpty( ref ) ) {
            throw new IllegalArgumentException( "annotation reference is empty or null" );
        }
        if ( ( ref.indexOf( ':' ) < 1 ) || ( ref.indexOf( ':' ) > ref.length() - 2 ) || ( ref.length() < 3 ) ) {
            throw new IllegalArgumentException( "illegal format for annotation reference: [" + ref + "]" );
        }
        _ref = ref;
        init();
    }

    public Annotation() {
        _ref = "";
        init();
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( !ForesterUtil.isEmpty( getRef() ) ? getRef() : getDesc() );
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        final Annotation ann = new Annotation( new String( getRef() ) );
        if ( getConfidence() != null ) {
            ann.setConfidence( ( Confidence ) getConfidence().copy() );
        }
        else {
            ann.setConfidence( null );
        }
        ann.setType( new String( getType() ) );
        ann.setDesc( new String( getDesc() ) );
        ann.setEvidence( new String( getEvidence() ) );
        ann.setSource( new String( getSource() ) );
        if ( getProperties() != null ) {
            ann.setProperties( ( PropertiesMap ) getProperties().copy() );
        }
        else {
            ann.setProperties( null );
        }
        if ( getUris() != null ) {
            ann.setUris( new ArrayList<Uri>() );
            for( final Uri uri : getUris() ) {
                if ( uri != null ) {
                    ann.getUris().add( uri );
                }
            }
        }
        return ann;
    }

    public Confidence getConfidence() {
        return _confidence;
    }

    public String getDesc() {
        return _desc;
    }

    public String getEvidence() {
        return _evidence;
    }

    public PropertiesMap getProperties() {
        return _properties;
    }

    public String getRef() {
        return _ref;
    }

    public String getSource() {
        return _source;
    }

    public String getType() {
        return _type;
    }

    private void init() {
        _desc = "";
        _type = "";
        _source = "";
        _evidence = "";
        _confidence = null;
        _properties = null;
        setUris( null );
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        final Annotation other = ( Annotation ) data;
        return getDesc().equalsIgnoreCase( other.getDesc() ) && getType().equals( other.getType() )
                && getSource().equals( other.getSource() ) && getRef().equals( other.getRef() );
    }

    public void setConfidence( final Confidence confidence ) {
        _confidence = confidence;
    }

    public void setDesc( final String desc ) {
        _desc = desc;
    }

    public void setEvidence( final String evidence ) {
        _evidence = evidence;
    }

    public void setProperties( final PropertiesMap property ) {
        _properties = property;
    }

    // public void setRef( final String ref ) {
    //     _ref = ref;
    // }
    public void setSource( final String source ) {
        _source = source;
    }

    public void setType( final String type ) {
        _type = type;
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( ( getConfidence() != null ) || ( getProperties() != null )
                || ( ( getUris() != null ) && !getUris().isEmpty() ) || !ForesterUtil.isEmpty( getDesc() ) ) {
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( indentation );
            PhylogenyDataUtil.appendOpen( writer,
                                          PhyloXmlMapping.ANNOTATION,
                                          PhyloXmlMapping.ANNOTATION_REF_ATTR,
                                          getRef(),
                                          PhyloXmlMapping.ANNOTATION_EVIDENCE_ATTR,
                                          getEvidence(),
                                          PhyloXmlMapping.ANNOTATION_TYPE_ATTR,
                                          getType(),
                                          PhyloXmlMapping.ANNOTATION_SOURCE_ATTR,
                                          getSource() );
            if ( !ForesterUtil.isEmpty( getDesc() ) ) {
                PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.ANNOTATION_DESC, getDesc(), indentation );
            }
            if ( getConfidence() != null ) {
                getConfidence().toPhyloXML( writer, level, indentation + PhylogenyWriter.PHYLO_XML_INTENDATION_BASE );
            }
            if ( getProperties() != null ) {
                getProperties().toPhyloXML( writer, level, indentation );
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
            PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.ANNOTATION );
        }
        else {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.ANNOTATION,
                                             PhyloXmlMapping.ANNOTATION_REF_ATTR,
                                             getRef(),
                                             PhyloXmlMapping.ANNOTATION_EVIDENCE_ATTR,
                                             getEvidence(),
                                             PhyloXmlMapping.ANNOTATION_TYPE_ATTR,
                                             getType(),
                                             PhyloXmlMapping.ANNOTATION_SOURCE_ATTR,
                                             getSource(),
                                             indentation );
        }
    }

    @Override
    public String toString() {
        return asText().toString();
    }

    @Override
    public void addUri( final Uri uri ) {
        if ( getUris() == null ) {
            setUris( new ArrayList<Uri>() );
        }
        getUris().add( uri );
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
    public void setUris( final List<Uri> uris ) {
        _uris = uris;
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
            return isEqual( ( Annotation ) o );
        }
    }

    @Override
    public int compareTo( final Annotation o ) {
        if ( equals( o ) ) {
            return 0;
        }
        if ( getRef().equals( o.getRef() ) ) {
            return getDesc().compareTo( o.getDesc() );
        }
        return getRef().compareTo( o.getRef() );
    }
}
