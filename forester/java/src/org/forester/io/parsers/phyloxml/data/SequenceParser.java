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

package org.forester.io.parsers.phyloxml.data;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Uri;

public class SequenceParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new SequenceParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private SequenceParser() {
    }

    @Override
    public Sequence parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        final Sequence sequence = new Sequence();
        if ( element.isHasAttribute( PhyloXmlMapping.SEQUENCE_TYPE ) ) {
            sequence.setType( element.getAttribute( PhyloXmlMapping.SEQUENCE_TYPE ) );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.SEQUENCE_SOURCE_ID ) ) {
            sequence.setSourceId( element.getAttribute( PhyloXmlMapping.SEQUENCE_SOURCE_ID ) );
        }
        for( int i = 0; i < element.getNumberOfChildElements(); ++i ) {
            final XmlElement child_element = element.getChildElement( i );
            if ( child_element.getQualifiedName().equals( PhyloXmlMapping.SEQUENCE_LOCATION ) ) {
                sequence.setLocation( child_element.getValueAsString() );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.SEQUENCE_NAME ) ) {
                sequence.setName( child_element.getValueAsString() );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.SEQUENCE_MOL_SEQ ) ) {
                if ( child_element.isHasAttribute( PhyloXmlMapping.SEQUENCE_MOL_SEQ_ALIGNED_ATTR ) ) {
                    sequence.setMolecularSequenceAligned( Boolean.parseBoolean( child_element
                            .getAttribute( PhyloXmlMapping.SEQUENCE_MOL_SEQ_ALIGNED_ATTR ) ) );
                }
                sequence.setMolecularSequence( child_element.getValueAsString() );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.ACCESSION ) ) {
                sequence.setAccession( ( Accession ) AccessionParser.getInstance().parse( child_element ) );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.SEQUENCE_SYMBOL ) ) {
                sequence.setSymbol( child_element.getValueAsString() );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.ANNOTATION ) ) {
                sequence.addAnnotation( ( Annotation ) AnnotationParser.getInstance().parse( child_element ) );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECURE ) ) {
                sequence.setDomainArchitecture( ( DomainArchitecture ) DomainArchitectureParser.getInstance()
                        .parse( child_element ) );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.URI ) ) {
                sequence.addUri( ( Uri ) UriParser.getInstance().parse( child_element ) );
            }
        }
        return sequence;
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
