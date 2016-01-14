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

import java.util.HashMap;
import java.util.Map;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlHandler;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.SequenceRelation;

public class SequenceRelationParser implements PhylogenyDataPhyloXmlParser {

    private static final Map<Phylogeny, SequenceRelationParser> _instances = new HashMap<Phylogeny, SequenceRelationParser>();
    private Phylogeny                                           _phylogeny;

    private SequenceRelationParser() {
    }

    @Override
    public SequenceRelation parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        final SequenceRelation seqRelation = new SequenceRelation();
        if ( element.isHasAttribute( PhyloXmlMapping.SEQUENCE_RELATION_TYPE ) ) {
            final String sType = element.getAttribute( PhyloXmlMapping.SEQUENCE_RELATION_TYPE );
            seqRelation.setType( SequenceRelation.SEQUENCE_RELATION_TYPE.valueOf( sType ) );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.SEQUENCE_RELATION_ID_REF0 ) && ( _phylogeny != null ) ) {
            final Sequence ref = PhyloXmlHandler.getSequenceMapByIdForPhylogeny( _phylogeny )
                    .get( element.getAttribute( PhyloXmlMapping.SEQUENCE_RELATION_ID_REF0 ) );
            if ( ref != null ) {
                seqRelation.setRef0( ref );
            }
        }
        if ( element.isHasAttribute( PhyloXmlMapping.SEQUENCE_RELATION_ID_REF1 ) && ( _phylogeny != null ) ) {
            final Sequence ref = PhyloXmlHandler.getSequenceMapByIdForPhylogeny( _phylogeny )
                    .get( element.getAttribute( PhyloXmlMapping.SEQUENCE_RELATION_ID_REF1 ) );
            if ( ref != null ) {
                seqRelation.setRef1( ref );
            }
        }
        if ( element.isHasAttribute( PhyloXmlMapping.SEQUENCE_RELATION_DISTANCE ) ) {
            seqRelation
            .setDistance( Double.valueOf( element.getAttribute( PhyloXmlMapping.SEQUENCE_RELATION_DISTANCE ) ) );
        }
        for( int i = 0; i < element.getNumberOfChildElements(); ++i ) {
            final XmlElement child_element = element.getChildElement( i );
            if ( child_element.getQualifiedName().equals( PhyloXmlMapping.CONFIDENCE ) ) {
                seqRelation.setConfidence( ( Confidence ) ConfidenceParser.getInstance().parse( child_element ) );
            }
        }
        return seqRelation;
    }

    public static PhylogenyDataPhyloXmlParser getInstance( final Phylogeny phylogeny ) {
        SequenceRelationParser instance = _instances.get( phylogeny );
        if ( instance == null ) {
            instance = new SequenceRelationParser();
            instance._phylogeny = phylogeny;
            _instances.put( phylogeny, instance );
        }
        return instance;
    }
}
