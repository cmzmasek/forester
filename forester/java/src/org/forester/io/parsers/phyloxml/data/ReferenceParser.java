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

package org.forester.io.parsers.phyloxml.data;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.Reference;
import org.forester.util.ForesterUtil;

public class ReferenceParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new ReferenceParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private ReferenceParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        String desc = "";
        String doi = "";
        if ( element.isHasAttribute( PhyloXmlMapping.REFERENCE_DOI_ATTR ) ) {
            doi = element.getAttribute( PhyloXmlMapping.REFERENCE_DOI_ATTR );
        }
        for( int i = 0; i < element.getNumberOfChildElements(); ++i ) {
            final XmlElement child_element = element.getChildElement( i );
            if ( child_element.getQualifiedName().equals( PhyloXmlMapping.REFERENCE_DESC ) ) {
                desc = child_element.getValueAsString();
                break;
            }
        }
        if ( !ForesterUtil.isEmpty( doi ) ) {
            return new Reference( desc, doi );
        }
        else {
            return new Reference( desc );
        }
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
