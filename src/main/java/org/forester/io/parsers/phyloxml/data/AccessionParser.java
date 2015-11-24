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
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.PhylogenyData;

public class AccessionParser implements PhylogenyDataPhyloXmlParser {

    private final static PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new AccessionParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private AccessionParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        if ( element.isHasAttribute( PhyloXmlMapping.ACCESSION_SOURCE_ATTR )
                && element.isHasAttribute( PhyloXmlMapping.ACCESSION_COMMENT_ATTR ) ) {
            return new Accession( element.getValueAsString(),
                                  element.getAttribute( PhyloXmlMapping.ACCESSION_SOURCE_ATTR ),
                                  element.getAttribute( PhyloXmlMapping.ACCESSION_COMMENT_ATTR ) );
        }
        else if ( element.isHasAttribute( PhyloXmlMapping.ACCESSION_SOURCE_ATTR ) ) {
            return new Accession( element.getValueAsString(),
                                  element.getAttribute( PhyloXmlMapping.ACCESSION_SOURCE_ATTR ) );
        }
        else if ( element.isHasAttribute( PhyloXmlMapping.ACCESSION_COMMENT_ATTR ) ) {
            return new Accession( element.getValueAsString(),
                                  "?",
                                  element.getAttribute( PhyloXmlMapping.ACCESSION_COMMENT_ATTR ) );
        }
        else {
            return new Accession( element.getValueAsString(), "?" );
        }
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
