// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

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
