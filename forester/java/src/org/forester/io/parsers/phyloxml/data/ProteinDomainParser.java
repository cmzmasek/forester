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
import org.forester.phylogeny.data.ProteinDomain;

public class ProteinDomainParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new ProteinDomainParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private ProteinDomainParser() {
    }

    @Override
    public ProteinDomain parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        String name = "";
        int f = -1;
        int t = -1;
        double conf = ProteinDomain.CONFIDENCE_DEFAULT;
        String id = ProteinDomain.IDENTIFIER_DEFAULT;
        try {
            f = Integer
                    .parseInt( element.getAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_FROM ) );
            t = Integer.parseInt( element.getAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_TO ) );
            conf = Double.parseDouble( element
                                       .getAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_CONFIDENCE ) );
            if ( element.isHasAttribute( PhyloXmlMapping.IDENTIFIER ) ) {
                id = element.getAttribute( PhyloXmlMapping.IDENTIFIER );
            }
        }
        catch ( final Exception e ) {
            throw new PhyloXmlDataFormatException( "failed to parse element [" + element + "]: " + e.getMessage() );
        }
        name = element.getValueAsString();
        if ( f == -1 ) {
            throw new PhyloXmlDataFormatException( "\"" + PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_FROM
                    + "\" attribute not set in: " + element );
        }
        if ( t == -1 ) {
            throw new PhyloXmlDataFormatException( "\"" + PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_TO
                    + "\" attribute not set in: " + element );
        }
        if ( conf < 0 ) {
            throw new PhyloXmlDataFormatException( "\""
                    + PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_CONFIDENCE
                    + "\" attribute either negative or not set in: " + element );
        }
        return new ProteinDomain( name, f, t, id, conf );
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
