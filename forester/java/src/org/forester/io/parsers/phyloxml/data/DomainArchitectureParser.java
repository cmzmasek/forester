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
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.ProteinDomain;

public class DomainArchitectureParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new DomainArchitectureParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private DomainArchitectureParser() {
    }

    @Override
    public DomainArchitecture parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        final DomainArchitecture architecure = new DomainArchitecture();
        if ( !element.isHasAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_LENGTH ) ) {
            throw new PhyloXmlDataFormatException( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_LENGTH
                                                   + " attribute is required for domain architecture" );
        }
        final String lenght_str = element.getAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_LENGTH );
        try {
            architecure.setTotalLength( Integer.parseInt( lenght_str ) );
        }
        catch ( final NumberFormatException e ) {
            throw new PhyloXmlDataFormatException( "could not extract domain architecture length from [" + lenght_str
                                                   + "]: " + e.getMessage() );
        }
        for( int i = 0; i < element.getNumberOfChildElements(); ++i ) {
            final XmlElement child_element = element.getChildElement( i );
            if ( child_element.getQualifiedName().equals( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_DOMAIN ) ) {
                architecure.addDomain( ( ProteinDomain ) ProteinDomainParser.getInstance().parse( child_element ) );
            }
        }
        return architecure;
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
