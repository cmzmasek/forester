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
        final double conf = confidence( element );
        String id = ProteinDomain.IDENTIFIER_DEFAULT;
        try {
            f = Integer
                    .parseInt( element.getAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_FROM ) );
            t = Integer.parseInt( element.getAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_TO ) );
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
        return new ProteinDomain( name, f, t, id, conf );
    }

    /**
     * The domain's confidence (an E-value, say), or {@link ProteinDomain#CONFIDENCE_DEFAULT} when it has none. The
     * attribute is OPTIONAL in phyloXML, and the writer saves a domain without one as the default, -1.0 -- so a
     * missing, unreadable or negative value means "no confidence", never a malformed domain. (Until 0.11.195 it
     * meant the latter: every domain saved without a confidence was dropped when its file was opened again.)
     */
    private static double confidence( final XmlElement element ) {
        if ( !element.isHasAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_CONFIDENCE ) ) {
            return ProteinDomain.CONFIDENCE_DEFAULT;
        }
        try {
            final double c = Double.parseDouble( element
                    .getAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_CONFIDENCE ).trim() );
            return ( Double.isNaN( c ) || ( c < 0 ) ) ? ProteinDomain.CONFIDENCE_DEFAULT : c;
        }
        catch ( final NumberFormatException e ) {
            return ProteinDomain.CONFIDENCE_DEFAULT;
        }
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
