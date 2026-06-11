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
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.PhylogenyData;

public class ConfidenceParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new ConfidenceParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private ConfidenceParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        final Confidence confidence = new Confidence();
        confidence.setValue( element.getValueAsDouble() );
        if ( element.isHasAttribute( PhyloXmlMapping.CONFIDENCE_TYPE_ATTR ) ) {
            confidence.setType( element.getAttribute( PhyloXmlMapping.CONFIDENCE_TYPE_ATTR ) );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.CONFIDENCE_SD_ATTR ) ) {
            try {
                confidence.setStandardDeviation( Double.parseDouble( element
                                                                     .getAttribute( PhyloXmlMapping.CONFIDENCE_SD_ATTR ) ) );
            }
            catch ( final NumberFormatException ex ) {
                throw new PhyloXmlDataFormatException( "attempt to parse ["
                        + element.getAttribute( PhyloXmlMapping.CONFIDENCE_SD_ATTR + "] into double" ) );
            }
        }
        return confidence;
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
