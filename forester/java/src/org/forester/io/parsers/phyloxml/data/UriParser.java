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

import java.net.URI;
import java.net.URISyntaxException;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.Uri;

public class UriParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new UriParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private UriParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        String type = "";
        String desc = "";
        URI uri = null;
        try {
            uri = new URI( element.getValueAsString() );
        }
        catch ( final URISyntaxException e ) {
            throw new PhyloXmlDataFormatException( "ill formatted Uri: " + element.getValueAsString() );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.URI_DESC_ATTR ) ) {
            desc = element.getAttribute( PhyloXmlMapping.URI_DESC_ATTR );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.TYPE_ATTR ) ) {
            type = element.getAttribute( PhyloXmlMapping.TYPE_ATTR );
        }
        return new Uri( uri, desc, type );
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
