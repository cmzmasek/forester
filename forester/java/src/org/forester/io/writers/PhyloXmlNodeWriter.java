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

package org.forester.io.writers;

import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.util.ForesterUtil;

public class PhyloXmlNodeWriter {

    public static void toPhyloXml( final Writer w, final PhylogenyNode node, final int level, final String indentation )
            throws IOException {
        String ind = "";
        if ( ( indentation != null ) && ( indentation.length() > 0 ) ) {
            ind = indentation + PhylogenyWriter.PHYLO_XML_INTENDATION_BASE;
        }
        if ( !ForesterUtil.isEmpty( node.getName() ) ) {
            PhylogenyDataUtil.appendElement( w, PhyloXmlMapping.NODE_NAME, node.getName(), indentation );
        }
        if ( node.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
            PhylogenyDataUtil.appendElement( w, PhyloXmlMapping.BRANCH_LENGTH, String.valueOf( ForesterUtil.round( node
                                                                                                                   .getDistanceToParent(), PhyloXmlUtil.ROUNDING_DIGITS_FOR_PHYLOXML_DOUBLE_OUTPUT ) ), indentation );
        }
        if ( node.getBranchData() != null ) {
            node.getBranchData().toPhyloXML( w, level, ind );
        }
        if ( node.getNodeData() != null ) {
            node.getNodeData().toPhyloXML( w, level, ind );
        }
    }
}
