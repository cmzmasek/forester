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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;

/*
 * Interface for data for annotating a Phylogeny.
 */
public interface PhylogenyData {

    public StringBuffer asSimpleText();

    public StringBuffer asText();

    /**
     * Creates a new PhylogenyData object with identical values as this
     * PhylogenyData.
     * This ~should~ return a deep copy, but not there yet.
     *
     *
     * @return a ~deep~ copy of this PhylogenyData
     */
    public PhylogenyData copy();

    /**
     * Compares this PhylogenyData to PhylogenyData data. In general, this
     * should return true if and only if all fiels are exactly identical.
     *
     * @param PhylogenyData
     *            the PhylogenyData to compare to
     * @return in general, true if and only if all fiels are exactly identical,
     *         false otherwise
     */
    public boolean isEqual( final PhylogenyData data );

    public StringBuffer toNHX();

    /**
     *  Writes a phyloXML representation of this phylogeny data.
     *
     * @param writer
     * @param level
     * @param indentation
     * @throws IOException
     */
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException;
}