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

package org.forester.phylogeny.factories;

import java.io.IOException;

import org.forester.phylogeny.Phylogeny;

/*
 * Interface for Phylogeny factories.
 *
 * @author Christian M. Zmasek
 */
public interface PhylogenyFactory {

    /**
     * This must create a Phylogeny from source (e.g. an XML file, an alignment,
     * pairwise distances) by using creator (e.g. an XML file parser, an
     * algorithm implementation).
     *
     * @param source
     *            a source to create a Phylogeny from
     * @param creator
     *            a means to create a Phylogeny
     * @return a Phylogeny[] based on argument source
     * @throws IOException
     */
    public Phylogeny[] create( Object source, Object creator ) throws IOException;
}
