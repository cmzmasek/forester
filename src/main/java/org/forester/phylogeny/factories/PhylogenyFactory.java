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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

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
