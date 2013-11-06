// $Id:
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
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

package org.forester.ws.seqdb;

import java.util.SortedSet;

import org.forester.go.GoTerm;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;

public interface SequenceDatabaseEntry {

    public String getAccession();

    public String getGeneName();

    public SortedSet<GoTerm> getGoTerms();

    public SortedSet<Annotation> getAnnotations();

    public String getProvider();

    public String getSequenceName();

    public String getSequenceSymbol();

    public String getTaxonomyIdentifier();

    public String getTaxonomyScientificName();

    public boolean isEmpty();

    public SortedSet<Accession> getCrossReferences();

    public String getMap();

    public String getChromosome();
}