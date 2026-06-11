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

package org.forester.ws.seqdb;

import java.util.List;
import java.util.SortedSet;

import org.forester.go.GoTerm;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.sequence.MolecularSequence;

public interface SequenceDatabaseEntry {

    public String getAccession();

    public String getPrimaryAccession();

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

    public MolecularSequence getMolecularSequence();

    public List<String> getTaxonomicLineage();

    public String getStrain();

    public String getViralHost();

    public String getViralCountry();

    public String getViralIsolate();

    public String getCollectionDate();

    public String getViralIsolationSource();

    public String getViralSegment();
}