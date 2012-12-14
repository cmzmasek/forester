
package org.forester.sdi;

import java.util.List;
import java.util.Set;
import java.util.SortedSet;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.sdi.SDIutil.TaxonomyComparisonBase;

public interface GSDII {

    public abstract int getSpeciationsSum();

    public abstract Set<PhylogenyNode> getMappedExternalSpeciesTreeNodes();

    public abstract SortedSet<String> getReMappedScientificNamesFromGeneTree();

    public abstract List<PhylogenyNode> getStrippedExternalGeneTreeNodes();

    public abstract List<PhylogenyNode> getStrippedSpeciesTreeNodes();

    public abstract TaxonomyComparisonBase getTaxCompBase();
}