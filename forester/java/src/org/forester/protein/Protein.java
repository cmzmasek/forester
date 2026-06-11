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

package org.forester.protein;

import java.util.List;

import org.forester.species.Species;

public interface Protein {

    public void addProteinDomain( final Domain protein_domain );

    /**
     * If in_nc_order is set to true, this should return true only and only if
     * the order in List 'domains' and this protein (as determined by the start positions
     * of the domains of this proteins, _not_ by their index) are the same
     * (interspersing, 'other', domains in this are ignored).
     * If in_nc_order is set to false, this should return true only and only if
     * this contains all domains listed in 'domains' (order and count do not matter).
     *
     * @param domains a list of domain ids in a certain order.
     * @param in_nc_order to consider order
     * @return
     */
    public boolean contains( final List<String> domains, final boolean in_nc_order );

    public String getAccession();

    public String getDescription();

    public String getName();

    public int getNumberOfProteinDomains();

    public Domain getProteinDomain( final int index );

    public int getProteinDomainCount( final String domain_id );

    public List<Domain> getProteinDomains();

    public List<Domain> getProteinDomains( final String domain_id );

    public ProteinId getProteinId();

    public int getLength();

    public Species getSpecies();

    public List<Domain> getDomainsSortedByPosition();

    public String toDomainArchitectureString( final String separator, final double ie_cutoff );
    
}