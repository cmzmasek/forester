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

package org.forester.archaeopteryx;

import java.util.ArrayList;
import java.util.List;

import org.forester.analysis.InternalTaxonomyAudit;
import org.forester.analysis.InternalTaxonomyAudit.Finding;

/** The findings of an {@link InternalTaxonomyAudit} as a {@link TableReport}: one row per internal annotation --
 *  the node, the annotation as it is, the verdict, the evidence, and the proposed correction. */
final class InternalTaxonomyAuditReport implements TableReport {

    static final String[]       COLUMNS = { "Node", "Annotation", "Verdict", "Evidence", "Proposed" };
    private final List<Finding> _findings;

    InternalTaxonomyAuditReport( final List<Finding> findings ) {
        _findings = findings;
    }

    List<Finding> getFindings() {
        return _findings;
    }

    @Override
    public String getSummary() {
        return InternalTaxonomyAudit.summary( _findings );
    }

    @Override
    public String[] getColumns() {
        return COLUMNS.clone();
    }

    @Override
    public List<String[]> getCells() {
        final List<String[]> out = new ArrayList<String[]>();
        for( final Finding f : _findings ) {
            out.add( new String[] { CladeAnnotationReport.nodeLabel( f.getNode() ), f.getCurrent(),
                    f.getVerdict().getLabel(), ( f.getDetail() == null ) ? "" : f.getDetail(), f.getProposed() } );
        }
        return out;
    }

    @Override
    public String toTsv() {
        return TableReport.tsvOf( COLUMNS, getCells() );
    }

    /** How many findings carry a proposed correction. */
    int proposalCount() {
        int n = 0;
        for( final Finding f : _findings ) {
            if ( f.hasProposal() ) {
                n++;
            }
        }
        return n;
    }
}
