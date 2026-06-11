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

package org.forester.pccx;

public class ExternalNodeBasedCoverageMethodOptions implements CoverageCalculationOptions {

    final private String _scoring_method;

    /**
     * This constructor sets the class name for the scoring method e.g.
     * "org.forester.tools.modeling.BranchCountingBasedScoringMethod"
     *
     * @param scoring_method
     *            class name for the scoring method
     */
    public ExternalNodeBasedCoverageMethodOptions( final String scoring_method ) {
        _scoring_method = scoring_method;
    }

    @Override
    public String asString() {
        final StringBuffer sb = new StringBuffer();
        sb.append( "scoring method: " );
        BranchCountingBasedScoringMethod scoring_method;
        try {
            scoring_method = ( BranchCountingBasedScoringMethod ) ( Class.forName( getScoringMethod() ) ).getDeclaredConstructor().newInstance();
        }
        catch ( final Exception e ) {
            sb.append( "?" );
            return sb.toString();
        }
        sb.append( scoring_method.getDesciption() );
        return sb.toString();
    }

    public String getScoringMethod() {
        return _scoring_method;
    }
}
