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

package org.forester.surfacing;

final class MappingResults {

    private String _description;
    private int    _sum_of_failures;
    private int    _sum_of_successes;

    public String getDescription() {
        return _description;
    }

    public int getSumOfFailures() {
        return _sum_of_failures;
    }

    public int getSumOfSuccesses() {
        return _sum_of_successes;
    }

    public void setDescription( final String description ) {
        _description = description;
    }

    public void setSumOfFailures( final int sum_of_failures ) {
        _sum_of_failures = sum_of_failures;
    }

    public void setSumOfSuccesses( final int sum_of_successes ) {
        _sum_of_successes = sum_of_successes;
    }
}
