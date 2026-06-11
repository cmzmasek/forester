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

package org.forester.sdi;

final class GSDIsummaryResult {

    private int _speciation_or_duplication_events_sum;
    private int _speciations_sum;
    private int _duplications_sum;

    GSDIsummaryResult() {
        _speciation_or_duplication_events_sum = 0;
        _speciations_sum = 0;
        _duplications_sum = 0;
    }

    final int getDuplicationsSum() {
        return _duplications_sum;
    }

    final int getSpeciationOrDuplicationEventsSum() {
        return _speciation_or_duplication_events_sum;
    }

    final int getSpeciationsSum() {
        return _speciations_sum;
    }

    final void increaseDuplicationsSum() {
        ++_duplications_sum;
    }

    final void increaseSpeciationOrDuplicationEventsSum() {
        ++_speciation_or_duplication_events_sum;
    }

    final void increaseSpeciationsSum() {
        ++_speciations_sum;
    }
}
