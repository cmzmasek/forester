
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
