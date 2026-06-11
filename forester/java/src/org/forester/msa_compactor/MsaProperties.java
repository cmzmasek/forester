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

package org.forester.msa_compactor;

import org.forester.msa.Msa;
import org.forester.msa.MsaMethods;

public final class MsaProperties {

    final private double _entropy21;
    final private double _entropy7;
    final private double _gap_ratio;
    final private int    _length;
    final private int    _number_of_sequences;
    final private double _avg_number_of_gaps;
    final private String _removed_seq;

    public MsaProperties( final int number_of_sequences,
                          final int length,
                          final double gap_ratio,
                          final double entropy7,
                          final double entropy21,
                          final double avg_number_of_gaps,
                          final String removed_seq ) {
        _number_of_sequences = number_of_sequences;
        _length = length;
        _gap_ratio = gap_ratio;
        _entropy7 = entropy7;
        _entropy21 = entropy21;
        _avg_number_of_gaps = avg_number_of_gaps;
        _removed_seq = removed_seq;
    }

    public MsaProperties( final Msa msa, final String removed_seq, final boolean calculate_normalized_shannon_entropy ) {
        _number_of_sequences = msa.getNumberOfSequences();
        _length = msa.getLength();
        _gap_ratio = MsaMethods.calcGapRatio( msa );
        _removed_seq = removed_seq;
        _avg_number_of_gaps = MsaMethods.calcNumberOfGapsStats( msa ).arithmeticMean();
        if ( calculate_normalized_shannon_entropy ) {
            _entropy7 = MsaMethods.calcNormalizedShannonsEntropy( 7, msa );
            _entropy21 = MsaMethods.calcNormalizedShannonsEntropy( 21, msa );
        }
        else {
            _entropy7 = -1;
            _entropy21 = -1;
        }
    }

    public final double getEntropy21() {
        return _entropy21;
    }

    public final double getEntropy7() {
        return _entropy7;
    }

    public final double getGapRatio() {
        return _gap_ratio;
    }

    public final double getAvgNumberOfGaps() {
        return _avg_number_of_gaps;
    }
    
    public final int getLength() {
        return _length;
    }

    public final int getNumberOfSequences() {
        return _number_of_sequences;
    }

    public final String getRemovedSeq() {
        return _removed_seq;
    }
}
