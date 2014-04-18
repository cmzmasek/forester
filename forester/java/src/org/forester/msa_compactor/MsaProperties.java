// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2014 Christian M. Zmasek
// Copyright (C) 2014 Sanford-Burnham Medical Research Institute
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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.msa_compactor;

import org.forester.msa.Msa;
import org.forester.msa.MsaMethods;

public final class MsaProperties {

    final private double _average_identity_ratio;
    final private double _gap_ratio;
    final private int    _length;
    final private int    _number_of_sequences;

    public MsaProperties( final Msa msa ) {
        _number_of_sequences = msa.getNumberOfSequences();
        _length = msa.getLength();
        _gap_ratio = MsaMethods.calcGapRatio( msa );
        _average_identity_ratio = MsaMethods.calculateIdentityRatio( 0, msa.getLength() - 1, msa ).arithmeticMean();
    }

    public MsaProperties( final int number_of_sequences,
                          final int length,
                          final double gap_ratio,
                          final double average_identity_ratio ) {
        _number_of_sequences = number_of_sequences;
        _length = length;
        _gap_ratio = gap_ratio;
        _average_identity_ratio = average_identity_ratio;
    }

    public final double getAverageIdentityRatio() {
        return _average_identity_ratio;
    }

    public final double getGapRatio() {
        return _gap_ratio;
    }

    public final int getLength() {
        return _length;
    }

    public final int getNumberOfSequences() {
        return _number_of_sequences;
    }
}
