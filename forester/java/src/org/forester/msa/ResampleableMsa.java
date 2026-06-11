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

package org.forester.msa;

import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;

public final class ResampleableMsa extends BasicMsa {

    private int[] _resampled_column_positions = null;

    public ResampleableMsa( final BasicMsa msa ) {
        super( msa );
    }

    @Override
    final public char getResidueAt( final int row, final int col ) {
        if ( _resampled_column_positions != null ) {
            return super.getResidueAt( row, _resampled_column_positions[ col ] );
        }
        return super.getResidueAt( row, col );
    }

    final public void resample( final int[] resampled_column_positions ) {
        if ( resampled_column_positions.length != getLength() ) {
            throw new IllegalArgumentException( "illegal attempt to use " + resampled_column_positions.length
                                                + " resampled column positions on msa of length " + getLength() );
        }
        _resampled_column_positions = resampled_column_positions;
    }

    @Override
    final public void setResidueAt( final int row, final int col, final char residue ) {
        throw new NoSuchMethodError( "illegal attempt to set residue in resampleable msa" );
    }

    @Override
    public MolecularSequence getSequence( final int row ) {
        return new BasicSequence( getIdentifier( row ), getSequenceAsString( row ).toString(), getType() );
    }
}
