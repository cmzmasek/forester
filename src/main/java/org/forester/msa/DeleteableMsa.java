// / $Id:
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

package org.forester.msa;

import java.util.List;

import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;

public final class DeleteableMsa extends BasicMsa {

    private int _length                 = 0;
    private int _mapped_col_positions[] = null;
    private int _mapped_row_positions[] = null;
    private int _seqs                   = 0;

    private DeleteableMsa( final BasicMsa msa ) {
        super( msa );
        _mapped_col_positions = new int[ msa.getLength() ];
        _mapped_row_positions = new int[ msa.getNumberOfSequences() ];
        for( int i = 0; i < _mapped_col_positions.length; ++i ) {
            _mapped_col_positions[ i ] = i;
        }
        for( int i = 0; i < _mapped_row_positions.length; ++i ) {
            _mapped_row_positions[ i ] = i;
        }
        _length = msa.getLength();
        _seqs = msa.getNumberOfSequences();
    }

    public final double[] calcGappiness() {
        final int length = getLength();
        final double gappiness[] = new double[ length ];
        final int seqs = getNumberOfSequences();
        for( int row = 0; row < seqs; ++row ) {
            for( int col = 0; col < length; ++col ) {
            }
        }
        return gappiness;
    }

    public static int calcGapSumPerColumn( final Msa msa, final int col ) {
        int gap_rows = 0;
        for( int j = 0; j < msa.getNumberOfSequences(); ++j ) {
            if ( msa.isGapAt( j, col ) ) {
                gap_rows++;
            }
        }
        return gap_rows;
    }

    public short determineMaxIdLength() {
        short max = 0;
        for( int row = 0; row < getNumberOfSequences(); ++row ) {
            final short l = ( short ) getIdentifier( row ).length();
            if ( l > max ) {
                max = l;
            }
        }
        return max;
    }

    final public void deleteGapColumns( final double max_allowed_gap_ratio ) {
        if ( ( max_allowed_gap_ratio < 0 ) || ( max_allowed_gap_ratio > 1 ) ) {
            throw new IllegalArgumentException( "max allowed gap ration is out of range: " + max_allowed_gap_ratio );
        }
        for( int col = getLength() - 1; col >= 0; --col ) {
            final boolean delete = ( ( double ) MsaMethods.calcGapSumPerColumn( this, col ) / getNumberOfSequences() ) > max_allowed_gap_ratio;
            if ( delete ) {
                deleteColumn( col );
            }
        }
    }

    final public void deleteGapOnlyColumns() {
        for( int col = getLength() - 1; col >= 0; --col ) {
            if ( isAllGap( col ) ) {
                deleteColumn( col );
            }
        }
    }

    final public MolecularSequence deleteRow( final String id, final boolean return_removed_seq ) {
        int row = -1;
        for( int r = 0; r < getNumberOfSequences(); ++r ) {
            if ( getIdentifier( r ).equals( id ) ) {
                row = r;
                break;
            }
        }
        if ( row < 0 ) {
            throw new IllegalArgumentException( "id [" + id + "] not found" );
        }
        MolecularSequence s = null;
        StringBuilder sb = null;
        if ( return_removed_seq ) {
            s = getSequence( row );
            final char[] x = s.getMolecularSequence();
            sb = new StringBuilder( x.length );
            for( final char element : x ) {
                if ( element != MolecularSequence.GAP ) {
                    sb.append( element );
                }
            }
        }
        deleteRow( row );
        if ( return_removed_seq ) {
            return new BasicSequence( new String( s.getIdentifier() ), sb.toString(), s.getType() );
        }
        else {
            return null;
        }
    }

    @Override
    final public String getIdentifier( final int row ) {
        checkRow( row );
        return super.getIdentifier( _mapped_row_positions[ row ] );
    }

    @Override
    final public int getLength() {
        return _length;
    }

    @Override
    final public int getNumberOfSequences() {
        return _seqs;
    }

    @Override
    final public char getResidueAt( final int row, final int col ) {
        checkRow( row );
        checkColumn( col );
        return super.getResidueAt( _mapped_row_positions[ row ], _mapped_col_positions[ col ] );
    }

    @Override
    public MolecularSequence getSequence( final int row ) {
        checkRow( row );
        return new BasicSequence( getIdentifier( row ), getSequenceAsString( row ).toString(), getType() );
    }

    final public boolean isAllGap( final int col ) {
        final int m_col = _mapped_col_positions[ col ];
        for( int j = 0; j < getNumberOfSequences(); ++j ) {
            if ( super.getResidueAt( _mapped_row_positions[ j ], m_col ) != MolecularSequence.GAP ) {
                return false;
            }
        }
        return true;
    }

    @Override
    final public void setIdentifier( final int row, final String id ) {
        checkRow( row );
        super.setIdentifier( _mapped_row_positions[ row ], id );
    }

    @Override
    final public void setResidueAt( final int row, final int col, final char residue ) {
        checkRow( row );
        checkColumn( col );
        super.setResidueAt( _mapped_row_positions[ row ], _mapped_col_positions[ col ], residue );
    }

    final private void checkColumn( final int col ) {
        if ( ( col >= _length ) || ( col < 0 ) ) {
            throw new IllegalArgumentException( "column " + col + " is out of range" );
        }
    }

    final private void checkRow( final int row ) {
        if ( ( row >= _seqs ) || ( row < 0 ) ) {
            throw new IllegalArgumentException( "row " + row + " is out of range" );
        }
    }

    final private void deleteColumn( final int col ) {
        checkColumn( col );
        for( int c = col; c < ( _length - 1 ); ++c ) {
            _mapped_col_positions[ c ] = _mapped_col_positions[ c + 1 ];
        }
        --_length;
    }

    final private void deleteRow( final int row ) {
        checkRow( row );
        for( int r = row; r < ( _seqs - 1 ); ++r ) {
            _mapped_row_positions[ r ] = _mapped_row_positions[ r + 1 ];
        }
        --_seqs;
    }

    public final static DeleteableMsa createInstance( final List<MolecularSequence> seqs ) {
        return new DeleteableMsa( ( BasicMsa ) BasicMsa.createInstance( seqs ) );
    }

    public final static DeleteableMsa createInstance( final Msa msa ) {
        return new DeleteableMsa( ( BasicMsa ) msa );
    }
}
