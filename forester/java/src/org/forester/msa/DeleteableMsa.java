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

import java.util.HashMap;

import org.forester.util.ForesterUtil;

public final class DeleteableMsa extends BasicMsa {

    private int                      _length                 = 0;
    private int                      _mapped_col_positions[] = null;
    private int                      _mapped_row_positions[] = null;
    private int                      _seqs                   = 0;
    private HashMap<String, Integer> _seq_id_to_row_map      = null;

    public DeleteableMsa( final BasicMsa msa ) {
        super( msa );
        _mapped_col_positions = new int[ msa.getLength() ];
        _mapped_row_positions = new int[ msa.getNumberOfSequences() ];
        for( int i = 0; i < _mapped_col_positions.length; ++i ) {
            _mapped_col_positions[ i ] = i;
        }
        for( int i = 0; i < _mapped_row_positions.length; ++i ) {
            _mapped_row_positions[ i ] = i;
        }
        _seq_id_to_row_map = new HashMap<String, Integer>();
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            _seq_id_to_row_map.put( msa.getIdentifier( row ), row );
        }
        _length = msa.getLength();
        _seqs = msa.getNumberOfSequences();
    }
    
    
    @Override
    public char[] getSequenceAsArray( final int row ) {
        return  super.getSequenceAsArray( _mapped_row_positions[ row ] );
    }

    public void deleteColumn( final int col ) {
        if ( col >= _length || col < 0 ) {
            throw new IllegalArgumentException( "column " + col + " is out of range" );
        }
        for( int c = col; c < _length - 1; ++c ) {
            _mapped_col_positions[ c ] = _mapped_col_positions[ c + 1 ];
        }
        --_length;
    }

    
  
    private void deleteRow( final int row ) {
        if ( row >= _seqs || row < 0 ) {
            throw new IllegalArgumentException( "row " + row + " is out of range" );
        }
        for( int r = row; r < _seqs - 1; ++r ) {
            _mapped_row_positions[ r ] = _mapped_row_positions[ r + 1 ];
        }
        --_seqs;
    }

    public void deleteRow( final String id ) {
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
        deleteRow( row );
    }

    @Override
    public String getIdentifier( final int row ) {
        return super.getIdentifier( _mapped_row_positions[ row ] );
    }

    @Override
    public int getLength() {
        return _length;
    }

    @Override
    public int getNumberOfSequences() {
        return _seqs;
    }

    
    @Override
    public char getResidueAt( final int row, final int col ) {
         return super.getResidueAt( _mapped_row_positions[ row ], _mapped_col_positions[ col ] );
    }

    @Override
    public void setIdentifier( final int row, final String id ) {
        super.setIdentifier( _mapped_row_positions[ row ], id );
    }

    @Override
    public void setResidueAt( final int row, final int col, final char residue ) {
        super.setResidueAt( _mapped_row_positions[ row ], _mapped_col_positions[ col ], residue );
    }
}
