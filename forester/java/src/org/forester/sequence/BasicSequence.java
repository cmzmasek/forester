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

package org.forester.sequence;

import org.forester.util.ForesterUtil;

public class BasicSequence implements MolecularSequence {

    private final char[] _mol_sequence;
    private String       _identifier;
    private final TYPE   _type;

    /**
     * Only use if you know what you are doing!
     *
     */
    public BasicSequence( final String identifier, final String mol_sequence, final TYPE type ) {
        check( identifier, mol_sequence );
        _mol_sequence = mol_sequence.toCharArray();
        _identifier = identifier;
        _type = type;
    }

    private static final void check( final String identifier, final String mol_sequence ) {
        if ( ForesterUtil.isEmpty( identifier ) ) {
            throw new IllegalArgumentException( "identifier of sequence cannot be empty" );
        }
        if ( ForesterUtil.isEmpty( mol_sequence ) ) {
            throw new IllegalArgumentException( "molecular sequence cannot be empty" );
        }
    }

    /**
     * Only use if you know what you are doing!
     *
     */
    public BasicSequence( final String identifier, final char[] mol_sequence, final TYPE type ) {
        if ( ForesterUtil.isEmpty( identifier ) ) {
            throw new IllegalArgumentException( "identifier of sequence cannot be empty" );
        }
        if ( ( mol_sequence == null ) || ( mol_sequence.length < 1 ) ) {
            throw new IllegalArgumentException( "molecular sequence cannot be empty" );
        }
        _mol_sequence = mol_sequence;
        _identifier = identifier;
        _type = type;
    }

    public void setIdentifier( final String id ) {
        _identifier = id;
    }

    @Override
    public String getIdentifier() {
        return _identifier;
    }

    @Override
    public int getLength() {
        return _mol_sequence.length;
    }

    @Override
    public char[] getMolecularSequence() {
        return _mol_sequence;
    }

    @Override
    public char getResidueAt( final int position ) {
        return _mol_sequence[ position ];
    }

    @Override
    public TYPE getType() {
        return _type;
    }

    @Override
    public int getNumberOfGapResidues() {
        int gaps = 0;
        for( final char element : _mol_sequence ) {
            if ( element == GAP ) {
                ++gaps;
            }
        }
        return gaps;
    }

    @Override
    public boolean equals( final Object obj ) {
        if ( obj == null ) {
            return false;
        }
        if ( obj.getClass() != getClass() ) {
            return false;
        }
        final MolecularSequence other = ( MolecularSequence ) obj;
        if ( getMolecularSequenceAsString().equals( other.getMolecularSequenceAsString() ) ) {
            return true;
        }
        return false;
    }

    @Override
    public int hashCode() {
        return getMolecularSequenceAsString().hashCode();
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append( _identifier.toString() );
        sb.append( ": " );
        sb.append( getMolecularSequenceAsString() );
        return sb.toString();
    }

    public static MolecularSequence copySequence( final MolecularSequence seq ) {
        final char[] s = new char[ seq.getMolecularSequence().length ];
        for( int i = 0; i < seq.getMolecularSequence().length; i++ ) {
            s[ i ] = seq.getMolecularSequence()[ i ];
        }
        return new BasicSequence( new String( seq.getIdentifier() ), s, seq.getType() );
    }

    public static MolecularSequence createSequence( final String identifier, final String mol_sequence ) {
        check( identifier, mol_sequence );
        final TYPE type = ForesterUtil.guessMolecularSequenceType( mol_sequence );
        final String re;
        final char repl;
        if ( type == TYPE.AA ) {
            re = AA_REGEXP;
            repl = UNSPECIFIED_AA;
        }
        else if ( type == TYPE.DNA ) {
            re = DNA_REGEXP;
            repl = UNSPECIFIED_NUC;
        }
        else if ( type == TYPE.RNA ) {
            re = RNA_REGEXP;
            repl = UNSPECIFIED_NUC;
        }
        else {
            throw new IllegalArgumentException( "could not determine sequence type for: " + mol_sequence);
        }
        return new BasicSequence( identifier, mol_sequence.toUpperCase().replaceAll( "\\.", GAP_STR )
                                  .replaceAll( re, Character.toString( repl ) ), type );
    }
    
    public static MolecularSequence createGeneralSequence( final String identifier, final String mol_sequence ) {
        check( identifier, mol_sequence );
        return new BasicSequence( identifier, mol_sequence.toUpperCase().replaceAll( "\\.", GAP_STR 
                                  ), TYPE.GENERAL );
    }
    
    public static MolecularSequence createAaSequence( final String identifier, final String mol_sequence ) {
        check( identifier, mol_sequence );
        return new BasicSequence( identifier, mol_sequence.toUpperCase().replaceAll( "\\.", GAP_STR )
                                  .replaceAll( AA_REGEXP, Character.toString( UNSPECIFIED_AA ) ), TYPE.AA );
    }

    public static MolecularSequence createDnaSequence( final String identifier, final String mol_sequence ) {
        check( identifier, mol_sequence );
        return new BasicSequence( identifier, mol_sequence.toUpperCase().replaceAll( "\\.", GAP_STR )
                                  .replaceAll( DNA_REGEXP, Character.toString( UNSPECIFIED_NUC ) ), TYPE.DNA );
    }

    public static MolecularSequence createRnaSequence( final String identifier, final String mol_sequence ) {
        check( identifier, mol_sequence );
        return new BasicSequence( identifier, mol_sequence.toUpperCase().replaceAll( "\\.", GAP_STR )
                                  .replaceAll( RNA_REGEXP, Character.toString( UNSPECIFIED_NUC ) ), TYPE.RNA );
    }

    @Override
    public String getMolecularSequenceAsString() {
        return new String( getMolecularSequence() );
    }

    @Override
    public boolean isGapAt( final int position ) {
        return getResidueAt( position ) == GAP;
    }
}
