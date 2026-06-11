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

import java.io.IOException;
import java.io.Writer;
import java.util.List;

import org.forester.sequence.MolecularSequence;
import org.forester.sequence.MolecularSequence.TYPE;

public interface Msa {

    public static enum MSA_FORMAT {
        FASTA, PHYLIP, NEXUS;
    }

    public String getIdentifier( int row );

    public void setIdentifier( int row, String identifier );

    public int getLength();

    public int getNumberOfSequences();

    public char getResidueAt( int row, int col );

    public boolean isGapAt( int row, int col );

    public List<Character> getColumnAt( int col );

    public MolecularSequence getSequence( final String id );

    public MolecularSequence getSequence( final int row );

    public List<MolecularSequence> asSequenceList();

    public StringBuffer getSequenceAsString( int row );

    public abstract TYPE getType();

    public void setResidueAt( final int row, final int col, final char residue );

    public void write( Writer w, MSA_FORMAT format ) throws IOException;
}
