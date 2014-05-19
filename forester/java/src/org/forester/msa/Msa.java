// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.msa;

import java.io.IOException;
import java.io.Writer;
import java.util.List;

import org.forester.sequence.Sequence;
import org.forester.sequence.Sequence.TYPE;

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

    public Sequence getSequence( final String id );

    public Sequence getSequence( final int row );

    public List<Sequence> asSequenceList();

    public StringBuffer getSequenceAsString( int row );

    public abstract TYPE getType();

    public void setResidueAt( final int row, final int col, final char residue );

    public void write( Writer w, MSA_FORMAT format ) throws IOException;
}
