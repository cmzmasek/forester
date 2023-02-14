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
//
//
// "java -Xmx1024m -cp path\to\forester.jar org.forester.application.fasta_split
//
//

package org.forester.application;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.forester.evoinference.distance.PairwiseDistanceCalculator;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.msa.DeleteableMsa;
import org.forester.util.ForesterUtil;

public final class msa2pwd {

    public static void main( final String args[] ) {
        try {
            final File infile = new File( args[ 0 ] );
            final File outfile = new File( args[ 1 ] );
            DeleteableMsa msa = null;
            final FileInputStream is = new FileInputStream( infile );
            if ( FastaParser.isLikelyFasta( infile ) ) {
                msa = DeleteableMsa.createInstance( FastaParser.parseMsa( is ) );
            }
            else {
                msa = DeleteableMsa.createInstance( GeneralMsaParser.parseMsa( is ) );
            }
            final BasicSymmetricalDistanceMatrix pwd = PairwiseDistanceCalculator.calcFractionalSimilarities( msa );
           
            BufferedWriter w = ForesterUtil.createBufferedWriter( outfile ); 
            pwd.write( w);
            w.close();
        }
        catch ( final FileNotFoundException e ) {
            e.printStackTrace();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }
}
