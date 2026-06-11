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

package org.forester.util;

import java.io.BufferedWriter;
import java.io.IOException;

public final class EasyWriter extends BufferedWriter {

    private final static String LINE_SEPARATOR = ForesterUtil.LINE_SEPARATOR;

    public EasyWriter( final BufferedWriter out ) {
        super( out );
    }

    public void println( final String s ) throws IOException {
        write( s );
        write( LINE_SEPARATOR );
    }
    

    public void println() throws IOException {
        write( LINE_SEPARATOR );
    }

    public void print( final String s ) throws IOException {
        write( s );
    }
 
}
