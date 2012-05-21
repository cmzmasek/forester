// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2012 Christian M. Zmasek
// Copyright (C) 2008-2012 Sanford Burnham Medical Research Institute
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
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx.tools;

import java.text.SimpleDateFormat;
import java.util.Calendar;

final public class ProcessRunning {

    private static int   count = 0;
    final private int    _id;
    final private String _name;
    final private String _start;

    public  int getId() {
        return _id;
    }
    
    public  String getName() {
        return _name;
    }
    
    public String getStart() {
        return _start;
    }
    
    public String toString() {
        return getName() + " [id=" + getId() + "] [start=" + getStart() + "]";
    }
    
    synchronized static ProcessRunning createInstance( final String name ) {
        final Calendar cal = Calendar.getInstance();
        final SimpleDateFormat sdf = new SimpleDateFormat( "HH:mm:ss" );
        return new ProcessRunning( count++, name, sdf.format( cal.getTime() ) );
    }
    
    private ProcessRunning( final int id, final String name, final String start ) {
        _id = id;
        _name = name;
        _start = start;  
    }

}
