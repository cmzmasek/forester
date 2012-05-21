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

import java.util.ArrayList;
import java.util.List;

public class ProcessPool {

    final static private ArrayList<ProcessRunning> _processes = new ArrayList<ProcessRunning>();

    public synchronized ProcessRunning getProcessByIndex( final int i ) {
        return getProcesses().get( i );
    }

    public synchronized int size() {
        return getProcesses().size();
    }

    public synchronized ProcessRunning getProcessById( final int id ) {
        for( final ProcessRunning p : getProcesses() ) {
            if ( p.getId() == id ) {
                return p;
            }
        }
        return null;
    }

    public synchronized int addProcess( final String name ) {
        final ProcessRunning p = ProcessRunning.createInstance( name );
        final int id = p.getId();
        if ( getProcessById( id ) != null ) {
            throw new IllegalStateException( "process with id " + id + "already exists" );
        }
        getProcesses().add( p );
        return id;
    }

    public synchronized boolean removeProcess( final int id ) {
        final int i = getProcessIndexById( id );
        if ( i >= 0 ) {
            getProcesses().remove( i );
            return true;
        }
        return false;
    }

    private synchronized int getProcessIndexById( final int id ) {
        final ProcessRunning p = getProcessById( id );
        if ( p != null ) {
            return p.getId();
        }
        return -1;
    }

    private synchronized List<ProcessRunning> getProcesses() {
        return _processes;
    }
}
