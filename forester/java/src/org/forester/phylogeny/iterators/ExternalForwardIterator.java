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
// WWW: www.phylosoft.org/forester

package org.forester.phylogeny.iterators;

import java.util.NoSuchElementException;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

public class ExternalForwardIterator implements PhylogenyNodeIterator {

    private PhylogenyNode       _current_node;
    private final PhylogenyNode _last_ext_node;
    private final PhylogenyNode _first_ext_node;

    /**
     * Constructor for ExternalForwardIterator.
     * 
     * @param tree
     *            the tree on which to iterate over all external nodes.
     */
    public ExternalForwardIterator( final Phylogeny phylogeny ) throws IllegalArgumentException {
        if ( phylogeny.isEmpty() ) {
            throw new IllegalArgumentException( "attempt to use ExternalForwardIterator on an empty phylogeny" );
        }
        PhylogenyNode n = phylogeny.getRoot();
        while ( !n.isExternal() ) {
            n = n.getLastChildNode();
        }
        _last_ext_node = n;
        _first_ext_node = phylogeny.getFirstExternalNode();
        reset();
    }

    private PhylogenyNode getCurrentNode() {
        return _current_node;
    }

    private PhylogenyNode getFirstExtNode() {
        return _first_ext_node;
    }

    private PhylogenyNode getLastExtNode() {
        return _last_ext_node;
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.util.Iterator#hasNext()
     */
    @Override
    public boolean hasNext() {
        return getCurrentNode() != null;
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.util.Iterator#next()
     */
    @Override
    public PhylogenyNode next() throws NoSuchElementException {
        if ( !hasNext() ) {
            throw new NoSuchElementException( "attempt to call \"next()\" on iterator which has no more next elements" );
        }
        final PhylogenyNode n = getCurrentNode();
        if ( n == getLastExtNode() ) {
            setCurrentNode( null );
        }
        else {
            setCurrentNode( n.getNextExternalNode() );
        }
        return n;
    }

    /**
     * Not supported.
     * 
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * DOCUMENT ME!
     */
    @Override
    public void reset() {
        setCurrentNode( getFirstExtNode() );
    }

    private void setCurrentNode( final PhylogenyNode current_node ) {
        _current_node = current_node;
    }
} // end of class ExternalForwardIterator
