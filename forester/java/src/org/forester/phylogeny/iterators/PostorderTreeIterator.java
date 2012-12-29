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

package org.forester.phylogeny.iterators;

import java.util.NoSuchElementException;
import java.util.Stack;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/*
 * *
 */
public final class PostorderTreeIterator implements PhylogenyNodeIterator {

    final private Phylogeny                   _tree;
    final private PhylogenyNode               _root;
    private boolean                           _has_next;
    final private Stack<PostOrderStackObject> _stack;

    /**
     * @param t
     *            Phylogeny for which a Iterator is to be constructed.
     */
    public PostorderTreeIterator( final Phylogeny tree ) throws IllegalArgumentException {
        if ( tree.isEmpty() ) {
            throw new IllegalArgumentException( "Attempt to use PostorderTreeIterator on an empty phylogeny." );
        }
        _tree = tree;
        _root = getTree().getRoot();
        _stack = new Stack<PostOrderStackObject>();
        reset();
    }

    final private PhylogenyNode getRoot() {
        return _root;
    }

    final private Stack<PostOrderStackObject> getStack() {
        return _stack;
    }

    final private Phylogeny getTree() {
        return _tree;
    }

    @Override
    final public boolean hasNext() {
        return _has_next;
    }

    /**
     * Advances the Iterator by one.
     */
    @Override
    final public PhylogenyNode next() throws NoSuchElementException {
        if ( !hasNext() ) {
            throw new NoSuchElementException( "Attempt to call \"next()\" on iterator which has no more next elements." );
        }
        while ( true ) {
            final PostOrderStackObject si = getStack().pop();
            final PhylogenyNode node = si.getNode();
            final int phase = si.getPhase();
            if ( phase > node.getNumberOfDescendants() ) {
                setHasNext( node != getRoot() );
                return node;
            }
            else {
                getStack().push( new PostOrderStackObject( node, ( phase + 1 ) ) );
                if ( node.isInternal() ) {
                    getStack().push( new PostOrderStackObject( node.getChildNode( phase - 1 ), 1 ) );
                }
            }
        }
    }

    @Override
    final public void remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    final public void reset() {
        setHasNext( true );
        getStack().clear();
        getStack().push( new PostOrderStackObject( getTree().getRoot(), 1 ) );
    }

    final private void setHasNext( final boolean has_next ) {
        _has_next = has_next;
    }
}
