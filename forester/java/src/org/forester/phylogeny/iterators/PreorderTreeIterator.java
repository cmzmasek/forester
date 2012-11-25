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
import java.util.Stack;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

public final class PreorderTreeIterator implements PhylogenyNodeIterator {

    final private Phylogeny            _tree;
    final private Stack<PhylogenyNode> _stack;

    /**
     * @param tree
     *            Phylogeny for which a Iterator is to be constructed.
     */
    public PreorderTreeIterator( final Phylogeny tree ) throws IllegalArgumentException {
        if ( tree.isEmpty() ) {
            throw new IllegalArgumentException( "Attempt to use PreorderTreeIterator on empty tree." );
        }
        _stack = new Stack<PhylogenyNode>();
        _tree = tree;
        reset();
    }

    public PreorderTreeIterator( final PhylogenyNode node ) throws IllegalArgumentException {
        _stack = new Stack<PhylogenyNode>();
        _tree = null;
        reset( node );
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.util.Iterator#hasNext()
     */
    @Override
    public final boolean hasNext() {
        return !_stack.isEmpty();
    }

    /**
     * Advances the Iterator by one.
     */
    @Override
    public final PhylogenyNode next() throws NoSuchElementException {
        if ( !hasNext() ) {
            throw new NoSuchElementException( "Attempt to call \"next()\" on iterator which has no more next elements." );
        }
        final PhylogenyNode node = _stack.pop();
        if ( !node.isExternal() ) {
            for( int i = node.getNumberOfDescendants() - 1; i >= 0; --i ) {
                _stack.push( node.getChildNode( i ) );
            }
        }
        return node;
    }

    /**
     * Not supported.
     * 
     */
    @Override
    public final void remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    public final void reset() {
        _stack.clear();
        _stack.push( _tree.getRoot() );
    }

    private final void reset( final PhylogenyNode node ) {
        _stack.clear();
        _stack.push( node );
    }
}
