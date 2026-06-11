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

package org.forester.datastructures;

import java.util.LinkedList;
import java.util.NoSuchElementException;

/*
 * A simple Queue data structure. Created: 10/23/2005 by Christian M. Zmasek.
 * Last modified: 10/23/2005 by Christian M. Zmasek.
 *
 * @author Christian M. Zmasek
 *
 * @version 1.000
 */
public class Queue {

    // Instance variables
    // ------------------
    private final LinkedList<Object> _data;

    // Constructor
    // -----------
    /**
     * This created a new, empty Queue object.
     */
    public Queue() {
        _data = new LinkedList<Object>();
    }

    /**
     * Removes all elements from this queue.
     */
    public void clear() {
        getData().clear();
    }

    /**
     * Dequeues one element from this queue.
     *
     * @return the dequeued object
     * @throws NoSuchElementException
     *             if this queue is empty
     */
    public Object dequeue() throws NoSuchElementException {
        if ( isEmpty() ) {
            throw new NoSuchElementException( "Attempt to dequeue from an empty Queue." );
        }
        return getData().removeFirst();
    }

    // Public methods
    // --------------
    /**
     * Adds Object element to thisqueue.
     *
     * @param element
     *            the Object to be enqueued
     */
    public void enqueue( final Object element ) {
        getData().add( element );
    }

    // Private methods
    // ---------------
    /**
     * Returns the LinkedList upon which this queue is based.
     *
     * @return the LinkedList upon which this queue is based
     */
    private LinkedList<Object> getData() {
        return _data;
    }

    /**
     * Returns whether or not this queue is empty.
     *
     * @return true if this queue is empty, false otherwise
     */
    public boolean isEmpty() {
        return getData().isEmpty();
    }
} // end of class Queue.
