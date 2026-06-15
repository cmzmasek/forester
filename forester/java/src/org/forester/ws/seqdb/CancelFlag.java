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

package org.forester.ws.seqdb;

/**
 * A cooperative cancellation signal a long-running resolution polls between items so the user can stop
 * it promptly (the GUI's background process sets it; see {@code SequenceTaxonomyResolver}).
 */
public interface CancelFlag {

    /** A flag that is never cancelled (for non-cancellable callers and tests). */
    CancelFlag NEVER = new CancelFlag() {

        @Override
        public boolean isCancelled() {
            return false;
        }
    };

    boolean isCancelled();
}
