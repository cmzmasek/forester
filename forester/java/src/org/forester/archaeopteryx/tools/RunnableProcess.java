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

package org.forester.archaeopteryx.tools;

import org.forester.archaeopteryx.AptxConstants;
import org.forester.archaeopteryx.MainFrame;
import org.forester.util.ForesterUtil;
import org.forester.ws.seqdb.CancelFlag;

/**
 * A background task registered with the {@link ProcessPool} (so it shows in the process menu while it
 * runs). It is a {@link CancelFlag}: a long-running task should poll {@link #isCancelled()} and stop
 * promptly when the user requests cancellation (the process menu calls {@link #requestCancel()}).
 */
public abstract class RunnableProcess implements Runnable, CancelFlag {

    long             _process_id;
    private volatile boolean _cancelled = false;

    @Override
    public boolean isCancelled() {
        return _cancelled;
    }

    /** Requests cooperative cancellation; a running task polling {@link #isCancelled()} stops at its next check. */
    public void requestCancel() {
        _cancelled = true;
    }

    long getProcessId() {
        return _process_id;
    }

    void setProcessId( final long process_id ) {
        _process_id = process_id;
    }

    public void start( final MainFrame mf, final String name ) {
        setProcessId( mf.getProcessPool().addProcess( name, this ) );
        mf.updateProcessMenu();
    }

    public void end( final MainFrame mf ) {
        final boolean removed = mf.getProcessPool().removeProcess( getProcessId() );
        if ( !removed ) {
            ForesterUtil.printWarningMessage( AptxConstants.PRG_NAME, "could not remove process " + getProcessId()
                                              + " from process pool" );
        }
        mf.updateProcessMenu();
    }
}
