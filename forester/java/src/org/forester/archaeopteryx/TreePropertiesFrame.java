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

package org.forester.archaeopteryx;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.Timer;

/**
 * The per-tab "Tree Properties" window (View menu, Cmd-I, tab double-click): a {@link TreePropertiesForm} in the
 * shared {@link EditorFrame} chrome. Non-modal and one per tree panel: opening it again brings the existing window
 * to the front. The panel marks it {@link #markStale() stale} whenever its tree changes (an edit, an undo, a
 * redo); the window then re-reads the tree once, shortly afterwards (edits made in quick succession are
 * coalesced), keeping any unwritten edits of its own.
 */
final class TreePropertiesFrame extends EditorFrame {

    private static final long        serialVersionUID = 1L;
    /** How long after the last change the window waits before re-reading the tree. */
    static final int                 REFRESH_DELAY_MS = 300;
    private final TreePanel          _tree_panel;
    private final TreePropertiesForm _form;
    private final Timer              _refresh;

    TreePropertiesFrame( final TreePanel tp ) {
        this( new TreePropertiesForm( tp ), tp );
    }

    private TreePropertiesFrame( final TreePropertiesForm form, final TreePanel tp ) {
        super( form, titleFor( form ), "This tree has property changes that have not been written to the tree." );
        _tree_panel = tp;
        _form = form;
        _refresh = new Timer( REFRESH_DELAY_MS, e -> rebindNow() );
        _refresh.setRepeats( false );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowActivated( final WindowEvent e ) {
                if ( _refresh.isRunning() ) { // a pending refresh: do it now, the user is looking
                    _refresh.stop();
                    rebindNow();
                }
            }
        } );
        sizeAndPlace( tp, 0, 48, 34 );
        setVisible( true );
    }

    private static String titleFor( final TreePropertiesForm form ) {
        return "Tree Properties: " + form.titleText();
    }

    TreePropertiesForm getForm() {
        return _form;
    }

    /** The tree changed (or was replaced): re-read it soon. Cheap to call often. */
    void markStale() {
        _refresh.restart();
    }

    /** For tests: whether a refresh is pending. */
    boolean isRefreshPendingForTest() {
        return _refresh.isRunning();
    }

    /** Re-reads the tree right away (also what the timer does). */
    void rebindNow() {
        _refresh.stop();
        _form.rebind();
        refreshState();
    }

    @Override
    protected void refreshState() {
        super.refreshState();
        if ( _form != null ) {
            final String t = titleFor( _form );
            if ( !t.equals( getTitle().replaceFirst( "^• ", "" ) ) ) {
                setBaseTitle( t );
            }
        }
    }

    @Override
    protected void onClosed() {
        _refresh.stop();
        if ( _tree_panel != null ) {
            _tree_panel.treePropertiesFrameClosed( this );
        }
    }
}
