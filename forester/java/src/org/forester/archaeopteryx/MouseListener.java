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

import java.awt.Point;
import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;

/*
 * @author Christian Zmasek
 */
final class MouseListener extends MouseAdapter implements MouseMotionListener {

    private final TreePanel _treepanel;
    private boolean         _being_dragged = false;
    private boolean         _dragging_legend = false;
    // Whether the current legend gesture has moved the legend past LEGEND_DRAG_THRESHOLD. Set in mouseDragged,
    // consulted (and cleared) by the trailing mouseClicked so a real drag is not ALSO read as a chip click.
    private boolean         _legend_moved  = false;
    private final Point     _click_point   = new Point();
    private final Point     _legend_press_point = new Point();
    // Pixels the pointer must travel from the press point before a press on the legend counts as a drag rather
    // than a click; below this, a tiny press+jitter on a chip stays a pure click and does not nudge the legend.
    private static final int LEGEND_DRAG_THRESHOLD = 5;
    // A press on an annotation column's header (rectangular) or ring (circular) that may become a drag to move the
    // column. The same click-vs-drag rule as the legend: below the threshold the gesture stays a click (which toggles
    // the column's legend); past it, it is a drag, and _column_moved swallows the trailing click.
    private boolean         _dragging_column = false;
    private boolean         _column_moved    = false;
    private final Point     _column_press_point = new Point();

    /**
     * Constructor.
     */
    MouseListener( final TreePanel tp ) {
        _treepanel = tp;
    }

    /**
     * Mouse clicked.
     */
    @Override
    public void mouseClicked( final MouseEvent e ) {
        if ( _column_moved ) {
            // only the tail of a column drag: the column has already moved, so this must not ALSO toggle its legend
            _column_moved = false;
            return;
        }
        if ( _legend_moved ) {
            // this "click" is only the tail of a real legend drag (press + move + release); the legend was
            // already repositioned, so do not ALSO treat it as a recolor/toggle click. See mouseDragged.
            _legend_moved = false;
            return;
        }
        // Legends are drawn color/rank -> size -> pie (pie on top), and startLegendDrag gives the top-most drag
        // priority in an overlap, so the click dispatch must hit-test in the same reverse order -- else a click on
        // a visibly-top legend would fall through to the one underneath (recolor/reset the wrong legend).
        if ( _treepanel.isOnDomainLegend( e ) ) {
            _treepanel.handleDomainLegendClick( e ); // double-click to reset the domain legend's position
            return; // a click on the legend is not a node action
        }
        if ( _treepanel.isOnInternalTaxaKey( e ) ) {
            _treepanel.handleInternalTaxaKeyClick( e ); // double-click to reset the internal-taxonomy key's position
            return; // a click on the legend is not a node action
        }
        if ( _treepanel.isOnAncestralPieLegend( e ) ) {
            _treepanel.handleAncestralPieLegendClick( e ); // double-click to reset the pie legend's position
            return; // a click on the legend is not a node action
        }
        if ( _treepanel.isOnSizeLegend( e ) ) {
            _treepanel.handleSizeLegendClick( e ); // double-click to reset the size legend's position
            return; // a click on the legend is not a node action
        }
        if ( _treepanel.isOnPropertyLegend( e ) ) {
            _treepanel.handleLegendClick( e ); // recolor a value row, or double-click to reset position
            return;
        }
        if ( _treepanel.handleAnnotationHeaderClick( e ) ) {
            return; // a click on a column header toggles that column's color legend
        }
        _click_point.setLocation( e.getX(), e.getY() );
        _treepanel.mouseClicked( e );
    }

    @Override
    public void mouseDragged( final MouseEvent e ) {
        if ( _dragging_legend ) {
            // sub-threshold jitter: keep treating the gesture as a click-in-progress and leave the legend put,
            // so a click on a chip does not also nudge the legend. Once the pointer travels far enough the
            // gesture becomes a drag (and _legend_moved swallows any trailing click).
            if ( !_legend_moved && ( Math.abs( e.getX() - _legend_press_point.x ) < LEGEND_DRAG_THRESHOLD )
                    && ( Math.abs( e.getY() - _legend_press_point.y ) < LEGEND_DRAG_THRESHOLD ) ) {
                return;
            }
            _legend_moved = true;
            _treepanel.dragLegend( e );
            return;
        }
        if ( _dragging_column ) {
            // below the threshold it is still a click on the header; a press on a header never pans the view
            if ( !_column_moved && ( Math.abs( e.getX() - _column_press_point.x ) < LEGEND_DRAG_THRESHOLD )
                    && ( Math.abs( e.getY() - _column_press_point.y ) < LEGEND_DRAG_THRESHOLD ) ) {
                return;
            }
            _column_moved = true;
            _treepanel.dragAnnotationColumn( e.getX(), e.getY() );
            return;
        }
        if ( ( e.getModifiersEx() == InputEvent.BUTTON1_DOWN_MASK )
                || ( e.getModifiersEx() == InputEvent.BUTTON3_DOWN_MASK ) ) {
            if ( !_treepanel.inOvRectangle( e ) ) {
                if ( !_being_dragged ) {
                    _being_dragged = true;
                    _treepanel.setLastMouseDragPointX( e.getX() );
                    _treepanel.setLastMouseDragPointY( e.getY() );
                }
                _treepanel.mouseDragInBrowserPanel( e );
            }
            else {
                if ( !_being_dragged ) {
                    _being_dragged = true;
                    _treepanel.setLastMouseDragPointX( e.getX() );
                    _treepanel.setLastMouseDragPointY( e.getY() );
                }
                _treepanel.mouseDragInOvRectangle( e );
            }
        }
    }

    @Override
    public void mouseMoved( final MouseEvent e ) {
        if ( _treepanel.isOnAnyLegend( e ) || ( _treepanel.annotationColumnGrabbedAt( e.getX(), e.getY() ) >= 0 ) ) {
            _treepanel.setCursor( TreePanel.MOVE_CURSOR ); // hint that the legend, or the column, can be dragged
            _treepanel.clearHoverPreview(); // don't leave a select/deselect preview on the tree behind the legend
            _treepanel.hideNodeDataPopup(); // ...nor the rollover popup (this early return skips its normal hide)
            return;
        }
        _treepanel.mouseMoved( e );
    }

    @Override
    public void mouseExited( final MouseEvent e ) {
        _treepanel.clearHoverPreview(); // don't leave a hover preview when the pointer leaves the panel
        // The rollover popup can be a heavyweight NATIVE window: leaving the canvas -- onto the control panel,
        // the menus, the scrollbars, or out of the window entirely -- must take it down, or it stays floating
        // on the desktop (no further mouse-move over the canvas will ever hide it).
        _treepanel.hideNodeDataPopup();
    }

    @Override
    public void mousePressed( final MouseEvent e ) {
        // a click is about to act (select / edit dialog / context menu): don't leave the rollover popup
        // floating over whatever comes up
        _treepanel.hideNodeDataPopup();
        _legend_moved = false; // a fresh gesture: no legend drag has happened yet
        _column_moved = false; // ... nor any column drag
        if ( ( e.getModifiersEx() == InputEvent.BUTTON1_DOWN_MASK ) && _treepanel.isOnAnyLegend( e ) ) {
            _dragging_legend = true;
            _legend_press_point.setLocation( e.getX(), e.getY() );
            _treepanel.startLegendDrag( e ); // decides which legend (color/rank/column or size) this drag moves
        }
        else if ( e.getModifiersEx() == InputEvent.BUTTON1_DOWN_MASK ) {
            // a legend drawn over a header wins (legends are on top); otherwise a header / ring can be dragged
            final int col = _treepanel.annotationColumnGrabbedAt( e.getX(), e.getY() );
            if ( col >= 0 ) {
                _dragging_column = true;
                _column_press_point.setLocation( e.getX(), e.getY() );
                _treepanel.startAnnotationColumnDrag( col );
            }
        }
        // A pan/overview drag is NOT started here: it begins on the first mouseDragged (which initializes
        // _being_dragged + the drag anchor), so a plain press stays a click. (A former press-time drag block here
        // set the anchor at the press point and called the pan handler with a zero delta -- no visible scroll --
        // which mouseDragged already does on the first move; removing it just anchors the pan one event later.)
    }

    /** Test hook: whether a browser-panel / overview pan drag is currently in progress. */
    boolean isDraggingForTest() {
        return _being_dragged;
    }

    @Override
    public void mouseReleased( final MouseEvent e ) {
        if ( _dragging_legend ) {
            _dragging_legend = false;
            _treepanel.endLegendDrag();
            return;
        }
        if ( _dragging_column ) {
            _dragging_column = false;
            // drop only a real drag; a sub-threshold press ends as a click, which the trailing mouseClicked handles
            _treepanel.endAnnotationColumnDrag( _column_moved );
            return;
        }
        if ( _being_dragged ) {
            _being_dragged = false;
        }
        _treepanel.mouseReleasedInBrowserPanel( e );
    }
}