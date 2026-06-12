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
    private final Point     _click_point   = new Point();

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
        if ( _treepanel.isOnPropertyLegend( e ) ) {
            _treepanel.handleLegendClick( e ); // recolor a value row, or double-click to reset position
            return; // a click on the legend is not a node action
        }
        _click_point.setLocation( e.getX(), e.getY() );
        _treepanel.mouseClicked( e );
    }

    @Override
    public void mouseDragged( final MouseEvent e ) {
        if ( _dragging_legend ) {
            _treepanel.dragLegend( e );
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
        if ( _treepanel.isOnPropertyLegend( e ) ) {
            _treepanel.setCursor( TreePanel.MOVE_CURSOR ); // hint that the legend can be dragged
            return;
        }
        _treepanel.mouseMoved( e );
    }

    @Override
    public void mousePressed( final MouseEvent e ) {
        if ( ( e.getModifiersEx() == InputEvent.BUTTON1_DOWN_MASK ) && _treepanel.isOnPropertyLegend( e ) ) {
            _dragging_legend = true;
            _treepanel.startLegendDrag( e );
            return;
        }
        //TODO is this a good idea? It is certainly not NEEDED.
        if ( e.getModifiersEx() == InputEvent.BUTTON1_DOWN_MASK ) {
            if ( !_being_dragged ) {
                _being_dragged = true;
                _treepanel.setLastMouseDragPointX( e.getX() );
                _treepanel.setLastMouseDragPointY( e.getY() );
            }
            if ( !_treepanel.inOvRectangle( e ) ) {
                _treepanel.mouseDragInBrowserPanel( e );
            }
            else {
                _treepanel.mouseDragInOvRectangle( e );
            }
        }
    }

    @Override
    public void mouseReleased( final MouseEvent e ) {
        if ( _dragging_legend ) {
            _dragging_legend = false;
            _treepanel.endLegendDrag();
            return;
        }
        if ( _being_dragged ) {
            _being_dragged = false;
        }
        _treepanel.mouseReleasedInBrowserPanel( e );
    }
}