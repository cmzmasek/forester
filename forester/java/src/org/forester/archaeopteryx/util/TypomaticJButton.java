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

package org.forester.archaeopteryx.util;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.Action;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.Timer;

public final class TypomaticJButton extends JButton implements ActionListener, MouseListener {

    private static final long serialVersionUID = 7435606927739361980L;
    private boolean           pressed          = false;
    private boolean           repeat_enabled   = true;
    private Timer             timer            = null;
    private int               initial_delay    = 600;
    private int               delay            = 200;
    private int               modifiers        = 0;

    public TypomaticJButton() {
        super();
        init();
    }

    public TypomaticJButton( final Action a ) {
        super( a );
        init();
    }

    public TypomaticJButton( final Icon icon ) {
        super( icon );
        init();
    }

    public TypomaticJButton( final String text ) {
        super( text );
        init();
    }

    public TypomaticJButton( final String text, final Icon icon ) {
        super( text, icon );
        init();
    }

    @Override
    final public void actionPerformed( final ActionEvent ae ) {
        if ( ae.getSource() == timer ) {
            final ActionEvent event = new ActionEvent( this,
                                                       ActionEvent.ACTION_PERFORMED,
                                                       super.getActionCommand(),
                                                       modifiers );
            super.fireActionPerformed( event );
        }
    }

    final public int getDelay() {
        return delay;
    }

    final public int getInitialDelay() {
        return initial_delay;
    }

    final private void init() {
        addMouseListener( this );
        timer = new Timer( delay, this );
        timer.setRepeats( true );
    }

    final public boolean isRepeatEnabled() {
        return repeat_enabled;
    }

    @Override
    final public void mouseClicked( final MouseEvent me ) {
        if ( me.getSource() == this ) {
            pressed = false;
            if ( timer.isRunning() ) {
                timer.stop();
            }
        }
    }

    @Override
    final public void mouseEntered( final MouseEvent e ) {
        if ( ( e.getSource() == this ) && isEnabled() && isRepeatEnabled() ) {
            if ( pressed && !timer.isRunning() ) {
                modifiers = e.getModifiersEx();
                timer.setInitialDelay( delay );
                timer.start();
            }
        }
    }

    @Override
    final public void mouseExited( final MouseEvent e ) {
        if ( e.getSource() == this ) {
            if ( timer.isRunning() ) {
                timer.stop();
            }
        }
    }

    @Override
    final public void mousePressed( final MouseEvent e ) {
        if ( ( e.getSource() == this ) && isEnabled() && isRepeatEnabled() ) {
            pressed = true;
            if ( !timer.isRunning() ) {
                modifiers = e.getModifiersEx();
                timer.setInitialDelay( initial_delay );
                timer.start();
            }
        }
    }

    @Override
    final public void mouseReleased( final MouseEvent e ) {
        if ( e.getSource() == this ) {
            pressed = false;
            if ( timer.isRunning() ) {
                timer.stop();
            }
        }
    }

    final public void setDelay( final int d ) {
        delay = d;
    }

    @Override
    final public void setEnabled( final boolean e ) {
        if ( e != super.isEnabled() ) {
            pressed = false;
            if ( timer.isRunning() ) {
                timer.stop();
            }
        }
        super.setEnabled( e );
    }

    final public void setInitialDelay( final int d ) {
        initial_delay = d;
    }

    final public void setRepeatEnabled( final boolean e ) {
        if ( !e ) {
            pressed = false;
            if ( timer.isRunning() ) {
                timer.stop();
            }
        }
        repeat_enabled = e;
    }
}