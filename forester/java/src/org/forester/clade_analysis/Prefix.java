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

package org.forester.clade_analysis;

import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.regex.Pattern;

public final class Prefix {

    private final static DecimalFormat df = new DecimalFormat( "0.0###" );
    private final String               _prefix;
    private final BigDecimal           _confidence;
    private final String               _separator;
    private final String               _first;

    public Prefix( final String prefix, final String confidence, final String separator ) {
        _prefix = prefix;
        _confidence = new BigDecimal( confidence );
        _separator = separator;
        if ( _prefix.indexOf( _separator ) < 0 ) {
            _first = _prefix;
        }
        else {
            _first = _prefix.substring( 0, _prefix.indexOf( _separator ) );
        }
    }

    public Prefix( final String prefix, final double confidence, final String separator ) {
        _prefix = prefix;
        _confidence = new BigDecimal( confidence );
        _separator = separator;
        if ( _prefix.indexOf( _separator ) < 0 ) {
            _first = _prefix;
        }
        else {
            _first = _prefix.substring( 0, _prefix.indexOf( _separator ) );
        }
    }

    public String getPrefix() {
        return _prefix;
    }
    
    public String getPrefixRemoveSeparator() {
        return _prefix.replaceAll( Pattern.quote( _separator ), "" );
    }

    public String getPrefixFirstElement() {
        return _first;
    }

    public double getConfidence() {
        return _confidence.doubleValue();
    }

    @Override
    public String toString() {
        return getPrefix() + ": " + df.format( getConfidence() );
    }
    
    public String toStringRemoveSeparator() {
        return getPrefixRemoveSeparator() + ": " + df.format( getConfidence() );
    }
}
