
package org.forester.clade_analysis;

import java.math.BigDecimal;

final class Prefix {

    private final String _prefix;
    private final BigDecimal _confidence;
    private final String _separator;
    private final String _first;

    Prefix( final String prefix, final String confidence, final String separator ) {
        _prefix = prefix;
        _confidence = new BigDecimal( confidence);
        _separator = separator ;
        if ( _prefix.indexOf( _separator ) < 0) {
            _first = _prefix;
        }
        else {
        _first = _prefix.substring( 0, _prefix.indexOf(_separator ) );
        }
    }
    
    Prefix( final String prefix, final double confidence , final String separator) {
        _prefix = prefix;
        _confidence = new BigDecimal( confidence);
        _separator = separator ;
        if ( _prefix.indexOf( _separator ) < 0) {
            _first = _prefix;
        }
        else {
            _first = _prefix.substring( 0, _prefix.indexOf(_separator ) );
            }
    }

    String getPrefix() {
        return _prefix;
    }
    
    String getPrefixFirstElement() {
       return _first;
    }
    double getConfidence() {
        return _confidence.doubleValue();
    }

    @Override
    public String toString() {
     return getPrefix() + ": " + getConfidence();
    }
}
