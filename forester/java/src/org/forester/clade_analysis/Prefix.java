package org.forester.clade_analysis;


final class Prefix {
    final String _prefix;
    final double _confidence;
    
    Prefix( final String prefix, final double confidence ) {
        _prefix = prefix;
        _confidence = confidence;
    }

    
   String getPrefix() {
        return _prefix;
    }

    
    double getConfidence() {
        return _confidence;
    }
    
    
    
}
