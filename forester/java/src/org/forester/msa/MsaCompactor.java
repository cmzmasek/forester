package org.forester.msa;


public class MsaCompactor {
    
    final private Msa _msa;
    
    public MsaCompactor( Msa msa ) {
        
        _msa = msa;
    }
    
    
    
    private calc
    
    private double[] calcGappiness() {
        final double gappiness[] = new double[ _msa.getLength() ];
        final int seqs =     _msa.getNumberOfSequences();
        for( int i = 0; i < gappiness.length; ++i ) {
            gappiness[ i ] = ( double ) MsaMethods.calcGapSumPerColumn( _msa, i )  / _msa.getNumberOfSequences();
            
        }
        return gappiness;
        
    }
    
}
