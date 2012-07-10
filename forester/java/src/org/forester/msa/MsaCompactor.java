package org.forester.msa;

import org.forester.sequence.Sequence;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;


public class MsaCompactor {
    
    final private Msa _msa;
    
    public MsaCompactor( Msa msa ) {
        
        _msa = msa;
    }
    
    
    
    private DescriptiveStatistics[] calc() {
        final double gappiness[] = calcGappiness();
        final DescriptiveStatistics stats[] = new DescriptiveStatistics[ _msa.getNumberOfSequences() ];
        for ( int row = 0; row < _msa.getNumberOfSequences(); ++row ) {
            stats[ row ] = new BasicDescriptiveStatistics();
            for( int col = 0; col < _msa.getLength(); ++col ) {
                if ( _msa.getResidueAt( row, col ) != Sequence.GAP ) {
                    stats[ row ].addValue( gappiness[ col ] );
                    
                }
            }
        }
        return stats;
    }
    
    private double[] calcGappiness() {
        final double gappiness[] = new double[ _msa.getLength() ];
        final int seqs =  _msa.getNumberOfSequences();
        for( int i = 0; i < gappiness.length; ++i ) {
            gappiness[ i ] = ( double ) MsaMethods.calcGapSumPerColumn( _msa, i )  / _msa.getNumberOfSequences();
            
        }
        return gappiness;
        
    }
    
}
