// $Id:
// $
//
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BasicDescriptiveStatistics implements DescriptiveStatistics {

    private List<Double> _data;
    private double       _sum;
    private double       _min;
    private double       _max;
    private double       _sigma;
    private boolean      _recalc_sigma;
    private String       _desc;

    public BasicDescriptiveStatistics() {
        init();
    }

    public BasicDescriptiveStatistics( final String desc ) {
        init();
        setDescription( desc );
    }

    @Override
    public void addValue( final double d ) {
        _recalc_sigma = true;
        _sum += d;
        _data.add( new Double( d ) );
        if ( d < _min ) {
            _min = d;
        }
        if ( d > _max ) {
            _max = d;
        }
    }

    @Override
    public double arithmeticMean() {
        validate();
        return getSum() / getN();
    }

    @Override
    public String asSummary() {
        if ( getN() > 1 ) {
            return arithmeticMean() + DescriptiveStatistics.PLUS_MINUS + sampleStandardDeviation() + " [" + getMin()
                    + "..." + getMax() + "]";
        }
        else {
            return "" + arithmeticMean();
        }
    }

    @Override
    public double coefficientOfVariation() {
        validate();
        return ( sampleStandardDeviation() / arithmeticMean() );
    }

    @Override
    public List<Double> getData() {
        return _data;
    }

    @Override
    public double[] getDataAsDoubleArray() {
        validate();
        final double[] data_array = new double[ getN() ];
        for( int i = 0; i < getN(); ++i ) {
            data_array[ i ] = getValue( i );
        }
        return data_array;
    }

    @Override
    public String getDescription() {
        return _desc;
    }

    @Override
    public double getMax() {
        validate();
        return _max;
    }

    @Override
    public double getMin() {
        validate();
        return _min;
    }

    @Override
    public int getN() {
        return _data.size();
    }

    @Override
    public double getSum() {
        validate();
        return _sum;
    }

    @Override
    public String getSummaryAsString() {
        validate();
        final double mean = arithmeticMean();
        final double sd = sampleStandardDeviation();
        return "" + mean + ( ( char ) 177 ) + sd + " [" + getN() + "] [" + getMin() + "-" + getMax() + "]";
    }

    @Override
    public double getValue( final int index ) {
        validate();
        return ( ( ( _data.get( index ) ) ).doubleValue() );
    }

    @Override
    public double median() {
        validate();
        double median = 0.0;
        if ( getN() == 1 ) {
            median = getValue( 0 );
        }
        else {
            final int index = ( getN() / 2 );
            final double[] data_array = getDataAsDoubleArray();
            Arrays.sort( data_array );
            if ( ( ( data_array.length ) % 2 ) == 0 ) {
                // even number of data values
                median = ( data_array[ index - 1 ] + data_array[ index ] ) / 2.0;
            }
            else {
                median = data_array[ index ];
            }
        }
        return median;
    }

    @Override
    public double midrange() {
        validate();
        return ( _min + _max ) / 2.0;
    }

    @Override
    public double pearsonianSkewness() {
        validate();
        final double mean = arithmeticMean();
        final double median = median();
        final double sd = sampleStandardDeviation();
        return ( ( 3 * ( mean - median ) ) / sd );
    }

    @Override
    public double sampleStandardDeviation() {
        return Math.sqrt( sampleVariance() );
    }

    @Override
    public double sampleStandardUnit( final double value ) {
        validate();
        return BasicDescriptiveStatistics.sampleStandardUnit( value, arithmeticMean(), sampleStandardDeviation() );
    }

    @Override
    public double sampleVariance() {
        validate();
        if ( getN() < 2 ) {
            return 0;
        }
        return ( sumDeviations() / ( getN() - 1 ) );
    }

    @Override
    public void setDescription( final String desc ) {
        _desc = desc;
    }

    @Override
    public double standardErrorOfMean() {
        validate();
        return ( sampleStandardDeviation() / Math.sqrt( getN() ) );
    }

    @Override
    public double sumDeviations() {
        validate();
        if ( _recalc_sigma ) {
            _recalc_sigma = false;
            _sigma = 0.0;
            final double mean = arithmeticMean();
            for( int i = 0; i < getN(); ++i ) {
                _sigma += Math.pow( ( getValue( i ) - mean ), 2 );
            }
        }
        return _sigma;
    }

    @Override
    public String toString() {
        if ( getN() < 1 ) {
            return "empty data set statistics";
        }
        final StringBuffer sb = new StringBuffer();
        sb.append( "Descriptive statistics:" );
        sb.append( ForesterUtil.getLineSeparator() );
        sb.append( "n                       : " + getN() );
        if ( getN() > 1 ) {
            sb.append( ForesterUtil.getLineSeparator() );
            sb.append( "min                     : " + getMin() );
            sb.append( ForesterUtil.getLineSeparator() );
            sb.append( "max                     : " + getMax() );
            sb.append( ForesterUtil.getLineSeparator() );
            sb.append( "midrange                : " + midrange() );
            sb.append( ForesterUtil.getLineSeparator() );
            sb.append( "median                  : " + median() );
            sb.append( ForesterUtil.getLineSeparator() );
            sb.append( "mean                    : " + arithmeticMean() );
            sb.append( ForesterUtil.getLineSeparator() );
            sb.append( "sd                      : " + sampleStandardDeviation() );
            sb.append( ForesterUtil.getLineSeparator() );
            sb.append( "variance                : " + sampleVariance() );
            sb.append( ForesterUtil.getLineSeparator() );
            sb.append( "standard error of mean  : " + standardErrorOfMean() );
            sb.append( ForesterUtil.getLineSeparator() );
            sb.append( "coefficient of variation: " + coefficientOfVariation() );
            sb.append( ForesterUtil.getLineSeparator() );
            sb.append( "pearsonian skewness     : " + pearsonianSkewness() );
        }
        return sb.toString();
    }

    private void init() {
        _data = new ArrayList<Double>();
        _sum = 0.0;
        _min = Double.MAX_VALUE;
        _max = -Double.MAX_VALUE;
        _sigma = 0.0;
        _recalc_sigma = true;
        _desc = "";
    }

    private void validate() throws ArithmeticException {
        if ( getN() < 1 ) {
            throw new ArithmeticException( "attempt to get a result from empty data set statistics" );
        }
    }

    public static int[] performBinning( final double[] values,
                                        final double min,
                                        final double max,
                                        final int number_of_bins ) {
        if ( min >= max ) {
            throw new IllegalArgumentException( "min [" + min + "] is larger than or equal to max [" + max + "]" );
        }
        if ( number_of_bins < 3 ) {
            throw new IllegalArgumentException( "number of bins is smaller than 3" );
        }
        final int[] bins = new int[ number_of_bins ];
        final double binning_factor = number_of_bins / ( max - min );
        final int last_index = number_of_bins - 1;
        for( final double d : values ) {
            if ( !( ( d > max ) || ( d < min ) ) ) {
                final int bin = ( int ) ( ( d - min ) * binning_factor );
                if ( bin > last_index ) {
                    ++bins[ last_index ];
                }
                else {
                    ++bins[ bin ];
                }
            }
        }
        return bins;
    }

    /**
     * Computes the sample standard unit (z-score). Used to compute 'value' in
     * terms of standard units. Note that 'value', 'mean' and 'sd' must be all
     * from the same sample data.
     *
     * @param value
     *            a double in the sample for which
     * @param mean
     *            the mean of the sample.
     * @param sd
     *            The standard deviation of the sample.
     * @return 'value' in terms of standard units
     */
    public static double sampleStandardUnit( final double value, final double mean, final double sd ) {
        return ( value - mean ) / sd;
    }
}
