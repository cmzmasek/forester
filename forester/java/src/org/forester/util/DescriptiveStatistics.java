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

package org.forester.util;

import java.util.List;

public interface DescriptiveStatistics {

    public final static String PLUS_MINUS = "" + ( char ) 177;

    public abstract void addValue( final double d );

    public abstract double arithmeticMean();

    public abstract String asSummary();

    /**
     * Computes the coefficient of variation. Used to express standard deviation
     * independent of units of measure.
     *
     * @return
     */
    public abstract double coefficientOfVariation();

    public abstract double[] getDataAsDoubleArray();

    public abstract List<Double> getData();

    public abstract double getMax();

    public abstract double getMin();

    public abstract int getN();

    public abstract double getSum();

    public abstract String getSummaryAsString();

    public abstract double getValue( final int index );

    public abstract double median();

    public abstract double midrange();

    /**
     * Determines relationship between the mean and the median. This reflects
     * how the data differs from the normal bell shaped distribution.
     *
     * @return
     */
    public abstract double pearsonianSkewness();

    public abstract double sampleStandardDeviation();

    public abstract double sampleStandardUnit( final double value );

    public abstract double sampleVariance();

    public abstract double standardErrorOfMean();

    public abstract double sumDeviations();

    @Override
    public abstract String toString();

    public abstract void setDescription( final String desc );

    public abstract String getDescription();
}