
package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.math.BigDecimal;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public class Point implements PhylogenyData {

    private final String       _geodetic_datum;
    private final BigDecimal   _lat;
    private final BigDecimal   _long;
    private final BigDecimal   _alt;
    private final String       _alt_unit;
    public static final String UNKNOWN_GEODETIC_DATUM = "?";

    public Point() {
        this( UNKNOWN_GEODETIC_DATUM, null, null, null, "" );
    }

    public Point( final String geodetic_datum, final BigDecimal lat, final BigDecimal longitude ) {
        this( geodetic_datum, lat, longitude, null, "" );
    }

    public boolean isEmpty() {
        return ( _lat == null ) && ( _long == null ) && ( _alt == null );
    }

    public Point( final String geodetic_datum,
                  final BigDecimal lat,
                  final BigDecimal longitude,
                  final BigDecimal alt,
                  final String alt_unit ) {
        if ( ForesterUtil.isEmpty( geodetic_datum ) ) {
            throw new IllegalArgumentException( "illegal attempt to use empty geodetic datum" );
        }
        if ( ( alt != null ) && ForesterUtil.isEmpty( alt_unit ) ) {
            throw new IllegalArgumentException( "altitude must hava a unit" );
        }
        _geodetic_datum = geodetic_datum;
        _lat = lat;
        _long = longitude;
        _alt = alt;
        _alt_unit = alt_unit;
    }

    @Override
    public StringBuffer asSimpleText() {
        if ( isEmpty() ) {
            return new StringBuffer();
        }
        else if ( getAltitude() == null ) {
            return new StringBuffer( "[" + getLatitude().toPlainString() + ", " + getLongitude() + "]" );
        }
        else {
            return new StringBuffer( "[" + getLatitude().toPlainString() + ", " + getLongitude() + ", " + getAltitude()
                    + getAltiudeUnit() + "]" );
        }
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        return new Point( getGeodeticDatum(),
                          getLatitude() == null ? null : new BigDecimal( getLatitude().toPlainString() ),
                          getLongitude() == null ? null : new BigDecimal( getLongitude().toPlainString() ),
                          getAltitude() == null ? null : new BigDecimal( getAltitude().toPlainString() ),
                          getAltiudeUnit() );
    }

    public BigDecimal getAltitude() {
        return _alt;
    }

    public String getAltiudeUnit() {
        return _alt_unit;
    }

    public String getGeodeticDatum() {
        return _geodetic_datum;
    }

    public BigDecimal getLatitude() {
        return _lat;
    }

    public BigDecimal getLongitude() {
        return _long;
    }

    @Override
    public boolean isEqual( final PhylogenyData point ) {
        throw new UnsupportedOperationException();
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( isEmpty() ) {
            return;
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        if ( getAltitude() != null ) {
            PhylogenyDataUtil.appendOpen( writer,
                                          PhyloXmlMapping.POINT,
                                          PhyloXmlMapping.POINT_GEODETIC_DATUM,
                                          getGeodeticDatum(),
                                          PhyloXmlMapping.POINT_ALTITUDE_UNIT_ATTR,
                                          getAltiudeUnit() );
        }
        else {
            PhylogenyDataUtil.appendOpen( writer,
                                          PhyloXmlMapping.POINT,
                                          PhyloXmlMapping.POINT_GEODETIC_DATUM,
                                          getGeodeticDatum() );
        }
        PhylogenyDataUtil.appendElement( writer,
                                         PhyloXmlMapping.POINT_LATITUDE,
                                         getLatitude().toPlainString(),
                                         indentation );
        PhylogenyDataUtil.appendElement( writer,
                                         PhyloXmlMapping.POINT_LONGITUDE,
                                         getLongitude().toPlainString(),
                                         indentation );
        if ( getAltitude() != null ) {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.POINT_ALTITUDE,
                                             getAltitude().toPlainString(),
                                             indentation );
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.POINT );
    }

    @Override
    public String toString() {
        return asSimpleText().toString();
    }

    static public final boolean isSeemsEmpty( final Point p ) {
        return ( ( ( p.getAltitude() == null ) || ( p.getAltitude().compareTo( BigDecimal.ZERO ) <= 0 ) )
                && ( ( p.getLongitude() == null ) || ( p.getLongitude().compareTo( BigDecimal.ZERO ) <= 0 ) )
                && ( ( p.getLatitude() == null ) || ( p.getLatitude().compareTo( BigDecimal.ZERO ) <= 0 ) )
                && ( ForesterUtil.isEmpty( p.getGeodeticDatum() ) || p.getGeodeticDatum()
                        .equalsIgnoreCase( UNKNOWN_GEODETIC_DATUM ) ) && ( ForesterUtil.isEmpty( p.getAltiudeUnit() ) || p
                .getAltiudeUnit().equalsIgnoreCase( "?" ) ) );
    }
}
