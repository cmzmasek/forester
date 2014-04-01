
package org.forester.phylogeny.data;

import java.awt.Color;
import java.awt.Font;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.util.ForesterUtil;

public final class NodeVisualData implements PhylogenyData {

    public static final String SIZE_REF             = "aptx_visualiation:node_size";
    public static final String SIZE_TYPE            = "xsd:decimal";
    private static final byte  DEFAULT_FONT_SIZE    = -1;
    private static final int   DEFAULT_SIZE         = -1;
    private static final int   DEFAULT_TRANSPARANCY = -1;
    private Color              _border_color;
    private Color              _fill_color;
    private NodeFill           _fill_type;
    private String             _font_name;
    private Color              _font_color;
    private byte               _font_size;
    private FontType           _font_style;
    private NodeShape          _shape;
    private float              _size;
    private float              _transparancy;
    private Font               _font;

    public NodeVisualData() {
        init();
    }

    public NodeVisualData( final String font_name,
                           final FontType font_style,
                           final byte font_size,
                           final Color font_color,
                           final NodeShape shape,
                           final NodeFill fill_type,
                           final Color border_color,
                           final Color fill_color,
                           final float size,
                           final float transparancy ) {
        setFontName( font_name );
        setFontStyle( font_style );
        setFontSize( font_size );
        setFontColor( font_color );
        setShape( shape );
        setFillType( fill_type );
        setBorderColor( border_color );
        setFillColor( fill_color );
        setSize( size );
        setTransparancy( transparancy );
    }

    @Override
    public final StringBuffer asSimpleText() {
        return asText();
    }

    @Override
    public final StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        return sb;
    }

    @Override
    public final PhylogenyData copy() {
        return new NodeVisualData( !ForesterUtil.isEmpty( getFontName() ) ? new String( getFontName() ) : null,
                                   getFontStyle(),
                                   getFontSize(),
                                   getFontColor() != null ? new Color( getFontColor().getRed(), getFontColor()
                                           .getGreen(), getFontColor().getBlue() ) : null,
                                   getShape(),
                                   getFillType(),
                                   getBorderColor() != null ? new Color( getBorderColor().getRed(), getBorderColor()
                                           .getGreen(), getBorderColor().getBlue() ) : null,
                                   getFillColor() != null ? new Color( getFillColor().getRed(), getFillColor()
                                           .getGreen(), getFillColor().getBlue() ) : null,
                                   getSize(),
                                   getTransparancy() );
    }

    public final Color getBorderColor() {
        return _border_color;
    }

    public final Color getFillColor() {
        return _fill_color;
    }

    public final NodeFill getFillType() {
        return _fill_type;
    }

    public final String getFontName() {
        return _font_name;
    }

    public final Color getFontColor() {
        return _font_color;
    }

    public final byte getFontSize() {
        return _font_size;
    }

    public final FontType getFontStyle() {
        return _font_style;
    }

    public final NodeShape getShape() {
        return _shape;
    }

    public final float getSize() {
        return _size;
    }

    public final float getTransparancy() {
        return _transparancy;
    }

    public final boolean isEmpty() {
        return ( ForesterUtil.isEmpty( getFontName() ) && ( getFontStyle() == FontType.PLAIN )
                && ( getFontSize() == DEFAULT_FONT_SIZE ) && ( getFontColor() == null )
                && ( getShape() == NodeShape.DEFAULT ) && ( getFillType() == NodeFill.DEFAULT )
                && ( getBorderColor() == null ) && ( getFillColor() == null ) && ( getSize() == DEFAULT_SIZE ) && ( getTransparancy() == DEFAULT_TRANSPARANCY ) );
    }

    @Override
    public final boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    public final void setBorderColor( final Color border_color ) {
        _border_color = border_color;
    }

    public final void setFillColor( final Color fill_color ) {
        _fill_color = fill_color;
    }

    public final void setFillType( final NodeFill fill_type ) {
        _fill_type = fill_type;
    }

    public final void setFontName( final String font_name ) {
        if ( !ForesterUtil.isEmpty( font_name ) ) {
            _font_name = font_name;
        }
        else {
            _font_name = null;
        }
        _font = null;
    }

    public final void setFontColor( final Color font_color ) {
        _font_color = font_color;
    }

    public final void setFontSize( final int font_size ) {
        if ( ( font_size != DEFAULT_FONT_SIZE ) && ( font_size < 0 ) ) {
            throw new IllegalArgumentException( "negative font size: " + font_size );
        }
        else if ( font_size > Byte.MAX_VALUE ) {
            throw new IllegalArgumentException( "font size is too large: " + font_size );
        }
        _font_size = ( byte ) font_size;
        _font = null;
    }

    public final void setFontStyle( final FontType font_style ) {
        _font_style = font_style;
        _font = null;
    }

    public final void setFontStyle( final int font_style ) {
        if ( ( font_style == ( Font.BOLD + Font.ITALIC ) ) ) {
            setFontStyle( FontType.BOLD_ITALIC );
        }
        else if ( font_style == Font.ITALIC ) {
            setFontStyle( FontType.ITALIC );
        }
        else if ( font_style == Font.BOLD ) {
            setFontStyle( FontType.BOLD );
        }
        else {
            setFontStyle( FontType.PLAIN );
        }
    }

    public final void setShape( final NodeShape shape ) {
        _shape = shape;
    }

    public final void setSize( final float size ) {
        _size = size;
    }

    public final void setTransparancy( final float transparancy ) {
        _transparancy = transparancy;
    }

    public final int getFontStyleInt() {
        if ( getFontStyle() == FontType.BOLD ) {
            return Font.BOLD;
        }
        else if ( getFontStyle() == FontType.ITALIC ) {
            return Font.ITALIC;
        }
        else if ( getFontStyle() == FontType.BOLD_ITALIC ) {
            return Font.BOLD + Font.ITALIC;
        }
        return Font.PLAIN;
    }

    public final Font getFont() {
        if ( _font != null ) {
            return _font;
        }
        else if ( !ForesterUtil.isEmpty( getFontName() ) ) {
            _font = new Font( getFontName(), getFontStyleInt(), getFontSize() );
            return _font;
        }
        return null;
    }

    @Override
    public final StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public final void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        throw new UnsupportedOperationException();
    }

    @Override
    public final String toString() {
        return asText().toString();
    }

    private final void init() {
        setFontName( null );
        setFontStyle( FontType.PLAIN );
        setFontSize( DEFAULT_FONT_SIZE );
        setFontColor( null );
        setShape( NodeShape.DEFAULT );
        setFillType( NodeFill.DEFAULT );
        setBorderColor( null );
        setFillColor( null );
        setSize( DEFAULT_SIZE );
        setTransparancy( DEFAULT_TRANSPARANCY );
        _font = null;
    }

    private final List<Property> toProperties() {
        final List<Property> properties = new ArrayList<Property>();
        properties.add( new Property( SIZE_REF, String.valueOf( getSize() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        properties.add( new Property( SIZE_REF, String.valueOf( getShape() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        properties.add( new Property( SIZE_REF, String.valueOf( getFillType() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        properties.add( new Property( SIZE_REF, String.valueOf( getTransparancy() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        properties.add( new Property( SIZE_REF, String.valueOf( getFillColor() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        properties.add( new Property( SIZE_REF, String.valueOf( getBorderColor() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        return properties;
    }

    public enum FontType {
        BOLD, BOLD_ITALIC, ITALIC, PLAIN
    }

    public enum NodeFill {
        DEFAULT, GRADIENT, NONE, SOLID
    }

    public enum NodeShape {
        CIRCLE, DEFAULT, RECTANGLE
    }
}
