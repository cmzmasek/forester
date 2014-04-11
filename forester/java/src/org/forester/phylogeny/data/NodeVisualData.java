
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

    public static final String APTX_VISUALIZATION_REF = "style:";
    public static final int    DEFAULT_SIZE           = -1;
    public static final String FONT_COLOR_REF         = APTX_VISUALIZATION_REF + "font_color";
    public static final String FONT_COLOR_TYPE        = "xsd:token";
    public static final String FONT_REF               = APTX_VISUALIZATION_REF + "font";
    public static final String FONT_SIZE_REF          = APTX_VISUALIZATION_REF + "font_size";
    public static final String FONT_SIZE_TYPE         = "xsd:unsignedByte";
    public static final String FONT_STYLE_BOLD        = "bold";
    public static final String FONT_STYLE_BOLD_ITALIC = "bold_italic";
    public static final String FONT_STYLE_ITALIC      = "italic";
    public static final String FONT_STYLE_PLAIN       = "plain";
    public static final String FONT_STYLE_REF         = APTX_VISUALIZATION_REF + "font_style";
    public static final String FONT_STYLE_TYPE        = "xsd:token";
    public static final String FONT_TYPE              = "xsd:token";
    public static final String NODE_COLOR_REF         = APTX_VISUALIZATION_REF + "node_color";
    public static final String NODE_COLOR_TYPE        = "xsd:token";
    public static final String NODE_FILL_GRADIENT     = "gradient";
    public static final String NODE_FILL_NONE         = "none";
    public static final String NODE_FILL_SOLID        = "solid";
    public static final String NODE_FILL_TYPE_REF     = APTX_VISUALIZATION_REF + "node_fill_type";
    public static final String NODE_FILL_TYPE_TYPE    = "xsd:token";
    public static final String NODE_SHAPE_CIRCLE      = "circle";
    public static final String NODE_SHAPE_RECTANGLE   = "rectangle";
    public static final String NODE_SHAPE_REF         = APTX_VISUALIZATION_REF + "node_shape";
    public static final String NODE_SHAPE_TYPE        = "xsd:token";
    public static final String NODE_SIZE_REF          = APTX_VISUALIZATION_REF + "node_size";
    public static final String NODE_SIZE_TYPE         = "xsd:float";
    public static final String NODE_TRANSPARENCY_REF  = APTX_VISUALIZATION_REF + "node_transparency";
    public static final String NODE_TRANSPARENCY_TYPE = "xsd:float";
    private static final byte  DEFAULT_FONT_SIZE      = -1;
    private static final int   DEFAULT_TRANSPARENCY   = -1;
    private NodeFill           _fill_type;
    private Font               _font;
    private Color              _font_color;
    private String             _font_name;
    private byte               _font_size;
    private FontType           _font_style;
    private Color              _node_color;
    private NodeShape          _shape;
    private float              _size;
    private float              _transparency;

    public NodeVisualData() {
        init();
    }

    public NodeVisualData( final String font_name,
                           final FontType font_style,
                           final byte font_size,
                           final Color font_color,
                           final NodeShape shape,
                           final NodeFill fill_type,
                           final Color node_color,
                           final float size,
                           final float transparency ) {
        setFontName( font_name );
        setFontStyle( font_style );
        setFontSize( font_size );
        setFontColor( font_color );
        setShape( shape );
        setFillType( fill_type );
        setNodeColor( node_color );
        setSize( size );
        setTransparency( transparency );
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
                                   getNodeColor() != null ? new Color( getNodeColor().getRed(), getNodeColor()
                                           .getGreen(), getNodeColor().getBlue() ) : null,
                                   getSize(),
                                   getTransparency() );
    }

    public final NodeFill getFillType() {
        return _fill_type;
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

    public final Color getFontColor() {
        return _font_color;
    }

    public final String getFontName() {
        return _font_name;
    }

    public final byte getFontSize() {
        return _font_size;
    }

    public final FontType getFontStyle() {
        return _font_style;
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

    public final Color getNodeColor() {
        return _node_color;
    }

    public final NodeShape getShape() {
        return _shape;
    }

    public final float getSize() {
        return _size;
    }

    public final float getTransparency() {
        return _transparency;
    }

    public final boolean isEmpty() {
        return ( ForesterUtil.isEmpty( getFontName() ) && ( getFontStyle() == FontType.PLAIN )
                && ( getFontSize() == DEFAULT_FONT_SIZE ) && ( getFontColor() == null )
                && ( getShape() == NodeShape.DEFAULT ) && ( getFillType() == NodeFill.DEFAULT )
                && ( getNodeColor() == null ) && ( getSize() == DEFAULT_SIZE ) && ( getTransparency() == DEFAULT_TRANSPARENCY ) );
    }

    @Override
    public final boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    public void parseProperty( final Property prop ) {
        if ( prop.getRef().equals( FONT_REF ) ) {
            setFontName( prop.getValue().trim() );
        }
        else if ( prop.getRef().equals( FONT_SIZE_REF ) ) {
            int s = -1;
            try {
                s = Integer.parseInt( prop.getValue() );
            }
            catch ( final NumberFormatException e ) {
                return;
            }
            if ( ( s >= 0 ) && ( s < Byte.MAX_VALUE ) ) {
                setFontSize( s );
            }
        }
        else if ( prop.getRef().equals( FONT_STYLE_REF ) ) {
            setFontStyle( prop.getValue() );
        }
        else if ( prop.getRef().equals( FONT_COLOR_REF ) ) {
            try {
                setFontColor( Color.decode( prop.getValue() ) );
            }
            catch ( final NumberFormatException e ) {
                return;
            }
        }
        else if ( prop.getRef().equals( NODE_SIZE_REF ) ) {
            float s = -1.0f;
            try {
                s = Float.parseFloat( prop.getValue() );
            }
            catch ( final NumberFormatException e ) {
                return;
            }
            if ( s >= 0 ) {
                setSize( s );
            }
        }
        else if ( prop.getRef().equals( NODE_COLOR_REF ) ) {
            try {
                setNodeColor( Color.decode( prop.getValue() ) );
            }
            catch ( final NumberFormatException e ) {
                return;
            }
        }
        else if ( prop.getRef().equals( NODE_SHAPE_REF ) ) {
            setShape( prop.getValue() );
        }
        else if ( prop.getRef().equals( NODE_FILL_TYPE_REF ) ) {
            setFillType( prop.getValue() );
        }
    }

    public final void setFillType( final NodeFill fill_type ) {
        _fill_type = fill_type;
    }

    public final void setFillType( final String fill ) {
        if ( fill.equalsIgnoreCase( NODE_FILL_NONE ) ) {
            setFillType( NodeFill.NONE );
        }
        else if ( fill.equalsIgnoreCase( NODE_FILL_SOLID ) ) {
            setFillType( NodeFill.SOLID );
        }
        else if ( fill.equalsIgnoreCase( NODE_FILL_GRADIENT ) ) {
            setFillType( NodeFill.GRADIENT );
        }
        else {
            setFillType( NodeFill.DEFAULT );
        }
    }

    public final void setFontColor( final Color font_color ) {
        _font_color = font_color;
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

    public final void setFontStyle( final String font_style ) {
        if ( font_style.equalsIgnoreCase( FONT_STYLE_BOLD ) ) {
            setFontStyle( FontType.BOLD );
        }
        else if ( font_style.equalsIgnoreCase( FONT_STYLE_ITALIC ) ) {
            setFontStyle( FontType.ITALIC );
        }
        else if ( font_style.equalsIgnoreCase( FONT_STYLE_BOLD_ITALIC ) ) {
            setFontStyle( FontType.BOLD_ITALIC );
        }
        else {
            setFontStyle( FontType.PLAIN );
        }
    }

    public final void setNodeColor( final Color node_color ) {
        _node_color = node_color;
    }

    public final void setShape( final NodeShape shape ) {
        _shape = shape;
    }

    public final void setShape( final String shape ) {
        if ( shape.equalsIgnoreCase( NODE_SHAPE_CIRCLE ) ) {
            setShape( NodeShape.CIRCLE );
        }
        else if ( shape.equalsIgnoreCase( NODE_SHAPE_RECTANGLE ) ) {
            setShape( NodeShape.RECTANGLE );
        }
        else {
            setShape( NodeShape.DEFAULT );
        }
    }

    public final void setSize( final float size ) {
        if ( ( size != DEFAULT_SIZE ) && ( size < 0 ) ) {
            throw new IllegalArgumentException( "negative size: " + size );
        }
        _size = size;
    }

    public final void setTransparency( final float transparency ) {
        _transparency = transparency;
    }

    @Override
    public final StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public final void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        for( final Property p : toProperties() ) {
            p.toPhyloXML( writer, level, indentation );
        }
    }

    @Override
    public final String toString() {
        return asText().toString();
    }

    private String colorToHex( final Color c ) {
        return String.format( "#%02x%02x%02x", c.getRed(), c.getGreen(), c.getBlue() );
    }

    private final void init() {
        setFontName( null );
        setFontStyle( FontType.PLAIN );
        setFontSize( DEFAULT_FONT_SIZE );
        setFontColor( null );
        setShape( NodeShape.DEFAULT );
        setFillType( NodeFill.DEFAULT );
        setNodeColor( null );
        setSize( DEFAULT_SIZE );
        setTransparency( DEFAULT_TRANSPARENCY );
        _font = null;
    }

    private final List<Property> toProperties() {
        final List<Property> properties = new ArrayList<Property>();
        if ( !ForesterUtil.isEmpty( getFontName() ) ) {
            properties.add( new Property( FONT_REF, getFontName(), "", FONT_TYPE, AppliesTo.NODE ) );
        }
        if ( getFontSize() != DEFAULT_FONT_SIZE ) {
            properties.add( new Property( FONT_SIZE_REF,
                                          String.valueOf( getFontSize() ),
                                          "",
                                          FONT_SIZE_TYPE,
                                          AppliesTo.NODE ) );
        }
        if ( getFontStyle() != FontType.PLAIN ) {
            String font_style = "";
            if ( getFontStyle() == FontType.ITALIC ) {
                font_style = FONT_STYLE_ITALIC;
            }
            else if ( getFontStyle() == FontType.BOLD ) {
                font_style = FONT_STYLE_BOLD;
            }
            else if ( getFontStyle() == FontType.BOLD_ITALIC ) {
                font_style = FONT_STYLE_BOLD_ITALIC;
            }
            else {
                throw new RuntimeException( "unknown font style" + getShape() );
            }
            properties.add( new Property( FONT_STYLE_REF, font_style, "", FONT_STYLE_TYPE, AppliesTo.NODE ) );
        }
        if ( getFontColor() != null ) {
            properties.add( new Property( FONT_COLOR_REF,
                                          colorToHex( getFontColor() ),
                                          "",
                                          FONT_COLOR_TYPE,
                                          AppliesTo.NODE ) );
        }
        if ( getShape() != NodeShape.DEFAULT ) {
            String shape = null;
            if ( getShape() == NodeShape.RECTANGLE ) {
                shape = NODE_SHAPE_RECTANGLE;
            }
            else if ( getShape() == NodeShape.CIRCLE ) {
                shape = NODE_SHAPE_CIRCLE;
            }
            else {
                throw new RuntimeException( "unknown node shape" + getShape() );
            }
            properties.add( new Property( NODE_SHAPE_REF, shape, "", NODE_SHAPE_TYPE, AppliesTo.NODE ) );
        }
        if ( getSize() != DEFAULT_SIZE ) {
            properties.add( new Property( NODE_SIZE_REF,
                                          String.valueOf( getSize() ),
                                          "",
                                          NODE_SIZE_TYPE,
                                          AppliesTo.NODE ) );
        }
        if ( getNodeColor() != null ) {
            properties.add( new Property( NODE_COLOR_REF,
                                          colorToHex( getNodeColor() ),
                                          "",
                                          NODE_COLOR_TYPE,
                                          AppliesTo.NODE ) );
        }
        if ( getFillType() != NodeFill.DEFAULT ) {
            String fill = null;
            if ( getFillType() == NodeFill.GRADIENT ) {
                fill = NODE_FILL_GRADIENT;
            }
            else if ( getFillType() == NodeFill.NONE ) {
                fill = NODE_FILL_NONE;
            }
            else if ( getFillType() == NodeFill.SOLID ) {
                fill = NODE_FILL_SOLID;
            }
            else {
                throw new RuntimeException( "unknown fill type " + getFillType() );
            }
            properties.add( new Property( NODE_FILL_TYPE_REF, fill, "", NODE_FILL_TYPE_TYPE, AppliesTo.NODE ) );
        }
        if ( getTransparency() != DEFAULT_TRANSPARENCY ) {
            properties.add( new Property( NODE_TRANSPARENCY_REF,
                                          String.valueOf( getTransparency() ),
                                          "",
                                          NODE_TRANSPARENCY_TYPE,
                                          AppliesTo.NODE ) );
        }
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
