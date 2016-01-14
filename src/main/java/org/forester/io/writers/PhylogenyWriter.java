// $Id:
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

package org.forester.io.writers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

import org.forester.io.parsers.nexus.NexusConstants;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PostOrderStackObject;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public final class PhylogenyWriter {

    public final static boolean         INDENT_PHYLOXML_DEAFULT         = true;
    public final static String          PHYLO_XML_INTENDATION_BASE      = "  ";
    public final static String          PHYLO_XML_VERSION_ENCODING_LINE = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
    public final static String          PHYLO_XML_NAMESPACE_LINE        = "<phyloxml xmlns:xsi=\""
            + ForesterConstants.XML_SCHEMA_INSTANCE
            + "\" xsi:schemaLocation=\""
            + ForesterConstants.PHYLO_XML_LOCATION
            + " "
            + ForesterConstants.PHYLO_XML_LOCATION
            + "/"
            + ForesterConstants.PHYLO_XML_VERSION
            + "/" + ForesterConstants.PHYLO_XML_XSD
            + "\" " + "xmlns=\""
            + ForesterConstants.PHYLO_XML_LOCATION
            + "\">";
    public final static String          PHYLO_XML_END                   = "</phyloxml>";
    private boolean                     _saw_comma;
    private StringBuffer                _buffer;
    private Writer                      _writer;
    private PhylogenyNode               _root;
    private boolean                     _has_next;
    private Stack<PostOrderStackObject> _stack;
    private boolean                     _nh_write_distance_to_parent;
    NH_CONVERSION_SUPPORT_VALUE_STYLE   _nh_conversion_support_style;
    private boolean                     _indent_phyloxml;
    private int                         _node_level;
    private int                         _phyloxml_level;
    private FORMAT                      _format;

    public PhylogenyWriter() {
        setIndentPhyloxml( INDENT_PHYLOXML_DEAFULT );
        setNhConversionSupportStyle( NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
    }

    private void appendPhylogenyLevelPhyloXml( final Writer writer, final Phylogeny tree ) throws IOException {
        final String indentation = new String();
        if ( !ForesterUtil.isEmpty( tree.getName() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.PHYLOGENY_NAME, tree.getName(), indentation );
        }
        if ( tree.getIdentifier() != null ) {
            if ( ForesterUtil.isEmpty( tree.getIdentifier().getProvider() ) ) {
                PhylogenyDataUtil.appendElement( writer,
                                                 PhyloXmlMapping.IDENTIFIER,
                                                 tree.getIdentifier().getValue(),
                                                 indentation );
            }
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.IDENTIFIER,
                                             tree.getIdentifier().getValue(),
                                             PhyloXmlMapping.IDENTIFIER_PROVIDER_ATTR,
                                             tree.getIdentifier().getProvider(),
                                             indentation );
        }
        if ( !ForesterUtil.isEmpty( tree.getDescription() ) ) {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.PHYLOGENY_DESCRIPTION,
                                             tree.getDescription(),
                                             indentation );
        }
        if ( tree.getConfidence() != null ) {
            if ( ForesterUtil.isEmpty( tree.getConfidence().getType() ) ) {
                PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.CONFIDENCE, tree.getConfidence().getValue()
                                                 + "", indentation );
            }
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.CONFIDENCE,
                                             tree.getConfidence().getValue() + "",
                                             PhyloXmlMapping.CONFIDENCE_TYPE_ATTR,
                                             tree.getConfidence().getType(),
                                             indentation );
        }
    }

    private StringBuffer createIndentation() {
        if ( !isIndentPhyloxml() ) {
            return null;
        }
        final StringBuffer sb = new StringBuffer( getNodeLevel() * 2 );
        for( int i = 0; i < getNodeLevel(); ++i ) {
            sb.append( PhylogenyWriter.PHYLO_XML_INTENDATION_BASE );
        }
        return sb;
    }

    private void decreaseNodeLevel() {
        --_node_level;
    }

    private StringBuffer getBuffer() {
        return _buffer;
    }

    private int getNodeLevel() {
        return _node_level;
    }

    private StringBuffer getOutput( final Phylogeny tree ) throws IOException {
        if ( getOutputFormt() == FORMAT.PHYLO_XML ) {
            throw new RuntimeException( "method inappropriately called" );
        }
        if ( tree != null ) {
            reset( tree );
            while ( isHasNext() ) {
                next();
            }
            if ( getOutputFormt() == FORMAT.NH ) {
                getBuffer().append( ';' );
            }
            return getBuffer();
        }
        else {
            return new StringBuffer( 0 );
        }
    }

    private FORMAT getOutputFormt() {
        return _format;
    }

    private int getPhyloXmlLevel() {
        return _phyloxml_level;
    }

    private PhylogenyNode getRoot() {
        return _root;
    }

    private Stack<PostOrderStackObject> getStack() {
        return _stack;
    }

    private Writer getWriter() {
        return _writer;
    }

    private void increaseNodeLevel() {
        ++_node_level;
    }

    private boolean isHasNext() {
        return _has_next;
    }

    private boolean isIndentPhyloxml() {
        return _indent_phyloxml;
    }

    private boolean isSawComma() {
        return _saw_comma;
    }

    private boolean isWriteDistanceToParentInNH() {
        return _nh_write_distance_to_parent;
    }

    private void next() throws IOException {
        while ( true ) {
            final PostOrderStackObject si = getStack().pop();
            final PhylogenyNode node = si.getNode();
            final int phase = si.getPhase();
            if ( phase > node.getNumberOfDescendants() ) {
                setHasNext( node != getRoot() );
                if ( ( getOutputFormt() != FORMAT.PHYLO_XML ) || node.isExternal() ) {
                    if ( !node.isRoot() && node.isFirstChildNode() ) {
                        increaseNodeLevel();
                    }
                    if ( getOutputFormt() == FORMAT.PHYLO_XML ) {
                        writeNode( node, createIndentation() );
                    }
                    else {
                        writeNode( node, null );
                    }
                }
                if ( !node.isRoot() ) {
                    if ( !node.isLastChildNode() ) {
                        writeCladeSeparator();
                    }
                    else {
                        writeCloseClade();
                    }
                }
                return;
            }
            else {
                getStack().push( new PostOrderStackObject( node, ( phase + 1 ) ) );
                if ( node.isInternal() ) {
                    getStack().push( new PostOrderStackObject( node.getChildNode( phase - 1 ), 1 ) );
                    writeOpenClade( node );
                    if ( getOutputFormt() == FORMAT.PHYLO_XML ) {
                        if ( phase == 1 ) {
                            writeNode( node, createIndentation() );
                        }
                    }
                }
            }
        }
    }

    private void reset( final Phylogeny tree ) {
        setBuffer( new StringBuffer() );
        setWriter( null );
        setSawComma( false );
        setHasNext( true );
        setRoot( tree.getRoot() );
        setStack( new Stack<PostOrderStackObject>() );
        getStack().push( new PostOrderStackObject( tree.getRoot(), 1 ) );
        setNodeLevel( 1 );
    }

    private void reset( final Writer writer, final Phylogeny tree ) {
        setBuffer( null );
        setWriter( writer );
        setSawComma( false );
        setHasNext( true );
        setRoot( tree.getRoot() );
        setStack( new Stack<PostOrderStackObject>() );
        getStack().push( new PostOrderStackObject( tree.getRoot(), 1 ) );
        setNodeLevel( 1 );
    }

    private void setBuffer( final StringBuffer buffer ) {
        _buffer = buffer;
    }

    private void setHasNext( final boolean has_next ) {
        _has_next = has_next;
    }

    public void setIndentPhyloxml( final boolean indent_phyloxml ) {
        _indent_phyloxml = indent_phyloxml;
    }

    private void setNodeLevel( final int level ) {
        _node_level = level;
    }

    private void setOutputFormt( final FORMAT format ) {
        _format = format;
    }

    private void setPhyloXmlLevel( final int phyloxml_level ) {
        _phyloxml_level = phyloxml_level;
    }

    private void setRoot( final PhylogenyNode root ) {
        _root = root;
    }

    private void setSawComma( final boolean saw_comma ) {
        _saw_comma = saw_comma;
    }

    private void setStack( final Stack<PostOrderStackObject> stack ) {
        _stack = stack;
    }

    private void setWriteDistanceToParentInNH( final boolean nh_write_distance_to_parent ) {
        _nh_write_distance_to_parent = nh_write_distance_to_parent;
    }

    private void setWriter( final Writer writer ) {
        _writer = writer;
    }

    public void toNewHampshire( final List<Phylogeny> trees,
                                final boolean write_distance_to_parent,
                                final File out_file,
                                final String separator ) throws IOException {
        final Iterator<Phylogeny> it = trees.iterator();
        final StringBuffer sb = new StringBuffer();
        while ( it.hasNext() ) {
            sb.append( toNewHampshire( it.next(), write_distance_to_parent ) );
            sb.append( separator );
        }
        writeToFile( sb, out_file );
    }

    public StringBuffer toNewHampshire( final Phylogeny tree,
                                        final boolean nh_write_distance_to_parent,
                                        final NH_CONVERSION_SUPPORT_VALUE_STYLE svs ) throws IOException {
        setOutputFormt( FORMAT.NH );
        setNhConversionSupportStyle( svs );
        setWriteDistanceToParentInNH( nh_write_distance_to_parent );
        return getOutput( tree );
    }

    public StringBuffer toNewHampshire( final Phylogeny tree, final boolean nh_write_distance_to_parent )
            throws IOException {
        setOutputFormt( FORMAT.NH );
        setWriteDistanceToParentInNH( nh_write_distance_to_parent );
        return getOutput( tree );
    }

    public void toNewHampshire( final Phylogeny tree, final boolean write_distance_to_parent, final File out_file )
            throws IOException {
        writeToFile( toNewHampshire( tree, write_distance_to_parent ), out_file );
    }

    public void toNewHampshire( final Phylogeny tree,
                                final boolean write_distance_to_parent,
                                final NH_CONVERSION_SUPPORT_VALUE_STYLE svs,
                                final File out_file ) throws IOException {
        writeToFile( toNewHampshire( tree, write_distance_to_parent, svs ), out_file );
    }

    public void toNewHampshire( final Phylogeny[] trees,
                                final boolean write_distance_to_parent,
                                final File out_file,
                                final String separator ) throws IOException {
        final StringBuffer sb = new StringBuffer();
        for( final Phylogeny element : trees ) {
            sb.append( toNewHampshire( element, write_distance_to_parent ) );
            sb.append( separator );
        }
        writeToFile( sb, out_file );
    }

    public void toNewHampshireX( final List<Phylogeny> trees, final File out_file, final String separator )
            throws IOException {
        final Iterator<Phylogeny> it = trees.iterator();
        final StringBuffer sb = new StringBuffer();
        while ( it.hasNext() ) {
            sb.append( toNewHampshireX( it.next() ) );
            sb.append( separator );
        }
        writeToFile( sb, out_file );
    }

    public StringBuffer toNewHampshireX( final Phylogeny tree ) throws IOException {
        setOutputFormt( FORMAT.NHX );
        return getOutput( tree );
    }

    public void toNewHampshireX( final Phylogeny tree, final File out_file ) throws IOException {
        writeToFile( toNewHampshireX( tree ), out_file );
    }

    public void toNewHampshireX( final Phylogeny[] trees, final File out_file, final String separator )
            throws IOException {
        final StringBuffer sb = new StringBuffer();
        for( final Phylogeny element : trees ) {
            sb.append( toNewHampshireX( element ) );
            sb.append( separator );
        }
        writeToFile( sb, out_file );
    }

    public void toNexus( final File out_file, final Phylogeny tree, final NH_CONVERSION_SUPPORT_VALUE_STYLE svs )
            throws IOException {
        final Writer writer = new BufferedWriter( new PrintWriter( out_file ) );
        final List<Phylogeny> trees = new ArrayList<Phylogeny>( 1 );
        trees.add( tree );
        writeNexusStart( writer );
        writeNexusTaxaBlock( writer, tree );
        writeNexusTreesBlock( writer, trees, svs );
        writer.flush();
        writer.close();
    }

    public StringBuffer toNexus( final Phylogeny tree, final NH_CONVERSION_SUPPORT_VALUE_STYLE svs ) throws IOException {
        final StringWriter string_writer = new StringWriter();
        final Writer writer = new BufferedWriter( string_writer );
        final List<Phylogeny> trees = new ArrayList<Phylogeny>( 1 );
        trees.add( tree );
        writeNexusStart( writer );
        writeNexusTaxaBlock( writer, tree );
        writeNexusTreesBlock( writer, trees, svs );
        writer.flush();
        writer.close();
        return string_writer.getBuffer();
    }

    public void toPhyloXML( final File out_file,
                            final List<Phylogeny> trees,
                            final int phyloxml_level,
                            final String separator ) throws IOException {
        final Writer writer = new BufferedWriter( new PrintWriter( out_file ) );
        toPhyloXML( writer, trees, phyloxml_level, separator );
        writer.flush();
        writer.close();
    }

    public void toPhyloXML( final File out_file, final Phylogeny tree, final int phyloxml_level ) throws IOException {
        final Writer writer = new BufferedWriter( new PrintWriter( out_file ) );
        writePhyloXmlStart( writer );
        toPhyloXMLNoPhyloXmlSource( writer, tree, phyloxml_level );
        writePhyloXmlEnd( writer );
        writer.flush();
        writer.close();
    }

    public StringBuffer toPhyloXML( final Phylogeny tree, final int phyloxml_level ) throws IOException {
        final StringWriter string_writer = new StringWriter();
        final Writer writer = new BufferedWriter( string_writer );
        setPhyloXmlLevel( phyloxml_level );
        setOutputFormt( FORMAT.PHYLO_XML );
        writePhyloXmlStart( writer );
        writeOutput( writer, tree );
        writePhyloXmlEnd( writer );
        writer.flush();
        writer.close();
        return string_writer.getBuffer();
    }

    public void toPhyloXML( final Phylogeny[] trees,
                            final int phyloxml_level,
                            final File out_file,
                            final String separator ) throws IOException {
        final Writer writer = new BufferedWriter( new PrintWriter( out_file ) );
        toPhyloXML( writer, trees, phyloxml_level, separator );
        writer.flush();
        writer.close();
    }

    public void toPhyloXML( final Phylogeny phy, final int phyloxml_level, final File out_file ) throws IOException {
        final Writer writer = new BufferedWriter( new PrintWriter( out_file ) );
        toPhyloXML( writer, phy, phyloxml_level );
        writer.flush();
        writer.close();
    }

    public void toPhyloXML( final Writer writer,
                            final List<Phylogeny> trees,
                            final int phyloxml_level,
                            final String separator ) throws IOException {
        writePhyloXmlStart( writer );
        final Iterator<Phylogeny> it = trees.iterator();
        while ( it.hasNext() ) {
            toPhyloXMLNoPhyloXmlSource( writer, it.next(), phyloxml_level );
            writer.write( separator );
        }
        writePhyloXmlEnd( writer );
    }

    public void toPhyloXML( final Writer writer, final Phylogeny tree, final int phyloxml_level ) throws IOException {
        setPhyloXmlLevel( phyloxml_level );
        setOutputFormt( FORMAT.PHYLO_XML );
        writePhyloXmlStart( writer );
        writeOutput( writer, tree );
        writePhyloXmlEnd( writer );
    }

    public void toPhyloXML( final Writer writer,
                            final Phylogeny[] trees,
                            final int phyloxml_level,
                            final String separator ) throws IOException {
        writePhyloXmlStart( writer );
        for( final Phylogeny phylogeny : trees ) {
            toPhyloXMLNoPhyloXmlSource( writer, phylogeny, phyloxml_level );
            writer.write( separator );
        }
        writePhyloXmlEnd( writer );
    }

    private void toPhyloXMLNoPhyloXmlSource( final Writer writer, final Phylogeny tree, final int phyloxml_level )
            throws IOException {
        setPhyloXmlLevel( phyloxml_level );
        setOutputFormt( FORMAT.PHYLO_XML );
        writeOutput( writer, tree );
    }

    private void writeCladeSeparator() {
        setSawComma( true );
        if ( ( getOutputFormt() == FORMAT.NHX ) || ( getOutputFormt() == FORMAT.NH ) ) {
            getBuffer().append( "," );
        }
    }

    private void writeCloseClade() throws IOException {
        decreaseNodeLevel();
        if ( getOutputFormt() == FORMAT.PHYLO_XML ) {
            getWriter().write( ForesterUtil.LINE_SEPARATOR );
            if ( isIndentPhyloxml() ) {
                getWriter().write( createIndentation().toString() );
            }
            PhylogenyDataUtil.appendClose( getWriter(), PhyloXmlMapping.CLADE );
        }
        else if ( ( getOutputFormt() == FORMAT.NHX ) || ( getOutputFormt() == FORMAT.NH ) ) {
            getBuffer().append( ")" );
        }
    }

    private void writeNode( final PhylogenyNode node, final StringBuffer indentation ) throws IOException {
        if ( getOutputFormt() == FORMAT.PHYLO_XML ) {
            if ( node.isExternal() ) {
                getWriter().write( ForesterUtil.LINE_SEPARATOR );
                if ( indentation != null ) {
                    getWriter().write( indentation.toString() );
                }
                PhylogenyDataUtil.appendOpen( getWriter(), PhyloXmlMapping.CLADE );
            }
            PhyloXmlNodeWriter.toPhyloXml( getWriter(),
                                           node,
                                           getPhyloXmlLevel(),
                                           indentation != null ? indentation.toString() : "" );
            if ( node.isExternal() ) {
                getWriter().write( ForesterUtil.LINE_SEPARATOR );
                if ( indentation != null ) {
                    getWriter().write( indentation.toString() );
                }
                PhylogenyDataUtil.appendClose( getWriter(), PhyloXmlMapping.CLADE );
            }
        }
        else if ( getOutputFormt() == FORMAT.NHX ) {
            getBuffer().append( node.toNewHampshireX() );
        }
        else if ( getOutputFormt() == FORMAT.NH ) {
            getBuffer().append( node.toNewHampshire( isWriteDistanceToParentInNH(), getNhConversionSupportStyle() ) );
        }
    }

    private NH_CONVERSION_SUPPORT_VALUE_STYLE getNhConversionSupportStyle() {
        return _nh_conversion_support_style;
    }

    private void setNhConversionSupportStyle( final NH_CONVERSION_SUPPORT_VALUE_STYLE nh_conversion_support_style ) {
        _nh_conversion_support_style = nh_conversion_support_style;
    }

    private void writeOpenClade( final PhylogenyNode node ) throws IOException {
        if ( !isSawComma() ) {
            if ( !node.isRoot() && node.isFirstChildNode() ) {
                increaseNodeLevel();
            }
            if ( getOutputFormt() == FORMAT.PHYLO_XML ) {
                getWriter().write( ForesterUtil.LINE_SEPARATOR );
                if ( isIndentPhyloxml() ) {
                    getWriter().write( createIndentation().toString() );
                }
                if ( node.isCollapse() ) {
                    PhylogenyDataUtil.appendOpen( getWriter(),
                                                  PhyloXmlMapping.CLADE,
                                                  PhyloXmlMapping.NODE_COLLAPSE,
                            "true" );
                }
                else {
                    PhylogenyDataUtil.appendOpen( getWriter(), PhyloXmlMapping.CLADE );
                }
            }
            else if ( ( getOutputFormt() == FORMAT.NHX ) || ( getOutputFormt() == FORMAT.NH ) ) {
                getBuffer().append( "(" );
            }
        }
        setSawComma( false );
    }

    private void writeOutput( final Writer writer, final Phylogeny tree ) throws IOException {
        if ( getOutputFormt() != FORMAT.PHYLO_XML ) {
            throw new RuntimeException( "method inappropriately called" );
        }
        if ( tree != null ) {
            reset( writer, tree );
            String unit = "";
            String type = "";
            if ( !ForesterUtil.isEmpty( tree.getDistanceUnit() ) ) {
                unit = tree.getDistanceUnit();
            }
            if ( !ForesterUtil.isEmpty( tree.getType() ) ) {
                type = tree.getType();
            }
            PhylogenyDataUtil.appendOpen( writer,
                                          PhyloXmlMapping.PHYLOGENY,
                                          PhyloXmlMapping.PHYLOGENY_IS_ROOTED_ATTR,
                                          tree.isRooted() + "",
                                          PhyloXmlMapping.PHYLOGENY_BRANCHLENGTH_UNIT_ATTR,
                                          unit,
                                          PhyloXmlMapping.PHYLOGENY_TYPE_ATTR,
                                          type,
                                          PhyloXmlMapping.PHYLOGENY_IS_REROOTABLE_ATTR,
                                          tree.isRerootable() + "" );
            appendPhylogenyLevelPhyloXml( writer, tree );
            while ( isHasNext() ) {
                next();
            }
            writer.write( ForesterUtil.LINE_SEPARATOR );
            PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.PHYLOGENY );
        }
    }

    private void writeToFile( final StringBuffer sb, final File out_file ) throws IOException {
        if ( out_file.exists() ) {
            throw new IOException( "attempt to overwrite existing file \"" + out_file.getAbsolutePath() + "\"" );
        }
        final PrintWriter out = new PrintWriter( new FileWriter( out_file ), true );
        if ( getOutputFormt() == FORMAT.PHYLO_XML ) {
            out.print( PHYLO_XML_VERSION_ENCODING_LINE );
            out.print( ForesterUtil.LINE_SEPARATOR );
            out.print( PHYLO_XML_NAMESPACE_LINE );
            out.print( ForesterUtil.LINE_SEPARATOR );
        }
        out.print( sb );
        if ( getOutputFormt() == FORMAT.PHYLO_XML ) {
            out.print( ForesterUtil.LINE_SEPARATOR );
            out.print( PHYLO_XML_END );
        }
        out.flush();
        out.close();
    }

    public static PhylogenyWriter createPhylogenyWriter() {
        return new PhylogenyWriter();
    }

    private static void writeNexusStart( final Writer writer ) throws IOException {
        writer.write( NexusConstants.NEXUS );
        writer.write( ForesterUtil.LINE_SEPARATOR );
    }

    public static void writeNexusTaxaBlock( final Writer writer, final Phylogeny tree ) throws IOException {
        writer.write( NexusConstants.BEGIN_TAXA );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( " " );
        writer.write( NexusConstants.DIMENSIONS );
        writer.write( " " );
        writer.write( NexusConstants.NTAX );
        writer.write( "=" );
        writer.write( String.valueOf( tree.getNumberOfExternalNodes() ) );
        writer.write( ";" );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( " " );
        writer.write( NexusConstants.TAXLABELS );
        for( final PhylogenyNodeIterator it = tree.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            writer.write( " " );
            String data = "";
            if ( !ForesterUtil.isEmpty( node.getName() ) ) {
                data = node.getName();
            }
            else if ( node.getNodeData().isHasTaxonomy() ) {
                if ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
                    data = node.getNodeData().getTaxonomy().getTaxonomyCode();
                }
                else if ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getScientificName() ) ) {
                    data = node.getNodeData().getTaxonomy().getScientificName();
                }
                else if ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getCommonName() ) ) {
                    data = node.getNodeData().getTaxonomy().getCommonName();
                }
            }
            else if ( node.getNodeData().isHasSequence() ) {
                if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getName() ) ) {
                    data = node.getNodeData().getSequence().getName();
                }
                else if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getSymbol() ) ) {
                    data = node.getNodeData().getSequence().getSymbol();
                }
                else if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getGeneName() ) ) {
                    data = node.getNodeData().getSequence().getGeneName();
                }
            }
            writer.write( ForesterUtil.santitizeStringForNH( data ).toString() );
        }
        writer.write( ";" );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( NexusConstants.END );
        writer.write( ForesterUtil.LINE_SEPARATOR );
    }

    public static void writeNexusTreesBlock( final Writer writer,
                                             final List<Phylogeny> trees,
                                             final NH_CONVERSION_SUPPORT_VALUE_STYLE svs ) throws IOException {
        writer.write( NexusConstants.BEGIN_TREES );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        int i = 1;
        for( final Phylogeny phylogeny : trees ) {
            writer.write( " " );
            writer.write( NexusConstants.TREE );
            writer.write( " " );
            if ( !ForesterUtil.isEmpty( phylogeny.getName() ) ) {
                writer.write( "\'" );
                writer.write( phylogeny.getName() );
                writer.write( "\'" );
            }
            else {
                writer.write( "tree" );
                writer.write( String.valueOf( i ) );
            }
            writer.write( "=" );
            if ( phylogeny.isRooted() ) {
                writer.write( "[&R]" );
            }
            else {
                writer.write( "[&U]" );
            }
            writer.write( phylogeny.toNewHampshire( svs ) );
            writer.write( ForesterUtil.LINE_SEPARATOR );
            i++;
        }
        writer.write( NexusConstants.END );
        writer.write( ForesterUtil.LINE_SEPARATOR );
    }

    private static void writePhyloXmlEnd( final Writer writer ) throws IOException {
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( PhylogenyWriter.PHYLO_XML_END );
    }

    private static void writePhyloXmlStart( final Writer writer ) throws IOException {
        writer.write( PhylogenyWriter.PHYLO_XML_VERSION_ENCODING_LINE );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( PhylogenyWriter.PHYLO_XML_NAMESPACE_LINE );
        writer.write( ForesterUtil.LINE_SEPARATOR );
    }

    public static enum FORMAT {
        NH, NHX, PHYLO_XML, NEXUS;
    }
}
