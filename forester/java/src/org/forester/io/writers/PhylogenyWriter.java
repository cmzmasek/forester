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

package org.forester.io.writers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import org.forester.io.parsers.nexus.NexusConstants;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PostOrderStackObject;
import org.forester.sequence.MolecularSequence;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public final class PhylogenyWriter {

    public final static boolean         INDENT_PHYLOXML_DEAFULT         = true;
    public final static String          PHYLO_XML_INTENDATION_BASE      = "  ";
    public final static String          PHYLO_XML_VERSION_ENCODING_LINE = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
    /*
     * The phyloXML NAMESPACE is written; the xsi:schemaLocation HINT deliberately is NOT.
     *
     * schemaLocation is only a hint, and the phyloxml.org domain it used to point at lapsed and is now held by
     * someone else (2026-09-10) -- so a tool that honours the hint and fetches it would be fetching whatever that
     * party chooses to serve. Nothing here ever needed it: the validating parser pins the schema explicitly from
     * the copy bundled in the jar (PhyloXmlParser.createPhyloXmlParserXsdValidating -> JAXP_SCHEMA_SOURCE), and
     * the namespace URI is only ever string-compared, never dereferenced.
     *
     * The NAMESPACE itself must not change. It is an opaque identifier written into every phyloXML file in
     * existence; changing it would break the format and every other tool rather than protect anyone.
     */
    public final static String          PHYLO_XML_NAMESPACE_LINE        = "<phyloxml xmlns=\""
            + ForesterConstants.PHYLO_XML_LOCATION
            + "\">";
    public final static String          PHYLO_XML_END                   = "</phyloxml>";
    private boolean                     _saw_comma;
    private StringBuffer                _buffer;
    // "node1", "node2", ... for external nodes that nothing else names, by tip index. Built per tree so
    // that the New Hampshire string and the Nexus Taxa/Characters blocks use the SAME placeholder for the
    // same tip: a placeholder applied in only one of them would put the blocks right back into the
    // disagreement that nexusTaxonLabel exists to prevent.
    private Map<PhylogenyNode, String>  _tip_placeholders = new HashMap<PhylogenyNode, String>();
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

    /** Writes the tree-level {@code <property applies_to="phylogeny">} elements. Per the phyloXML schema the
     *  {@code property} children of {@code <phylogeny>} come AFTER {@code <clade>}, so this is called by
     *  {@link #writeOutput} once the clade body has been written, NOT from {@link #appendPhylogenyLevelPhyloXml}
     *  (which emits the name/id/description/confidence header BEFORE the clade). */
    private void appendPhylogenyLevelProperties( final Writer writer, final Phylogeny tree ) throws IOException {
        if ( tree.isHasProperties() ) {
            tree.getProperties().toPhyloXML( writer, getPhyloXmlLevel(), new String() );
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
        _tip_placeholders = placeholdersNeeded( tree );
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
        _tip_placeholders = placeholdersNeeded( tree );
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
        final Writer writer = new BufferedWriter( new PrintWriter( out_file, ForesterConstants.UTF_8 ) );
        final List<Phylogeny> trees = new ArrayList<Phylogeny>( 1 );
        trees.add( tree );
        writeNexusStart( writer );
        writeNexusTaxaBlock( writer, tree );
        writeNexusCharactersBlock( writer, tree );
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
        writeNexusCharactersBlock( writer, tree );
        writeNexusTreesBlock( writer, trees, svs );
        writer.flush();
        writer.close();
        return string_writer.getBuffer();
    }

    public void toPhyloXML( final File out_file,
                            final List<Phylogeny> trees,
                            final int phyloxml_level,
                            final String separator ) throws IOException {
        final Writer writer = new BufferedWriter( new PrintWriter( out_file, ForesterConstants.UTF_8 ) );
        toPhyloXML( writer, trees, phyloxml_level, separator );
        writer.flush();
        writer.close();
    }

    public void toPhyloXML( final File out_file, final Phylogeny tree, final int phyloxml_level ) throws IOException {
        final Writer writer = new BufferedWriter( new PrintWriter( out_file, ForesterConstants.UTF_8 ) );
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
            getBuffer().append( node.toNewHampshireX( _tip_placeholders.get( node ) ) );
        }
        else if ( getOutputFormt() == FORMAT.NH ) {
            getBuffer().append( node.toNewHampshire( isWriteDistanceToParentInNH(),
                                                     getNhConversionSupportStyle(),
                                                     false,
                                                     _tip_placeholders.get( node ) ) );
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
            appendPhylogenyLevelProperties( writer, tree ); // schema order: property* comes after clade
            writer.write( ForesterUtil.LINE_SEPARATOR );
            PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.PHYLOGENY );
        }
    }

    private void writeToFile( final StringBuffer sb, final File out_file ) throws IOException {
        if ( out_file.exists() ) {
            throw new IOException( "attempt to overwrite existing file \"" + out_file.getAbsolutePath() + "\"" );
        }
        final PrintWriter out = new PrintWriter( out_file, ForesterConstants.UTF_8 );
        out.print( sb );
        out.flush();
        out.close();
    }

    /** Only the New Hampshire and NHX paths read the placeholders, and building them walks every tip -- so
     *  phyloXML, which never consults the map, should not pay for it on a hundred-thousand-tip tree. */
    private Map<PhylogenyNode, String> placeholdersNeeded( final Phylogeny tree ) {
        return ( ( getOutputFormt() == FORMAT.NH ) || ( getOutputFormt() == FORMAT.NHX ) )
                ? tipPlaceholders( tree ) : new HashMap<PhylogenyNode, String>();
    }

    /**
     * A placeholder label for every external node, by tip index, in the order the tips appear in the tree
     * (which is the order they appear in a New Hampshire string and in iteratorExternalForward alike).
     * Only used for a node that nothing else names.
     */
    static Map<PhylogenyNode, String> tipPlaceholders( final Phylogeny tree ) {
        final Map<PhylogenyNode, String> m = new HashMap<PhylogenyNode, String>();
        if ( ( tree == null ) || tree.isEmpty() ) {
            return m;
        }
        // A tip may literally be called "node2". Minting that same token for a DIFFERENT tip gives two taxa
        // one label, which is illegal Nexus -- and the parser reads the repeat as an interleaved
        // continuation, so both tips come back carrying the two sequences concatenated. Collect the labels
        // the tree already produces and step over them.
        final java.util.Set<String> taken = new java.util.HashSet<String>();
        for( final PhylogenyNodeIterator it = tree.iteratorExternalForward(); it.hasNext(); ) {
            final String label = nexusTaxonLabel( it.next(), null );
            if ( label.length() > 0 ) {
                taken.add( label );
            }
        }
        int i = 1;
        int spare = tree.getNumberOfExternalNodes() + 1;
        for( final PhylogenyNodeIterator it = tree.iteratorExternalForward(); it.hasNext(); ) {
            String candidate = "node" + i++;
            while ( taken.contains( candidate ) ) {
                candidate = "node" + spare++;
            }
            taken.add( candidate );
            m.put( it.next(), candidate );
        }
        return m;
    }

    public static PhylogenyWriter createPhylogenyWriter() {
        return new PhylogenyWriter();
    }

    private static void writeNexusStart( final Writer writer ) throws IOException {
        writer.write( NexusConstants.NEXUS );
        writer.write( ForesterUtil.LINE_SEPARATOR );
    }

    public static void writeNexusTaxaBlock( final Writer writer, final Phylogeny tree ) throws IOException {
        final Map<PhylogenyNode, String> placeholders = tipPlaceholders( tree );
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
            writer.write( nexusTaxonLabel( node, placeholders.get( node ) ) );
        }
        writer.write( ";" );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( NexusConstants.END );
        writer.write( ForesterUtil.LINE_SEPARATOR );
    }

    /**
     * The label a node is written under in a Nexus file.
     *
     * Delegates to the node's own New Hampshire label, because the Trees block IS New Hampshire: any other rule
     * here produces a file whose TaxLabels and matrix rows name taxa that the tree does not contain. That is not
     * hypothetical -- this method used to carry its own name-first fallback chain, while
     * PhylogenyNode.toNewHampshire prefers a sequence ACCESSION when there is one, so every tree with accessions
     * was written with TaxLabels that did not match its own trees block.
     *
     * The style argument is irrelevant for an external node (support is null there, so neither
     * AS_INTERNAL_NODE_NAMES nor IN_SQUARE_BRACKETS adds anything), and the distance is switched off, so what
     * comes back is the bare sanitized label.
     */
    static String nexusTaxonLabel( final PhylogenyNode node, final String placeholder ) {
        return node.toNewHampshire( false, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE, false, placeholder );
    }

    private static String molecularSequenceOf( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasSequence() ) {
            return null;
        }
        final String seq = node.getNodeData().getSequence().getMolecularSequence();
        return ForesterUtil.isEmpty( seq ) ? null : seq;
    }

    /**
     * Writes the external nodes' molecular sequences as a Nexus Characters block, or nothing when the tree carries
     * none.
     *
     * A Characters block is used rather than a Data block because a Taxa block has already been written: the taxon
     * labels come from there, and NTax in a Characters block's Dimensions is illegal without NEWTAXA, so only NChar
     * is written. (BasicMsa writes a Data block instead, correctly, because it writes no Taxa block.)
     *
     * A Nexus matrix is rectangular, so this writes nothing when the sequences are of unequal length -- unaligned
     * sequences are not a character matrix, and padding them to a common width would state an alignment that does
     * not exist. A comment says so, rather than leaving the caller to wonder where the data went.
     *
     * Tips that carry no sequence get a row of the missing-data symbol, so the matrix covers every taxon in the
     * Taxa block. Only external nodes are written: a Nexus matrix is keyed on taxa, and an internal node is not one.
     */
    public static void writeNexusCharactersBlock( final Writer writer, final Phylogeny tree ) throws IOException {
        if ( ( tree == null ) || tree.isEmpty() ) {
            return;
        }
        final Map<PhylogenyNode, String> placeholders = tipPlaceholders( tree );
        final List<PhylogenyNode> tips = new ArrayList<PhylogenyNode>();
        final List<PhylogenyNode> with_seq = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator it = tree.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            tips.add( node );
            if ( molecularSequenceOf( node ) != null ) {
                with_seq.add( node );
            }
        }
        if ( with_seq.isEmpty() ) {
            return;
        }
        // A matrix is keyed on the taxon label, so two tips sharing one cannot be told apart: the parser
        // treats the second row as an interleaved continuation of the first and hands BOTH tips the two
        // sequences joined together. The Taxa and Trees blocks have always written such a tree (invalid
        // Nexus, but only cosmetically so); a matrix would make it corrupting, so it is not written.
        final java.util.Set<String> labels = new java.util.HashSet<String>();
        for( final PhylogenyNode node : tips ) {
            if ( !labels.add( nexusTaxonLabel( node, placeholders.get( node ) ) ) ) {
                writer.write( "[ Molecular sequences were not written: two or more tips share the taxon "
                        + "label " );
                writer.write( nexusTaxonLabel( node, placeholders.get( node ) ) );
                writer.write( ", and a character matrix keyed on an ambiguous label cannot be read back. ]" );
                writer.write( ForesterUtil.LINE_SEPARATOR );
                return;
            }
        }
        final int nchar = molecularSequenceOf( with_seq.get( 0 ) ).length();
        for( final PhylogenyNode node : with_seq ) {
            if ( molecularSequenceOf( node ).length() != nchar ) {
                writer.write( "[ Molecular sequences were not written: they are of unequal length (" );
                writer.write( String.valueOf( nchar ) );
                writer.write( " vs " );
                writer.write( String.valueOf( molecularSequenceOf( node ).length() ) );
                writer.write( "), so they are not an alignment and cannot form a Nexus character matrix. ]" );
                writer.write( ForesterUtil.LINE_SEPARATOR );
                return;
            }
        }
        // The datatype is a property of the whole matrix, so it is decided by ALL the sequences, not by the
        // first one that guesses non-null. guessMolecularSequenceType looks for residues that only protein
        // has, so a short protein made of nucleotide letters guesses DNA -- and a matrix wrongly declared DNA
        // is read back with every non-nucleotide residue replaced by N. Protein therefore wins any
        // disagreement: calling a nucleotide alignment protein keeps the residues readable, the reverse
        // destroys them.
        //
        // The strip of gap and missing symbols below CANNOT change the answer today, and is kept as a
        // guard rather than as working code: guessMolecularSequenceType tests membership of L/I/E/H/D/Q,
        // then T and U -- all letters -- while the stripped characters are -.?*, a disjoint set, and a
        // sequence that empties under the strip guesses null either way. Measured over 2,000,000 random
        // matrices with and without it: zero differences (the Archaeopteryx.js session derived the same
        // thing independently). It earns its place only if the guesser ever tests a character that is not
        // a letter, so a mutation that removes it is EXPECTED to survive -- that is not a gap in the tests.
        String type_str = "Protein";
        boolean saw_aa = false;
        boolean saw_nt = false;
        boolean rna = false;
        for( final PhylogenyNode node : with_seq ) {
            final String bare = molecularSequenceOf( node ).replaceAll( "[-.?*]", "" );
            if ( bare.length() < 1 ) {
                continue;
            }
            final MolecularSequence.TYPE t = ForesterUtil.guessMolecularSequenceType( bare );
            if ( t == MolecularSequence.TYPE.DNA ) {
                saw_nt = true;
            }
            else if ( t == MolecularSequence.TYPE.RNA ) {
                saw_nt = true;
                rna = true;
            }
            else if ( t != null ) {
                saw_aa = true;
            }
        }
        if ( saw_nt && !saw_aa ) {
            type_str = rna ? "RNA" : "DNA";
        }
        int max = 0;
        for( final PhylogenyNode node : tips ) {
            final int l = nexusTaxonLabel( node, placeholders.get( node ) ).length();
            if ( l > max ) {
                max = l;
            }
        }
        ++max;
        final StringBuilder missing = new StringBuilder( nchar );
        for( int i = 0; i < nchar; ++i ) {
            missing.append( '?' );
        }
        writer.write( NexusConstants.BEGIN_CHARACTERS );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( " " );
        writer.write( NexusConstants.DIMENSIONS );
        writer.write( " " );
        writer.write( NexusConstants.NCHAR );
        writer.write( "=" );
        writer.write( String.valueOf( nchar ) );
        writer.write( ";" );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( " " );
        writer.write( NexusConstants.FORMAT );
        writer.write( " " );
        writer.write( NexusConstants.DATATYPE );
        writer.write( "=" );
        writer.write( type_str );
        writer.write( " Interleave=No Gap=- Missing=?;" );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( " " );
        writer.write( NexusConstants.MATRIX );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        for( final PhylogenyNode node : tips ) {
            final String seq = molecularSequenceOf( node );
            writer.write( "  " );
            writer.write( ForesterUtil
                    .pad( nexusTaxonLabel( node, placeholders.get( node ) ), max, ' ', false ).toString() );
            writer.write( " " );
            writer.write( seq == null ? missing.toString() : seq );
            writer.write( ForesterUtil.LINE_SEPARATOR );
        }
        writer.write( " ;" );
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
                // The name is always written in single quotes, so an apostrophe INSIDE it must be doubled (the
                // Nexus escape) -- unescaped, "Seba's tree" would end the quoted name three characters in and
                // leave a stray quote that swallows the rest of the file for any reader.
                writer.write( "\'" );
                writer.write( phylogeny.getName().replace( "'", "''" ) );
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




