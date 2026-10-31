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

package org.forester.archaeopteryx;

import java.awt.event.InputEvent;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JComboBox;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;

import org.forester.archaeopteryx.NodeDataDraft.PropertyDraft;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;

/**
 * Headless tests for {@link NodeDataForm} (a plain {@code JPanel}, so it builds without a display): the widgets
 * reflect the node, typing makes the form dirty but touches nothing, invalid input outlines its field and blocks
 * the write, a valid write reaches the node exactly once, sequence cards / confidence rows / property rows come and
 * go, sections fold by content, and the read-only VIEW mode shows only what is there.
 */
public final class NodeDataFormTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeDataForm: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return editOpensClean() && dirtyValidateWrite() && sequenceCards() && confidenceRows()
                    && propertyRows() && sections() && viewMode() && viewLinks() && headerText();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean editOpensClean() throws Exception {
        final PhylogenyNode n = NodeDataDraftTest.richNode();
        final NodeDataForm f = new NodeDataForm( n, null, NodeDataForm.Mode.EDIT );
        return check( "editable", f.isEditable() )
                && eq( "widgets reflect the node", NodeDataDraft.from( n ), f.collect() )
                && check( "clean on open", !f.isDirty() ) && check( "valid on open", f.problems().isEmpty() )
                && check( "name is a text field", f.fieldForTest( NodeDataDraft.NAME ) instanceof JTextField )
                && check( "rank is a combo", f.fieldForTest( NodeDataDraft.TAX_RANK ) instanceof JComboBox )
                && check( "synonyms is an area", f.fieldForTest( NodeDataDraft.TAX_SYNONYMS ) instanceof JTextArea )
                && check( "mol seq is an area", f.fieldForTest( NodeDataDraft.sequenceKey( 0,
                                                                                            NodeDataDraft.SEQ_MOL_SEQ ) ) instanceof JTextArea )
                && check( "aligned is a checkbox", f.fieldForTest( NodeDataDraft.sequenceKey( 0,
                                                                                               NodeDataDraft.SEQ_ALIGNED ) ) instanceof JCheckBox )
                && check( "the page never scrolls sideways",
                          f.scrollPaneForTest().getHorizontalScrollBarPolicy() == JScrollPane.HORIZONTAL_SCROLLBAR_NEVER )
                && check( "a write with nothing to write is a harmless yes", f.write() );
    }

    private static boolean dirtyValidateWrite() throws Exception {
        final PhylogenyNode n = NodeDataDraftTest.richNode();
        final NodeDataForm f = new NodeDataForm( n, null, NodeDataForm.Mode.EDIT );
        final int[] changes = { 0 };
        f.addChangeListener( () -> changes[ 0 ]++ );
        f.setTextForTest( NodeDataDraft.NAME, "renamed" );
        boolean ok = check( "typing makes it dirty", f.isDirty() )
                && check( "listener fired", changes[ 0 ] > 0 )
                && eq( "typing does not reach the node", "BRCA1 clade", n.getName() );
        f.setTextForTest( NodeDataDraft.BRANCH_LENGTH, "abc" );
        ok = ok && eq( "problem reported", NodeDataDraft.BRANCH_LENGTH, f.problems().get( 0 ).key )
                && check( "field outlined", f.isOutlinedForTest( NodeDataDraft.BRANCH_LENGTH ) )
                && check( "name not outlined", !f.isOutlinedForTest( NodeDataDraft.NAME ) )
                && check( "write refused while invalid", !f.write() )
                && eq( "nothing written", "BRCA1 clade", n.getName() )
                && check( "still dirty", f.isDirty() );
        f.setTextForTest( NodeDataDraft.BRANCH_LENGTH, "0.75" );
        f.setTextForTest( NodeDataDraft.TAX_RANK, "Genus" ); // combo, any case
        f.setTextForTest( NodeDataDraft.sequenceKey( 0, NodeDataDraft.SEQ_ALIGNED ), "false" );
        ok = ok && check( "outline cleared once fixed", !f.isOutlinedForTest( NodeDataDraft.BRANCH_LENGTH ) )
                && check( "valid again", f.problems().isEmpty() ) && check( "write succeeds", f.write() )
                && eq( "name written", "renamed", n.getName() )
                && check( "branch length written", n.getDistanceToParent() == 0.75 )
                && eq( "rank written lower-case", "genus", n.getNodeData().getTaxonomy().getRank() )
                && check( "aligned flag written", !n.getNodeData().getSequence( 0 ).isMolecularSequenceAligned() )
                && check( "clean after write", !f.isDirty() )
                && eq( "baseline moved to what was written", f.collect(), f.baseline() )
                && check( "a second write is a no-op yes", f.write() );
        return ok;
    }

    private static boolean sequenceCards() throws Exception {
        final PhylogenyNode n = NodeDataDraftTest.richNode();
        final Sequence second = n.getNodeData().getSequence( 1 );
        final NodeDataForm f = new NodeDataForm( n, null, NodeDataForm.Mode.EDIT );
        boolean ok = eq( "two cards", 2, f.sequenceCardCountForTest() );
        f.addSequenceForTest();
        ok = ok && eq( "three cards", 3, f.sequenceCardCountForTest() )
                && check( "an empty new card is not dirty by itself", !f.isDirty() );
        f.setTextForTest( NodeDataDraft.sequenceKey( 2, NodeDataDraft.SEQ_NAME ), "third" );
        f.setTextForTest( NodeDataDraft.sequenceKey( 2, NodeDataDraft.SEQ_SYMBOL ), "bad symbol" );
        ok = ok && check( "dirty now", f.isDirty() )
                && check( "card field outlined",
                          f.isOutlinedForTest( NodeDataDraft.sequenceKey( 2, NodeDataDraft.SEQ_SYMBOL ) ) )
                && check( "refused", !f.write() );
        f.setTextForTest( NodeDataDraft.sequenceKey( 2, NodeDataDraft.SEQ_SYMBOL ), "THIRD" );
        f.setTextForTest( NodeDataDraft.sequenceKey( 2, NodeDataDraft.SEQ_TYPE ), "dna" );
        ok = ok && check( "written", f.write() ) && eq( "node has 3 sequences", 3,
                                                          n.getNodeData().getSequences().size() )
                && eq( "third name", "third", n.getNodeData().getSequence( 2 ).getName() )
                && eq( "third type", "dna", n.getNodeData().getSequence( 2 ).getType() );
        // editing the new card again must mutate the SAME object the first write created
        final Sequence created = n.getNodeData().getSequence( 2 );
        f.setTextForTest( NodeDataDraft.sequenceKey( 2, NodeDataDraft.SEQ_LOCATION ), "chr1" );
        ok = ok && check( "written again", f.write() )
                && check( "same object", n.getNodeData().getSequence( 2 ) == created )
                && eq( "location", "chr1", created.getLocation() );
        f.removeSequenceForTest( 0 );
        ok = ok && eq( "two cards after remove", 2, f.sequenceCardCountForTest() ) && check( "dirty", f.isDirty() )
                && check( "written", f.write() ) && eq( "node has 2", 2, n.getNodeData().getSequences().size() )
                && check( "the former second is now first, same object", n.getNodeData().getSequence( 0 ) == second );
        return ok;
    }

    private static boolean confidenceRows() throws Exception {
        final PhylogenyNode n = NodeDataDraftTest.richNode();
        final NodeDataForm f = new NodeDataForm( n, null, NodeDataForm.Mode.EDIT );
        boolean ok = eq( "two rows", 2, f.confidenceRowCountForTest() );
        f.addConfidenceForTest();
        f.setTextForTest( NodeDataDraft.confidenceKey( 2, NodeDataDraft.CONF_TYPE ), "posterior" );
        ok = ok && eq( "type without value is a problem", NodeDataDraft.confidenceKey( 2, NodeDataDraft.CONF_VALUE ),
                       f.problems().get( 0 ).key );
        f.setTextForTest( NodeDataDraft.confidenceKey( 2, NodeDataDraft.CONF_VALUE ), "0.5" );
        ok = ok && check( "written", f.write() ) && eq( "three confidences", 3,
                                                          n.getBranchData().getConfidences().size() )
                && eq( "type", "posterior", n.getBranchData().getConfidence( 2 ).getType() );
        f.removeConfidenceForTest( 0 );
        ok = ok && check( "written", f.write() ) && eq( "two left", 2, n.getBranchData().getConfidences().size() )
                && check( "the old second is first", n.getBranchData().getConfidence( 0 ).getValue() == 0.98 );
        return ok;
    }

    private static boolean propertyRows() throws Exception {
        final PhylogenyNode n = NodeDataDraftTest.richNode();
        final NodeDataForm f = new NodeDataForm( n, null, NodeDataForm.Mode.EDIT );
        boolean ok = eq( "two rows", 2, f.propertyTableForTest().getRowCount() )
                && check( "reference cell editable", f.propertyTableForTest().isCellEditable( 0, 0 ) );
        f.addPropertyForTest( new PropertyDraft( "noprefix", "7", "", "xsd:decimal", AppliesTo.NODE ) );
        ok = ok && check( "dirty", f.isDirty() )
                && eq( "problem on the new row's ref", NodeDataDraft.propertyKey( 2, NodeDataDraft.PROP_REF ),
                       f.problems().get( 0 ).key )
                && check( "cell marked", f.isPropertyCellMarkedForTest( 2, NodeDataDraft.PROP_REF ) )
                && check( "value cell not marked", !f.isPropertyCellMarkedForTest( 2, NodeDataDraft.PROP_VALUE ) )
                && check( "refused", !f.write() );
        // Edit through the table's CELL EDITOR, exactly as the UI does: JTable commits by calling setValueAt while
        // its editor is still installed, and the form's change handling must not stop that editor again (it did:
        // a StackOverflowError on the very first property edit).
        final javax.swing.JTable table = f.propertyTableForTest();
        if ( !table.editCellAt( 2, 0 ) ) {
            return check( "the reference cell must be editable through the table", false );
        }
        ( (javax.swing.JTextField) table.getEditorComponent() ).setText( "data:count" );
        table.getCellEditor().stopCellEditing(); // what Enter / focus-lost does
        ok = ok && check( "editor gone after commit", !table.isEditing() )
                && eq( "typed value reached the model", "data:count", table.getValueAt( 2, 0 ) );
        ok = ok && check( "mark cleared", !f.isPropertyCellMarkedForTest( 2, NodeDataDraft.PROP_REF ) )
                && check( "written", f.write() ) && eq( "three properties", 3, n.getNodeData().getProperties().size() )
                && eq( "value", "7", n.getNodeData().getProperties().getProperties( "data:count" ).get( 0 ).getValue() );
        f.propertyTableForTest().setValueAt( "8", 2, 1 );
        f.propertyTableForTest().setValueAt( AppliesTo.CLADE, 2, 4 );
        ok = ok && check( "written", f.write() )
                && eq( "edited value", "8", n.getNodeData().getProperties().getProperties( "data:count" ).get( 0 ).getValue() )
                && eq( "edited applies_to", AppliesTo.CLADE,
                       n.getNodeData().getProperties().getProperties( "data:count" ).get( 0 ).getAppliesTo() );
        f.removePropertyForTest( 0 );
        f.removePropertyForTest( 0 );
        f.removePropertyForTest( 0 );
        ok = ok && check( "written", f.write() ) && check( "no properties left", !n.getNodeData().isHasProperties() );
        return ok;
    }

    private static boolean sections() throws Exception {
        final NodeDataForm rich = new NodeDataForm( NodeDataDraftTest.richNode(), null, NodeDataForm.Mode.EDIT );
        boolean ok = true;
        for( final String s : new String[] { NodeDataDraft.SEC_BASIC, NodeDataDraft.SEC_TAXONOMY,
                NodeDataDraft.SEC_SEQUENCES, NodeDataDraft.SEC_EVENTS, NodeDataDraft.SEC_DATE,
                NodeDataDraft.SEC_DISTRIBUTION, NodeDataDraft.SEC_REFERENCE, NodeDataDraft.SEC_PROPERTIES } ) {
            ok = ok && check( "rich: section " + s + " present", rich.hasSectionForTest( s ) )
                    && check( "rich: section " + s + " open (it has data)", rich.isSectionExpandedForTest( s ) );
        }
        rich.toggleSectionForTest( NodeDataDraft.SEC_DATE );
        ok = ok && check( "toggle folds", !rich.isSectionExpandedForTest( NodeDataDraft.SEC_DATE ) );
        final PhylogenyNode bare = new PhylogenyNode();
        bare.setName( "tip" );
        final NodeDataForm empty = new NodeDataForm( bare, null, NodeDataForm.Mode.EDIT );
        ok = ok && check( "bare: Basic open", empty.isSectionExpandedForTest( NodeDataDraft.SEC_BASIC ) )
                && check( "bare: Taxonomy present but folded", empty.hasSectionForTest( NodeDataDraft.SEC_TAXONOMY )
                        && !empty.isSectionExpandedForTest( NodeDataDraft.SEC_TAXONOMY ) )
                && check( "bare: Sequences folded", !empty.isSectionExpandedForTest( NodeDataDraft.SEC_SEQUENCES ) )
                && check( "bare: Properties folded", !empty.isSectionExpandedForTest( NodeDataDraft.SEC_PROPERTIES ) )
                && check( "a tip has no Events section", !empty.hasSectionForTest( NodeDataDraft.SEC_EVENTS ) );
        return ok;
    }

    private static boolean viewMode() throws Exception {
        final PhylogenyNode n = NodeDataDraftTest.richNode();
        final NodeDataForm v = new NodeDataForm( n, null, NodeDataForm.Mode.VIEW );
        boolean ok = check( "not editable", !v.isEditable() ) && check( "never dirty", !v.isDirty() )
                && check( "write is refused", !v.write() )
                // a VIEW page is display, not a source: nothing is read back, so nothing can be "invalid" (a typed
                // confidence used to be re-read from its "95 (bootstrap)" text and flagged as not a number)
                && check( "view has no problems", v.problems().isEmpty() )
                && check( "view outlines nothing", !v.isOutlinedForTest( NodeDataDraft.confidenceKey( 0,
                                                                                                        NodeDataDraft.CONF_VALUE ) )
                        && !v.isOutlinedForTest( NodeDataDraft.NAME ) )
                && check( "view collect() is the node's draft", v.collect() == v.baseline() )
                && eq( "view draft equals the node", NodeDataDraft.from( n ), v.collect() )
                && check( "name shown read-only", ( v.fieldForTest( NodeDataDraft.NAME ) instanceof JTextField )
                        && !( (JTextField) v.fieldForTest( NodeDataDraft.NAME ) ).isEditable() )
                && eq( "name value", "BRCA1 clade", ( (JTextField) v.fieldForTest( NodeDataDraft.NAME ) ).getText() )
                && check( "rank shown as text, not a combo", v.fieldForTest( NodeDataDraft.TAX_RANK ) instanceof JTextField )
                && check( "properties table read-only", !v.propertyTableForTest().isCellEditable( 0, 0 ) )
                && check( "all data sections present", v.hasSectionForTest( NodeDataDraft.SEC_EVENTS )
                        && v.hasSectionForTest( NodeDataDraft.SEC_PROPERTIES ) );
        final PhylogenyNode bare = new PhylogenyNode();
        bare.setName( "tip" );
        final NodeDataForm vb = new NodeDataForm( bare, null, NodeDataForm.Mode.VIEW );
        ok = ok && check( "bare view: Basic only", vb.hasSectionForTest( NodeDataDraft.SEC_BASIC )
                && !vb.hasSectionForTest( NodeDataDraft.SEC_TAXONOMY ) && !vb.hasSectionForTest( NodeDataDraft.SEC_SEQUENCES )
                && !vb.hasSectionForTest( NodeDataDraft.SEC_DATE ) && !vb.hasSectionForTest( NodeDataDraft.SEC_PROPERTIES ) )
                && check( "bare view: empty fields hidden", vb.fieldForTest( NodeDataDraft.BRANCH_LENGTH ) == null );
        return ok;
    }

    /**
     * VIEW: every value that names a page is a link to it, in page order and nothing else is; a plain left click
     * opens it, a modified or double click does not; a metadata value that is a web address links too; the edit
     * page has no links; a rebind (undo/redo) rebuilds them.
     */
    private static boolean viewLinks() throws Exception {
        final PhylogenyNode n = NodeDataDraftTest.richNode();
        final Sequence s1 = n.getNodeData().getSequences().get( 0 );
        s1.addCrossReference( new Accession( "IPR002552", "interpro" ) );
        s1.addCrossReference( new Accession( "fig|1000.9.peg.1", "BV-BRC" ) ); // a source we do not know
        final DomainArchitecture da = new DomainArchitecture();
        da.addDomain( new ProteinDomain( "NB-ARC", 1, 10, "PF00931" ) );
        da.addDomain( new ProteinDomain( "WD40", 20, 60 ) ); // a name only: a search, as in the tree
        da.setTotalLength( 100 );
        s1.setDomainArchitecture( da );
        n.getNodeData().getSequences().get( 1 ).setAccession( new Accession( "fig|1000.9.peg.1", "BV-BRC" ) );
        final Sequence s3 = new Sequence(); // a shape only the declared type can place: the card passes the type on
        s3.setType( "rna" );
        s3.setAccession( new Accession( "XYZ12345678", "ncbi" ) );
        n.getNodeData().addSequence( s3 );
        n.getNodeData().getProperties().addProperty( new Property( "data:page", "https://example.org/sample/7", "",
                                                                   "xsd:anyURI", AppliesTo.NODE ) );
        final List<String> want = Arrays.asList( "Homo sapiens", NCBI_TAX + "9606", //
                                                 "HUMAN", "https://www.uniprot.org/taxonomy?query=mnemonic:HUMAN", //
                                                 "9606", NCBI_TAX + "9606", //
                                                 "https://example.org/tax/9606", "https://example.org/tax/9606", //
                                                 "P38398", "https://www.uniprot.org/uniprot/P38398", //
                                                 "https://example.org/seq/P38398", "https://example.org/seq/P38398", //
                                                 "interpro: IPR002552",
                                                 "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR002552/", //
                                                 "NB-ARC (1-10)", "https://www.ebi.ac.uk/interpro/entry/pfam/PF00931/", //
                                                 "WD40 (20-60)", "https://www.ebi.ac.uk/interpro/search/text/WD40/", //
                                                 "XYZ12345678", "https://www.ncbi.nlm.nih.gov/nuccore/XYZ12345678", //
                                                 "10.1126/science.7545954", "https://doi.org/10.1126/science.7545954" );
        final NodeDataForm v = new NodeDataForm( n, null, NodeDataForm.Mode.VIEW );
        final List<String> opened = new ArrayList<>();
        v.setLinkOpenerForTest( opened::add );
        boolean ok = eq( "the links, in page order", want, flatten( v.linksForTest() ) );
        final JComponent bv = v.fieldForTest( NodeDataDraft.sequenceKey( 1, NodeDataDraft.SEQ_ACC ) );
        ok &= check( "an accession of an unknown source is shown, as plain text",
                     ( bv instanceof JTextField ) && "fig|1000.9.peg.1".equals( ( (JTextField) bv ).getText() )
                             && ( bv.getClientProperty( FormWidgets.LINK_URL_PROPERTY ) == null ) );
        ok &= check( "a plain value is not a link", !v.clickLinkForTest( "Linnaeus, 1758", 0, 1 ) );
        ok &= check( "a plain left click opens", v.clickLinkForTest( "HUMAN", 0, 1 ) )
                && eq( "what it opened", Arrays.asList( "https://www.uniprot.org/taxonomy?query=mnemonic:HUMAN" ), opened );
        opened.clear();
        v.clickLinkForTest( "9606", InputEvent.META_DOWN_MASK, 1 );
        v.clickLinkForTest( "P38398", 0, 2 );
        ok &= eq( "a modified click and the second click of a double open nothing", new ArrayList<String>(), opened );
        ok &= eq( "a metadata value that is an address links", "https://example.org/sample/7", v.propertyLinkForTest( 2 ) )
                && eq( "a metadata value that is not, does not", null, v.propertyLinkForTest( 0 ) );
        v.clickPropertyValueForTest( 0, 0, 1 );
        ok &= eq( "a click on a plain metadata value opens nothing", new ArrayList<String>(), opened );
        v.clickPropertyValueForTest( 2, InputEvent.META_DOWN_MASK, 1 );
        v.clickPropertyValueForTest( 2, 0, 2 );
        ok &= eq( "a modified or double click on a linked metadata value opens nothing", new ArrayList<String>(),
                  opened );
        v.clickPropertyValueForTest( 2, 0, 1 );
        ok &= eq( "a click on a linked metadata value opens it", Arrays.asList( "https://example.org/sample/7" ), opened );
        ok &= linkSurface( n );
        // a rebind (undo/redo) to a node with OTHER data: the links are the new node's, not kept from the old page
        final PhylogenyNode n2 = NodeDataDraftTest.richNode();
        n2.getNodeData().getTaxonomy().setIdentifier( new org.forester.phylogeny.data.Identifier( "10090", "ncbi" ) );
        n2.getNodeData().getReference().setDoi( "10.1038/nature12373" );
        v.rebind( n2 );
        final List<String> after = flatten( v.linksForTest() );
        ok &= check( "a rebind to other data: the new id and DOI link: " + after,
                     after.contains( NCBI_TAX + "10090" ) && after.contains( "https://doi.org/10.1038/nature12373" )
                             && !after.contains( NCBI_TAX + "9606" )
                             && !after.contains( "https://doi.org/10.1126/science.7545954" ) );
        final NodeDataForm e = new NodeDataForm( n, null, NodeDataForm.Mode.EDIT );
        ok &= eq( "the edit page has no links", new ArrayList<String>(), flatten( e.linksForTest() ) )
                && eq( "nor in its metadata table", null, e.propertyLinkForTest( 2 ) );
        return ok;
    }

    /**
     * A link is its TEXT, not the row it is stretched across: the colour, the URL tooltip; a hand cursor and an
     * underline only over the text; a click beside the text opens nothing. The colour survives a look-and-feel
     * update (a theme switch). The metadata table: its columns cannot be dragged (the cell code reads view columns
     * as model columns), and its link cell is painted as one.
     */
    private static boolean linkSurface( final PhylogenyNode n ) {
        // a look and feel supplies its link colour as a UIResource -- the kind an update replaces; without FlatLaf
        // (a standalone run) there is none, so supply one the same way
        final boolean supplied = javax.swing.UIManager.getColor( "Component.linkColor" ) == null;
        if ( supplied ) {
            javax.swing.UIManager.put( "Component.linkColor", new javax.swing.plaf.ColorUIResource( 0x2675BF ) );
        }
        try {
            return linkSurface( new NodeDataForm( n, null, NodeDataForm.Mode.VIEW ) );
        }
        finally {
            if ( supplied ) {
                javax.swing.UIManager.put( "Component.linkColor", null );
            }
        }
    }

    private static boolean linkSurface( final NodeDataForm v ) {
        final List<String> opened = new ArrayList<>();
        v.setLinkOpenerForTest( opened::add );
        final javax.swing.text.JTextComponent link = v.linkForTest( "9606" );
        final java.awt.Point on = NodeDataForm.onTextForTest( link );
        final java.awt.Point beside = new java.awt.Point( link.getWidth() - 4, on.y );
        boolean ok = check( "precondition: the value is stretched well beyond its text (" + link.getWidth() + " px)",
                            FormWidgets.isOverText( link, on ) && !FormWidgets.isOverText( link, beside )
                                    && ( link.getWidth() > ( on.x + 100 ) ) );
        ok &= check( "precondition: below the text is not on it either",
                     !FormWidgets.isOverText( link, new java.awt.Point( on.x, link.getHeight() + 50 ) ) );
        NodeDataForm.mouseForTest( link, java.awt.event.MouseEvent.MOUSE_CLICKED, beside, 0, 1 );
        ok &= eq( "a click beside the text opens nothing", new ArrayList<String>(), opened );
        ok &= check( "the link colour and the address as tooltip",
                     link.getForeground().getRGB() == FormWidgets.linkColor().getRGB()
                             && ( NCBI_TAX + "9606" ).equals( link.getToolTipText() ) );
        final java.awt.Cursor base = link.getCursor();
        NodeDataForm.mouseForTest( link, java.awt.event.MouseEvent.MOUSE_MOVED, on, 0, 0 );
        ok &= check( "over the text: a hand, underlined",
                     ( link.getCursor().getType() == java.awt.Cursor.HAND_CURSOR )
                             && java.awt.font.TextAttribute.UNDERLINE_ON
                                     .equals( link.getFont().getAttributes().get( java.awt.font.TextAttribute.UNDERLINE ) ) );
        NodeDataForm.mouseForTest( link, java.awt.event.MouseEvent.MOUSE_MOVED, beside, 0, 0 );
        ok &= check( "beside it: neither", ( link.getCursor() == base ) && ( link.getFont().getAttributes()
                .get( java.awt.font.TextAttribute.UNDERLINE ) == null ) );
        final javax.swing.text.JTextComponent plain = (javax.swing.text.JTextComponent) v
                .fieldForTest( NodeDataDraft.TAX_COMMON_NAME );
        // a theme switch: the look and feel's link colour CHANGES (light 0x2675BF vs dark 0x579BF6), then every
        // component is updated; the link must take the NEW theme's colour, not keep the old one or lose it
        final int before = FormWidgets.linkColor().getRGB();
        javax.swing.UIManager.put( "Component.linkColor", new javax.swing.plaf.ColorUIResource( 0x579BF7 ) );
        try {
            javax.swing.SwingUtilities.updateComponentTreeUI( v );
            ok &= check( "precondition: the switch changed the link colour", FormWidgets.linkColor().getRGB() != before )
                    && check( "after a theme switch a link wears the new theme's link colour",
                              link.getForeground().getRGB() == FormWidgets.linkColor().getRGB() )
                    && check( "control: a plain value takes the text colour, not the link's",
                              plain.getForeground().getRGB() != FormWidgets.linkColor().getRGB() );
        }
        finally {
            javax.swing.UIManager.put( "Component.linkColor", null ); // the look and feel's own (or the caller's) again
        }
        final javax.swing.JTable t = v.propertyTableForTest();
        ok &= check( "the metadata columns cannot be dragged", !t.getTableHeader().getReorderingAllowed() );
        // one renderer component serves every cell: read each cell's look before preparing the next
        final java.awt.Component cell = t.prepareRenderer( t.getCellRenderer( 2, 1 ), 2, 1 );
        ok &= check( "the linked metadata value is painted as a link, with its address",
                     ( cell.getForeground().getRGB() == FormWidgets.linkColor().getRGB() )
                             && "https://example.org/sample/7".equals( ( (JComponent) cell ).getToolTipText() ) );
        final java.awt.Component ref = t.prepareRenderer( t.getCellRenderer( 2, 0 ), 2, 0 );
        ok &= check( "its reference is not", ( ref.getForeground().getRGB() != FormWidgets.linkColor().getRGB() )
                && ( ( (JComponent) ref ).getToolTipText() == null ) );
        return ok;
    }

    private static final String NCBI_TAX = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=";

    private static List<String> flatten( final List<String[]> pairs ) {
        final List<String> out = new ArrayList<>();
        for( final String[] p : pairs ) {
            out.add( p[ 0 ] );
            out.add( p[ 1 ] );
        }
        return out;
    }

    private static boolean headerText() throws Exception {
        final PhylogenyNode n = NodeDataDraftTest.richNode();
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( n );
        n.setDistanceToParent( 0.25 );
        boolean ok = eq( "label: name", "BRCA1 clade", NodeDataForm.nodeLabel( n ) )
                && eq( "subtitle", "Internal node · 2 children · 2 tips · depth 1 · 0.25 from root",
                       NodeDataForm.headerSubtitle( n ) )
                && eq( "root subtitle", "Root · 1 child · 2 tips", NodeDataForm.headerSubtitle( root ) );
        final PhylogenyNode tip = n.getChildNode( 0 );
        tip.setName( "" );
        ok = ok && eq( "tip subtitle", "External node · depth 2 · 0.25 from root", NodeDataForm.headerSubtitle( tip ) )
                && eq( "label falls back to the id", "node " + tip.getId(), NodeDataForm.nodeLabel( tip ) )
                && eq( "data label: the name", "BRCA1 clade", NodeDataForm.nodeDataLabel( n ) )
                // a node without a label of its own has no data label (never the session-only id)
                && eq( "no data label", null, NodeDataForm.nodeDataLabel( tip ) );
        n.setName( "" );
        ok = ok && eq( "label falls back to the scientific name", "Homo sapiens", NodeDataForm.nodeLabel( n ) )
                && eq( "data label falls back to the scientific name", "Homo sapiens", NodeDataForm.nodeDataLabel( n ) );
        n.getNodeData().getTaxonomy().setScientificName( "" );
        // root-free (a declared-unrooted tree shown unrooted): r <- (i <- (t1, t2), t3, t4)
        final PhylogenyNode r = new PhylogenyNode();
        final PhylogenyNode i = new PhylogenyNode();
        final PhylogenyNode t1 = new PhylogenyNode();
        i.addAsChild( t1 );
        i.addAsChild( new PhylogenyNode() );
        r.addAsChild( i );
        r.addAsChild( new PhylogenyNode() );
        r.addAsChild( new PhylogenyNode() );
        ok = ok && eq( "root-free internal", "Internal node · 3 neighbours · tips around 1 · 1 · 2",
                       NodeDataForm.headerSubtitle( i, true ) )
                && eq( "root-free stored root is not called Root", "Internal node · 3 neighbours · tips around 1 · 1 · 2",
                       NodeDataForm.headerSubtitle( r, true ) )
                && eq( "root-free tip", "External node", NodeDataForm.headerSubtitle( t1, true ) )
                && eq( "the one-argument form is the rooted view", NodeDataForm.headerSubtitle( i, false ),
                       NodeDataForm.headerSubtitle( i ) );
        return ok && eq( "then the taxonomy code", "HUMAN", NodeDataForm.nodeLabel( n ) );
    }

    private static boolean eq( final String what, final Object expected, final Object actual ) {
        if ( ( expected == null ) ? ( actual == null ) : expected.equals( actual ) ) {
            return true;
        }
        System.out.println( "  [NodeDataFormTest] " + what + ": expected <" + expected + "> but got <" + actual + ">" );
        return false;
    }

    private static boolean check( final String what, final boolean condition ) {
        if ( condition ) {
            return true;
        }
        System.out.println( "  [NodeDataFormTest] " + what );
        return false;
    }

    private NodeDataFormTest() {
        // not instantiable
    }
}
