// $Id:
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
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

package org.forester.ws.seqdb;

import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.go.GoTerm;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.sequence.MolecularSequence;
import org.forester.util.ForesterUtil;

public final class EbiDbEntry implements SequenceDatabaseEntry {
    
    private final static boolean  DEBUG = false;
    
    private final static Pattern LETTERS_PATTERN = Pattern.compile( "^[A-Z]+" );
    private final static Pattern chromosome_PATTERN = Pattern.compile( "\\s+/chromosome=\"(\\w+)\"" );
    private final static Pattern map_PATTERN = Pattern.compile( "\\s+/map=\"([\\w\\.]+)\"" );
    private final static Pattern gene_PATTERN = Pattern.compile( "\\s+/gene=\"(.+)\"" );
    private final static Pattern mim_PATTERN = Pattern.compile( "\\s+/db_xref=\"MIM:(\\d+)\"" );
    private final static Pattern taxon_PATTERN = Pattern.compile( "\\s+/db_xref=\"taxon:(\\d+)\"" );
    private final static Pattern interpro_PATTERN = Pattern.compile( "\\s+/db_xref=\"InterPro:([A-Z0-9]+)\"" );
    private final static Pattern uniprot_PATTERN = Pattern.compile( "\\s+/db_xref=\"UniProtKB/[A-Za-z-]*:(\\w+)\"" );
    private final static Pattern hgnc_PATTERN = Pattern.compile( "\\s+/db_xref=\"[A-Z:]*HGNC:(\\d+)\"" );
    private final static Pattern geneid_PATTERN = Pattern.compile( "\\s+/db_xref=\"GeneID:(\\d+)\"" );
    private final static Pattern pdb_PATTERN = Pattern.compile( "\\s+/db_xref=\"PDB:([A-Z0-9]+)\"" );
    private final static Pattern ec_PATTERN = Pattern.compile( "\\s+/EC_number=\"([\\.\\-\\d]+)\"" );
    private final static Pattern product_PATTERN = Pattern.compile( "\\s+/product=\"(\\w{1,10})\"" );
 
    private SortedSet<Annotation> _annotations;
    private String                _chromosome;
    private SortedSet<Accession>  _cross_references;
    private String                _de;
    private String                _gene_name;
    private String                _map;
    private String                _os;
    private String                _pa;
    private String                _provider;
    private String                _symbol;
    private String                _tax_id;
  
    private EbiDbEntry() {
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        throw new CloneNotSupportedException();
    }

    @Override
    public String getAccession() {
        return _pa;
    }

    @Override
    public SortedSet<Annotation> getAnnotations() {
        return _annotations;
    }

    @Override
    public String getChromosome() {
        return _chromosome;
    }

    @Override
    public SortedSet<Accession> getCrossReferences() {
        return _cross_references;
    }

    @Override
    public String getGeneName() {
        return _gene_name;
    }

    @Override
    public SortedSet<GoTerm> getGoTerms() {
        return null;
    }

    @Override
    public String getMap() {
        return _map;
    }

    @Override
    public String getProvider() {
        return _provider;
    }

    @Override
    public String getSequenceName() {
        return _de;
    }

    @Override
    public String getSequenceSymbol() {
        return _symbol;
    }

    @Override
    public String getTaxonomyIdentifier() {
        return _tax_id;
    }

    @Override
    public String getTaxonomyScientificName() {
        return _os;
    }

    @Override
    public boolean isEmpty() {
        return ( ForesterUtil.isEmpty( getAccession() ) && ForesterUtil.isEmpty( getSequenceName() )
                && ForesterUtil.isEmpty( getTaxonomyScientificName() )
                && ForesterUtil.isEmpty( getTaxonomyIdentifier() ) && ForesterUtil.isEmpty( getSequenceSymbol() ) );
    }


    @Override
    public MolecularSequence getMolecularSequence() {
        // TODO Auto-generated method stub
        return null;
    }
    private void addAnnotation( final Annotation annotation ) {
        if ( _annotations == null ) {
            _annotations = new TreeSet<Annotation>();
        }
        _annotations.add( annotation );
    }

    private void addCrossReference( final Accession accession ) {
        if ( _cross_references == null ) {
            _cross_references = new TreeSet<Accession>();
        }
        if ( DEBUG ) {
            System.out.println( "XREF ADDED: " + accession );
        }
        _cross_references.add( accession );
    }

    private void setAccession( final String pa ) {
        if ( _pa == null ) {
            _pa = pa;
        }
    }

    private void setChromosome( final String chromosome ) {
        _chromosome = chromosome;
    }

    private void setGeneName( final String gene_name ) {
        if ( _gene_name == null ) {
            _gene_name = gene_name;
        }
    }

    private void setMap( final String map ) {
        _map = map;
    }

    private void setSequenceName( final String rec_name ) {
        if ( _de == null ) {
            _de = rec_name;
        }
    }

    private void setSequenceSymbol( final String symbol ) {
        _symbol = symbol;
    }

    private void setTaxId( final String tax_id ) {
        if ( _tax_id == null ) {
            _tax_id = tax_id;
        }
    }

    private void setTaxonomyScientificName( final String os ) {
        if ( _os == null ) {
            _os = os;
        }
    }
   
    private static void append( final StringBuilder sb, final String s ) {
        if ( sb.length() > 0 ) {
            sb.append( " " );
        }
        sb.append( s.trim() );
    }
    
    public final static SequenceDatabaseEntry createInstance( final List<String> lines ) {
        
        final EbiDbEntry e = new EbiDbEntry();
        final StringBuilder def = new StringBuilder();
        boolean in_definition = false;
        boolean in_features = false;
        boolean in_source = false;
        boolean in_gene = false;
        boolean in_cds = false;
        boolean in_mrna = false;
        boolean in_protein = false;
        for( final String line : lines ) {
            if ( line.startsWith( "ACCESSION " ) ) {
                e.setAccession( SequenceDbWsTools.extractFrom( line, "ACCESSION" ) );
                in_definition = false;
            }
            else if ( line.startsWith( "ID " ) ) {
                e.setAccession( SequenceDbWsTools.extractFromTo( line, "ID", ";" ) );
                in_definition = false;
            }
            else if ( line.startsWith( "DEFINITION " ) || ( line.startsWith( "DE " ) ) ) {
                boolean definiton = false;
                if ( line.startsWith( "DEFINITION " ) ) {
                    definiton = true;
                }
                if ( line.indexOf( "[" ) > 0 ) {
                    if ( definiton ) {
                        append( def, ( SequenceDbWsTools.extractFromTo( line, "DEFINITION", "[" ) ) );
                    }
                    else {
                        append( def, ( SequenceDbWsTools.extractFromTo( line, "DE", "[" ) ) );
                    }
                }
                else if ( line.indexOf( "." ) > 0 ) {
                    if ( definiton ) {
                        append( def, ( SequenceDbWsTools.extractFromTo( line, "DEFINITION", "." ) ) );
                    }
                    else {
                        append( def, ( SequenceDbWsTools.extractFromTo( line, "DE", "." ) ) );
                    }
                }
                else {
                    if ( definiton ) {
                        append( def, ( SequenceDbWsTools.extractFrom( line, "DEFINITION" ) ) );
                    }
                    else {
                        append( def, ( SequenceDbWsTools.extractFrom( line, "DE" ) ) );
                    }
                }
                if ( definiton ) {
                    in_definition = true;
                }
            }
            else if ( line.startsWith( "  ORGANISM " ) ) {
                if ( line.indexOf( "(" ) > 0 ) {
                    e.setTaxonomyScientificName( SequenceDbWsTools.extractFromTo( line, "  ORGANISM", "(" ) );
                }
                else {
                    e.setTaxonomyScientificName( SequenceDbWsTools.extractFrom( line, "  ORGANISM" ) );
                }
            }
            else if ( line.startsWith( "OS " ) ) {
                if ( line.indexOf( "(" ) > 0 ) {
                    e.setTaxonomyScientificName( SequenceDbWsTools.extractFromTo( line, "OS", "(" ) );
                }
                else {
                    e.setTaxonomyScientificName( SequenceDbWsTools.extractFrom( line, "OS" ) );
                }
            }
            else if ( line.startsWith( " " ) && in_definition ) {
                def.append( " " );
                if ( line.indexOf( "[" ) > 0 ) {
                    def.append( SequenceDbWsTools.extractTo( line, "[" ) );
                }
                else if ( line.indexOf( "." ) > 0 ) {
                    def.append( SequenceDbWsTools.extractTo( line, "." ) );
                }
                else {
                    def.append( line.trim() );
                }
            }
            else {
                in_definition = false;
            }
            if ( !line.startsWith( "FT " ) && LETTERS_PATTERN.matcher( line ).find() ) {
                in_features = false;
                in_source = false;
                in_gene = false;
                in_cds = false;
                in_mrna = false;
                in_protein = false;
            }
            if ( line.startsWith( "FEATURES " ) || line.startsWith( "FT " ) ) {
                in_features = true;
            }
            if ( in_features && ( line.startsWith( "     source " ) || line.startsWith( "FT   source " ) ) ) {
                in_source = true;
                in_gene = false;
                in_cds = false;
                in_mrna = false;
                in_protein = false;
            }
            if ( in_features && ( line.startsWith( "     gene " ) || line.startsWith( "FT   gene " ) ) ) {
                in_source = false;
                in_gene = true;
                in_cds = false;
                in_mrna = false;
                in_protein = false;
            }
            if ( in_features && ( line.startsWith( "     CDS " ) || line.startsWith( "FT   CDS " ) ) ) {
                in_source = false;
                in_gene = false;
                in_cds = true;
                in_mrna = false;
                in_protein = false;
            }
            if ( in_features && ( line.startsWith( "     Protein " ) || line.startsWith( "FT   Protein " ) ) ) {
                in_source = false;
                in_gene = false;
                in_cds = false;
                in_mrna = false;
                in_protein = true;
            }
            if ( in_features && ( line.startsWith( "     mRNA " ) || line.startsWith( "FT   mRNA " ) ) ) {
                in_source = false;
                in_gene = false;
                in_cds = false;
                in_mrna = true;
                in_protein = false;
            }
            if ( in_source ) {
                final Matcher ti = taxon_PATTERN.matcher( line );
                if ( ti.find() ) {
                    e.setTaxId( ti.group( 1 ) );
                }
                final Matcher chr = chromosome_PATTERN.matcher( line );
                if ( chr.find() ) {
                    e.setChromosome( chr.group( 1 ) );
                }
                final Matcher map = map_PATTERN.matcher( line );
                if ( map.find() ) {
                    e.setMap( map.group( 1 ) );
                }
            }
            if ( in_cds || in_gene ) {
                final Matcher hgnc = hgnc_PATTERN.matcher( line );
                if ( hgnc.find() ) {
                    e.addCrossReference( new Accession( hgnc.group( 1 ), "hgnc" ) );
                }
                final Matcher geneid = geneid_PATTERN.matcher( line );
                if ( geneid.find() ) {
                    e.addCrossReference( new Accession( geneid.group( 1 ), "geneid" ) );
                }
            }
            if ( in_protein || in_cds || in_gene || in_mrna ) {
                final Matcher ec = ec_PATTERN.matcher( line );
                if ( ec.find() ) {
                    e.addAnnotation( new Annotation( "EC", ec.group( 1 ) ) );
                }
                final Matcher gene = gene_PATTERN.matcher( line );
                if ( gene.find() ) {
                    e.setGeneName( gene.group( 1 ) );
                }
                final Matcher uniprot = uniprot_PATTERN.matcher( line );
                if ( uniprot.find() ) {
                    e.addCrossReference( new Accession( uniprot.group( 1 ), "uniprot" ) );
                }
                final Matcher interpro = interpro_PATTERN.matcher( line );
                if ( interpro.find() ) {
                    e.addCrossReference( new Accession( interpro.group( 1 ), "interpro" ) );
                }
                final Matcher mim = mim_PATTERN.matcher( line );
                if ( mim.find() ) {
                    e.addCrossReference( new Accession( mim.group( 1 ), "mim" ) );
                }
                final Matcher product = product_PATTERN.matcher( line );
                if ( product.find() ) {
                    e.setSequenceSymbol( product.group( 1 ) );
                }
                final Matcher pdb = pdb_PATTERN.matcher( line );
                if ( pdb.find() ) {
                    e.addCrossReference( new Accession( pdb.group( 1 ), "pdb" ) );
                }
            }
        }
        if ( def.length() > 0 ) {
            e.setSequenceName( def.toString().trim() );
        }
        return e;
    }

}
