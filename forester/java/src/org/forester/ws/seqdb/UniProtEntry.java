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

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.go.BasicGoTerm;
import org.forester.go.GoNameSpace;
import org.forester.go.GoTerm;
import org.forester.util.ForesterUtil;

public final class UniProtEntry implements SequenceDatabaseEntry {

    public final static Pattern GO_PATTERN = Pattern.compile( "GO;\\s+(GO:\\d+);\\s+([PF]):([^;]+);" );
    private String              _ac;
    private String              _name;
    private String              _symbol;
    private String              _gene_name;
    private String              _os_scientific_name;
    private String              _tax_id;
    private List<GoTerm> _go_terms;

    private UniProtEntry() {
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        throw new CloneNotSupportedException();
    }

    public static SequenceDatabaseEntry createInstanceFromPlainText( final List<String> lines ) {
        final UniProtEntry e = new UniProtEntry();
        for( final String line : lines ) {
            //System.out.println( line );
            if ( line.startsWith( "AC" ) ) {
                e.setAc( DatabaseTools.extract( line, "AC", ";" ) );
            }
            else if ( line.startsWith( "DE" ) && ForesterUtil.isEmpty( e.getSequenceName() ) ) {
                if ( ( line.indexOf( "RecName:" ) > 0 ) && ( line.indexOf( "Full=" ) > 0 ) ) {
                    e.setSequenceName( DatabaseTools.extract( line, "Full=", ";" ) );
                }
                else if ( ( line.indexOf( "SubName:" ) > 0 ) && ( line.indexOf( "Full=" ) > 0 ) ) {
                    e.setSequenceName( DatabaseTools.extract( line, "Full=", ";" ) );
                }
            }
            else if ( line.startsWith( "DE" ) && ForesterUtil.isEmpty( e.getSequenceSymbol() ) ) {
                if ( line.indexOf( "Short=" ) > 0 ) {
                    e.setSequenceSymbol( DatabaseTools.extract( line, "Short=", ";" ) );
                }
            }
            else if ( line.startsWith( "GN" ) && ForesterUtil.isEmpty( e.getGeneName() ) ) {
                if ( line.indexOf( "Name=" ) > 0 ) {
                    e.setGeneName( DatabaseTools.extract( line, "Name=", ";" ) );
                }
            }
            else if ( line.startsWith( "DR" ) ) {
                if ( line.indexOf( "GO;" ) > 0 ) {
                    Matcher m = GO_PATTERN.matcher( line );
                    if ( m.find() ) {
                        String id = m.group( 1 );
                        String ns_str = m.group( 2 );
                        String desc = m.group( 3 );
                        String gns = GoNameSpace.BIOLOGICAL_PROCESS_STR;
                        if ( ns_str.equals( "F" ) ) { 
                            gns =  GoNameSpace.MOLECULAR_FUNCTION_STR;
                        }    
                        
                        System.out.println( "GO:" + id + " " + desc + " " + ns_str );
                      
                        e.addGoTerm( new BasicGoTerm( id, desc, gns, false ) ); 
                    }
                }
            }
            else if ( line.startsWith( "OS" ) ) {
                if ( line.indexOf( "(" ) > 0 ) {
                    e.setOsScientificName( DatabaseTools.extract( line, "OS", "(" ) );
                }
                else {
                    e.setOsScientificName( DatabaseTools.extract( line, "OS", "." ) );
                }
            }
            else if ( line.startsWith( "OX" ) ) {
                if ( line.indexOf( "NCBI_TaxID=" ) > 0 ) {
                    e.setTaxId( DatabaseTools.extract( line, "NCBI_TaxID=", ";" ) );
                }
            }
        }
        return e;
    }

    private void addGoTerm( BasicGoTerm g ) {
        if ( _go_terms == null ) {
            _go_terms = new ArrayList<GoTerm>();
        }
        _go_terms.add( g );
        
    }

    private void setSequenceSymbol( String symbol ) {
        _symbol = symbol;
    }

    @Override
    public String getAccession() {
        return _ac;
    }

    private void setAc( final String ac ) {
        if ( _ac == null ) {
            _ac = ac;
        }
    }

    @Override
    public String getSequenceName() {
        return _name;
    }

    private void setSequenceName( final String name ) {
        if ( _name == null ) {
            _name = name;
        }
    }

    @Override
    public String getTaxonomyScientificName() {
        return _os_scientific_name;
    }

    private void setOsScientificName( final String os_scientific_name ) {
        if ( _os_scientific_name == null ) {
            _os_scientific_name = os_scientific_name;
        }
    }

    @Override
    public String getTaxonomyIdentifier() {
        return _tax_id;
    }

    private void setTaxId( final String tax_id ) {
        if ( _tax_id == null ) {
            _tax_id = tax_id;
        }
    }

    private void setGeneName( final String gene_name ) {
        if ( _gene_name == null ) {
            _gene_name = gene_name;
        }
    }
    
    @Override
    public List<GoTerm> getGoTerms() {
        return _go_terms;
    }
    

    @Override
    public String getSequenceSymbol() {
        return _symbol;
    }

    @Override
    public boolean isEmpty() {
        return ( ForesterUtil.isEmpty( getAccession() ) && ForesterUtil.isEmpty( getSequenceName() )
                && ForesterUtil.isEmpty( getTaxonomyScientificName() ) && ForesterUtil.isEmpty( getSequenceSymbol() )
                && ForesterUtil.isEmpty( getGeneName() ) && ForesterUtil.isEmpty( getTaxonomyIdentifier() ) && ForesterUtil
                .isEmpty( getSequenceSymbol() ) && ( getGoTerms() == null || getGoTerms().isEmpty() ) );
    }

    @Override
    public String getProvider() {
        return "uniprot";
    }

    @Override
    public String getGeneName() {
        return _gene_name;
    }
}
