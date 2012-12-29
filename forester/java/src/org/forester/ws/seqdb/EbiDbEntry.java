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

import org.forester.util.ForesterUtil;

public final class EbiDbEntry implements SequenceDatabaseEntry {

    //http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/emb/AAR37336/
    private String _pa;
    private String _de;
    private String _os;
    private String _tax_id;
    private String _symbol;
    private String _provider;

    private EbiDbEntry() {
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        throw new CloneNotSupportedException();
    }

    public static SequenceDatabaseEntry createInstanceFromPlainTextForRefSeq( final List<String> lines ) {
        final EbiDbEntry e = new EbiDbEntry();
        for( final String line : lines ) {
            //  System.out.println( "-" + line );
            if ( line.startsWith( "ACCESSION" ) ) {
                e.setPA( DatabaseTools.extract( line, "ACCESSION" ) );
            }
            else if ( line.startsWith( "DEFINITION" ) ) {
                if ( line.indexOf( "[" ) > 0 ) {
                    e.setDe( DatabaseTools.extract( line, "DEFINITION", "[" ) );
                }
                else {
                    e.setDe( DatabaseTools.extract( line, "DEFINITION" ) );
                }
            }
            else if ( line.startsWith( "SOURCE" ) ) {
                if ( line.indexOf( "(" ) > 0 ) {
                    e.setOs( DatabaseTools.extract( line, "SOURCE", "(" ) );
                }
                else {
                    e.setOs( DatabaseTools.extract( line, "SOURCE" ) );
                }
            }
        }
        return e;
    }

    public static SequenceDatabaseEntry createInstanceFromPlainText( final List<String> lines ) {
        final EbiDbEntry e = new EbiDbEntry();
        for( final String line : lines ) {
            if ( line.startsWith( "PA" ) ) {
                e.setPA( DatabaseTools.extract( line, "PA" ) );
            }
            else if ( line.startsWith( "DE" ) ) {
                // if ( ( line.indexOf( "RecName:" ) > 0 ) && ( line.indexOf( "Full=" ) > 0 ) ) {
                e.setDe( DatabaseTools.extract( line, "DE" ) );
                //}
            }
            //  else if ( line.startsWith( "GN" ) ) {
            //      if ( ( line.indexOf( "Name=" ) > 0 ) ) {
            //          e.setSymbol( extract( line, "Name=", ";" ) );
            //      }
            //  }
            else if ( line.startsWith( "OS" ) ) {
                if ( line.indexOf( "(" ) > 0 ) {
                    e.setOs( DatabaseTools.extract( line, "OS", "(" ) );
                }
                else {
                    e.setOs( DatabaseTools.extract( line, "OS" ) );
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

    @Override
    public String getAccession() {
        return _pa;
    }

    private void setPA( final String pa ) {
        if ( _pa == null ) {
            _pa = pa;
        }
    }

    @Override
    public String getSequenceName() {
        return _de;
    }

    private void setDe( final String rec_name ) {
        if ( _de == null ) {
            _de = rec_name;
        }
    }

    @Override
    public String getTaxonomyScientificName() {
        return _os;
    }

    private void setOs( final String os ) {
        if ( _os == null ) {
            _os = os;
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

    @Override
    public String getSequenceSymbol() {
        return _symbol;
    }

    private void setSymbol( final String symbol ) {
        if ( _symbol == null ) {
            _symbol = symbol;
        }
    }

    @Override
    public boolean isEmpty() {
        return ( ForesterUtil.isEmpty( getAccession() ) && ForesterUtil.isEmpty( getSequenceName() )
                && ForesterUtil.isEmpty( getTaxonomyScientificName() )
                && ForesterUtil.isEmpty( getTaxonomyIdentifier() ) && ForesterUtil.isEmpty( getSequenceSymbol() ) );
    }

    @Override
    public String getProvider() {
        return _provider;
    }

    public void setProvider( final String provider ) {
        _provider = provider;
    }
}
