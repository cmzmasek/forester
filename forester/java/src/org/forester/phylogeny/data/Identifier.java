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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public final class Identifier implements PhylogenyData {

    final private String _value;
    final private String _provider;
    final private String _value_provider;

    public Identifier() {
        _value = "";
        _provider = "";
        _value_provider = "";
    }

    public Identifier( final String value ) {
        _value = value;
        _provider = "";
        _value_provider = value;
    }

    public Identifier( final String value, final String provider ) {
        _value = value;
        _provider = provider;
        if ( provider != null ) {
            _value_provider = value + provider;
        }
        else {
            _value_provider = value;
        }
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( getValue() );
    }

    @Override
    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        if ( !ForesterUtil.isEmpty( getProvider() ) ) {
            sb.append( "[" );
            sb.append( getProvider() );
            sb.append( "] " );
        }
        sb.append( getValue() );
        return sb;
    }

    public String getValuePlusProvider() {
        return _value_provider;
    }

    /**
     * Provider labels that all denote the SAME identifier namespace: the NCBI taxonomy. The UniProt taxonomy IS the
     * NCBI taxonomy -- a UniProt taxon id is an NCBI taxid -- and an id with no provider at all is, in practice,
     * one too. They arrive under different labels only because of where they were read from: "ncbi" from the NCBI
     * services and from an induced taxonomy tree, "uniprot" from UniProt records and from ids parsed out of
     * sequence labels.
     */
    private static final java.util.Set<String> NCBI_TAXONOMY_PROVIDERS = new java.util.HashSet<String>(
            java.util.Arrays.asList( "ncbi", "uniprot", "ncbitaxon", "ncbi_taxonomy", "uniprot.taxonomy" ) );

    /** The provider label reduced to the NAMESPACE it names, so equivalent labels compare equal (see
     *  {@link #NCBI_TAXONOMY_PROVIDERS}). Anything else is kept, lower-cased, so a genuinely different namespace
     *  (a GTDB id, say) still cannot be confused with an NCBI taxid. */
    public static String normalizedProvider( final String provider ) {
        if ( ( provider == null ) || provider.trim().isEmpty() ) {
            return "ncbi";
        }
        final String p = provider.trim().toLowerCase( java.util.Locale.ROOT );
        return NCBI_TAXONOMY_PROVIDERS.contains( p ) ? "ncbi" : p;
    }

    /**
     * The identifier as a comparison KEY: its value plus its {@link #normalizedProvider(String) normalized}
     * provider. Use this, not {@link #getValuePlusProvider()}, whenever two trees' identifiers are matched --
     * otherwise a gene tree annotated from UniProt can never map onto a species tree built from NCBI taxonomy,
     * even though every id is identical.
     */
    public String getValuePlusNormalizedProvider() {
        return _value + normalizedProvider( _provider );
    }

    @Override
    public PhylogenyData copy() {
        return new Identifier( getValue(), getProvider() );
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            return false;
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return isEqual( ( Identifier ) o );
        }
    }

    public String getProvider() {
        return _provider;
    }

    public String getValue() {
        return _value;
    }

    @Override
    public int hashCode() {
        return _value_provider.hashCode();
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        if ( this == data ) {
            return true;
        }
        if ( ( data == null ) || ( getValue() == null ) ) {
            return false;
        }
        final Identifier a = ( Identifier ) data;
        if ( ( getProvider() != null ) && ( a.getProvider() != null ) ) {
            return ( a.getValue().equals( getValue() ) && a.getProvider().equals( getProvider() ) );
        }
        return ( a.getValue().equals( getValue() ) );
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( !org.forester.util.ForesterUtil.isEmpty( getProvider() ) ) {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.IDENTIFIER,
                                             getValue(),
                                             PhyloXmlMapping.IDENTIFIER_PROVIDER_ATTR,
                                             getProvider(),
                                             indentation );
        }
        else {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.IDENTIFIER, getValue(), indentation );
        }
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
