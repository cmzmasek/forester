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
import java.util.ArrayList;
import java.util.List;

import org.forester.util.ForesterUtil;

public class PropertiesList implements PhylogenyData {

    // INSERTION ORDER is the contract: a node's properties stay in the order they were added, which for every
    // parser and importer is the order the SOURCE states -- the phyloXML document, the imported table's columns,
    // the NHX/JSON token order.
    //
    // This list used to re-sort itself by ref on every add. That alphabetized a node's properties on read and wrote
    // them back out alphabetized, so opening a phyloXML whose properties are grouped meaningfully (a gene
    // presence/absence matrix ordered core / accessory / mobile, say) and saving it silently rewrote the author's
    // column order -- and it disagreed with Archaeopteryx.js, whose reader and writer both preserve document order.
    // Order is data; it is kept. (The sort itself long predates the intent documented for it: it arrived in 0.9.6
    // beta with its comparator arguments swapped, sorting DESCENDING unnoticed for years, and 0.11.111 only
    // corrected the direction.)
    //
    // Nothing downstream needs alphabetical order: every surface that wants one sorts for itself -- the Annotation
    // Fields inventory and Tree Facts by display name, the Color-by candidates by PropertyColorScheme.VIS_ORDER --
    // and every other consumer looks a property up by ref.
    private final List<Property> _properties;

    public PropertiesList() {
        _properties = new ArrayList<Property>();
    }

    public int size() {
        return _properties.size();
    }

    /** Appends a property, keeping insertion (source) order; several properties may share one ref. */
    public void addProperty(final Property property) throws IllegalArgumentException {
        _properties.add(property);
    }

    @Override
    public StringBuffer asSimpleText() {
        final StringBuffer sb = new StringBuffer();
        boolean first = true;
        for (final Property p : getProperties()) {
            if (first) {
                first = false;
            } else {
                sb.append("\n");
            }
            sb.append(p.asText());
        }
        return sb;
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        final PropertiesList new_one = new PropertiesList();
        for (final Property r : getProperties()) {
            new_one.addProperty(r);
        }
        return new_one;
    }

    public List<Property> getProperties() {
        return _properties;
    }

    public List<Property> getPropertiesWithGivenReferencePrefix(final String ref_prefix)
            throws IllegalArgumentException {
        if (ForesterUtil.isEmpty(ref_prefix)) {
            throw new IllegalArgumentException("reference prefix is null or empty");
        }
        final String my_ref_prefix = new String(ref_prefix.trim());
        final List<Property> props = new ArrayList<Property>();
        for (final Property p : getProperties()) {
            if (p.getRef().startsWith(my_ref_prefix)) {
                props.add(p);
            }
        }
        return props;
    }

    public List<Property> getPropertiesWithGivenRef(final String ref)
            throws IllegalArgumentException {
        if (ForesterUtil.isEmpty(ref)) {
            throw new IllegalArgumentException("reference is null or empty");
        }
        final String my_ref = new String(ref.trim());
        final List<Property> props = new ArrayList<Property>();
        for (final Property p : getProperties()) {
            if (p.getRef().equals(my_ref)) {
                props.add(p);
            }
        }
        return props;
    }

    public List<Property> getProperties(final String ref) throws IllegalArgumentException {
        final List<Property> props = new ArrayList<Property>();
        for (final Property p : getProperties()) {
            if (p.getRef().equals(ref)) {
                props.add(p);
            }
        }
        return props;
    }

    @Override
    public boolean isEqual(final PhylogenyData data) {
        throw new UnsupportedOperationException();
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML(final Writer writer, final int level, final String indentation) throws IOException {
        if (getProperties() != null) {
            for (final Property p : getProperties()) {
                p.toPhyloXML(writer, level, indentation);
            }
        }
    }

    @Override
    public String toString() {
        return asSimpleText().toString();
    }
}
