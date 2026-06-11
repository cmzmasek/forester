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

package org.forester.ws.seqdb;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class UniprotData {

    final String _id;
    final String _entry;
    final String _entry_name;
    final String _status;
    final String _protein_names;
    final String _gene_names;
    final String _organism;
    final int _length;
    private final static Pattern FASTA = Pattern.compile(">.+\\|(.+)\\|.*?\\s+(.+?)\\s+OS=");

    public UniprotData(final String line, final String src, final boolean from_fasta) throws IOException {
        if (from_fasta) {
            final Matcher m = FASTA.matcher(line);
            if (m.find()) {
                _id = m.group(1);
                _protein_names = m.group(2);
            } else {
                if (line.startsWith(">")) {
                    System.out.println("Could not extract protein names from: " + line);
                }
                _id = "";
                _protein_names = "";
            }
            _entry = "";
            _status = "";
            _entry_name = "";
            _gene_names = "";
            _organism = "";
            _length = 0;
        } else {
            final String[] split_line = line.split("\t");
            if (split_line.length != 8) {
                throw new IOException("source: [" + src + "]: line has illegal format: \"" + line + "\" (expected 8 elements, got " + split_line.length + ")");
            }
            _id = split_line[0];
            _entry = split_line[1];
            _entry_name = split_line[2];
            _status = split_line[3];
            _protein_names = split_line[4];
            _gene_names = split_line[5];
            _organism = split_line[6];
            try {
                _length = Integer.parseInt(split_line[7]);
            } catch (final NumberFormatException e) {
                throw new IOException("source: [" + src + "]: could not parse length from " + split_line[8]);
            }
        }
    }

    public String getId() {
        return _id;
    }

    public String getEntry() {
        return _entry;
    }

    public String getEntryName() {
        return _entry_name;
    }

    public String getStatus() {
        return _status;
    }

    public String getProteinNames() {
        return _protein_names;
    }

    public String getGeneNames() {
        return _gene_names;
    }

    public String getOrganism() {
        return _organism;
    }

    public int getLength() {
        return _length;
    }
}