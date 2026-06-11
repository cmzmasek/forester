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

package org.forester.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class BasicTableParser {

    private final static String START_OF_COMMENT_LINE_DEFAULT = "#";

    private BasicTableParser() {
    }

    /**
     * Splits a (already trimmed) line into its column values. A space delimiter
     * denotes whitespace-aligned columns, so runs of spaces collapse into a single
     * separator. For structured delimiters (tab, comma, semicolon, ...) empty fields
     * are preserved, so an empty cell (two adjacent delimiters) becomes an empty column.
     */
    private static String[] splitLine(final String line, final char column_delimiter) {
        if (column_delimiter == ' ') {
            return line.split(" +", -1);
        }
        return line.split(column_delimiter + "", -1);
    }

    public static BasicTable<String> parse(final Object source, final char column_delimiter) throws IOException {
        return BasicTableParser.parse(source, column_delimiter, false, false, START_OF_COMMENT_LINE_DEFAULT, false)
                .get(0);
    }

    public static BasicTable<String> parse(final Object source,
                                           final char column_delimiter,
                                           final boolean use_first_separator_only,
                                           final boolean use_last_separator_only) throws IOException {
        return BasicTableParser.parse(source,
                column_delimiter,
                use_first_separator_only,
                use_last_separator_only,
                START_OF_COMMENT_LINE_DEFAULT,
                false).get(0);
    }

    public static List<BasicTable<String>> parse(final Object source,
                                                 final char column_delimiter,
                                                 final boolean use_first_separator_only,
                                                 final boolean use_last_separator_only,
                                                 final String start_of_comment_line,
                                                 final boolean tables_separated_by_single_string_line)
            throws IOException {
        if (use_first_separator_only && use_last_separator_only) {
            throw new IllegalArgumentException();
        }
        final BufferedReader reader = ForesterUtil.obtainReader(source);
        final List<BasicTable<String>> tables = new ArrayList<BasicTable<String>>();
        BasicTable<String> table = new BasicTable<String>();
        int row = 0;
        String line;
        boolean saw_first_table = false;
        final boolean use_start_of_comment_line = !(ForesterUtil.isEmpty(start_of_comment_line));
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (!ForesterUtil.isEmpty(line)
                    && (((line.charAt(0) == '"') && (line.charAt(line.length() - 1) == '"') && (ForesterUtil
                    .countChars(line, '"') == 2)) || ((line.charAt(0) == '\'')
                    && (line.charAt(line.length() - 1) == '\'') && (ForesterUtil.countChars(line, '\'') == 2)))) {
                line = line.substring(1, line.length() - 1).trim();
            }
            if (saw_first_table
                    && (ForesterUtil.isEmpty(line) || (tables_separated_by_single_string_line && (line
                    .indexOf(column_delimiter) < 0)))) {
                if (!table.isEmpty()) {
                    tables.add(table);
                }
                table = new BasicTable<String>();
                row = 0;
            } else if (!ForesterUtil.isEmpty(line)
                    && (!use_start_of_comment_line || !line.startsWith(start_of_comment_line))) {
                saw_first_table = true;
                if (use_last_separator_only) {
                    final String e[] = splitLine(line, column_delimiter);
                    final StringBuffer rest = new StringBuffer();
                    for (int i = 0; i < (e.length - 1); ++i) {
                        rest.append(e[i].trim());
                    }
                    table.setValue(0, row, rest.toString());
                    table.setValue(1, row, e[e.length - 1]);
                } else {
                    final String e[] = splitLine(line, column_delimiter);
                    if (e.length > 0) {
                        table.setValue(0, row, e[0].trim());
                    }
                    if (use_first_separator_only) {
                        final StringBuffer rest = new StringBuffer();
                        for (int i = 1; i < e.length; ++i) {
                            rest.append(e[i].trim());
                        }
                        table.setValue(1, row, rest.toString());
                    } else {
                        for (int i = 1; i < e.length; ++i) {
                            table.setValue(i, row, e[i].trim());
                        }
                    }
                }
                ++row;
            }
        }
        if (!table.isEmpty()) {
            tables.add(table);
        }
        reader.close();
        return tables;
    }
}
