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

import java.util.Locale;

/**
 * A snapshot of the persistent taxonomy cache for display in the Settings dialog: whether the on-disk
 * cache is usable, whether the user has it enabled, its file location, size, entry count, and the age
 * of its oldest entry. Immutable. Produced by {@code NcbiTaxonomyLineageService.getCacheStatus()}.
 *
 * <p>When {@link #isAvailable()} is false the cache could not be read/written (e.g. a read-only home
 * directory); Archaeopteryx still resolves taxa, just without remembering them between sessions, and
 * {@link #getUnavailableReason()} explains why.
 */
public final class TaxonomyCacheStatus {

    private static final long ONE_DAY_MS = 24L * 60L * 60L * 1000L;

    private final boolean     _available;
    private final boolean     _enabled;
    private final String      _path;
    private final long        _bytes;
    private final int         _entries;
    private final long        _oldest_epoch_ms;
    private final String      _unavailable_reason;

    public TaxonomyCacheStatus( final boolean available,
                                final boolean enabled,
                                final String path,
                                final long bytes,
                                final int entries,
                                final long oldest_epoch_ms,
                                final String unavailable_reason ) {
        _available = available;
        _enabled = enabled;
        _path = path;
        _bytes = bytes;
        _entries = entries;
        _oldest_epoch_ms = oldest_epoch_ms;
        _unavailable_reason = unavailable_reason;
    }

    /** True if the cache file/directory can be read and written. */
    public boolean isAvailable() {
        return _available;
    }

    /** True if the user has the persistent cache switched on. */
    public boolean isEnabled() {
        return _enabled;
    }

    /** The absolute path of the cache file (shown to the user even when it does not yet exist). */
    public String getPath() {
        return _path;
    }

    public long getBytes() {
        return _bytes;
    }

    /** Number of distinct, non-expired taxa currently cached on disk. */
    public int getEntries() {
        return _entries;
    }

    /** Epoch-millis of the oldest non-expired entry, or 0 if the cache is empty. */
    public long getOldestEpochMs() {
        return _oldest_epoch_ms;
    }

    /** Why the cache is unavailable, or {@code null} when it is available. */
    public String getUnavailableReason() {
        return _unavailable_reason;
    }

    /** A human-readable file size, e.g. "0 B", "42 KB", "1.4 MB". Pure. */
    public static String formatBytes( final long bytes ) {
        if ( bytes < 1024L ) {
            return bytes + " B";
        }
        if ( bytes < ( 1024L * 1024L ) ) {
            return ( bytes / 1024L ) + " KB";
        }
        return String.format( Locale.ROOT, "%.1f MB", bytes / ( 1024.0 * 1024.0 ) );
    }

    /**
     * A human-readable description of the oldest entry's age at {@code now_ms}, e.g. "Oldest entry:
     * today" / "Oldest entry: 1 day ago" / "Oldest entry: 12 days ago"; empty when the cache is empty
     * or the timestamp is in the future. Pure (time is passed in) so it is unit-testable.
     */
    public static String describeAge( final long oldest_epoch_ms, final long now_ms ) {
        if ( oldest_epoch_ms <= 0L ) {
            return "";
        }
        final long days = ( now_ms - oldest_epoch_ms ) / ONE_DAY_MS;
        if ( days <= 0L ) {
            return "Oldest entry: today";
        }
        if ( days == 1L ) {
            return "Oldest entry: 1 day ago";
        }
        return "Oldest entry: " + days + " days ago";
    }
}
