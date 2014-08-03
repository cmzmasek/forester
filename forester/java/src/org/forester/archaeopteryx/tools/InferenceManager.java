
package org.forester.archaeopteryx.tools;

import java.io.File;

import org.forester.archaeopteryx.Configuration;

public final class InferenceManager {

    private final static String DEFAULT_PATHS[] = { "C:\\Program Files\\", "C:\\Program Files (x86)\\", "/bin/",
            "/usr/local/bin/", "/usr/bin/"     };
    private final File          _path_to_local_mafft;
    private final File          _path_to_local_fastme;
    private final File          _path_to_local_raxml;

    public static InferenceManager createInstance( final Configuration c ) {
        return new InferenceManager( c.getPathToLocalMafft(), c.getPathToLocalFastme(), c.getPathToLocalRaxml() );
    }

    public boolean canDoMsa() {
        return ( getPathToLocalMafft() != null );
    }

    public File getPathToLocalMafft() {
        return _path_to_local_mafft;
    }

    public File getPathToLocalFastme() {
        return _path_to_local_fastme;
    }

    public File getPathToLocalRaxml() {
        return _path_to_local_raxml;
    }

    private final static File createLocalPath( final File path, final String name ) {
        if ( ( path != null ) && path.canExecute() && !path.isDirectory() ) {
            return path;
        }
        final File p1 = new File( name );
        if ( p1.canExecute() && !p1.isDirectory() ) {
            return p1;
        }
        for( final String path_str : DEFAULT_PATHS ) {
            try {
                final File p2 = new File( path_str + name );
                if ( p2.canExecute() && !p2.isDirectory() ) {
                    return p2;
                }
                final File p3 = new File( path_str + name + ".exe" );
                if ( p3.canExecute() && !p3.isDirectory() ) {
                    return p3;
                }
                final File p4 = new File( path_str + name + ".bat" );
                if ( p4.canExecute() && !p4.isDirectory() ) {
                    return p4;
                }
            }
            catch ( final Exception e ) {
            }
        }
        return null;
    }

    private InferenceManager( final File path_to_local_mafft,
                              final File path_to_local_fastme,
                              final File path_to_local_raxml ) {
        _path_to_local_mafft = createLocalPath( path_to_local_mafft, "mafft" );
        _path_to_local_fastme = createLocalPath( path_to_local_fastme, "fastme" );
        _path_to_local_raxml = createLocalPath( path_to_local_raxml, "raxml" );
    }
}
