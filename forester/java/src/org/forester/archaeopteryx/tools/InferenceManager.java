
package org.forester.archaeopteryx.tools;

import java.io.File;
import java.io.IOException;

import org.forester.archaeopteryx.Configuration;
import org.forester.util.ForesterUtil;

public final class InferenceManager {

    private final File _path_to_local_mafft;
    private final File _path_to_local_kalign;
    private final File _path_to_local_fastme;
    private final File _path_to_local_raxml;
    private final File _path_to_local_clustalo;

    public static InferenceManager createInstance( final Configuration c ) {
        return new InferenceManager( c.getpathToLocalMafft(),
                                     c.getPathToLocalKalign(),
                                     c.getPathToLocalFastme(),
                                     c.getPathToLocalRaxml(),
                                     c.getPathToLocalClustalOmega() );
    }

    public boolean canDoMsa() {
        return ( getPathToLocalMafft() != null ) || ( getPathToLocalKalign() != null )
                || ( getPathToLocalClustalo() != null );
    }

    public File getPathToLocalMafft() {
        return _path_to_local_mafft;
    }

    public File getPathToLocalKalign() {
        return _path_to_local_kalign;
    }

    public File getPathToLocalFastme() {
        return _path_to_local_fastme;
    }

    public File getPathToLocalRaxml() {
        return _path_to_local_raxml;
    }

    public File getPathToLocalClustalo() {
        return _path_to_local_clustalo;
    }

    private final static File createLocalPath( final File path ) {
        if ( path == null ) {
            return null;
        }
        try {
            if ( path.getCanonicalFile().canExecute() && !path.getCanonicalFile().isDirectory() ) {
                return new File( path.getCanonicalFile().toString() );
            }
        }
        catch ( final IOException e ) {
            return null;
        }
        return null;
    }

    private InferenceManager( final File path_to_local_mafft,
                              final File path_to_local_kalign,
                              final File path_to_local_fastme,
                              final File path_to_local_raxml,
                              final File path_to_local_clustalo ) {
        _path_to_local_mafft = createLocalPath( path_to_local_mafft ) != null ? createLocalPath( path_to_local_mafft )
                : createLocalPath( new File( "mafft" ) );
        _path_to_local_kalign = createLocalPath( path_to_local_kalign ) != null ? createLocalPath( path_to_local_kalign )
                : createLocalPath( new File( "kalign" ) );
        _path_to_local_fastme = createLocalPath( path_to_local_fastme ) != null ? createLocalPath( path_to_local_fastme )
                : createLocalPath( new File( "fastme" ) );
        _path_to_local_raxml = createLocalPath( path_to_local_raxml ) != null ? createLocalPath( path_to_local_raxml )
                : createLocalPath( new File( "raxml" ) );
        _path_to_local_clustalo = createLocalPath( path_to_local_clustalo ) != null ? createLocalPath( path_to_local_clustalo )
                : createLocalPath( new File( ForesterUtil.isWindowns() ? "clustalo.exe" : "clustalo" ) );
    }
}
