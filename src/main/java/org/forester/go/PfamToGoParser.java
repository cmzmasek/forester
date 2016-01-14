
package org.forester.go;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.util.ForesterUtil;

public class PfamToGoParser {

    // Pfam:PF00001 7tm_1 > GO:rhodopsin-like receptor activity ; GO:0001584
    private static final String  PFAM_TO_GO_FORMAT     = "Pfam:\\S+\\s+(\\S+)\\s*>\\s*GO:.+;\\s*(\\S+)";
    private static final Pattern PFAM_TO_GO_PATTERN    = Pattern.compile( PFAM_TO_GO_FORMAT );
    private static final String  PFAMACC_TO_GO_FORMAT  = "Pfam:(\\S+)\\s+\\S+\\s*>\\s*GO:.+;\\s*(\\S+)";
    private static final Pattern PFAMACC_TO_GO_PATTERN = Pattern.compile( PFAMACC_TO_GO_FORMAT );
    private final File           _input_file;
    private int                  _mapping_count;
    private boolean              _use_acc;

    public PfamToGoParser( final File input_file ) {
        _input_file = input_file;
        init();
    }

    private File getInputFile() {
        return _input_file;
    }

    public int getMappingCount() {
        return _mapping_count;
    }

    private void init() {
        setMappingCount( 0 );
        setUseAccessors( false );
    }

    public boolean isUseAccessors() {
        return _use_acc;
    }

    public List<PfamToGoMapping> parse() throws IOException {
        final String error = ForesterUtil.isReadableFile( getInputFile() );
        if ( !ForesterUtil.isEmpty( error ) ) {
            throw new IOException( error );
        }
        final BufferedReader br = new BufferedReader( new FileReader( getInputFile() ) );
        String line;
        final List<PfamToGoMapping> mappings = new ArrayList<PfamToGoMapping>();
        int line_number = 0;
        try {
            while ( ( line = br.readLine() ) != null ) {
                line_number++;
                line = line.trim();
                if ( ( line.length() > 0 ) && !line.startsWith( "!" ) ) {
                    Matcher m = null;
                    if ( isUseAccessors() ) {
                        m = PFAMACC_TO_GO_PATTERN.matcher( line );
                    }
                    else {
                        m = PFAM_TO_GO_PATTERN.matcher( line );
                    }
                    if ( !m.matches() ) {
                        throw new IOException( "unexpected format [\"" + line + "\"]" );
                    }
                    if ( m.groupCount() != 2 ) {
                        throw new IOException( "unexpected format [\"" + line + "\"]" );
                    }
                    final String pfam = m.group( 1 );
                    final String go = m.group( 2 );
                    if ( ForesterUtil.isEmpty( pfam ) || ForesterUtil.isEmpty( go ) ) {
                        throw new IOException( "unexpected format [\"" + line + "\"]" );
                    }
                    final PfamToGoMapping map = new PfamToGoMapping( pfam, new GoId( go ) );
                    ++_mapping_count;
                    mappings.add( map );
                }
            } // while ( ( line = br.readLine() ) != null )
        }
        catch ( final Exception e ) {
            throw new IOException( "parsing problem: " + e.getMessage() + " [at line " + line_number + "]" );
        }
        return mappings;
    }

    private void setMappingCount( final int mapping_count ) {
        _mapping_count = mapping_count;
    }

    public void setUseAccessors( final boolean use_ids ) {
        _use_acc = use_ids;
    }
}
