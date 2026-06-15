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

package org.forester.archaeopteryx;

import java.io.File;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.util.ForesterUtil;

public final class Archaeopteryx {

    public static MainFrame createApplication( final Phylogeny phylogeny ) {
        final Phylogeny[] phylogenies = new Phylogeny[ 1 ];
        phylogenies[ 0 ] = phylogeny;
        return createApplication( phylogenies );
    }

    public static MainFrame createApplication( final Phylogeny phylogeny, final Configuration config, final String title ) {
        final Phylogeny[] phylogenies = new Phylogeny[ 1 ];
        phylogenies[ 0 ] = phylogeny;
        return MainFrameApplication.createInstance( phylogenies, config, title );
    }

    public static MainFrame createApplication( final Phylogeny[] phylogenies ) {
        return MainFrameApplication.createInstance( phylogenies, new Configuration(), "" );
    }

    public static void main( final String args[] ) {
        // macOS integration: must be set before any AWT/Swing class is initialized.
        if ( ForesterUtil.isMac() ) {
            System.setProperty( "apple.awt.application.name", AptxConstants.PRG_NAME );
            System.setProperty( "apple.awt.application.appearance", "system" );
            // Render the menu bar inside the window so it is styled by FlatLaf; the native
            // macOS screen menu bar cannot be themed and would not match the rest of the UI.
            System.setProperty( "apple.laf.useScreenMenuBar", "false" );
        }
        // Register the bundled default figure font BEFORE the default font family is chosen (in
        // Configuration); must come after the macOS AWT properties above (this initializes AWT).
        FontResources.registerBundledFonts();
        // Configuration files are no longer supported: the -c option and its parser have been removed.
        if ( ( args.length > 0 ) && AptxUtil.isConfigFileOption( args[ 0 ] ) ) {
            ForesterUtil.fatalError( AptxConstants.PRG_NAME,
                                     "configuration files are no longer supported: the -c option has been removed. "
                                             + "Archaeopteryx now uses its built-in defaults; change settings at runtime via the Settings dialog." );
        }
        Phylogeny[] phylogenies = null;
        final Configuration conf = new Configuration();
        File f = null;
        try {
            int filename_index = 0;
            if ( args.length > 0 ) {
                if ( args[ 0 ].startsWith( "-open" ) ) {
                    filename_index += 1;
                }
                if ( args.length > filename_index ) {
                    f = new File( args[ filename_index ] );
                    final String err = ForesterUtil.isReadableFile( f );
                    if ( !ForesterUtil.isEmpty( err ) ) {
                        ForesterUtil.fatalError( AptxConstants.PRG_NAME, err );
                    }
                    boolean nhx_or_nexus = false;
                    final PhylogenyParser p = ParserUtils.createParserDependingOnFileType( f, conf
                                                                                           .isValidatePhyloXmlAgainstSchema() );
                    if ( p instanceof NHXParser ) {
                        nhx_or_nexus = true;
                        final NHXParser nhx = ( NHXParser ) p;
                        nhx.setReplaceUnderscores( conf.isReplaceUnderscoresInNhParsing() );
                        nhx.setIgnoreQuotes( false );
                        nhx.setTaxonomyExtraction( conf.getTaxonomyExtraction() );
                    }
                    else if ( p instanceof NexusPhylogeniesParser ) {
                        nhx_or_nexus = true;
                        final NexusPhylogeniesParser nex = ( NexusPhylogeniesParser ) p;
                        nex.setReplaceUnderscores( conf.isReplaceUnderscoresInNhParsing() );
                        nex.setIgnoreQuotes( false );
                    }
                    phylogenies = PhylogenyMethods.readPhylogenies( p, f );
                    if ( nhx_or_nexus && conf.isInternalNumberAreConfidenceForNhParsing() ) {
                        for( final Phylogeny phy : phylogenies ) {
                            PhylogenyMethods.transferInternalNodeNamesToConfidence( phy, "" );
                        }
                    }
                }
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( AptxConstants.PRG_NAME, "failed to start: " + e.getLocalizedMessage() );
        }
        String title = "";
        if ( f != null ) {
            title = f.getName();
        }
        File current_dir = null;
        if ( ( phylogenies != null ) && ( phylogenies.length > 0 ) ) {
            current_dir = new File( "." );
        }
        try {
            MainFrameApplication.createInstance( phylogenies, conf, title, current_dir );
        }
        catch ( final OutOfMemoryError e ) {
            AptxUtil.outOfMemoryError( e );
        }
        catch ( final Exception e ) {
            AptxUtil.unexpectedException( e );
        }
        catch ( final Error e ) {
            AptxUtil.unexpectedError( e );
        }
    }
}