// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2017 Christian M. Zmasek
// Copyright (C) 2017 J. Craig Venter Institute
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
// Contact: phyloxml @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.forester.clade_analysis.AnalysisMulti;
import org.forester.clade_analysis.AnalysisSingle;
import org.forester.clade_analysis.Prefix;
import org.forester.clade_analysis.ResultMulti;
import org.forester.clade_analysis.ResultSingle;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public final class cladinator {

    final static private String        PRG_NAME      = "cladinator";
    final static private String        PRG_VERSION   = "0.100";
    final static private String        PRG_DATE      = "170823";
    final static private String        PRG_DESC      = "clades within clades -- analysis of pplacer type outputs";
    final static private String        E_MAIL        = "phyloxml@gmail.com";
    final static private String        WWW           = "https://sites.google.com/site/cmzmasek/home/software/forester";
    final static private String        HELP_OPTION_1 = "help";
    final static private String        HELP_OPTION_2 = "h";
    final static private String        SEP_OPTION    = "s";
    private final static DecimalFormat df2           = new DecimalFormat( "0.0#" );

    public static void main( final String args[] ) {
        try {
            ForesterUtil.printProgramInformation( PRG_NAME,
                                                  PRG_DESC,
                                                  PRG_VERSION,
                                                  PRG_DATE,
                                                  E_MAIL,
                                                  WWW,
                                                  ForesterUtil.getForesterLibraryInformation() );
            CommandLineArguments cla = null;
            try {
                cla = new CommandLineArguments( args );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
            }
            if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) ) {
                System.out.println();
                print_help();
                System.exit( 0 );
            }
            else if ( ( args.length != 2 && args.length != 3 ) ) {
                System.out.println();
                System.out.println( "Wrong number of arguments." );
                System.out.println();
                print_help();
                System.exit( -1 );
            }
            final List<String> allowed_options = new ArrayList<String>();
            allowed_options.add( SEP_OPTION );
            final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
            if ( dissallowed_options.length() > 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
            }
            final String separator;
            if ( cla.isOptionSet( SEP_OPTION ) ) {
                separator = cla.getOptionValue( SEP_OPTION );
            }
            else {
                separator = null;
            }
            final File intreefile = cla.getFile( 0 );
            final String query = cla.getName( 1 );
            System.out.println( "Input tree: " + intreefile );
            System.out.println( "Query     : " + query );
            if ( !ForesterUtil.isEmpty( separator ) ) {
                System.out.println( "Separator : " + separator );
            }
            else {
                System.out.println( "Separator : none" );
            }
            Phylogeny p = null;
            try {
                final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
                final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intreefile, true );
                p = factory.create( intreefile, pp )[ 0 ];
            }
            catch ( final IOException e ) {
                System.out.println( "\nCould not read \"" + intreefile + "\" [" + e.getMessage() + "]\n" );
                System.exit( -1 );
            }
          
            final ResultMulti res = AnalysisMulti.execute( p, query, separator, 0.5 );
            
            System.out.println();
            System.out.println( "Result:" );
            System.out.println( "Query                        : " + query );
            
            ///////////////////
            
         
         
            System.out.println( "Collapsed:" );
          
              for( final Prefix prefix : res.getCollapsedMultiHitPrefixes() ) {
                  System.out.println( prefix );
              }
              if ( _has_specifics ) {
                 
                  System.out.println( "Specifics:" );
                 
                  for( final Prefix prefix : _cleaned_spec ) {
                      System.out.println( prefix );
                     
                  }
                  
                  System.out.println( "Collapsed With Specifics:" );
                 
                  for( final Prefix prefix : _collapsed ) {
                      System.out.println( prefix );
                      
                      for( final Prefix spec : _cleaned_spec ) {
                          if ( spec.getPrefix().startsWith( prefix.getPrefix() ) ) {
                              System.out.println( "    " + spec );
                             
                          }
                      }
                  }
              }
              if ( !ForesterUtil.isEmpty( _all_down ) ) {
                  
                  System.out.println( "Collapsed Down:" );
                  
                  for( final Prefix prefix : _collapsed_down ) {
                      System.out.println( prefix );
                      
                  }
              
              }
              if ( !ForesterUtil.isEmpty( _all_up ) ) {
                  
           
                  System.out.println( "Collapsed Up:" );
                 
                  for( final Prefix prefix : _collapsed_up ) {
                      System.out.println( prefix );
                     
                  }
             
              }
            
            ///////////////////
            
            
            System.out.print( "Greatest Common Prefix       : " + res.getGreatestCommonPrefix() );
            if ( !ForesterUtil.isEmpty( res.getGreatestCommonPrefix() )
                    && !ForesterUtil.isEmpty( res.getGreatestCommonCladeSubtreeConfidence() ) ) {
                System.out.println( "\t(" + res.getGreatestCommonCladeSubtreeConfidence() + ")" );
            }
            else {
                System.out.println();
            }
            System.out.print( "Greatest Common Prefix Up    : " + res.getGreatestCommonPrefixUp() );
            if ( !ForesterUtil.isEmpty( res.getGreatestCommonPrefixUp() )
                    && !ForesterUtil.isEmpty( res.getGreatestCommonCladeUpSubtreeConfidence() ) ) {
                System.out.println( "\t(" + res.getGreatestCommonCladeUpSubtreeConfidence() + ")" );
            }
            else {
                System.out.println();
            }
            System.out.print( "Greatest Common Prefix Down  : " + res.getGreatestCommonPrefixDown() );
            if ( !ForesterUtil.isEmpty( res.getGreatestCommonPrefixDown() )
                    && !ForesterUtil.isEmpty( res.getGreatestCommonCladeDownSubtreeConfidence() ) ) {
                System.out.println( "\t(" + res.getGreatestCommonCladeDownSubtreeConfidence() + ")" );
            }
            else {
                System.out.println();
            }
            System.out.println( "Least Encompassing Clade size: " + res.getLeastEncompassingCladeSize()
                    + " external nodes" );
            final double lec_ratio = ( 100.0 * res.getLeastEncompassingCladeSize() ) / res.getTreeSize();
            System.out.println( "Least Encompassing Clade size: " + df2.format( lec_ratio ) + "%" );
            System.out.println( "Total tree size              : " + res.getTreeSize() + " external nodes" );
            if ( res.getWarnings().size() > 0 ) {
                System.out.println( "Warnings:" );
                for( final String s : res.getWarnings() ) {
                    System.out.println( s );
                }
            }
            System.out.println();
        }
        catch ( final IllegalArgumentException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            ForesterUtil.fatalError( PRG_NAME, "Unexpected errror!" );
        }
    }

    private final static void print_help() {
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + " [options] <gene tree file> <query>" );
        System.out.println();
        System.out.println( " options:" );
        System.out.println( "  -" + SEP_OPTION + "=<separator>: the separator to be used" );
        System.out.println();
        System.out.println( "Example:" );
        System.out.println();
        System.out.println( " " + PRG_NAME + " -s=. my_tree.xml A.1.1.1" );
        System.out.println();
    }
}
