/*
# bootstrap_cz
# ------------
# Copyright (C) 1999-2002 Washington University School of Medicine
# and Howard Hughes Medical Institute
# All rights reserved
#
# Author: Christian M. Zmasek 
#         zmasek@genetics.wustl.edu
#         http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Created: 06/06/01
#
# Last modified: 01/27/02
#
# Purpose:
# Bootstrap resamples an alignment in PHYLIP sequential format <bootstraps> times.
# Bootstrapping is not done randomly but according to a BSP (bootstrap 
# positions) file.
# The BSP file can be created with the Perl program "bootstrap_cz.pl"
# in mode 0.
# This prgram has the same functionality as "bootstrap_cz.pl" in mode 1.
# Sequence names are normalized to LENGTH_OF_NAME characters.
# The output alignment is in PHYLIP's sequential or interleaved format.
# (These two are the same in this case, since all the seqs will be one
# line in length (no returns in seq).)
# 
# Usage: bootstrap_cz <bootstraps> <alignment infile> <BSP file> <outfile>
#        [number of processors]
*/



#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


#define LENGTH_OF_NAME 26


static char **names,       /* This stores the sequence names */
            **sequences;   /* This stores the sequences */
static int  number_of_seqs,       
            number_of_colm;
               

void readInAlignmnet( const char * );
void bootstrapAccordingToBSPfile( int, const char *, const char * );
void checkForMemAllocFailure( void * );
int  fileExists( const char *);
void errorInCommandLine();





/* Reads the seqs and seq-names from inalignment     */
/* into **sequences and **sequences.                 */
/* Inalignment must be in PHYLIP sequential format.  */
/* Last modified: 06/25/01                           */
void readInAlignment( const char *inalignment ) {
 
    FILE *inalignment_fp              = NULL;
    char *str                         = NULL;
    int  max_length                   = 0; 
    register char c                   = ' ';   
    register int  i                   = 0,
                  ii                  = 0,
                  z                   = 0,
                  seq                 = 0;

    number_of_seqs = 0;
    number_of_colm = 0;
         
    inalignment_fp = fopen( inalignment, "r" );
    if ( inalignment_fp == NULL ) {
        printf( "\nbootstrap_cz: Error: Could not open alignment file for reading.\n" );
        exit( -1 );
    }

    if ( fscanf( inalignment_fp, "%d", &number_of_seqs ) != 1 ) {
        printf( "\nbootstrap_cz: Error: Could not read in number of seqs.\n" );
        exit( -1 );
    }
    if ( fscanf( inalignment_fp, "%d", &number_of_colm ) != 1 ) {
        printf( "\nbootstrap_cz: Error: Could not read in number of columns.\n" );
        exit( -1 );
    }

    names = malloc( number_of_seqs * sizeof( char *) );
    checkForMemAllocFailure( names );
    for ( i = 0; i < number_of_seqs; ++i ) {
        names[ i ] = malloc( LENGTH_OF_NAME * sizeof( char ) );
        checkForMemAllocFailure( names[ i ] );
    }

    sequences = malloc( number_of_seqs * sizeof( char * ) );
    checkForMemAllocFailure( sequences );
    for ( i = 0; i < number_of_seqs; ++i ) {
        sequences[ i ] = malloc( number_of_colm * sizeof( char ) );
        checkForMemAllocFailure( sequences[ i ] );
    }

    max_length = ( 30 * LENGTH_OF_NAME ) + number_of_colm;

    str = malloc( max_length * sizeof( char * ) );
    checkForMemAllocFailure( str );

    while ( fgets( str, max_length, inalignment_fp ) != NULL ) {

        if ( !isspace( str[ 0 ] ) != 0 ) {

            i = 0;
            while ( str[ i ] != ' ' ) {
                names[ seq ][ i ] = str[ i ];
                i++;
            }

            ii = i;
            while ( ii < LENGTH_OF_NAME ) {
                names[ seq ][ ii ] = ' ';
                ii++;
            }
           
            z = 0;
            
            while ( str[ i ] != '\n' && str[ i ] != '\r' && str[ i ] != '\0' ) {
                c = str[ i ];
                if ( c != ' ' ) {
                    if ( isupper( c ) != 0 || c == '-' ) {
                        sequences[ seq ][ z++ ] = c;
                    }
                    else {
                        printf( "\nbootstrap_cz: Error: Sequence must be represented by uppercase letters A-Z and \"-\" only.\n" );
                        exit( -1 );
                    }
                }
                i++;
                if ( z > number_of_colm ) {
                    printf( "\nbootstrap_cz: Error: line in \"%s\" contains more than %d columns.\n",
                            inalignment, number_of_colm );
                    exit( -1 );
                }
            }
            if ( z != number_of_colm ) {
                printf( "\nbootstrap_cz: Error: line in \"%s\" contains a incorrect number of columns.\n",
                        inalignment );
                exit( -1 );
            }

            seq++;
  
            if ( seq > number_of_seqs ) {
                printf( "\nbootstrap_cz: Error: \"%s\" contains more than %d seqs.\n",
                         inalignment, number_of_seqs );
                exit( -1 );
            }          
        }
        
       
    } /* while ( fgets ) */

    if ( seq != number_of_seqs ) {
        printf( "\nbootstrap_cz: Error: \"%s\" contains a incorrect number of seqs.\n",
                inalignment );
        exit( -1 );
    }  

    fclose( inalignment_fp ); 

    return;

} /* readInAlignment */



/* Rearrenges the aa in sequences according to    */
/* the bsp (bootstrap positions) file bsp_file.   */
/* Writes the results to outfile                  */
/* Last modified: 06/07/01                        */
void bootstrapAccordingToBSPfile( int bootstraps,
                                  const char *bsp_file,
                                  const char *outfile ) {
  
    FILE          *bsp_file_fp = NULL,
                  *outfile_fp  = NULL;
    int           *positions   = NULL,
                  p            = 0;
    register int  boot         = 0,
                  seq          = 0, 
                  i            = 0; 

    positions = malloc( number_of_colm * sizeof( int ) );
    checkForMemAllocFailure( positions );

         
    bsp_file_fp = fopen( bsp_file, "r" );
    if ( bsp_file_fp == NULL ) {
        printf( "\nbootstrap_cz: Error: could not open file \"%s\" for reading.\n", 
                bsp_file  );
        exit( -1 );
    }

    outfile_fp = fopen( outfile, "w" );
    if ( outfile_fp == NULL ) {
        printf( "\nbootstrap_cz: Error: could not open file \"%s\" for writing.\n",
                outfile );
        exit( -1 );
    }

    for ( boot = 0; boot < bootstraps; ++boot ) {
        
        for ( i = 0; i < number_of_colm; ++i ) {
            if ( fscanf( bsp_file_fp, "%d", &p ) != 1 ) { 
                printf( "\nbootstrap_cz: Error: file \"%s\" does not correspond to alignment.\n",
                        bsp_file );
                exit( -1 );
            }
            positions[ i ] = p;
        }

        fprintf( outfile_fp, " %d  %d\n", number_of_seqs, number_of_colm );
        for ( seq = 0; seq < number_of_seqs; ++seq ) {
            for ( i = 0; i < LENGTH_OF_NAME; ++i ) {
                fprintf( outfile_fp, "%c", names[ seq ][ i ] );   
            }
            for ( i = 0; i < number_of_colm; ++i ) {
                fprintf( outfile_fp, "%c", sequences[ seq ][ positions[ i ] ] );  
            }
            fprintf( outfile_fp, "\n" );
        }
    }

    /* Now, the bsp file must not contain any more numbers */
    if ( fscanf( bsp_file_fp, "%d", &p ) == 1 ) { 
        printf( "\nbootstrap_cz: Error: file \"%s\" does not correspond to alignment (too long).\n",
                bsp_file );
        printf( ">%d<\n", p );
        printf( "number of seqs=%d\n", number_of_seqs );
        exit( -1 );
    }

    fclose( bsp_file_fp ); 
    fclose( outfile_fp ); 
    free( positions );
    return;

} /* bootstrapAccordingToBSPfile */



/* Rearrenges the aa in sequences according to    */
/* the bsp (bootstrap positions) file bsp_file.   */
/* Writes the results to outfile                  */
/* Last modified: 01/25/02                        */
void bootstrapAccordingToBSPfileP( int bootstraps,
                                   int processors,
                                   const char *bsp_file,
                                   const char *outfile ) {
  
    FILE          *bsp_file_fp  = NULL,
                  *outfile_fp   = NULL;
    int           *positions    = NULL,
                  p             = 0;
    char          *outfile_     = NULL;
    register int  boot          = 0,
                  seq           = 0, 
                  i             = 0,
                  j             = 0,
                  z             = 0,
                  flag          = 0;
    int           block_size    = 0,
                  larger_blocks = 0;

    block_size    = ( int ) bootstraps / processors;
    larger_blocks = bootstraps - ( block_size * processors ); /* number of blocks which have a size of
                                                                 block_size + 1 */

    positions = malloc( number_of_colm * sizeof( int ) );
    checkForMemAllocFailure( positions );

    outfile_ = malloc( ( strlen( outfile ) + 20 ) * sizeof( char ) );
    checkForMemAllocFailure( outfile_ );

    bsp_file_fp = fopen( bsp_file, "r" );
    if ( bsp_file_fp == NULL ) {
        printf( "\nbootstrap_cz: Error: could not open file \"%s\" for reading.\n", 
                bsp_file  );
        exit( -1 );
    }

    j    = -1;
    flag = 1;
    z    = 0;

    for ( boot = 0; boot < bootstraps; ++boot ) {
        
        for ( i = 0; i < number_of_colm; ++i ) {
            if ( fscanf( bsp_file_fp, "%d", &p ) != 1 ) { 
                printf( "\nbootstrap_cz: Error: file \"%s\" does not correspond to alignment.\n",
                        bsp_file );
                exit( -1 );
            }
            positions[ i ] = p;
        }

        j++;

        if ( larger_blocks > 0 ) {
            if ( j >= block_size + 1 ) {
                flag = 1;
                j    = 0;
                larger_blocks--;
            }
        }
        else if ( j >= block_size ) {
            flag = 1;
            j    = 0;
        }

        if ( flag == 1 ) {
            if ( boot > 0 ) {
                fclose( outfile_fp );
            }
            sprintf( outfile_, "%s%d", outfile, z++ );
            if ( fileExists( outfile_ ) == 1 ) {
                printf( "\nbootstrap_cz: Error: outfile \"%s\" already exists.\n",
                        outfile_ );
                exit( -1 );
            }
            outfile_fp = fopen( outfile_, "w" );
            if ( outfile_fp == NULL ) {
                printf( "\nbootstrap_cz: Error: could not open file \"%s\" for writing.\n",
                        outfile_ );
                exit( -1 );
            }
            flag = 0;
        }

        fprintf( outfile_fp, " %d  %d\n", number_of_seqs, number_of_colm );
        for ( seq = 0; seq < number_of_seqs; ++seq ) {
            for ( i = 0; i < LENGTH_OF_NAME; ++i ) {
                fprintf( outfile_fp, "%c", names[ seq ][ i ] );   
            }
            for ( i = 0; i < number_of_colm; ++i ) {
                fprintf( outfile_fp, "%c", sequences[ seq ][ positions[ i ] ] );  
            }
            fprintf( outfile_fp, "\n" );
        }
    }

    /* Now, the bsp file must not contain any more numbers */
    if ( fscanf( bsp_file_fp, "%d", &p ) == 1 ) { 
        printf( "\nbootstrap_cz: Error: file \"%s\" does not correspond to alignment (too long).\n",
                bsp_file );
        printf( ">%d<\n", p );
        printf( "number of seqs=%d\n", number_of_seqs );
        exit( -1 );
    }

    fclose( bsp_file_fp ); 
    fclose( outfile_fp );

    free( positions );
    free( outfile_ );

    return;

} /* bootstrapAccordingToBSPfileP */




/* Exits if *p is NULL.        */
/* Last modified: 06/06/01     */
void checkForMemAllocFailure( void *p ) {
    if ( p == NULL ) {
        printf( "\nbootstrap_cz: Memory allocation failed.\n" );
        exit( -1 );
    }
    else {
        return;
    }
} /* checkForMemAllocFailure */



/* Returns 1 if filename can be opened. */
/* Returns 0 otherwise.                 */
/* Last modified: 06/07/01              */
int fileExists( const char *filename ) {
    FILE *fp = NULL;
    if ( ( fp = fopen( filename, "r" ) ) != NULL ) {
        fclose( fp );
        return 1;
    }
    else {
        return 0;
    }
} /* fileExists */



void errorInCommandLine() {
    printf( "\n" );
    printf( " bootstrap_cz  version 3.000\n" );
    printf( " ---------------------------\n\n" );
    printf( " Purpose:\n" );
    printf( " Bootstrap resamples an alignment in PHYLIP sequential format <bootstraps> times.\n" );
    printf( " Bootstrapping is not done randomly but according to a BSP (bootstrap\n" );
    printf( " positions) file.\n" );
    printf( " The BSP file can be created with the Perl program \"bootstrap_cz.pl\"\n" );
    printf( " in mode 0.\n" );
    printf( " This prgram has the same functionality as \"bootstrap_cz.pl\" in mode 1.\n" );
    printf( " Sequence names are normalized to LENGTH_OF_NAME characters.\n" );
    printf( " The output alignment is in PHYLIP's sequential or interleaved format.\n" );
    printf( " (These two are the same in this case, since all the seqs will be one\n" );
    printf( " line in length (no returns in seq).)\n\n" );
    printf( " Usage: bootstrap_cz <bootstraps> <alignment infile> <BSP file> <outfile>\n" );
    printf( "        [number of processors]\n\n" );
} /* errorInCommandLine */



int main( int argc, char *argv[] ) {
   
    char *inalign   = NULL,
         *bsp_file  = NULL,
         *outfile   = NULL;
    int  bootstraps = 0,
         processors = 0;
         

    if ( argc != 5 && argc != 6 ) {
        errorInCommandLine();
        exit( -1 );
    }
   
    bootstraps = atoi( argv[ 1 ] );
    inalign    = argv[ 2 ];
    bsp_file   = argv[ 3 ];
    outfile    = argv[ 4 ];

    if ( bootstraps < 1 ) {
        errorInCommandLine();
        exit( -1 );
    }

    if ( argc == 6 ) {
        processors = atoi( argv[ 5 ] );
        if ( processors < 1 ) {
            errorInCommandLine();
            exit( -1 );
        }
        if ( processors > bootstraps ) {
            processors = bootstraps;
        }
    }
    
    if ( argc == 5 && fileExists( outfile ) == 1 ) {
        printf( "\nbootstrap_cz: Error: outfile \"%s\" already exists.\n",
                outfile );
        exit( -1 );
    }
   
    readInAlignment( inalign );

    if ( argc == 5 ) { 
        bootstrapAccordingToBSPfile( bootstraps,
                                     bsp_file,
                                     outfile );
    }
    else {
        bootstrapAccordingToBSPfileP( bootstraps,
                                      processors,
                                      bsp_file,
                                      outfile );
    }
    
    return 0;

} /* main */
