#!/usr/bin/perl -w

my $in  = $ARGV[ 0 ];
my $out = $ARGV[ 1 ];

if ( -e $out ) {
    print "File $out already exists.\n";
    exit( -1 );
} 

if ( !-e $in ) {
    print "File $in does not exist.\n";
    exit( -1 );
} 

open( IN, $in ) ;
open ( OUT, ">>$out"  ) || die ( "Could not open file $out for writing!\n" );

while ( my $line = <IN> ) {
    my $newline = &proc_line( $line );
    if ( length( $newline ) > 0 ) {
        print OUT $newline;
    }
}


close( OUT ) or die( "can't close $out: $!" );
close( IN ) or die( "can't close $in: $!" );

sub proc_line {
   my $line = shift;
   
   
   if ( $line =~ /^#/ ) {
       return "";
   }
   if ( $line =~ /^Predicted coding sequence\(s\):/ ) {
       return "";
   }
   elsif ( $line =~ /^>.*_aa\s*$/ ) {
       return "";
   }
   elsif ( $line =~ /^>/ ) {
       return $line;
   }
   elsif ( $line !~ /[a-z]/ ) {
       return "";
   }
   else {
       return $line;;
   } 
}
