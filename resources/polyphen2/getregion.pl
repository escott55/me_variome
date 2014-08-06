#!/usr/bin/perl -w

use strict;
use Bio::DB::BigWig 'binMean';
#use Mysql;
use List::MoreUtils qw/ uniq /;
use Getopt::Long;

my %ALLWIGS;

######################################################################
# openWig
######################################################################
sub openWig{
    my $chrom = $_[0];
    my $tree = $_[1];
    my $dir = "/home/escott/resources/phastcons";

    my $filename;
    if( $chrom =~ /chr/ ) {
        #print "Chrom:",$chrom,"\n";
        $filename = $dir."/".$tree."bigwig/".$chrom.".phastCons46way.bigwig";
    }elsif( $chrom ){
        #print "Chrom:",$chrom,"\n";
        $filename = $dir."/".$tree."bigwig/chr".$chrom.".phastCons46way.bigwig";
    } else {
        $filename = $dir."/".$tree."bigwig/chr12.phastCons46way.bigwig";
    }

    #my $wig  = Bio::DB::BigWig->new(-bigwig=>'ExampleData/dpy-27-variable.bw' )
    my $wig  = Bio::DB::BigWig->new(-bigwig=>$filename );
    return $wig
} # End openWig

######################################################################
# getRange
######################################################################
sub getRange{
    my $chrom = $_[0];
    my $start = $_[1];
    my $end = $_[2];
    my $tree = $_[3];

    my $wig = openWig($chrom, $tree);
    # Fetch individual intervals
    # fetch the individual data points in the wig file over a region of interest
    my @points;
    if( $chrom =~ /chr/ ){
        @points = $wig->features(-seq_id=>$chrom,-start=>$start,-end=>$end);
    } else {
        @points = $wig->features(-seq_id=>"chr".$chrom,-start=>$start,-end=>$end);
    }
    my @allvals;
    for my $p (@points) {
        my $start = $p->start;
        my $end   = $p->end;
        my $val   = $p->score;
        #print "$start..$end : $val\n";
        #push(@allvals,$start."\t".sprintf("%.3f",$val));
        push(@allvals,$start."\t".$val);
    }
    return @allvals;
} #END getRange

######################################################################
# getBase
######################################################################
sub getBase{
    my $chrom = $_[0];
    my $point = $_[1];
    my $tree = $_[2];

    my $wig;
    #print "Chrom ",$chrom," End: ",$point,"\n";
    if( exists $ALLWIGS{$chrom} ){
        $wig = $ALLWIGS{$chrom};
    } else {
        $wig = openWig($chrom,$tree);
        $ALLWIGS{$chrom} = $wig;
    }
   
    # Get specific point
    my @points;
    if( $chrom =~ /chr/ ){
        @points = $wig->features(-seq_id=>$chrom,-start=>$point,-end=>$point);
    }else{
        @points = $wig->features(-seq_id=>"chr".$chrom,-start=>$point,-end=>$point);
    }

    my $val;
    my $start;
    my $end;
    for my $p (@points) {
        $start = $p->start;
        $end   = $p->end;
        $val   = $p->score;
        #print "$start..$end : $val\n";
        #print $val."\n"
    }
    return $start."\t".$val;
    #return $start.":".sprintf("%.3f",$val)
} #END getBase

######################################################################
# runFile
######################################################################
sub runFile{
    my $chrom;
    my $filename = "overlap.txt";

    open FILE, $filename or die $1;
    while( <FILE> ){
        chomp;
        my @atts = split("\t");
        if( scalar(@atts) == 0 ){
            next;
        }
        if ($atts[2] =~ /\D/){
            next;
        }elsif ( $atts[2] == 23 ) {
            $chrom = "X";
        }elsif ( $atts[2] == 24 ) {
            $chrom = "Y";
        } else {
            $chrom = $atts[2];
        }
        my $val = getBase($chrom, $atts[3]);
        #print $chrom,":",$atts[3],":",$atts[5]," - ",$val,"\n";
        print join("\t",@atts),"\t",$val,"\n";
    } 
} # End runFile

#######################################################################
## updateDatabase
#######################################################################
#sub updateDatabase {
    #my $host = "localhost";
    #my $database = "gleeson2";
    #my $user = "root";
    #my $pw = "(umulus88";

    #$connect = Mysql->connect( $host, $database, $user, $pw);
    #$connect->selectdb($database);
    #my $variantquery = "Select id, chrom, pos, ref, mut FROM Variants WHERE vertPhastCons IS NULL";
    ##limit 0,30";
    #my $execute = $connect->query($variantquery);
    #while( my @results = $execute->fetchrow() ){
        #my $val = getBase($results[1], $results[2]);
        #if( defined($val) ){
            #print join(", ", @results),"\t",$val,"\n";
            #my $updatequery = "Update Variants SET vertPhastCons = ".$val." WHERE id= ".$results[0]." "; 
            #my $updateex = $connect->query($updatequery);
            ##print $updateex,"\n";
        #}
    #}
#}
## END updateDatabase

######################################################################
# printHelp
######################################################################
sub printHelp {
    my $programname=$_[0];
    print "Program:".$programname,"\n";
    exit 1;
}

######################################################################
# Main
######################################################################
{

    my $verbose = '';
    #my $chrom = "chr1";
    #my $start = 45567859;
    #my $end = 45567889;
    my $chrom;
    my $start;
    my $end;
    my $tree = "vert"; # vert, primates, placentalMammals
    my $reverse = '';
    my $debug = '';
    GetOptions( "verbose!"=>\$verbose,
            "chrom=s"=>\$chrom,
            "tree:s",\$tree,
            "start=i"=> \$start,
            "end:i"=> \$end,
            "debug!"=> \$debug, 
            "reverse!"=> \$reverse ); 

    if( not $chrom or not $start ){
        print "Error: parameters have not been set correctly";
        printHelp($0);
    }

    # Fix Chromosome
    if( $chrom !~ /chr/ ){
        $chrom = "chr".$chrom;
    }

    if( $debug ){
        print "Verbose: ",$verbose,"\n";
        print "Chrom: ",$chrom,"\n";
        print "tree: ",$tree,"\n";
        print "start: ",$start,"\n";
        if( $end ){
            print "end: ",$end,"\n";
        }
        print "reverse: ",$reverse,"\n";
    }

    if( not $end ){
        print "Get Base\n" if($debug);
        my $val = getBase( $chrom, $start, $tree );
        #print sprintf("%.3f",$val),"\n";
        print $val,"\n";
    } else {
        print "Get Range\n" if($debug);
        my @range = getRange($chrom, $start, $end, $tree);
        if( $reverse ){
            print "Rev:",$reverse,"\n" if($debug);
            @range = reverse(@range)
        }
        # Print Values
        for my $val (@range){
            print $val,"\n";
        }
    }

    #runFile( );
    #updateDatabase();
} # END MAIN

######################################################################
######################################################################
    # Get User Input
    #my $chrom = shift;
    #my $position = shift;

    #if ( not $position ){
        #$position = 45567859;
        #print "Test position: ".$position."\n"; 
    #}
    #if( not $chrom ){
        #$chrom = 12;
    #}
 
