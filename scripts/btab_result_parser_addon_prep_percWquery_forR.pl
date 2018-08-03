#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Assumption is that you have run the btab_result_parser addon (queryclustering_v2).
#One of the outputs of that script is a file that has printed the percentage of each genome cluster per query.
#File should have "_clusteredgenome_percentWquery.txt" on the end.
#v1
#   1) Convert that matrix to something that R will understand.


# - - - - - O P T I O N S  - - - - - -
my %options=();
getopts("i:m:n:g:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = x_clusteredgenome_percentWquery.txt output file from btab_result_parser_addon_queryclustering_v2.pl script.\n";
	print "-m = Metadata information for colouring queries by pathway or other common category. (Format: query".'\t'."name) (optional)\n";
        print "-n = Query name mapping file for new query names (Format: old".'\t'."new) (optional)\n";
        print "-g = Genome cluster name mapping file for new genome cluster names (Format: old".'\t'."new) (optional)\n";
        print "-h = This help message\n\n";
        print 'Note: \t'." = tab character\n\n";
	die;
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my %PercentDict;
my %QueryName;
my %GenomeName;
my %Colors;

open(IN, "<$options{i}") or die "\n\nFile $options{i} does not exist or was not given. Try -h for the help file.\n\n";
my @percdata = <IN>; close(IN);
my $firstline = shift(@percdata);
chomp($firstline);
my @split_header = split('\t', $firstline);

if ($options{n})
    {   open(IN2, "<$options{n}") or die "\n\nFile $options{n} does not exist or was not given. Try -h for the help file.\n\n";
        my @queryname = <IN2>; close(IN2);
        foreach my $q (@queryname)
            {   chomp($q);
                my @splitq = split('\t', $q);
                $QueryName{$splitq[0]}{'newname'} = $splitq[1];
            }
    }

if ($options{g})
    {   open(IN3, "<$options{g}") or die "\n\nFile $options{g} does not exist or was not given. Try -h for the help file.\n\n";
        my @genomename = <IN3>; close(IN3);
        foreach my $g (@genomename)
            {   chomp($g);
                my @splitg = split('\t', $g);
                $GenomeName{$splitg[0]}{'newname'} = $splitg[1];
            }
    }
    
if ($options{m})
    {  open(IN4, "<$options{m}") or die "\n\nFile $options{m} does not exist or was not given. Try -h for the help file.\n\n";
        my @colorcat = <IN4>; close(IN4);
        foreach my $c (@colorcat)
            {   chomp($c);
                my @splitc = split('\t', $c);
                $Colors{$splitc[0]}{'newcat'} = $splitc[1];
            }
    }

my $unid = 1;
foreach my $i (@percdata)
    {   chomp($i);
        my @split_data = split('\t', $i);
        foreach my $j (1..$#split_data)
            {   $PercentDict{$unid}{'query'} = $split_data[0];
                $PercentDict{$unid}{'genome'} = $split_header[$j];
                $PercentDict{$unid}{'value'} = $split_data[$j];
                if ($options{m})
                    {   $PercentDict{$unid}{'category'} = $Colors{$split_data[0]}{'newcat'};}
                if ($options{n})
                    {   $PercentDict{$unid}{'query'} = $QueryName{$split_data[0]}{'newname'};}
                if ($options{g})
                    {   $PercentDict{$unid}{'genome'} = $GenomeName{$split_header[$j]}{'newname'};}
                $unid += 1;
            }
    }
my $final_key = $unid - 1;


print ",Query_Gene,";
if ($options{m})
    {   print "Category,";}
print "Genome_Cluster,Value\n";

foreach my $i (1..$final_key)
    {   print "$i,$PercentDict{$i}{'query'},";
        if ($options{m})
            {   print "$PercentDict{$i}{'category'},";}
        print "$PercentDict{$i}{'genome'},$PercentDict{$i}{'value'}\n";
    }
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -