#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Btab format from blast run is the input.
#v1
#   1) Choose the best hit for a gene from a single genome that hits multiple queries and filter other hits.

# - - - - - O P T I O N S  - - - - - -
my %options=();
getopts("i:pmhe:d:", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Output file from blast (btab format)\n";
        print "-p = Base filtering decisions on percent identity. (Default = e-value).\n";
        print "-m = Base filtering decisions on e-value magnitude only (i.e. 10^-30). (Default 1E-30 better hit than 2E-30)\n";
        print "-e = Evalue cut off for filtering\n";
        print "-d = Percent identity cut off for filtering\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my %BtabHits;

open(IN, "<$options{i}") or die "\n\nFile $options{i} does not exist or was not given. Try -h for the help file.\n\n";
my @btabdata = <IN>; close(IN);
my $unid = 100001;
foreach my $line (@btabdata)
    {   chomp($line);
        my @split_line = split('\t', $line);
        $BtabHits{$unid}{'query'} = $split_line[0];
        $BtabHits{$unid}{'genome_hit'} = $split_line[1];
        $BtabHits{$unid}{'pident'} = $split_line[2];
        $BtabHits{$unid}{'length'} = $split_line[3];
        $BtabHits{$unid}{'mismatch'} = $split_line[4];
        $BtabHits{$unid}{'gapopen'} = $split_line[5];
        $BtabHits{$unid}{'qstart'} = $split_line[6];
        $BtabHits{$unid}{'qend'} = $split_line[7];
        $BtabHits{$unid}{'sstart'} = $split_line[8];
        $BtabHits{$unid}{'send'} = $split_line[9];
        $BtabHits{$unid}{'evalue'} = $split_line[10];
        $BtabHits{$unid}{'bitscore'} = $split_line[11];
        $unid += 1; 
    }

foreach my $i (sort keys %BtabHits)
    {   my $evalue = $BtabHits{$i}{'evalue'};
        my $pids = $BtabHits{$i}{'pident'};
        if ($options{e})
            {   if ($evalue >= $options{e})
                    {   delete $BtabHits{$i};
                    }
            }
        if ($options{d})
            {   if ($pids <=  $options{d})
                    {   delete $BtabHits{$i};   
                    }
            }
    }

my %Multiquery;
foreach my $i (sort keys %BtabHits)
    {   foreach my $j (sort keys %BtabHits)
            {   unless ($i == $j)
                    {if ($BtabHits{$i}{'genome_hit'} eq $BtabHits{$j}{'genome_hit'})
                        {   $Multiquery{$BtabHits{$i}{'genome_hit'}}{'btab_key'} .= $i."~";
                            $Multiquery{$BtabHits{$i}{'genome_hit'}}{'btab_key'} .= $j."~"; 
                        }
                    }
            }
    }

foreach my $i (sort keys %Multiquery)
    {   my $keysOI = $Multiquery{$i}{'btab_key'};
        my @keysOI_split = split('~', $keysOI);
        my @uniq_keys = uniq @keysOI_split;
        my $new_keys = join('~', @uniq_keys);
        $Multiquery{$i}{'btab_key'} = $new_keys;
    }

foreach my $i (sort keys %Multiquery)
    {   my $btab_keys_of_interest = $Multiquery{$i}{'btab_key'};
        my @split_btab = split('~', $btab_keys_of_interest);
        my $choice_number = 0;
        my $choice_list = "0;";
        foreach my $j (1..$#split_btab)
            {   if ($options{p})
                    {   if ($BtabHits{$split_btab[$j]}{'pident'} > $BtabHits{$split_btab[$choice_number]}{'pident'})
                            {   $choice_number = $j;
                                $choice_list = "$j;";
                            }
                        elsif ($BtabHits{$split_btab[$j]}{'pident'} == $BtabHits{$split_btab[$choice_number]}{'pident'})
                            {   $choice_list .= "$j;";
                            }
                    }
                else
                    {   if ($options{m})
                            {   my $choice_mag = &ROUNDEXP($BtabHits{$split_btab[$choice_number]}{'evalue'});
                                my $test_mag = &ROUNDEXP($BtabHits{$split_btab[$j]}{'evalue'});
                                if ($test_mag < $choice_mag)
                                    {   $choice_number = $j;
                                        $choice_list = "$j;";
                                    }
                                elsif ($test_mag == $choice_mag)
                                    {   $choice_list .= "$j;";
                                    }
                            }
                        else
                            {   if ($BtabHits{$split_btab[$j]}{'evalue'} < $BtabHits{$split_btab[$choice_number]}{'evalue'})
                                    {   $choice_number = $j;
                                        $choice_list = "$j;";
                                    }
                                elsif ($BtabHits{$split_btab[$j]}{'evalue'} == $BtabHits{$split_btab[$choice_number]}{'evalue'})
                                    {   $choice_list .= "$j;";
                                    }
                            }
                    }
            }
        my @choice_split = split(';', $choice_list);
        my @yes_please;
        foreach my $k (@choice_split)
            { push(@yes_please, $split_btab[$k]);
            }
        my @no_thank_you;
        foreach my $b (@split_btab)
            {   my $keepq = "FALSE";
                foreach my $f (@yes_please)
                    {   if ($b eq $f)
                            {   $keepq = "TRUE";
                            }
                    }
                if ($keepq eq "FALSE")
                    {   push(@no_thank_you, $b);
                    }
            }
        foreach my $deleteme (@no_thank_you)
            {   delete $BtabHits{$deleteme};
            }
    }
    
foreach my $i (sort keys %BtabHits)
    {   print "$BtabHits{$i}{'query'}\t";
        print "$BtabHits{$i}{'genome_hit'}\t";
        print "$BtabHits{$i}{'pident'}\t";
        print "$BtabHits{$i}{'length'}\t";
        print "$BtabHits{$i}{'mismatch'}\t";
        print "$BtabHits{$i}{'gapopen'}\t";
        print "$BtabHits{$i}{'qstart'}\t";
        print "$BtabHits{$i}{'qend'}\t";
        print "$BtabHits{$i}{'sstart'}\t";
        print "$BtabHits{$i}{'send'}\t";
        print "$BtabHits{$i}{'evalue'}\t";
        print "$BtabHits{$i}{'bitscore'}\n";
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub ROUNDEXP
{	my $num = $_[0];
        if ($num == 0)
            {   my $good_num = $num - 1000000000000;
                return $good_num;
            }
        else
            {   my $lognum = &LOG10($num);
                my $roundnum;
                if ($lognum >= 0)
                    {   $roundnum = int($lognum + 0.5);}
                if ($lognum < 0)
                    {   $roundnum = int($lognum - 0.5);}
                return $roundnum;
            }
}

sub LOG10
{   my $n = $_[0];
    return log($n)/log(10);   
}


    
# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -