#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#v1
#Take a btab (blast, tab-delimited) file and parse the results to a hit/nohit table display.
#A note on assumptions about the query and blastdbs:
#   1) The query can be anything, but it would behoove you to have the name up till the first space be unique and informative.
#   2) For the blastdb, the assumption for the header, is PC_#######_#####_genomename. The first set of numbers is the anvio protein cluster id, the second is the unique protein id.
#v2 features
#   1) Files for genome order and query order.
#   2) Tell us what PC-bin the query is in.
#v3 features
#   1) Make the grepping for PC-bin optional.


# - - - - - O P T I O N S  - - - - - -
my %options=();
getopts("i:ae:p:hg:q:f:b", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Tab-delimited (outformat 6) BLAST output (btab).\n";
        print "-b = Add column, connecting queries with PC-bins (slows down program significantly).\n";
        print "-f = THE file exported from Anvio. (required only if -b is called)\n";
        print "-g = Genome order file (one per line) (optional)\n";
        print "-q = Query order file (one per line) (optional)\n";
        print "-a = Print all hits from the btab file.\n";
        print "-e = E-value cutoff for reporting a hit.\n";
        print "-p = Percent identity (PID) cutoff for reporting a hit.\n";
	print "-h = This help message\n\n";
	die;
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my %BlastHits;
#my %Sequences;

open(IN, "<$options{i}") or die "\n\nFile $options{i} does not exist or was not given. Try -h for the help file.\n\n";
my @btabdata = <IN>; close(IN);
my $unid = 10000000000001;
my @query;
my @genomes;
foreach my $line (@btabdata)
    {       #my $counting = $unid - 10000000000000;
            #print STDERR "On number $counting of $#btabdata\n";
            chomp($line);
            my @data = split('\t', $line);
            $BlastHits{$unid}{'query'} = $data[0];
            $BlastHits{$unid}{'db_entry'} = $data[1];
            $BlastHits{$unid}{'pident'} = $data[2];
            $BlastHits{$unid}{'length'} = $data[3];
            $BlastHits{$unid}{'mismatch'} = $data[4];
            $BlastHits{$unid}{'gapopen'} = $data[5];
            $BlastHits{$unid}{'qstart'} = $data[6];
            $BlastHits{$unid}{'qend'} = $data[7];
            $BlastHits{$unid}{'sstart'} = $data[8];
            $BlastHits{$unid}{'send'} = $data[9];
            $BlastHits{$unid}{'evalue'} = $data[10];
            $BlastHits{$unid}{'bitscore'} = $data[11];
            my $parse_me = $data[1];
            $parse_me =~ m/(PC_[0-9]+)_([0-9]+)_(.+)/;
            $BlastHits{$unid}{'protein_cluster'} = $1;
            $BlastHits{$unid}{'uniq_protein_id'} = $2;
            $BlastHits{$unid}{'genome_name'} = $3;
            if ($options{b})
                {   my $seduse = $BlastHits{$unid}{'uniq_protein_id'} + 1; 
                    my $command = "cat $options{f} | sed '".$seduse."q;d'";
                    my $testing_sed = `$command`;
                    my @parse_results = split('\t', $testing_sed);
                    if ($parse_results[0] == $BlastHits{$unid}{'uniq_protein_id'} && $parse_results[1] eq $BlastHits{$unid}{'protein_cluster'})
                        {   $BlastHits{$unid}{'PC-bin'} = $parse_results[2];
                        }
                    else {print STDERR "Possible problems: No PC-bin match for unique_protein_id $BlastHits{$unid}{'uniq_protein_id'}\n";}
                    #print STDERR "$BlastHits{$unid}{'PC-bin'}\n";
                }
            unshift(@genomes, $3);
            push(@query, $data[0]);
            $unid += 1;
    }
my @uniq_query = uniq @query;
my @test_query = @uniq_query;
my @uniq_genomes = uniq @genomes;

if ($options{g})
    {   open(GIN, "<$options{g}") or die "\n\nFile $options{g} does not exist or was not given. Try -h for the help file.\n\n";
        my @genome_order = <GIN>; close(GIN);
        my @tt;
        foreach my $p (@genome_order)
            {   chomp($p);
                unless ($p eq "")
                    {   push(@tt, $p);
                    }
            }
        @uniq_genomes = @tt;
    }

if ($options{q})
    {   open(QIN, "<$options{q}") or die "\n\nFile $options{q} does not exist or was not given. Try -h for the help file.\n\n";
        my @quer_order = <QIN>; close(QIN);
        my @ff;
        foreach my $p (@quer_order)
            {   chomp($p);
                unless ($p eq "")
                    {   push(@ff, $p);
                    }
            }
        @uniq_query = @ff;
        foreach my $f (@uniq_query)
                {       my $exists = "FALSE";
                        foreach my $g (sort keys %BlastHits)
                                {       if ($BlastHits{$g}{'query'} eq $f)
                                                {       $exists = "TRUE";     
                                                }
                                }
                        if ($exists eq "FALSE")
                                {       $unid += 1;
                                        $BlastHits{$unid}{'query'} = $f;
                                        $BlastHits{$unid}{'genome_name'} = "";
                                        $BlastHits{$unid}{'uniq_protein_id'} = "";
                                        $BlastHits{$unid}{'protein_cluster'} = "";
                                        $BlastHits{$unid}{'db_entry'} = "";
                                        $BlastHits{$unid}{'pident'} = "";
                                        $BlastHits{$unid}{'evalue'} = "";
                                        if ($options{b}) {      $BlastHits{$unid}{'PC-bin'} = "";}
                                }
                        
                        
                }
        foreach my $lets (@test_query)
                {       my $oops = "FALSE";
                        foreach my $sort (@uniq_query)
                                {       if ($lets eq $sort)
                                                {       $oops = "TRUE";       
                                                }
                                }
                        if ($oops eq "FALSE")
                                {       print STDERR "This query from the blast results is not in your query mapping file: $lets\n";
                                }
                }
    }

print "Query\t";
foreach my $i (0..$#uniq_genomes)
    {   print "$uniq_genomes[$i]\t";   
    }
if ($options{b})
    {   print "PC-bin\n";
    }
else {  print "\n";
     }

my %Query;  #Dictionary to store the results from a query
foreach my $i (sort keys %BlastHits)
    {   my $insert = $BlastHits{$i}{'query'};
        $Query{$insert}{'genome_name'} .= $BlastHits{$i}{'genome_name'}."~";
        $Query{$insert}{'uniq_protein_id'} .= $BlastHits{$i}{'uniq_protein_id'}."~";
        $Query{$insert}{'protein_cluster'} .= $BlastHits{$i}{'protein_cluster'}."~";
        $Query{$insert}{'db_entry'} .= $BlastHits{$i}{'db_entry'}."~";
        $Query{$insert}{'pident'} .= $BlastHits{$i}{'pident'}."~";
        $Query{$insert}{'evalue'} .= $BlastHits{$i}{'evalue'}."~";
        if ($options{b}) {  $Query{$insert}{'PC-bin'} .= $BlastHits{$i}{'PC-bin'}."~";}
        $Query{$insert}{'count'} += 1;
    }

if ($options{e})
    {   foreach my $i (0..$#uniq_query)
            {   print "$uniq_query[$i]\t";
                my $querygenomes = $Query{$uniq_query[$i]}{'genome_name'};
                my @split_genomes = split('~', $querygenomes);
                my $queryevalues = $Query{$uniq_query[$i]}{'evalue'};
                my @split_evalues = split('~', $queryevalues);
                my $queryuniqpid = $Query{$uniq_query[$i]}{'uniq_protein_id'};
                my @split_uniqpid = split('~', $queryuniqpid);
                my $queryPC = $Query{$uniq_query[$i]}{'protein_cluster'};
                my @split_PC = split('~', $queryPC);
                my $querypercident = $Query{$uniq_query[$i]}{'pident'};
                my @split_percident = split('~', $querypercident);
                #my $querypcbin = $Query{$i}{'pc_bin'};
                #my @split_pcbin = split('~', $querypcbin);
                foreach my $j (0..$#uniq_genomes)
                    {   foreach my $k (0..$#split_genomes)
                            {   if ($uniq_genomes[$j] eq $split_genomes[$k])
                                    {   if ($split_evalues[$k] <= $options{e})
                                            {   print "$split_PC[$k]_$split_uniqpid[$k] ($split_evalues[$k], $split_percident[$k]);";
                                            }
                                    }
                            }
                        print " \t";
                    }
                if ($options{b})
                    {   my $pcbin = $Query{$uniq_query[$i]}{'PC-bin'};
                        my @split_pcbin = split('~', $pcbin);
                        my @unique_pcbin = uniq @split_pcbin;
                        my $unieuqpcbin = join(';', @unique_pcbin);
                        print "$unieuqpcbin\n";
                    }
                else {print "\n";}
            }
    }


if ($options{p})
    {   foreach my $i (0..$#uniq_query)
            {   print "$uniq_query[$i]\t";
                my $querygenomes = $Query{$uniq_query[$i]}{'genome_name'};
                my @split_genomes = split('~', $querygenomes);
                my $queryevalues = $Query{$uniq_query[$i]}{'evalue'};
                my @split_evalues = split('~', $queryevalues);
                my $queryuniqpid = $Query{$uniq_query[$i]}{'uniq_protein_id'};
                my @split_uniqpid = split('~', $queryuniqpid);
                my $queryPC = $Query{$uniq_query[$i]}{'protein_cluster'};
                my @split_PC = split('~', $queryPC);
                my $querypercident = $Query{$uniq_query[$i]}{'pident'};
                my @split_percident = split('~', $querypercident);
                #my $querypcbin = $Query{$i}{'pc_bin'};
                #my @split_pcbin = split('~', $querypcbin);
                foreach my $j (0..$#uniq_genomes)
                    {   foreach my $k (0..$#split_genomes)
                            {   if ($uniq_genomes[$j] eq $split_genomes[$k])
                                    {   if ($split_percident[$k] >= $options{p})
                                            {   print "$split_PC[$k]_$split_uniqpid[$k] ($split_evalues[$k], $split_percident[$k]);";
                                            }
                                    }
                            }
                        print " \t";
                    }
                if ($options{b})
                    {   my $pcbin = $Query{$uniq_query[$i]}{'PC-bin'};
                        my @split_pcbin = split('~', $pcbin);
                        my @unique_pcbin = uniq @split_pcbin;
                        my $unieuqpcbin = join(';', @unique_pcbin);
                        print "$unieuqpcbin\n";
                    }
                else {print "\n";}
            }
    }

if ($options{a})
    {   foreach my $i (0..$#uniq_query)
            {   print "$uniq_query[$i]\t";
                my $querygenomes = $Query{$uniq_query[$i]}{'genome_name'};
                my @split_genomes = split('~', $querygenomes);
                my $queryevalues = $Query{$uniq_query[$i]}{'evalue'};
                my @split_evalues = split('~', $queryevalues);
                my $queryuniqpid = $Query{$uniq_query[$i]}{'uniq_protein_id'};
                my @split_uniqpid = split('~', $queryuniqpid);
                my $queryPC = $Query{$uniq_query[$i]}{'protein_cluster'};
                my @split_PC = split('~', $queryPC);
                my $querypercident = $Query{$uniq_query[$i]}{'pident'};
                my @split_percident = split('~', $querypercident);
                #my $querypcbin = $Query{$i}{'pc_bin'};
                #my @split_pcbin = split('~', $querypcbin);
                foreach my $j (0..$#uniq_genomes)
                    {   foreach my $k (0..$#split_genomes)
                            {   if ($uniq_genomes[$j] eq $split_genomes[$k])
                                    {   print "$split_PC[$k]_$split_uniqpid[$k] ($split_evalues[$k], $split_percident[$k]);";
                                    }
                            }
                        print " \t";
                    }
                if ($options{b})
                    {   my $pcbin = $Query{$uniq_query[$i]}{'PC-bin'};
                        my @split_pcbin = split('~', $pcbin);
                        my @unique_pcbin = uniq @split_pcbin;
                        my $unieuqpcbin = join(';', @unique_pcbin);
                        print "$unieuqpcbin\n";
                    }
                else {print "\n";}
            }
    }
    
if ($options{e} && $options{p})
    {   foreach my $i (0..$#uniq_query)
            {   print "$uniq_query[$i]\t";
                my $querygenomes = $Query{$uniq_query[$i]}{'genome_name'};
                my @split_genomes = split('~', $querygenomes);
                my $queryevalues = $Query{$uniq_query[$i]}{'evalue'};
                my @split_evalues = split('~', $queryevalues);
                my $queryuniqpid = $Query{$uniq_query[$i]}{'uniq_protein_id'};
                my @split_uniqpid = split('~', $queryuniqpid);
                my $queryPC = $Query{$uniq_query[$i]}{'protein_cluster'};
                my @split_PC = split('~', $queryPC);
                my $querypercident = $Query{$uniq_query[$i]}{'pident'};
                my @split_percident = split('~', $querypercident);
                #my $querypcbin = $Query{$i}{'pc_bin'};
                #my @split_pcbin = split('~', $querypcbin);
                foreach my $j (0..$#uniq_genomes)
                    {   foreach my $k (0..$#split_genomes)
                            {   if ($uniq_genomes[$j] eq $split_genomes[$k])
                                    {   if ($split_percident[$k] >= $options{p} && $split_evalues[$k] <= $options{e})
                                            {   print "$split_PC[$k]_$split_uniqpid[$k] ($split_evalues[$k], $split_percident[$k]);";
                                            }
                                    }
                            }
                        print " \t";
                    }
                if ($options{b})
                    {   my $pcbin = $Query{$uniq_query[$i]}{'PC-bin'};
                        my @split_pcbin = split('~', $pcbin);
                        my @unique_pcbin = uniq @split_pcbin;
                        my $unieuqpcbin = join(';', @unique_pcbin);
                        print "$unieuqpcbin\n";
                    }
                else {print "\n";}
            }
    }


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -