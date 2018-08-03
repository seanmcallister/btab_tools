#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Assumption is that you have run the btab_result_parser program. The input to this is the output from that.
#v1
#   1) Hierarchical clustering of queries
#v2
#   2) Output an additional file showing the percentage of genomes in a cluster that have the query of interest.

# - - - - - O P T I O N S  - - - - - -
my %options=();
getopts("i:g:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Output file from btab_result_parser_v3.pl script. Doesn't support PC-Bin column.\n";
	print "-g = Genome clustering information (i.e. ZOTU cluster groups). Two rows: genome_name".'\t'."cluster_name\n";
        print "-h = This help message\n\n";
	die;
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my %Genomes;
my %Clusters;

open(IN3, "<$options{g}") or die "\n\nNADA $options{g} you FOOL!!!\n\n";
my @sampleclusters = <IN3>; close(IN3);
foreach my $j (@sampleclusters)
    {   chomp($j);
        my @clustersplit = split('\t', $j);
        $Genomes{$clustersplit[0]}{'cluster'} = $clustersplit[1];
        $Clusters{$clustersplit[1]}{'genomes'} .= $clustersplit[0].";";
    }

open(IN, "<$options{i}") or die "\n\nFile $options{i} does not exist or was not given. Try -h for the help file.\n\n";
my @btabdata = <IN>; close(IN);
my $firstline = shift(@btabdata);
chomp($firstline);

my @split_header = split('\t', $firstline);
my @genome_order = @split_header;
#my $pcbin_test = "FALSE";
#if ($split_header[$#split_header] eq "PC-bin")
#    {   pop(@genome_order);
#        pop(@split_header);
#        $firstline = join("\t", @split_header);
#        $pcbin_test = "TRUE";
#    }
shift(@genome_order);
foreach my $i (sort keys %Clusters)
    {   my @coord;
        my $gens = $Clusters{$i}{'genomes'};
        my @Gens = split(';', $gens);
        #print "$i\n";
        #print "<<@Gens>>\n";
        foreach my $j (@Gens)
            {   foreach my $k (1..$#split_header)
                    {   #print "$split_header[$k]\n";
                    #    print "$j\n";
                        if ($split_header[$k] eq $j)
                            {   push(@coord, $k);
                            }
                    }
            }
        my $coordinates = join(';', @coord);
        $Clusters{$i}{'coordinates'} = $coordinates;
    }

open(OUT2, ">".$options{i}."_clusteredgenome_countdata.txt");
print OUT2 "Query\t";
foreach my $i (sort keys %Clusters)
    {   print OUT2 "$i\t";
    }
print OUT2 "\n";


open(OUT3, ">".$options{i}."_clusteredgenome_percentWquery.txt");
print OUT3 "Query\t";
foreach my $i (sort keys %Clusters)
    {   print OUT3 "$i\t";
    }
print OUT3 "\n";

open(OUT, ">".$options{i}."_countdata.txt");
print OUT "$firstline\n";



foreach my $line (@btabdata)
    {       chomp($line);
            my @splitline = split('\t', $line);
            my $query = shift(@splitline);
            #if ($pcbin_test eq "TRUE")
            #    {   pop(@splitline);
            #    }
            print OUT "$query\t";
            foreach my $i (@splitline)
                {   my $count = $i =~ tr/;/;/;
                    print OUT "$count\t";
                }
            print OUT "\n";
    }

close(OUT);

foreach my $line (@btabdata)
    {   chomp($line);
        my @splitline = split('\t', $line);
        my $query = $splitline[0];
        print OUT2 "$query\t";
        foreach my $i (sort keys %Clusters)
            {   my $coords = $Clusters{$i}{'coordinates'};
                my @Coords = split(';', $coords);
                my $count = 0;
                foreach my $j (@Coords)
                    {   my $value = $splitline[$j];
                        $count += $value =~ tr/;/;/;
                    }
                print OUT2 "$count\t";
            }
        print OUT2 "\n";
    }
close(OUT2);


foreach my $line (@btabdata)
    {   chomp($line);
        my @splitline = split('\t', $line);
        my $query = $splitline[0];
        print OUT3 "$query\t";
        foreach my $i (sort keys %Clusters)
            {   my $coords = $Clusters{$i}{'coordinates'};
                my @Coords = split(';', $coords);
                my $totalgenomesincluster = $#Coords + 1;
                #print "$i\n";
                #print "$totalgenomesincluster\n";
                my $genonmesWhitcount = 0;
                foreach my $j (@Coords)
                    {   my $value = $splitline[$j];
                        if ($value =~ m/;/)
                            {   $genonmesWhitcount += 1;
                            }
                    }
                my $perc1 = $genonmesWhitcount / $totalgenomesincluster;
                my $percentgenomesincluster = $perc1 * 100;
                my $roundedperc = &ROUND($percentgenomesincluster,2);
                print OUT3 "$roundedperc\t";
            }
        print OUT3 "\n";
    }
close(OUT3);


my $wd = `pwd`;
chomp ($wd);
open(RSCRIPT, ">rcode_count_hclust.R");
print RSCRIPT "#! /usr/bin/env Rscript\n";
print RSCRIPT "setwd(\"$wd\")\n\n";
print RSCRIPT "count <- read.table(\"$options{i}_countdata.txt\", header = TRUE, sep = \"\\t\", row.names=1)\n";
print RSCRIPT "count_dist <- dist(count, method = \"euclidean\", diag = FALSE, upper = FALSE)\n";
print RSCRIPT "complete <- hclust(count_dist, method = \"complete\")\n\n";
print RSCRIPT "complete_dendro <- as.dendrogram(complete)\n";
print RSCRIPT "fileConn<-file(\"leaforder_allgenomes_completelinkage.txt\")\n";
print RSCRIPT "writeLines(c(labels(complete_dendro)), fileConn)\n";
print RSCRIPT "close(fileConn)\n\n";
print RSCRIPT "averagelinkage <- hclust(count_dist, method = \"average\")\n\n";
print RSCRIPT "averagelinkage_dendro <- as.dendrogram(averagelinkage)\n";
print RSCRIPT "fileConn<-file(\"leaforder_allgenomes_averagelinkage.txt\")\n";
print RSCRIPT "writeLines(c(labels(averagelinkage_dendro)), fileConn)\n";
print RSCRIPT "close(fileConn)\n\n";
print RSCRIPT "pdf(\"hierarchical_clustering_allgenomes_completelinkage.pdf\", width = 9, height = 6)\n";
print RSCRIPT "plot_complete <- plot(complete, check = TRUE, cex = 0.25, main = \"Blast Query Hierarchical Clustering (All genomes, Complete)\")\n";
print RSCRIPT "dev.off()\n\n";
print RSCRIPT "pdf(\"hierarchical_clustering_allgenomes_averagelinkage.pdf\", width = 9, height = 6)\n";
print RSCRIPT "plot_average <- plot(averagelinkage, check = TRUE, cex = 0.25, main = \"Blast Query Hierarchical Clustering (All genomes, Average)\")\n";
print RSCRIPT "dev.off()\n\n";
#Set up counts by genome cluster (AKA ZOTU)
print RSCRIPT "clustered_genomes <- read.table(\"$options{i}_clusteredgenome_countdata.txt\", header = TRUE, sep = \"\\t\", row.names=1)\n";
print RSCRIPT "count_dist_clu <- dist(clustered_genomes, method = \"euclidean\", diag = FALSE, upper = FALSE)\n";
print RSCRIPT "complete_clu <- hclust(count_dist_clu, method = \"complete\")\n\n";
print RSCRIPT "complete_clu_dendro <- as.dendrogram(complete_clu)\n";
print RSCRIPT "fileConn<-file(\"leaforder_clusteredgenomes_completelinkage.txt\")\n";
print RSCRIPT "writeLines(c(labels(complete_clu_dendro)), fileConn)\n";
print RSCRIPT "close(fileConn)\n\n";
print RSCRIPT "averagelinkage_clu <- hclust(count_dist_clu, method = \"average\")\n\n";
print RSCRIPT "averagelinkage_clu_dendro <- as.dendrogram(averagelinkage_clu)\n";
print RSCRIPT "fileConn<-file(\"leaforder_clusteredgenomes_averagelinkage.txt\")\n";
print RSCRIPT "writeLines(c(labels(averagelinkage_clu_dendro)), fileConn)\n";
print RSCRIPT "close(fileConn)\n\n";
print RSCRIPT "pdf(\"hierarchical_clustering_clusteredgenomes_completelinkage.pdf\", width = 9, height = 6)\n";
print RSCRIPT "plot_complete_clu <- plot(complete_clu, check = TRUE, cex = 0.25, main = \"Blast Query Hierarchical Clustering (Clustered genomes, Complete)\")\n";
print RSCRIPT "dev.off()\n\n";
print RSCRIPT "pdf(\"hierarchical_clustering_clusteredgenomes_averagelinkage.pdf\", width = 9, height = 6)\n";
print RSCRIPT "plot_average_clu <- plot(averagelinkage_clu, check = TRUE, cex = 0.25, main = \"Blast Query Hierarchical Clustering (Clustered genomes, Average)\")\n";
print RSCRIPT "dev.off()\n\n";
close(RSCRIPT);

system("Rscript rcode_count_hclust.R");
    
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sub ROUND
{	my $num = shift(@_);
	my $decimals = shift(@_);
	my $roundnum = int(($num * 10**$decimals) + 0.5)/(10**$decimals);
	return $roundnum;
}    

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -