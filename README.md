# BTAB (Blast TAB-delimited) tools

A collection of tools for filtering and parsing the blast btab format (outfmt 6). When you call the program use ```-h``` to see the help information.

Tools:

1. ```btab_besthitfilter.pl```
2. ```btab_result_parser_v3.pl```
3. ```btab_result_parser_addon_queryclustering_v2.pl```
4. ```btab_result_parser_addon_prep_percWquery_forR.pl```

I developed these scripts for use with some Anvio protein cluster files I was using, so some are specifically designed to parse contig names from those files (mentioned in the script header when you ```less``` the script itself).

The generic input is from a blastp run, with the default ```-outfmt 6``` settings.

It is helpful to run these tools in order. To run Tool 3, you must have the output from Tool 2. To run Tool 4, you must have the output from Tool 3. 

Tool 1: This is probably the most widely useable script. The only goal is to filter out the best hit of a subject (i.e. contig/gene from a genome) that hits multiple queries (i.e. multiple genes of interest). I have a database w/ proximal and distal versions of cbb3, and wanted a fast way to filter the blast results given that genes that were cbb3 would hit both of the cbb3 queries. The options allow the user to filter based on strict e-value, magnitude e-value, or percent identity (pid). The user can also select the maximum e-value/pid to report.

Tool 2 takes a btab output and parses it into a table with query on the y-axis and genome/bin on the x-axis. The assumption for headers is: ```PC_#######_#####_genomename```, where the first set of numbers is the Anvio protein cluster id, and the second is the unique protein id. Output is very useful for humans to read! Much better than the btab format itself.

Tool 3 takes the output from Tool 2, and a file showing clustering of genomes (see help information). The output is a series of pdfs showing the hierarchical clustering of queries based on all genomes (separated) and genome clusters. This creates and runs an R script, so you should have R installed. Additionally, you get an output table showing the percentage of each genome cluster that contains each query.

Tool 4 converts the table from Tool 3 (```*clusteredgenome_percentWquery.txt```) into something that can be imported into R to generate figures.
