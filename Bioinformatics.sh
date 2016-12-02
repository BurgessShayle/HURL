#!/bin/bash
# tell terminal waht it's getting fed

# fasta file header: isolate sample name in first field of header line, use program "sed" stream editor, s= substitution mode
sed 's/::/ /' YostF.fna > YostF_name.fna


# reformat fasta IDs for qiime using "deal with texas" script HP
python deal_with_texas.py -f YostF_name.fna -o YostF_relabeled.fna
# writing python first tells it to use py.

source activate qiime

## Add few lines to subset samples of interest
filter_fasta.py -f YostF_relabeled.fna -o TripSamples.fna --sample_id_fp SamplesTrips.txt

# Start at Trim primer sequences using (RC below)


# First, subset only ITSD sequences using awk (google it),save it as new file
##awk '/\-ITSD/ {print; getline; print}' YostF_relabeled.fna >ITSD.fasta


# Trim primer sequences using cutadapt
# Trim forward primers (allow error rate of 15% (3 indels/mismatches)) and keep only trimmed sequences
cutadapt -g GTGAATTGCAGAACTCCGTG -e 0.15 TripSamples.fna --trimmed-only -o trimF.fasta
    # cutadapt is the call, --trimmed-only, on/off <-discares all read that don't contain adapter
    # -o is the output file, give it a name 


# Trim forward primers again (may be multiple internal primer sequences), but do not discard sequences that do not contain primer
cutadapt -g GTGAATTGCAGAACTCCGTG -e 0.15 trimF.fasta -o trimF2.fasta


# Trim reverse primers using cutadapt and keep only trimmed sequences
cutadapt -a AAGCATATAAGTAAGCGGAGG -e 0.15 trimF2.fasta --trimmed-only -o trimF2_trimR.fasta

# Trim reverse primers again for any remaining internally and remove (-m) sequences shorter than 250 bp
cutadapt -a AAGCATATAAGTAAGCGGAGG -e 0.15 trimF2_trimR.fasta -m 250 -o combined_seqs_trimmed.fasta



# After al trimmed: Cluster by (use whichever script makes sense)



# Cluster entire dataset at 97% similarity, s=similarity threshold, i=input, o=output directory (new folder)
pick_otus.py -i combined_seqs_trimmed.fasta -s 0.97 --optimal_uclust -o otus_97/

# Remove singletons before picking rep set and assigning taxonomy (to save time)
awk '$3 ~ /./ {print}' otus_97/combined_seqs_trimmed_otus.txt > otus_97/nosingles_otus.txt
# ^ look in field 3 and see if you see anything. ie first field is denovo, second field is the first OTU, third field would be either empty or not...
# to check new number of OTU, do >wc which tells you number of lines, aka number of OTUs left

# Get rep set of 97% OTUs
pick_rep_set.py -i otus_97/nosingles_otus.txt \
-f combined_seqs_trimmed.fasta \
-o otus_97/97_otus_rep_set.fasta
# to view: <less 97_otus_rep_set.fasta


# Make OTU table
make_otu_table.py -i otus_97/nosingles_otus.txt -o otus_97/97_otus.biom
rm -f  otus_97/97_otus.tsv  # delete old OTU .tsv if present bc you can't overwrite the file
biom convert -i otus_97/97_otus.biom -o otus_97/97_otus.tsv --to-tsv

### Faux blasting
R --vanilla <run_nw.R --args otus_97/97_otus_rep_set.fasta ITS2db_trimmed_derep.fasta




# Cluster entire dataset at 100% similarity, s=similarity threshold, i=input, o=output directory (new folder)
pick_otus.py -i combined_seqs_trimmed.fasta -s 1.0 --optimal_uclust -o otus_100/

#FIX THIS STUFF BELOW!!!!!!!!!!!!!!!!!
# Remove singletons before picking rep set and assigning taxonomy (to save time)
awk '$3 ~ /./ {print}' otus_100/combined_seqs_trimmed_otus.txt > otus_100/nosingles_otus.txt
# ^ look in field 3 and see if you see anything. ie first field is denovo, second field is the first OTU, third field would be either empty or not...
# to check new number of OTU, do >wc which tells you number of lines, aka number of OTUs left

# Get rep set of 100% OTUs
pick_rep_set.py -i otus_100/nosingles_otus.txt \
-f combined_seqs_trimmed.fasta \
-o otus_100/100_otus_rep_set.fasta
# to view: <less 100_otus_rep_set.fasta


# Make OTU table
make_otu_table.py -i otus_100/nosingles_otus.txt -o otus_100/100_otus.biom
rm -f  otus_100/100_otus.tsv  # delete old OTU .tsv if present bc you can't overwrite the file
biom convert -i otus_100/100_otus.biom -o otus_100/100_otus.tsv --to-tsv

### Faux blasting
R --vanilla <run_nw.R --args otus_100/100_otus_rep_set.fasta ITS2db_trimmed_derep.fasta


## 97% within sample
~/SymITS2/otus_97_bysample.sh combined_seqs_trimmed.fasta otus_97_bysample
								# data                     # where to save it
								# mergers all 100% similar ones that were already grouped by 97%
#global taxonomy
R --vanilla <~/SymITS2/run_nw.R --args otus_97_bysample/all_rep_set_rep_set.fasta ITS2db_trimmed_derep.fasta


