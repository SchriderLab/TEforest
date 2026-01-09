# this shows an example command line used to run the mcclintock metapipeline in the paper
# reads were processed with fastp and downsampled to a target coverage prior to this analysis 
# using the TEforest pipeline. 

python mcclintock.py -p 32 -r <input_data_dir>/ISO1_GCF_000001215.4_Release_6_prepped.fasta \
    -c <input_data_dir>/consensusTEs.fasta -g <input_data_dir>/ISO1.gff3 -t <input_data_dir>/taxonomy.tsv \
    -1 <fq_dir>/AKA-017_GIM-024_1.fq -2 <fq_dir>/AKA-017_GIM-024_2.fq \ -o <output_dir>
    -a <input_data_dir>/consensusTEs.fasta -m popoolationte,popoolationte2,temp,temp2,teflon,retroseq,tepid