# Cisco-Linkage-Map
repository of various scripts used to generate the cisco linkage map

# Preparing Input Data for lepmap and rqtl2
## haploids
1. get the final populations.haplotype.vcf file from stacks
2. remove the consensus SNPs with removeConsensus.pl
3. remove the Ns: sed -E '/^#/! s/\t[^\t]*N[^\t]*/\t-/g' populations.haplotypes.tsv
4. run it through ryans pipeline: https://github.com/rwaples/ml-psv
5. grab the outputs and convert from mst to lepmap with haploid_mst_to_lepmap.py
6. transpose the output using the lepmap script using transpose_file.py
7. use combineLinkageFiles.pl to combine families into a one file
8. find duplicate loci using duplcate_renaming_chrom1.py and duplcate_renaming_chrom2.py (or more depending on how many families you used (I needed to make three of these scripts) and use duplicate_loci_renmaning_female_concenses_LOD_11_CHR_POS.txt as an example set up file 
9. can't do rqtl on haploids so run lepmap
10. LOD increment until it stabilizes and you can use graphing_diff_lod.R to look at it
11. view your map with Female_linkage_map_graphing_script.R
12. view your duplicated links with circle_plot.R (this was last minuet and poorly coded, sorry)

here is the script that was used to visualize the rank orders of all the salmonids with genomic resources with duplicate regions

## diploids
1. get final population.hap.vcf file from stacks
2. make sure you have all the input files to run grad_project_code_commented.py
   * this script asks if you have everything in the middle so you can start and stop depending on where you are at for intermediate file making with other programs
   * here is a summary of what this script does:
   * recodes genotypes
   * filters out unmappable loci
   * builds lepmap file
   * transposes file for new lepmap format 
   * builds rqtl file 
3. use combineLinkageFiles.pl to combine families into a one file
4. run lepmap
5. view your map with male_linkage_map_graphing_script.R
6. run rqtl with rqtl2_script_normalR.R or https://kbroman.org/qtl2/



