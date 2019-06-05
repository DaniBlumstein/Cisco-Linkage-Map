# Cisco-Linkage-Map
repository if various scripts and things used to make the linkage map

## Goals
To outline how to generate plots of similarities between tetrasomically pairing regions of salmonid genomes using any high-quality assembly.

## Files
  * stacks output

# Preparing Input Data for lepmap and rqtl2
## haploids

1. get the final populations.haplotype.vcf file from stacks
2. remove the Ns
3. run it trough ryans pipeline
4. grab the outputs
5. transpose the output using the lepmap script
6. use garretts script to combine everthing into a giant file
7. add the shit duplicate renaming scripts?
8. can't do rqtl on haploids so run lepmap
9. LOD incremint until it stabilizes and you can use the graphing script to look at it
9. graph with the graphing script
8. align with blast and look at the outputs ect
9. circle plot?


## diploids
1. get final population.hap.vcf file from stacks
2. make sure you have all the input files
3. run it through the crazy script thing that does all
4. transpose the lepmap with the lepmap thing
5. use garretts script to combine everthing into a giant file
6. run lepmap
7. run rqtl


