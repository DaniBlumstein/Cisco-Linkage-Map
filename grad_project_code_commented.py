# script parse vcf file get alleles per individual
import os, string, sys, csv
os.getcwd()
os.chdir("Z:/Geog441/dblum940/Grad_project/")
family = raw_input("Enter family name here ")

# open vcf file read all lines into an array, close file
raw_vcf_file = open(family + "_populations.haps.vcf", "r")
raw_vcf_array = raw_vcf_file.readlines()
raw_vcf_file.close()

# open file you are going to write new GENOTYPE CALL results to
out_file = open(family + "_parsed_genotype_calls.txt", "w")
# for each line in vcf file
for i in raw_vcf_array:
    # this is trying to recognize the header line
    if i.startswith("#CHROM"):
        # hard code column names you want
            out_file.write("Chrom" + "\t")
            header_line = i.rstrip().split("\t")
            # grabbing individual names
            z = 0
            for j in header_line:
                z = z + 1
                #there are nine headers before the data 
                if z > 9 and z < len(header_line)+1:
                    out_file.write(j + "\t")
            out_file.write("\n")
            # at this point we should have header line

    # use the if not # to go to the data
    elif "#" not in i:
        # splits each thing into their respective cells by tabs
        split_ind_line = i.rstrip().split("\t")
        # outputs locus name, you're going to need to output stuff from other columns too
        out_file.write(split_ind_line[0] + "\t")
        # iterates through individuals
        z = 0
        for j in split_ind_line:
            z = z + 1
            if z > 9:
                # split the genotype cell
                gen_data = j.split("/t")
                zygosity = gen_data[0]
                if zygosity == "0/0":
                    genotype_to_print = "A/A"
                elif zygosity == "0/1":
                    genotype_to_print = "A/B"
                elif zygosity == "1/1":
                    genotype_to_print = "B/B"
                elif zygosity == "0/2":
                    genotype_to_print = "A/C"
                elif zygosity == "1/2":
                    genotype_to_print = "B/C"
                elif zygosity =="2/2":
                    genotype_to_print = "C/C"
                else:
                    genotype_to_print = "0"
                out_file.write(genotype_to_print + "\t")
        out_file.write("\n")
out_file.close()

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

#logic, need to read in diploid genotypes, create matrix of possible genotypes (punnent square)
#and get the unique values from that matrix
#then match progeny genotypes to that matrix to determine which allele came from mom or dad
#if offspring has an impossible genotype code as missing data, also need to keep track of proportion of missing data, missing data coded as dashes
#"rescue" loci where one parent is missing, we will only do this in loci with 2 alleles for simplicity (only losing 3 loci from HX13 this way)
#to rescue we need to calculate the proportion of AA, BB and AB in the offspring, take the homo genotype with the highest frequency and if the freq
#of the top homo and the het are > .4 we can infer the missing genotype i.e if parent with genotype is AB then other missing should be the 
#common homo and if parent with genotype is homo then missing should be het.  

#assume C and D loci will not be included with new genotyping script so take them out
#script only works on 1 family at a time, will write another script to combine families

####filters
#remove loci with < X prop missing data
#make really low...? or .25
prop_missing_cutoff = raw_input("Remove loci with < X proportion of missing data." + "\n" + "Suggested is 0.25" + "\n" + "Type your X here: ")
#remove loci with < X prop het, all mappable loci should have prop het >= 25%
min_prop_het_cutoff = raw_input("Remove loci with < X proportion of heterozygotes" + "\n" + "Note: all mappable loci should have a proportion of heterozygotes >= 25%" + "\n" + "Suggested is 0.1" + "\n" + "Type your X here: ")
#remove loci with > X prop het, all mappable loci should have prop het <= 50%
max_prop_het_cutoff = raw_input("Remove loci with > X prop heterozygote" + "\n" + "Note: all mappable loci should have prop a proportion of heterozygotes <= 50%" + "\n" + "Suggested is 0.75" + "\n" + "Type your X here: ")

prop_missing_cutoff = float(prop_missing_cutoff)
min_prop_het_cutoff = float(min_prop_het_cutoff)
max_prop_het_cutoff = float(max_prop_het_cutoff)

#read in file, make each line an arrary entry
diploid_genotypes_file = open(family + "_parsed_genotype_calls.txt","r")
diploid_genotypes_array = diploid_genotypes_file.readlines()
diploid_genotypes_file.close()

#test_locus 
test_locus=diploid_genotypes_array[49].rstrip().split("\t")

#step 1: for each locus do punnet square of parents and get uniques
header_line = diploid_genotypes_array[0]
header_array = header_line.rstrip().split("\t")
num_inds = len(header_array)-1
ind_names_list = header_array[0:num_inds+1]

def get_unique_alleles(genotypes):
    unique_alleles =[]
    for i in genotypes:
        #if missing don't put into alleles list
        if "-" in i:
            continue
        #if het split and check to see if one or both alleles in list
        #if homo and not in unique_alleles list put in list
        elif i not in unique_alleles:
            unique_alleles.append(i)
    #sort alleles alphabetically
    unique_alleles_alpha=sorted(unique_alleles)
    return(unique_alleles_alpha)

def recode_gens_as_nums(genotype):
    new_gen = ""
    if genotype == "A/A":
        new_gen="1 0 0 0 0 0 0 0 0 0"
    if genotype == "A/B":
        new_gen="0 1 0 0 0 0 0 0 0 0"
    if genotype == "B/B":
        new_gen="0 0 0 0 1 0 0 0 0 0"
    return(new_gen)        

#functions reads in locus + genotypes from parents and offspring in array
def code_gens(locus_array):
    #make a punnett square with all possible genotypes
    #split parental genotypes
    split_fem_gen=locus_array[2].split("/")
    split_male_gen=locus_array[1].split("/")
    fem_gen1=split_fem_gen[0]
    fem_gen2=split_fem_gen[1]
    male_gen1=split_male_gen[0]
    male_gen2=split_male_gen[1]
    possible_genotypes_fwd = [fem_gen1 + "/" + male_gen1, fem_gen1 + "/" + male_gen2, fem_gen2 + "/" + male_gen1, fem_gen2 + "/"+male_gen2]
    #need to reverse gens to make sure none are missed
    possible_gens_rev = [possible_genotypes_fwd[0][::-1], possible_genotypes_fwd[1][::-1], possible_genotypes_fwd[2][::-1], possible_genotypes_fwd[3][::-1]]
    possible_genotypes = possible_genotypes_fwd+possible_gens_rev

    #get unique alleles
    unique_alleles = []
    for i in possible_genotypes:
            #if missing don't put into alleles list
        if "0" in i:
            continue
        elif i not in unique_alleles:
            unique_alleles.append(i)

    #iterate over offspring genotypes add to male and female if possible, if not possible make missing and append, if already missing append
    #initialize arrays with locus name
    updated_gens = [locus_array[0]]
    updated_gens.append(recode_gens_as_nums(locus_array[1]))
    updated_gens.append(recode_gens_as_nums(locus_array[2]))
    
    for i in range(3, num_inds + 1):
        #if missing append 0 0
        if locus_array[i] == "0":
            updated_gens.append("0 0 0 0 0 0 0 0 0 0")
            continue
        if locus_array[i] not in unique_alleles:
            updated_gens.append("0 0 0 0 0 0 0 0 0 0") 
            continue
        #if genotype is found add to output
        for j in unique_alleles:
            if locus_array[i] == j:
                updated_gens.append(recode_gens_as_nums(locus_array[i]))
    #returning updated genotypes
    return(updated_gens)    

#
####function to calculate proportion of different genotypes in offspring
def get_gen_props(current_locus_input):
    missing_gens = 0
    AA_gens = 0
    BB_gens = 0
    AB_gens = 0
    for i in range(3, num_inds):
        if current_locus_input[i] == "0":
            missing_gens+= 1 
        elif current_locus_input[i] == "A/A":
            AA_gens+= 1
        elif current_locus_input[i] == "B/B":
            BB_gens+= 1
        elif current_locus_input[i] == "A/B":
            AB_gens+= 1
    total_gens=AA_gens+BB_gens+AB_gens
    total_gens_with_missing = AA_gens + BB_gens + AB_gens + missing_gens
    prop_AA = 0
    prop_BB = 0
    prop_AB = 0
    prop_missing = 0
    if total_gens > int(num_inds) * 0.5:
        prop_AA = round(float(AA_gens)/float(total_gens), 3)
        prop_BB = round(float(BB_gens)/float(total_gens), 3)
        prop_AB = round(float(AB_gens)/float(total_gens), 3)
        prop_missing = round(float(missing_gens)/float(total_gens_with_missing), 3)

    gen_props = [prop_AA, prop_AB, prop_BB, prop_missing]
    return(gen_props)

#iterate through loci if one or both of genotypes missing skip, will handle later
#if it is mappable, call function to store locus
good_locus_dict = {}
for i in diploid_genotypes_array:
    current_locus = i.rstrip().split("\t")
    current_locus_name = current_locus[0]
    #skip first row
    if "CHR" in current_locus[0]:
        continue
    female_gen = current_locus[1]
    male_gen = current_locus[2]
    #if "C" or "D" in current_locus skip
    more_than_2_alleles = "FALSE"
    for i in current_locus:
        if "C" in i or "D" in i:
            more_than_2_alleles = "TRUE"
            break
    if(more_than_2_alleles == "TRUE"):
        continue
    #if neither parent has a het not mappable so skip
    if "A/B" not in female_gen and "A/B" not in male_gen:
        #print(current_locus)
        continue
    #skip if missing gens for a parent, will handle later
    if "0" in female_gen or "0" in male_gen:
        continue
    #look at proportion of missing data and prop het, if prop missing < .8 and 
    #prop het > .1 and prop het < .8 it is a good locus and should be mappable
    prop_gens_current_locus = get_gen_props(current_locus)
    prop_missing_current_locus = prop_gens_current_locus[3]
    prop_het_current_locus = prop_gens_current_locus[1]
    if prop_missing_current_locus < prop_missing_cutoff and prop_het_current_locus > min_prop_het_cutoff and prop_het_current_locus < max_prop_het_cutoff:  
        recoded_gens=code_gens(current_locus)
        #final check to make sure once genotypes are recoded nothing went wrong
        num_missing_updated = 0
        for j in recoded_gens:
            if "0 0 0 0 0 0 0 0 0 0" in j:
              num_missing_updated = num_missing_updated + 1
        prop_missing_updated = round(float(num_missing_updated)/float(num_inds), 3)
        if prop_missing_updated < prop_missing_cutoff:
            good_locus_dict[current_locus[0]] = code_gens(current_locus)  

##############build lepmap output file
lepmap_out_file = open(family + "_lepmap.linkage","w")

list_of_inds = ind_names_list = header_array[1:num_inds + 1]

#this is sorted by not necessarily in proper numeric order
list_of_snps_in_order = sorted(good_locus_dict.keys())
#write header line
lepmap_out_file.write("CHR\tCHR\tCHR\tCHR\tCHR\tCHR\t")
for i in list_of_snps_in_order:
    lepmap_out_file.write(i + "\t")
lepmap_out_file.write("\n")

#write full file
male=list_of_inds[0]
female=list_of_inds[1]

#write male line
lepmap_out_file.write(family + "\t"+male+"\t"+"0"+"\t"+"0"+"\t"+"2"+"\t"+"0"+"\t")
for i in list_of_snps_in_order:
    current_gens = good_locus_dict[i]
    lepmap_out_file.write(current_gens[1]+"\t")
lepmap_out_file.write("\n")

#write female line
lepmap_out_file.write(family+"\t"+female+"\t"+"0"+"\t"+"0"+"\t"+"1"+"\t"+"0"+"\t")
for i in list_of_snps_in_order:
    current_gens = good_locus_dict[i]
    lepmap_out_file.write(current_gens[2]+"\t")
lepmap_out_file.write("\n")

#write offspring genotypes
for i in range(2,num_inds):
    ind_name = list_of_inds[i]
    lepmap_out_file.write(family+"\t" + ind_name + "\t" + male + "\t" + female + "\t" + "0" + "\t" + "0" + "\t")
    for j in list_of_snps_in_order:
        current_gens = good_locus_dict[j]
        #print(current_gens)
        lepmap_out_file.write(current_gens[i + 1] + "\t")
        #print(current_gens)
    lepmap_out_file.write("\n")
        
lepmap_out_file.close()

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

with open(family + "_lepmap.linkage") as infile, open(family + "_lepmap.csv",'w') as outfile: 
    for line in infile: 
        outfile.write(line.replace('\t',','))
        
from itertools import izip
lepmap_csv = zip(*csv.reader(open(family + "_lepmap.csv","rb")))
genotypes_file1 = csv.writer(open(family + "_genotypes_final1.csv", "wb")).writerows(lepmap_csv)


genotypes_file2 = family + "_genotypes_final1.csv"
genotypes_file3 = family + "_genotypes_final1.txt"

with open(genotypes_file3, "w") as my_output_file:
    with open(genotypes_file2, "r") as my_input_file:
        [ my_output_file.write("\t".join(row)+'\n') for row in csv.reader(my_input_file)]
    my_output_file.close()

genotypes_file3 = open(family + "_genotypes_final1.txt", "r")
genotypes_file4 = open(family + "_genotypes_final.txt", "w")

raw_vcf_array = genotypes_file3.readlines()
genotypes_file3.close()

i = 0
for line in raw_vcf_array:
    i = i + 1
    if i == 2:
        genotypes_file4.write(line)
    elif i == 7:
        genotypes_file4.write(line)
    elif i > 7:
        genotypes_file4.write(line)

genotypes_file4.close()

replacements = {"1 0 0 0 0 0 0 0 0 0" : "AA", "0 1 0 0 0 0 0 0 0 0": "AB", "0 0 0 0 1 0 0 0 0 0" : "BB", "0 0 0 0 0 0 0 0 0 0" : "-"}
lines = []
with open(family + "_genotypes_final.txt") as infile:
    for line in infile:
        for src, target in replacements.iteritems():
            line = line.replace(src, target)
        lines.append(line)
with open(family + "_genotypes_final.txt", 'w') as outfile:
    for line in lines:
        outfile.write(line)

os.remove(family + "_genotypes_final1.txt")
os.remove(family + "_genotypes_final1.csv")
os.remove(family + "_lepmap.csv")
os.remove(family + "_parsed_genotype_calls.txt")

#--------------------------------------------------------------------------------------------------------------------------------------------------

def get_unique_from_list(input_list):
    unique_values = []
    for i in input_list:
        #if missing don't put into alleles list
        if i not in unique_values:
            unique_values.append(i)
    #sort alleles alphabetically
    unique_values_alpha = sorted(unique_values)
    return(unique_values_alpha)

user_answer = raw_input("I have add the data to the folder so for the purpose of this project the answer is yes" + "\n" + "Have you run the external lepmap program yet, yes or no: ")
if user_answer != "yes":
    print "Please go run the external lepmap program."
    sys.exit()
if user_answer == "yes":

user_answer = raw_input("I have add the data to the folder so for the purpose of this project the answer is yes" + "\n" + "Have you added the phenotypic data to the folder, yes or no: ")
if user_answer != "yes":
    print "Please add the phenotypic data to the folder."
    sys.exit()
if user_answer == "yes":

    #read in diplod genotypes (slightly modified output of recode_haplotype_file_as_diploid_gens.py)
    #make sure ind names match phenotypes file
    genotypes_file = open(family + "_genotypes_final.txt","r")
    genotypes_array = genotypes_file.readlines()
    genotypes_file.close()

    #read in phenotypes file, file is a tab delimited text file with a header row 
    #of form individual, phenotype 1, phenotype 2, ... etc and should
    #take as many phenotypes as desired, make sure individual IDs match diploid genotypes
    phenotypes_file = open(family + "_phenotypes_for_qtl.txt", "r")
    phenotypes_array = phenotypes_file.readlines()
    phenotypes_file.close()

    #read in map file, file is a tab delimited text file with a header row "Tag", "lg","cM" with
    #markers in map order, make sure to order in excel if necessary
    #run unique on map to get rid of duplicate tag names
    map_file = open("linkage_map_for_rqtl.txt","r")
    map_array = map_file.readlines()
    map_file.close()

    #open ouput file in csv format for rqtl
    rqtl_out_file = open(family + "_rqtl_input.csv", "w")


    #basically make a structure that makes a dictionary with the key each individual and the value a list 
    #containing all the pertinant information about that individual, start with phenotypes, then go through marker info
    #but only for markers on the map

    #first subset the genotypes file by the markers that are on the map, need list of markers on the map, dict of marker on map info
    linkage_map_dict = {}
    markers_on_map = []
    for i in map_array:
        if "Tag" in i:
            continue
        else:
            map_line = i.rstrip().split("\t")
            map_locus = map_line[0]
            markers_on_map.append(map_locus)
            lg = map_line[1]
            cM = map_line[2]
            map_list = [lg,cM]
            linkage_map_dict[map_locus] = map_list       
        
        
    #subset genotypes file by markers on map
    #these markers are not in order, make a dictionary of key=marker name value=list of genotypes for all individuals
    z = 0
    loci_and_genotype_dict = {}
    list_of_genotyped_inds = []
    for i in genotypes_array:
        z+= 1
        genotype_line = i.rstrip().split("\t")
        if z == 1:
            list_of_genotyped_inds = (genotype_line[1:len(genotype_line)])
        #if marker is on map keep
        elif genotype_line[0] in markers_on_map:
            locus_name = genotype_line[0]
            progeny_genotypes = genotype_line[1:len(genotype_line)]
            loci_and_genotype_dict[locus_name] = progeny_genotypes
            
    num_inds = len(list_of_genotyped_inds)   

    #start making dictionary for each individual, add phenotypes to dict first, get good inds from phenotypes file
    master_dict_by_ind = {}
    list_of_phenotypes = []
    list_of_phenotyped_inds = []
    z = 0
    for i in phenotypes_array:
        z+= 1
        phenotype_line = i.rstrip().split("\t")
        num_phenotypes = len(phenotype_line) - 1
        phenotyped_ind = phenotype_line[0]
        if z == 1:
            list_of_phenotypes = phenotype_line[1:num_phenotypes + 1]
        else:
            master_dict_by_ind[phenotyped_ind ] = phenotype_line[1:num_phenotypes + 1]
            list_of_phenotyped_inds.append(phenotyped_ind)
        
    #now work on genotypes, need to get them in map order (easy part because can look up with dict), and
    #need to subset by inds that were genotyped using a which type function that I write
    indices_of_phenotyped_inds_in_genotypes = []
    for i in list_of_phenotyped_inds:
        z = 0
        for j in list_of_genotyped_inds:
            if i == j:
                indices_of_phenotyped_inds_in_genotypes.append(z)
            z+= 1
                  
    #add genotypes onto master dictionary
    for i in markers_on_map:
        if i in loci_and_genotype_dict.keys():
            current_locus = loci_and_genotype_dict[i]
            for j in indices_of_phenotyped_inds_in_genotypes:
                current_ind_index = j
                current_ind = list_of_genotyped_inds[current_ind_index]
                current_genotype=current_locus[current_ind_index]
                #append new genotype
                master_dict_by_ind[current_ind].append(current_genotype)

    #build output file
    #build header line 1, count ind id line and num of phenotypes, these will be blank in rows 2 and 3
    blank_spaces = 1
    rqtl_out_file.write("Ind_ID" + ",")
    for i in list_of_phenotypes:
        blank_spaces+=1
        rqtl_out_file.write(i + ",")
    for i in markers_on_map:
        if i in loci_and_genotype_dict.keys():
            rqtl_out_file.write(i + ",")
    rqtl_out_file.write("\n")

    #build header line 2
    for i in range(0,blank_spaces):
        rqtl_out_file.write(",")
    for i in markers_on_map:
        if i in loci_and_genotype_dict.keys():
            map_info = linkage_map_dict[i]
            rqtl_out_file.write(map_info[0] + ",")
    rqtl_out_file.write("\n")

    #build header line 3
    for i in range(0,blank_spaces):
        rqtl_out_file.write(",")
    for i in markers_on_map:
        if i in loci_and_genotype_dict.keys():
            map_info = linkage_map_dict[i]
            rqtl_out_file.write(map_info[1] + ",")
    rqtl_out_file.write("\n")

    #build main file
    for i in list_of_phenotyped_inds:
        rqtl_out_file.write(i + ",")
        ind_info = master_dict_by_ind[i]
        for i in ind_info:
            rqtl_out_file.write(i + ",")
        rqtl_out_file.write("\n")
rqtl_out_file.close()
