import collections, os, string, sys, bisect, numpy 
from collections import Counter
from collections import OrderedDict


os.getcwd()
os.chdir("F:/consensus_files_for_lepmap/")


chr_file = open("duplicate_loci_renmaning_female_concenses_LOD_11_CHR_POS.txt", "r")
sys.stdout = open("chrom2_data.txt", "w")
line = chr_file.readline()
split = line.rstrip().split("\t")

families = [split[2], split[4], split[6], split[8], split[10], split[12]]
x1_46 = split[2]
x2_46 = split[4]
x1_47 = split[6]
x2_47 = split[8]
x1_95 = split[10]
x2_95 = split[12]

temp = sys.stdout

NA = "#N/A"


while True:
    line = chr_file.readline()
    if not line: break
    line_split = line.rstrip().split("\t")
    chrom = []

    #grabbing families with a SNP on a chrom and dropping out the ones that have no
    #data for the SNP
    if line_split[1] != NA:
        chrom1 = line_split[1]
        chrom.append(chrom1)
    if line_split[3] != NA:
        chrom2 = line_split[3]
        chrom.append(chrom2)
    if line_split[5] != NA:
        chrom3 = line_split[5]
        chrom.append(chrom3)
    if line_split[7] != NA:
        chrom4 = line_split[7]
        chrom.append(chrom4)
    if line_split[9] != NA:
        chrom5 = line_split[9]
        chrom.append(chrom5)
    if line_split[11] != NA:
        chrom6 = line_split[11]
        chrom.append(chrom6)

    #putting the unique choms for that snp into a list
    unique_chrom = []
    used_chrom = [unique_chrom.append(x) for x in chrom if x not in unique_chrom]

    v1 = []
    v2 = []
    v3 = []



    #using the the unique list of choms that the snp shows up on, this will make a list of the
    #pos that the snp occurs at for the chrom
    try:
        t1 = unique_chrom[0]
        t2 = unique_chrom[1]
        t3 = unique_chrom[2]
        if line_split[1] == t1:
            pos1 = line_split[2]
            v1.append(pos1)
        elif line_split[1] == t2:
            pos1 = line_split[2]
            v2.append(pos1)
        elif line_split[1] == t3:
            pos1 = line_split[2]
            v3.append(pos1)
    except IndexError:
        try:
            t2 = unique_chrom[1]
            if line_split[1] == t2:
               pos1 = line_split[2]
               v2.append(pos1)
        except IndexError:
            try:
                t3 = unique_chrom[2]
                if line_split[1] == t3:
                    pos1 = line_split[2]
                    v3.append(pos1)
            except IndexError:
                pass

    try:
        t1 = unique_chrom[0]
        t2 = unique_chrom[1]
        t3 = unique_chrom[2]
        if line_split[3] == t1:
            pos2 = line_split[4]
            v1.append(pos2)
        elif line_split[3] == t2:
            pos2 = line_split[4]
            v2.append(pos2)
        elif line_split[3] == t3:
            pos2 = line_split[4]
            v3.append(pos2)
    except IndexError:
        try:
            t2 = unique_chrom[1]
            if line_split[3] == t2:
                pos2 = line_split[4]
                v2.append(pos2)
        except IndexError:
            try:
                t3 = unique_chrom[2]
                if line_split[3] == t3:
                    pos2 = line_split[4]
                    v3.append(pos2)
            except IndexError:
                pass
            
    try:
        t1 = unique_chrom[0]
        t2 = unique_chrom[1]
        t3 = unique_chrom[2]
        if line_split[5] == t1:
            pos3 = line_split[6]
            v1.append(pos3)
        elif line_split[5] == t2:
            pos3 = line_split[6]
            v2.append(pos3)
        elif line_split[5] == t3:
            pos3 = line_split[6]
            v3.append(pos3)
    except IndexError:
        try:
            t2 = unique_chrom[1]
            if line_split[5] == t2:
                pos3 = line_split[6]
                v2.append(pos3)
        except IndexError:
            try:
                t3 = unique_chrom[2]
                if line_split[5] == t3:
                    pos3 = line_split[6]
                    v3.append(pos3)
            except IndexError:
                pass
    
    try:
        t1 = unique_chrom[0]
        t2 = unique_chrom[1]
        t3 = unique_chrom[2]
        if line_split[7] == t1:
            pos4 = line_split[8]
            v1.append(pos4)
        elif line_split[7] == t2:
            pos4 = line_split[8]
            v2.append(pos4)
        elif line_split[7] == t3:
            pos4 = line_split[8]
            v3.append(pos4)
    except IndexError:
        try:
            t2 = unique_chrom[1]
            if line_split[7] == t2:
                pos4 = line_split[8]
                v2.append(pos4)
        except IndexError:
            try:
                t3 = unique_chrom[2]
                if line_split[7] == t3:
                    pos4 = line_split[8]
                    v3.append(pos4)
            except IndexError:
                pass

    try:
        t1 = unique_chrom[0]
        t2 = unique_chrom[1]
        t3 = unique_chrom[2]
        if line_split[9] == t1:
            pos5 = line_split[10]
            v1.append(pos5)
        elif line_split[9] == t2:
            pos5 = line_split[10]
            v2.append(pos5)
        elif line_split[9] == t3:
            pos5 = line_split[10]
            v3.append(pos5)
    except IndexError:
        try:
            t2 = unique_chrom[1]
            if line_split[9] == t2:
                pos5 = line_split[10]
                v2.append(pos5)
        except IndexError:
            try:
                t3 = unique_chrom[2]
                if line_split[9] == t3:
                    pos5 = line_split[10]
                    v3.append(pos5)
            except IndexError:
                pass
        
    try:
        t1 = unique_chrom[0]
        t2 = unique_chrom[1]
        t3 = unique_chrom[2]
        if line_split[11] == t1:
            pos6 = line_split[12]
            v1.append(pos6)
        elif line_split[11] == t2:
            pos6 = line_split[12]
            v2.append(pos6)
        elif line_split[11] == t3:
            pos6 = line_split[12]
            v3.append(pos6)
    except IndexError:
        try:
            t2 = unique_chrom[1]
            if line_split[11] == t2:
                pos6 = line_split[12]
                v2.append(pos6)
        except IndexError:
            try:
                t3 = unique_chrom[2]
                if line_split[11] == t3:
                    pos6 = line_split[12]
                    v3.append(pos6)
            except IndexError:
                pass


    #these need to be real numbers so we can bin the positions that are within a certain range of each other
    ########################################################################################################
    #this must get changed!
    ########################################################################################################
    v2 = map(float, v2)
    #print v2

    #for this first pass I used 20 cM bins (feels a little lenient but it can always be changed)
    score_ranges = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360]
    binning = []
    for a in v2:
        binning.append(bisect.bisect_right(score_ranges, a))
    #instead of the index i might be able to get the pos?
        
    #score_ranges.insert(bisect(score_ranges, a), a)
    #the binning works by taking the values in v1 and changeing them to the repsective bin it falls in (i.e if
    #v1 is in the 0-19 bin it will be given a 1, if v1 falls in 20-39 it will be giving a value of 2, ect)
    #then we can call the unique vlaues in the binned list and get how many true versons we have on each chrom
    unique_bin = []
    used_bin = [unique_bin.append(x) for x in binning if x not in unique_bin]
        
    if line_split[2] != NA:
        line_split[2] = float(line_split[2])
    if line_split[4] != NA:
        line_split[4] = float(line_split[4])
    if line_split[6] != NA:
        line_split[6] = float(line_split[6])
    if line_split[8] != NA:
        line_split[8] = float(line_split[8])
    if line_split[10] != NA:
        line_split[10] = float(line_split[10])
    if line_split[12] != NA:
       line_split[12] = float(line_split[12])


    ########################################################################################################
    #this must get changed!
    ########################################################################################################
    dicts = collections.OrderedDict()
       
    if line_split[2] in v2:
        chrom1 = line_split[2]
        dicts[x1_46] = chrom1
    if line_split[4] in v2:
        chrom2 = line_split[4]
        dicts[x2_46] = chrom2
    if line_split[6] in v2:
        chrom3 = line_split[6]
        dicts[x1_47] = chrom3
    if line_split[8] in v2:
        chrom4 = line_split[8]
        dicts[x2_47] = chrom4
    if line_split[10] in v2:
        chrom5 = line_split[10]
        dicts[x1_95] = chrom5
    if line_split[12] in v2:
        chrom6 = line_split[12]
        dicts[x2_95] = chrom6

    
    binning_fam = {}
    for value, key in dicts.items():
        binning_fam[(bisect.bisect_right(score_ranges, key))] = value
    print binning_fam


sys.stdout.close()
sys.stdout = temp

