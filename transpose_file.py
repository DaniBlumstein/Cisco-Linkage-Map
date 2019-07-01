#Part four: this part will transpose the lepmap file so it is in the right orientation for rqtl. this will also remove the lepmap specific lines so we are left with
#only the marker names and the individual genotypes at each marker

#this line opens the lepmap output file and re writes it as a csv file
with open(family + "_lepmap.linkage") as infile, open(family + "_lepmap.csv",'w') as outfile: 
    for line in infile: 
        outfile.write(line.replace('\t',','))
        
from itertools import izip
#then the csv is read and a new csv is transposed and written
#there are other ways to do this without converting  to a csv as I've shown above but I wanted to use less lines and the csv package allows for that
lepmap_csv = zip(*csv.reader(open(family + "_lepmap.csv","rb")))
genotypes_file1 = csv.writer(open(family + "_genotypes_final1.csv", "wb")).writerows(lepmap_csv)

#the file is then reconverted into a text file because text files are easier to manipulate in python
with open(family + "_genotypes_final1.csv") as infile, open(family + "_genotypes_final1.txt",'w') as outfile: 
    for line in infile: 
        outfile.write(line.replace(',','\t'))