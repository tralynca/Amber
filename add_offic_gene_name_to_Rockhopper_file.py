# script_name: add_official_gene_name_to_Rockhopper_list.py
# Author: Tracey Calvert-Joshua - tralynca@gmail.com - Oct 2021

## USE ##
# script to reorganize the operon prediction file from Rockhopper_list
# to look similar to that of COSMO
# We added  column containing the official gene names of the Operon
# because Rockhopper uses the HGNC symbols

##############################################################################
import os.path, glob, itertools

print("\nPlease add the path to the directory file produced by Rockhopper and then press ENTER\n e.g. /home/amber/Desktop/Input_files/Input_files/Rockhopper_control.txt \n\n")
Rockhopper_list_path = input()

# Rockhopper_list_path = "/home/tracey/Desktop/Amber/Input_files/Rockhopper/Rockhopper_control.txt"


out_file_name  = os.path.basename(Rockhopper_list_path).split(".txt")[0] +  "_replaced.txt"
Heading  = "Operon" + "\t" + "Start" + "\t" + "End" + "\t" + "Length" + "\t" + "Strand" + "\t" + "Number_of_Genes" + "\t" + "Genes"  + "\n"
print(out_file_name)
out_file = open(out_file_name, "w")
out_file.write(Heading)

##################################################################################
# Format of file:
# Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
# 190..255	+	66	AAC73112.1	thrL	b0001	-	-	thr operon leader peptide
################################################################################

# open the .ptt file from which we will retrieve the official gene name
# store the lines to a list

ptt_file_list = []
open_ptt_file = open("/home/tracey/Desktop/Amber/GCA_000005845.2_ASM584v2_genomic.gb.ptt")
for ptt_line in itertools.islice(open_ptt_file, 4 , None): # start at line 4; stop = None
    ptt_line = ptt_line.strip()
    ptt_file_list.append(ptt_line)
# print(ptt_file_list)

##################################################################################
# Format of file:
# Operon_name First_gene Last_gene start_coordinate end_coordinate strand no_of_genes
# b3196 - b3206 yrbG npr 3340275 3348238 + 11
##################################################################################

#operon the validated ioperon list from RegulonDB and save the lines to a list

Rockhopper_list = []
open_Rockhopper_list = open(Rockhopper_list_path)
for line in open_Rockhopper_list:
    line = line.strip()
    if line.startswith("Start"):
        continue
    Rockhopper_list.append(line)
# print(Rockhopper_list[0])
    # break
#
gene_name_to_symbol = dict()
for ptt_element in ptt_file_list:
    gene_name_ptt = ptt_element.strip().split("\t")[4] # this is the name e.g. rpoB
    gene_symbol_ptt = ptt_element.strip().split("\t")[5] # this is the Rv symbol name
    gene_name_to_symbol[gene_name_ptt] = gene_symbol_ptt

# print(gene_name_to_symbol)
    # break
#
for list in Rockhopper_list: # format: 337	5020	+	3	thrA, thrB, thrC
    item = list.strip().split("\t")
    strand = item[2]
    start, end, operon_start, operon_end = item[0], item[1], item[4].split(",")[0], item[4].split(",")[-1]
    no_of_genes =  item[3]
    length = int(end) -int(start) + 1
    # print(operon_start, operon_end, strand, no_of_genes, length) #genes_HGNC)
    # break
    mix_genes = item[4] # this is a group of genes making op operon, sep by ","
    mix_gene = mix_genes.strip().split(",")
    # print(mix_genes,mix_gene)
    # break
    for i, gene in enumerate(mix_gene):
        gene = gene.strip()
        if gene in gene_name_to_symbol:
            mix_gene[i] = gene_name_to_symbol[gene] #.replace(HGNC_gene, official_gene_symbol)
        else:
            mix_gene[i] = gene

    if strand == "+" or strand == "-":
        new_operon_name = f"{mix_gene[0]} - {mix_gene[-1]}"
        out_line = f"{new_operon_name}\t{start}\t{end}\t{length}\t{strand}\t{mix_genes}\t{no_of_genes}\n"
        # print(out_line) # similar to COSMO, except that the coverage is replaced with mix_genes

        out_file.write(out_line)


open_ptt_file .close()
open_Rockhopper_list.close()

print("\nThe output file called: ", out_file_name, "\nWas successfully printed to:", os.getcwd() + "\n")
