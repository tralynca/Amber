# script_name: add_official_gene_name_to_EVO_list.py

# Author: Tracey Calvert-Joshua - tralynca@gmail.com - Oct 2021

##############################################################################
import os.path, glob, itertools

EVO_list_path = "/home/tracey/Desktop/Amber/OperonSet_Exp_strong.txt"
out_file_name  = os.path.basename(EVO_list_path).split(".txt")[0] +  "_with_OFFICIAL_gene_names.txt"
Heading  = "Operon_name" + " " + "First_gene" + " " + "Last_gene" + " " + "start_coordinate" + " " + "end_coordinate" + " " + "strand" + " " + "no_of_genes" + "\n"

# print(out_file_name)
out_file = open(out_file_name, "w")
out_file.write(Heading)


ptt_file_list = []
open_ptt_file = open("/home/tracey/Desktop/Amber/GCA_000005845.2_ASM584v2_genomic.gb.ptt")
for ptt_line in itertools.islice(open_ptt_file, 4 , None): # start at line 4; stop = None
    ptt_line = ptt_line.strip()
    ptt_file_list.append(ptt_line)
# print(ptt_file_list)

EVO_list = []
open_EVO_list = open(EVO_list_path)
for line in open_EVO_list:
    line = line.strip()
    EVO_list.append(line)
# print(EVO_list)
#
gene_name_to_symbol = dict()
for ptt_element in ptt_file_list:
    gene_name_ptt = ptt_element.strip().split("\t")[4] # this is the name e.g. rpoB
    gene_symbol_ptt = ptt_element.strip().split("\t")[5] # this is the Rv symbol name
    gene_name_to_symbol[gene_name_ptt] = gene_symbol_ptt

# print(gene_name_to_symbol)

for list in EVO_list:
    # print(list)
    item = list.strip().split("\t")
    strand = item[3]
    if strand =="forward":
        strand = "+"
    elif strand == "reverse":
        strand = "-"
    start, end, operon_start, operon_end, no_of_genes = item[1], item[2], item[5].split(",")[0], item[5].split(",")[-1], item[4]
    # print(operon_start, operon_end, strand, no_of_genes, genes_HGNC)
    mix_genes = item[5] # this is a group of genes making op operon, sep by ","
    mix_gene = mix_genes.strip().split(",")
    # print(mix_genes, mix_gene)

    for i, gene in enumerate(mix_gene):
        gene = gene.strip()
        if gene in gene_name_to_symbol:
            mix_gene[i] = gene_name_to_symbol[gene] #.replace(gene, gene_symbol_gb)
            # print(mix_gene)
        else:
            mix_gene[i] = gene
        # print(mix_gene)
    if strand == "+":
        new_operon_name = f"{mix_gene[0]} - {mix_gene[-1]}"
        out_line = f"{new_operon_name} {operon_start} {operon_end} {start} {end} {strand} {no_of_genes}\n"
        # print(out_line)
        out_file.write(out_line)
    elif strand == "-":
    # print(new_operon_name)
        new_operon_name = f"{mix_gene[-1]} - {mix_gene[0]}"
        out_line = f"{new_operon_name} {operon_end} {operon_start} {start} {end} {strand} {no_of_genes}\n"
        # print(out_line)
        out_file.write(out_line)

open_ptt_file .close()
open_EVO_list.close()

print("\nThe output file called: ", out_file_name, "\nWas successfully printed to:", os.getcwd() + "\n")
