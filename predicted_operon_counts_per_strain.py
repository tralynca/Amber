
## script name: predicted_operon_counts_per_strain.py # WAS predicted_operons_counts_per_lineage.py

#######################################################################################

# script that goes through all COSMO predicted files

## OUTPUT
# Output file with extention "_calls_&_cov_all_lineages.txt" : The Prediction (TP, FP, FN) **PER OPERON** for all isolates/strains
# Output file with extention "_counts_per_isolate" : Prediction **PER STRAIN ** with count of how many TPs, FPs and FNs were predicted per strain
    # -The file also outputs precision, recall and F1 scores

# Requirements:
# 1. The csv files produced by COSMO containing the operon predictions
# 2. The operon list file containing all Mtb operons and their positions and genes

########################################################################################

import string
import os.path
import glob
import math



dir_pred_operon_files = input("\nPlease add the path to the directory containing all the csv files produced by COSMO and then press ENTER\n e.g. /home/Amber/Desktop/Input_files/COSMO/COSMO_*f6.csv* \n\nThe '*.csv' part, tells this Python script that all of the input files have different names, but they all end in '.csv'\n\n")
# dir_pred_operon_files = "/home/tracey/Desktop/tracey/Input_files/COSMO/COSMO_*f6.csv"

valid_operon_list = input("\n\nPlease add the path to the file with your verified operons from RegulonDB, and then press ENTER\n e.g. /home/Amber/Desktop/OperonSet_Exp_strong_with_OFFICIAL_gene_names.csv \n\n")
# valid_operon_list = "/home/tracey/Desktop/Amber/OperonSet_Exp_strong_with_OFFICIAL_gene_names.csv"
pred_operons_path = glob.glob(dir_pred_operon_files)

print("#######################################")
for path in pred_operons_path:
    print("reading in file: ", path)
# print(pred_operons_path)

# create output file names
outfile_name = dir_pred_operon_files.split("/")[-1] + "_calls_&_cov_all_lineages.txt" #.split(".")[1]
outfile_name_counts = dir_pred_operon_files.split("/")[-1] +  "_counts_per_isolate.txt"

outfile = open(outfile_name, "w")
outfile_counts = open(outfile_name_counts, "w")

##############################################################
# add the lines from each file with VALIDATED operons to a list
#############################################################
operon_lists = []

header = "Call" + "\t" + "Operon" + "\t" + "Predicted_operon" + "\t" + "Operon_length" + "\t" + "Pred_operon_length" + "\t" + "Coverage" + "\t" + "Strain" + "\n\n"
outfile.write(header)

################################################################################
# Format Validated operon file

# Operon_name	First_gene	Last_gene	start_coordinate	end_

################################################################################

# operon the file containing the verified operons:
# Format of file above
open_operon_files = open(valid_operon_list)
for operon_line in open_operon_files:
    if not operon_line.startswith("Operon"):
        operon_line = operon_line.strip().split(",")
        operon_lists.append(operon_line) # Format = [['b3196 - b3206', 'yrbG', 'npr', '3340275', '3348238', '+', '11']]

###############################################################
# add the lines from each file with PREDICTED operons to a list
##############################################################

Total_TPs, Total_FNs, Total_FN_counterpart, Total_FPs, Total_FP_counterpart = [],[],[],[],[]

for pred_operon_files in pred_operons_path:

    TP_list, TP_operon_counterpart = [],[]
    FP_list, FP_operon_counterpart = [],[]
    FN_list, FN_operon_counterpart = [],[]
    FN_single_gene_list, FN_single_gene_counterpart = [], []
    pred_operon_list = []

################################################################
    # Open the COSMO prediction files:
    # Format: Operon,start,end,length,Strand,coverage,Number_of_Genes,Genes, cov_each_gene, cov_each_IGR
    #        b0002 - b0005	337	5530	5193	+	69.3	4
	#						b0002	75.17
	#						b0003	89.43	157.5

################################################################

    pred_operon_filename  = os.path.basename(pred_operon_files)#.split(".")[0]
    strain = pred_operon_filename.rstrip(".csv")
    pred_operons = open(pred_operon_files)
    for row in pred_operons:
        row = row.strip().split(",")
        pred_operon_list.append(row) # Format = ['b4396', '4634441', '4635310', '869', '-', '605.12', '1']
# print(pred_operon_list)

    # the validated operon list
    for operon_column in operon_lists:
        operon_name = operon_column[0]
        # print (operon_name)
        operon = operon_column[3] + " - " + operon_column[4]
        # print (operon)
        operon_len = operon_column[-1]
        operon_start = int(operon_column[3])
        operon_end = int(operon_column[4])

# ###########################################################
# # Loop through the list containing the real operons and
# # check if they were predicted by algorithm
# ###########################################################

        for elem in pred_operon_list: # Format of 'pred_operon_list' = ['b4396', '4634441', '4635310', '869', '-', '605.12', '1']
            if elem[0]:
                pred_operon_name = elem[0]
                pred_operon = elem[1] + " - " + elem[2] #coordinates
                pred_operon_start = int(elem[1])
                pred_operon_end = int(elem[2])
                pred_operon_len = (pred_operon_end - pred_operon_start +1)
                pred_operon_single = elem[0]
                pred_operon_single_start = int(elem[1])
                pred_operon_single_end = int(elem[2])
                pred_operon_len = ((int(pred_operon_end) - int(pred_operon_start)) +1)
            if elem[4]:
                strand = elem[4]
            if elem[5]:
                coverage = elem[5]
            if elem[6]:
                pred_gene_count = elem[6]
                # print(pred_gene_count)
                # break

# ##################################################################
# # if the real operon matches the predicted operon, print TPs
    # ##################################################################

                if operon == pred_operon:
                    # print(operon, pred_operon, strain, operon_len, pred_gene_count, coverage)
                    TP_list.append(operon_name)
                    TP_operon_counterpart.append(operon_name) # not really necessary, but just keeping with format
                    Total_TPs.append(operon_name)
                    output_TP = "TP" + "\t" + operon_name + "\t" + pred_operon_name + "\t" + str(operon_len) + "\t" + str(pred_gene_count) + "\t" + coverage + "\t" + strain + "\n"
                    outfile.write (output_TP)
                    # print(output_TP)
#
    # ################################################################################################
    #         # if the opredicted operon is longer than the validated operon, print FPs
    #         # the entire operon must be contained in it through
    #         # the predicted operon must contain all the genes of the real operon
    #         # AND have extra genes up/dowstream of the real operon
#     # #################################################################################################
                elif operon != pred_operon\
                and elem [0] \
                and not "EBG" in elem[0] \
                and pred_operon_end > pred_operon_start:

                    if not coverage == "NOT EXPRESSED" and int(pred_gene_count) > 1 :
                        if (pred_operon_start < operon_start and pred_operon_end > operon_end) \
                        or (pred_operon_start < operon_start and operon_end == pred_operon_end) \
                        or (pred_operon_start == operon_start and pred_operon_end > operon_end):
                            FP_list.append(pred_operon_name)
                            FP_operon_counterpart.append(operon_name) # keep track of the verified operons
                            Total_FPs.append(pred_operon_name) # keep track of the total no of FPs
                            Total_FP_counterpart.append(operon_name)
                            output_FP = "FP" + "\t" + operon_name + "\t" + pred_operon_name + "\t" + str(operon_len) + "\t" + str(pred_gene_count) + "\t" +  coverage + "\t" + strain + "\n"
                            outfile.write (output_FP)
                            Total_FPs.append(operon_name)
                            # print(output_FP)
                            # break
#############################################################################################
              # if the predicted operon is shorter than it should be,
              # but it's still longer than 1 gene, print FNs
#############################################################################################
#
                        elif (pred_operon_start > operon_start \
                        and pred_operon_end < operon_end) \
                        or (pred_operon_start > operon_start \
                        and operon_end == pred_operon_end) \
                        or (pred_operon_start == operon_start and pred_operon_end < operon_end) :

                            FN_list.append(pred_operon_name)
                            Total_FNs.append(pred_operon_name)
                            FN_operon_counterpart.append(operon_name)
                            Total_FN_counterpart.append(operon_name)
                            output_FN = "FN" + "\t" + operon_name + "\t" + pred_operon_name + "\t" + str(operon_len) + "\t" + str(pred_gene_count) + "\t" +  coverage + "\t" + strain + "\n"
                            outfile.write (output_FN)

#############################################################################################
              # if the predicted operon is a single gene
              # that is EXPRESSED, print FNs,
              # but state that it's a single gene
#############################################################################################

                    elif elem [0] and not "-" in elem[0]:
                        if not coverage == "NOT EXPRESSED"\
                        and pred_operon_single_start == operon_start \
                        or pred_operon_single_start == operon_end \
                        or pred_operon_single_end  == operon_start \
                        or pred_operon_single_end  == operon_end :
                            pred_operon_gene_len_single = "1"

                            FN_single_gene_list.append(pred_operon_single)
                            FN_list.append(pred_operon_single)
                            Total_FNs.append(pred_operon_single)
                            FN_single_gene_counterpart.append(operon_name)
                            FN_operon_counterpart.append(operon_name)
                            Total_FN_counterpart.append(operon_name)
                            output_FN_single_gene = "FN_single_gene" + "\t" + operon_name + "\t" + pred_operon_single + "\t" + str(operon_len) + "\t" + pred_operon_gene_len_single + "\t" +  coverage + "\t" + strain + "\n"
                            outfile.write (output_FN_single_gene)
                            # print(output_FN_single_gene )
#############################################################################################
              # if the predicted operon is a single gene
              # that is UNEXPRESSED, print FNs,
              # but state that it's a single UNEXPRESSED gene
#############################################################################################

                        elif coverage == "NOT EXPRESSED" \
                        and pred_operon_single_start == operon_start \
                            or pred_operon_single_start == operon_end \
                            or pred_operon_single_end  == operon_start \
                            or pred_operon_single_end  == operon_end:
                                pred_operon_single_unexp = elem[0]
                                pred_operon_gene_len_single_unexp = "0"
                                # pred_operon_single_start = int(elem[1])
                                # pred_operon_single_end = int(elem[2])
                                FN_single_gene_list.append(pred_operon_single_unexp)
                                FN_list.append(pred_operon_single_unexp)
                                Total_FNs.append(pred_operon_single_unexp)
                                FN_single_gene_counterpart.append(operon_name)
                                FN_operon_counterpart.append(operon_name)
                                Total_FN_counterpart.append(operon_name)
                                output_FN_unexpressed = "FN_NOT_EXPRESSED" + "\t" + operon_name + "\t" + pred_operon_single_unexp + "\t" + str(operon_len) + "\t" + pred_operon_gene_len_single_unexp + "\t" +  coverage + "\t" + strain + "\n"
                                outfile.write (output_FN_unexpressed)
                            # outfile.write("\n\n")
                                # print(output_FN_unexpressed)
        #
# ################################################################################
# # write the exact operon that was predicted for the strain to a files
# # together with the call: eg.g TP, FP, FN
# ################################################################################

    header_total_TPs = f"Predictions made by COSMO per isolate {strain} \n"
    strain_title = strain + "\n\n"
    outfile_counts.write(strain_title)

    outfile_counts.write(header_total_TPs)
    for TP in set(Total_TPs):
        out_TP = TP + "\t" + "TP" + "\t" + strain + "\n"
        outfile_counts.write(out_TP)
    outfile_counts.write("\n")


    for FP in set(Total_FPs):
        out_FP = FP + "\t" + "FP" + "\t" + strain + "\n"
        outfile_counts.write(out_FP)
    outfile_counts.write("\n")

    for FN in set(Total_FNs):
        out_FN = FN + "\t" + "FN" + "\t" + strain + "\n"
        outfile_counts.write(out_FN)
    outfile_counts.write("\n\n\n")
#
##################################################################
# store only unique calls by converting a list to a set
#################################################################

    set_TPs = set(TP_list)

    set_FPs = set(FP_list)
    set_FP_counterpart = set(FP_operon_counterpart)

    set_FNs = set(FN_list)
    set_FN_counterpart = set(FN_operon_counterpart)

#####################################################################
# Calculate the Algorithm's performance metrics
# and write the counts of the calls to a file
# also write the performance metrics to the file
#####################################################################

    correct_FN_count = len(operon_lists) - (len(set_TPs) + len(set_FPs))
    # print(correct_FN_count)
    Total_correct  = math.ceil((len(set_TPs)/len(operon_lists))*100)
    PPV = math.ceil((len(set_TPs)/(len(set_TPs) + len(set_FP_counterpart)))*100)
    recall = math.ceil((len(set_TPs)/(len(set_TPs) + correct_FN_count))*100)
    F1 = math.ceil((2*((PPV * recall/(PPV + recall)))))
    # print(Total_correct,PPV, recall, F1)
    # print(set_TPs, len(set_TPs))

    TP_counts = "TP count"  + " = " + str(len(set_TPs))  + "\t" + strain +  "\t" + "\n"
    FP_counts =  "FP count"  + " = " + str(len(set_FP_counterpart)) + "\t" + strain + "\t" +"\n"
    FN_counts = "FN count"  + " = " + str(len(set_FN_counterpart)) + "\t" + strain + "\t" + "\n"
    corrected_FN_count = "corrected FN count"  + " = " + str(correct_FN_count) + "\t" + strain + "\t" + "\n"

    Total_correct_perc = "Total correct"  + " = " + str(Total_correct) + "%"  + "\t" + strain +  "\t" + "\n"
    PPV_output =  "PPV/Precision"  + " = " + str(PPV) + "%" + "\t" + strain +  "\t" + "\n"
    recall_output = "Sensitivity/Recall"  + " = " + str(recall)  + "%" + "\t" + strain +  "\t" + "\n"
    F1_output = "F1 score"  + " = " + str(F1)  + "%" + "\t" + strain +  "\t" + "\n"


    # write evrything to files
    header_raw_counts = "Raw counts:\n"
    header_total_calls = "Summary of predictions made by COSMO per Isolate/strain:\n"
    outfile_counts.write(header_raw_counts)
    outfile_counts.write(header_total_calls)
    outfile_counts.write(TP_counts)
    outfile_counts.write(FP_counts)
    outfile_counts.write(FN_counts)
    outfile_counts.write(corrected_FN_count)

    outfile_counts.write("\n\n")
    header_stats = "Statistics:\n"
    outfile_counts.write(header_stats)
    outfile_counts.write(Total_correct_perc)
    outfile_counts.write(PPV_output)
    outfile_counts.write(recall_output)
    outfile_counts.write(F1_output )
    end_line = "##################################################################################\n\n"
    outfile_counts.write("\n\n")
    outfile_counts.write(end_line)

print("\nThe file:", outfile_name, "has been successfully printed", "to the directory:", os.getcwd() + "\n")
print("The file:", outfile_name_counts, "has been successfully printed", "to the directory:", os.getcwd() + "\n")

outfile.close()
outfile_counts.close()
