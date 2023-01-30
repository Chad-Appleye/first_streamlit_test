#import python modules
import argparse
import json
from pandas import DataFrame

#define positional and optional arguments and their helpful help messages
parser = argparse.ArgumentParser(description = """Generate LaTex results CSV files from a validation simulation performance results JSON object""")
parser.add_argument("prefix", type = str, help = "outfile prefix") #required, positional
parser.add_argument("summary_json_file_path", type = str, help = "path to validation simulation performance results JSON object") #required, positional
parser.add_argument("organism_class", type = str, help = "Organism class. Choices = [Bacteria, Fungi, Parasite, Virus]", choices = ["Bacteria", "Fungi", "Parasite", "Virus"]) #required, positional
parser.add_argument("nucleic_acid_type", type = str, help = "Nucleic acid type. Choices = [Rna or Dna]", choices = ["Rna", "Dna"]) #required, positional
parser.add_argument("--outfile_dir", type = str, help = "output file directory. default = current directory", default = ".") #optional
args = parser.parse_args()

#set outfile names
print_outfile_handle = open(args.outfile_dir + "/" + args.prefix + args.organism_class + args.nucleic_acid_type + "_val_sim_latex_tables.log", "w")
variables_outfile_handle = open(args.outfile_dir + "/" + args.prefix + args.organism_class + args.nucleic_acid_type + "Variables.csv", "w")
details_outfile_handle = open(args.outfile_dir + "/" + args.prefix + args.organism_class + args.nucleic_acid_type + "DetailsData.csv", "w")
fn_outfile = args.outfile_dir + "/" + args.prefix + args.organism_class + args.nucleic_acid_type + "FNegResultsData.csv" #file created only if fn results
major_fp_outfile = args.outfile_dir + "/" + args.prefix + args.organism_class + args.nucleic_acid_type + "MajorFPosResultsData.csv" #file created only if major fp results
minor_fp_outfile = args.outfile_dir + "/" + args.prefix + args.organism_class + args.nucleic_acid_type + "MinorFPosResultsData.csv" #file created only if minor fp results

#define functions
def format_latex_percentage(n):
    #convert percentage to float if percentage is no 0 % or 100 %
    if(n == 0 or n == 100):
        return int(n)
    return float(n)

def latex_fn_row_merge_mask(fn_nested_list):
    #initialize for masked fn results
    masked_fn_nested_list = []
    #convert to dataframe
    fn_df = DataFrame(fn_nested_list)
    #store boolean series of rows to apply mask operation to
    rows_to_mask = fn_df.duplicated(subset = 0, keep = "last") #FN Organism is first Column, keep last
    #for each row and index
    for i, row in enumerate(rows_to_mask):
        #if encounter a row to mask for FN Organism
        if(row):
            masked_fn_nested_list.append(["", "", "", "", "", fn_df.iloc[i,5]]) #include only the FNresult (column index 5)
        else:
            masked_fn_nested_list.append(fn_df.iloc[i].tolist())
    return masked_fn_nested_list

def latex_fp_row_merge_mask(fp_nested_list):
    #initialze for masked fp results
    masked_fp_nested_list = []
    #convert to dataframe
    fp_df = DataFrame(fp_nested_list)
    #store boolean series of rows to apply mask operation to
    rows_to_mask1 = fp_df.duplicated(subset = 0, keep = "last") #FP Organism is first Column, keep last
    rows_to_mask2 = fp_df.duplicated(subset = 6, keep = "last") #Simulated Organism is 6th Column, keep last
    #for each row and index
    for i, row1 in enumerate(rows_to_mask1):
        row2 = rows_to_mask2[i]
        #if encounter a row to mask for both FP Organism and Simulated Organism
        if(row1 and row2):
            masked_fp_nested_list.append(["", "", "", "", "", fp_df.iloc[i,5], "", "", fp_df.iloc[i,8]])
        #elif mask only FP Organism
        elif(row1):
            masked_fp_nested_list.append(["", "", "", "", ""] + fp_df.iloc[i,5:9].tolist())
        else:
            masked_fp_nested_list.append(fp_df.iloc[i].tolist())
    return masked_fp_nested_list

def write_nested_list_to_csv_file(list_of_lists, outfile_handle):
    for nested_list in list_of_lists:
        for item in nested_list:
            outfile_handle.write(str(item) + ",")
        outfile_handle.write("\n")

########################################################################################################################################################################

#get summary json dict
summary_json_dict = json.load(open(args.summary_json_file_path))

#store class_type, namespace, nucleic_acid of interest based on user-specified arguments
nucleic_acid_of_interest = "NA"
if(args.organism_class == "Bacteria"):
    class_type_of_interest = "bacterial"
    if(args.nucleic_acid_type == "Rna"):
        namespace_of_interest = "16S"
    elif(args.nucleic_acid_type == "Dna"):
        namespace_of_interest = "bg"
elif(args.organism_class == "Fungi"):
    class_type_of_interest = "fungal"
    if(args.nucleic_acid_type == "Rna"):
        namespace_of_interest = "18S"
    elif(args.nucleic_acid_type == "Dna"):
        namespace_of_interest = "fg"
elif(args.organism_class == "Parasite"):
    class_type_of_interest = "parasite"
    if(args.nucleic_acid_type == "Rna"):
        namespace_of_interest = "18S"
    elif(args.nucleic_acid_type == "Dna"):
        namespace_of_interest = "pg"
elif(args.organism_class == "Virus"):
    class_type_of_interest = "viral"
    namespace_of_interest = "viral"
    if(args.nucleic_acid_type == "Rna"):
        nucleic_acid_of_interest = "RNA"
    elif(args.nucleic_acid_type == "Dna"):
        nucleic_acid_of_interest = "DNA"

#initialize variables and latex nested list headers
var_n_org, var_n_pos, var_n_tpos, var_n_fneg_org, var_n_no_fneg_org, var_n_neg, var_n_tneg, var_n_fpos_org, var_n_no_fpos_org = [0]*9
#variables
variables_list_of_lists = []
variables_header_list_of_lists = [[args.nucleic_acid_type + "Org",
                                   args.nucleic_acid_type + "Accuracy",
                                   args.nucleic_acid_type + "PosCount",
                                   args.nucleic_acid_type + "TPosCount",
                                   args.nucleic_acid_type + "FNegOrg",
                                   args.nucleic_acid_type + "NoFNegOrg",
                                   args.nucleic_acid_type + "NoFNegOrgPercent",
                                   args.nucleic_acid_type + "Specificity",
                                   args.nucleic_acid_type + "NegCount",
                                   args.nucleic_acid_type + "TNegCount",
                                   args.nucleic_acid_type + "FPosOrg",
                                   args.nucleic_acid_type + "NoFPosOrg",
                                   args.nucleic_acid_type + "NoFPosOrgPercent"]]
#details                          
details_list_of_lists = []
details_header_list_of_lists = [["Organism",
                                 args.nucleic_acid_type.upper() + "_SensitiveThreshold",
                                 args.nucleic_acid_type.upper() + "_SpecificThreshold",
                                 args.nucleic_acid_type.upper() + "_TP",
                                 args.nucleic_acid_type.upper() + "_P",
                                 args.nucleic_acid_type.upper() + "_A",
                                 args.nucleic_acid_type.upper() + "_TN",
                                 args.nucleic_acid_type.upper() + "_N",
                                 args.nucleic_acid_type.upper() + "_S"]]
#fn
fn_list_of_lists = []
fn_header_list_of_lists = [["Organism",
                            "FN",
                            "P",
                            "FNrate",
                            "SpecificThreshold",
                            "FNresult"]]
#major fp                    
major_fp_list_of_lists = []
major_fp_header_list_of_lists = [["Organism",
                                  "FP",
                                  "N",
                                  "FPrate",
                                  "SpecificThreshold",
                                  "FPresult",
                                  "Sim_Organism",
                                  "Sim_SpecificThreshold",
                                  "Sim_result"]]
#minor fp
minor_fp_list_of_lists = []
minor_fp_header_list_of_lists = [["Organism",
                                  "FP",
                                  "N",
                                  "FPrate",
                                  "SpecificThreshold",
                                  "FPresult",
                                  "Sim_Organism",
                                  "Sim_SpecificThreshold",
                                  "Sim_result"]]

#loop through json results dicts and populate variables and latex nested lists
for rep_id in summary_json_dict.keys():
    #if has a meta dict                         ********** REMOVE IN FUTURE WHEN THIS CONDITIONAL ISN'T NEEDED ************
    if(summary_json_dict[rep_id].get("meta")):
        #store rep_id meta info
        class_type = summary_json_dict[rep_id]["meta"]["class_type"]
        namespace = summary_json_dict[rep_id]["meta"]["namespace"]
        nucleic_acid = summary_json_dict[rep_id]["meta"]["nucleic_acid"]
        specific_threshold = summary_json_dict[rep_id]["meta"]["specific_threshold"]
        sensitive_threshold = summary_json_dict[rep_id]["meta"]["sensitive_threshold"]
        #if class_type, namespace, nucleic_acid of interest
        if(class_type == class_type_of_interest and namespace == namespace_of_interest and nucleic_acid == nucleic_acid_of_interest):
            #store within-namespace performance results
            reporting_name = summary_json_dict[rep_id][namespace]["name"]
            TP_list = summary_json_dict[rep_id][namespace]["TP"]
            FN_list = summary_json_dict[rep_id][namespace]["FN"]
            TN_count = summary_json_dict[rep_id][namespace]["TN"]
            FP_list = summary_json_dict[rep_id][namespace]["FP"]
            #update variables
            var_n_org += 1
            var_n_pos += len(TP_list + FN_list)
            var_n_tpos += len(TP_list)
            var_n_neg += (TN_count + len(FP_list))
            var_n_tneg += TN_count
            if(len(FN_list) != 0):
                var_n_fneg_org += 1
            else:
                var_n_no_fneg_org += 1
            if(len(FP_list) != 0):
                var_n_fpos_org += 1
            else:
                var_n_no_fpos_org += 1
            #update details
            details_list_of_lists.append([reporting_name, 
                                          sensitive_threshold, 
                                          specific_threshold, 
                                          len(TP_list),
                                          len(TP_list + FN_list),
                                          format_latex_percentage(100*len(TP_list)/len(TP_list + FN_list)),
                                          TN_count,
                                          TN_count + len(FP_list),
                                          format_latex_percentage(100*TN_count/(TN_count+ len(FP_list)))])
            #update fn
            for i in range(0, len(FN_list)):
                fn_list_of_lists.append([reporting_name,
                                         len(FN_list),
                                         len(TP_list + FN_list),
                                         format_latex_percentage(100*len(FN_list)/len(TP_list + FN_list)),
                                         specific_threshold,
                                         format_latex_percentage(FN_list[i])])
            #count major vs minor fp (pass 1)
            FP_major_count = 0
            FP_minor_count = 0
            for i in range(0, len(FP_list)):
                #if major fp
                if(FP_list[i]["coverage"] > FP_list[i]["expected_coverage"]):
                    FP_major_count += 1
                else:
                    FP_minor_count += 1
            #update fp (pass 2)
            for i in range(0, len(FP_list)):
                #if major fp
                if(FP_list[i]["coverage"] > FP_list[i]["expected_coverage"]):
                    major_fp_list_of_lists.append([reporting_name,
                                                   FP_major_count,
                                                   TN_count + len(FP_list),
                                                   format_latex_percentage(100*FP_major_count/(TN_count + len(FP_list))),
                                                   specific_threshold,
                                                   FP_list[i]["coverage"],
                                                   FP_list[i]["expected_name"],
                                                   FP_list[i]["expected_specific_threshold"],
                                                   FP_list[i]["expected_coverage"]])    
                else:
                    minor_fp_list_of_lists.append([reporting_name,
                                                   FP_major_count,
                                                   TN_count + len(FP_list),
                                                   format_latex_percentage(100*FP_major_count/(TN_count + len(FP_list))),
                                                   specific_threshold,
                                                   FP_list[i]["coverage"],
                                                   FP_list[i]["expected_name"],
                                                   FP_list[i]["expected_specific_threshold"],
                                                   FP_list[i]["expected_coverage"]])

#write variables latex table
variables_list_of_lists.append([var_n_org,
                                format_latex_percentage(100*var_n_tpos/var_n_pos),
                                var_n_pos,
                                var_n_tpos,
                                var_n_fneg_org,
                                var_n_no_fneg_org,
                                format_latex_percentage(100*var_n_no_fneg_org/var_n_org),
                                format_latex_percentage(100*var_n_tneg/var_n_neg),
                                var_n_neg,
                                var_n_tneg,
                                var_n_fpos_org,
                                var_n_no_fpos_org,
                                format_latex_percentage(100*var_n_no_fpos_org/var_n_org)])
write_nested_list_to_csv_file(variables_header_list_of_lists + variables_list_of_lists, variables_outfile_handle)
#write details latex table
details_list_of_lists.sort()
write_nested_list_to_csv_file(details_header_list_of_lists + details_list_of_lists, details_outfile_handle)                        
#write fn latex table (IF fn results)
if(len(fn_list_of_lists) > 0):
    fn_list_of_lists.sort(key = lambda x: x[0]) #sort by FN organism
    masked_fn_list_of_lists = latex_fn_row_merge_mask(fn_list_of_lists)
    fn_outfile_handle = open(fn_outfile, "w")
    write_nested_list_to_csv_file(fn_header_list_of_lists + masked_fn_list_of_lists, fn_outfile_handle) 
#write major fp latex table (IF major fp results)
if(len(major_fp_list_of_lists) > 0):
    major_fp_list_of_lists.sort(key = lambda x: x[6]) #sort by simulated organism then by FP organism
    major_fp_list_of_lists.sort(key = lambda x: x[0])
    masked_major_fp_list_of_lists = latex_fp_row_merge_mask(major_fp_list_of_lists)
    major_fp_outfile_handle = open(major_fp_outfile, "w")
    write_nested_list_to_csv_file(major_fp_header_list_of_lists + masked_major_fp_list_of_lists, major_fp_outfile_handle)
#write minor fp latex table (IF minor fp results)
if(len(minor_fp_list_of_lists) > 0):
    minor_fp_list_of_lists.sort(key = lambda x: x[6]) #sort by simulated organism then by FP organism
    minor_fp_list_of_lists.sort(key = lambda x: x[0])
    masked_minor_fp_list_of_lists = latex_fp_row_merge_mask(minor_fp_list_of_lists)
    minor_fp_outfile_handle = open(minor_fp_outfile, "w")
    write_nested_list_to_csv_file(minor_fp_header_list_of_lists + masked_minor_fp_list_of_lists, minor_fp_outfile_handle)   

#write to log file
print("For", args.organism_class, "and", args.nucleic_acid_type, ":", file = print_outfile_handle)
print("The class type to summarize performance results for is:", class_type_of_interest, file = print_outfile_handle)
print("The namespace to summarize performance results for is:", namespace_of_interest, "\n", file = print_outfile_handle)

#close outfile handles
print_outfile_handle.close()
variables_outfile_handle.close()
details_outfile_handle.close()
