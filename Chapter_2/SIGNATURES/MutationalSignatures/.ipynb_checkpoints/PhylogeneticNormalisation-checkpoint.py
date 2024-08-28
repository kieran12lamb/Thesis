# script.py
import argparse
import pandas as pd
import numpy as np
import math
import multiprocessing
import compiled_functions as cf
from Bio import SeqIO
import json
from Bio import AlignIO
import warnings
warnings.filterwarnings('ignore')
from datetime import datetime

def position_checker(positions):
    # print(positions)
    inner_mutations = mutations[mutations.POS.isin(positions)]
    mutation_array = []
    # Loop through the each position
    inner_mutations["HasParent"] = False
    # print(inner_mutations)
    for pos in inner_mutations.POS.unique():
        # Extract rows for a single position
        pos_mutations = inner_mutations[inner_mutations.POS == pos]
        # Loop through each mutation at the position
        for alt_nuc in pos_mutations.ALT.unique():
            pos_alt_mutations = pos_mutations[pos_mutations.ALT == alt_nuc]
            for index in pos_alt_mutations.index:
                row = pos_alt_mutations.loc[index]
                # Split the lineage of the current row
                lin = row.lineage.split(".")
                # Loop up through the lineage levels i.e B.1.1.7->B.1.1->B.1->B
                for j in range(len(lin)-1,0,-1):
                    parent_lineage = '.'.join(lin[0:j])
                    # print(".".join(lin),parent_lineage)
                    # If a parent lineage is also at this position 
                    if parent_lineage in pos_alt_mutations.lineage.values:
                        pos_mutations.loc[index,"HasParent"] = True
        mutation_array.append(pos_mutations[pos_mutations.HasParent ==False])
    mutation_array = pd.concat(mutation_array,axis=0)
    return mutation_array

def inverted_position_checker(positions):
    # print(positions)
    inner_mutations = mutations.loc[mutations.POS.isin(positions)]
    mutation_array = []
    # Loop through the each position
    pos_time =  datetime.now()
    for pos in list(inner_mutations.POS.unique()):
        # Extract rows for a single position
        pos_mutations = inner_mutations.loc[inner_mutations.POS == pos]
        pos_mutations["Keep"] = True
        # If mutation is unique, keep it
        # keep = [pd.DataFrame()]
        # Loop through each mutation at the position
        for lineage in list(pos_mutations.lineage.unique()):
            # print(lineage)
            # print(pos_mutations.lineage.unique())
            #Get lineage rows that are not the currrent lineage, and are below in the hierarchy
            pos_lineage_mutations = pos_mutations.loc[(pos_mutations.lineage != lineage) & (pos_mutations.lineage.str.contains(lineage) == True)]
            # print(pos_lineage_mutations.lineage.unique())
            #Get all mutant nucleotides from current lineage
            mutant_nucleotides = pos_mutations.loc[pos_mutations.lineage == lineage].ALT.unique()
            #Get direct descendants of current lineage
            lineage_depths = [len(split_lineage.split(".")) for split_lineage in pos_lineage_mutations.lineage]
            pos_lineage_mutations.loc[:,"depths"] = lineage_depths
            direct_descendants = pos_lineage_mutations.loc[pos_lineage_mutations["depths"]<len(lineage.split("."))+1]
            for mutant in mutant_nucleotides:
                #Check if any mutations are in the directly descendant lineages
                keep_mutations = False
                
                #If a lineage has no descendants keep its mutations
                if len(direct_descendants.lineage.unique()) == 0:
                    keep_mutations = True
                #If a lineage has decendants only keep mutations that are inherited
                else:
                    for descendant_lineage in list(direct_descendants.lineage.unique()):
                        #If mutant is in any direct descendant of current lineage 
                        if mutant == str(pseudo_references[descendant_lineage])[pos]:
                            #Set Keep value to true if mutation becomes a direct descendant 
                            keep_mutations = True
                            break
                #Identify mutations to get rid of
                if keep_mutations == False:
                    # #Check if those that arent, the proportion of sequences that have the mutation
                    # proportions = mutations_with_duplicates[(mutations_with_duplicates.lineage == lineage) & (mutations_with_duplicates.ALT == mutant)]
                    # total_sequences_with_mutant = len(proportions.loc[proportions["ALT"] ==mutant])
                    # #If there are more than 10 sequences with the mutation, keep the mutation
                    # if total_sequences_with_mutant>=10:
                    #     pos_mutations['Keep'].loc[(pos_mutations['lineage'] == lineage) & (pos_mutations['ALT'] == mutant)] = True
                    # else:
                    pos_mutations['Keep'].loc[(pos_mutations['lineage'] == lineage) & (pos_mutations['ALT'] == mutant)] = False
        print(f"Removed {len(pos_mutations[pos_mutations['Keep'] == False])}/{len(pos_mutations[pos_mutations['Keep'] == False])+len(pos_mutations[pos_mutations['Keep'] == True])}")
        mutation_array.append(pos_mutations[pos_mutations['Keep'] == True])
        # print(f'POS {pos} Finished in Start Time:{datetime.now()-pos_time}')
    mutation_array = pd.concat(mutation_array,axis=0)
    return mutation_array


def position_ratio_checker(inner_mutations):
    # Loop through the each position
    inner_mutations["Remove"] = False
    inner_mutations.index = inner_mutations[["REF","ALT","POS","lineage"]].astype(str).sum(axis=1)
    for pos in inner_mutations.POS.unique():
        # Extract rows for a single position
        pos_mutations = inner_mutations[inner_mutations.POS == pos]
        # Loop through each mutation at the position
        for alt_nuc in pos_mutations.ALT.unique():
            pos_alt_mutations = pos_mutations[pos_mutations.ALT == alt_nuc]
            for index in pos_alt_mutations.index:
                row = pos_alt_mutations.loc[index]
                # Split the lineage of the current row
                lin = row.lineage.split(".")
                current_lineage = row.lineage
                current_lineage_alias = convert_lineage_name(row.lineage,lineage_aliases)
                #Get most recent parent (Should never really loop except for sequences where parents are not in dataset)
                for j in range(len(lin)-1,0,-1):
                    parent_lineage = '.'.join(lin[0:j])
                    parent_lineage_alias = convert_lineage_name(parent_lineage,lineage_aliases)
                    #If lineage has sequences with the mutation
                    if parent_lineage in pos_alt_mutations.lineage.values:
                        nucleotides = pseudo_reference_counts[pseudo_reference_counts.lineage ==current_lineage_alias][str(pos)]
                        nucleotides = nucleotides/nucleotides.sum()
                        nucleotides.index = pseudo_reference_counts[pseudo_reference_counts.lineage ==current_lineage_alias].nucleotide
                        parent_nucleotides = pseudo_reference_counts[pseudo_reference_counts.lineage ==parent_lineage_alias][str(pos)]
                        parent_nucleotides = parent_nucleotides/parent_nucleotides.sum()
                        parent_nucleotides.index = pseudo_reference_counts[pseudo_reference_counts.lineage ==parent_lineage_alias].nucleotide
                        #Check that lineages have sequences in the dataset i.e B.1.617 has no sequences 
    #                     print(alt_nuc,pos_mutations.loc[index].name)
                        if nucleotides.size != 0 and parent_nucleotides.size != 0:
                            if nucleotides.loc[alt_nuc]>=0.5 and parent_nucleotides.loc[alt_nuc] >=0.5:
                                inner_mutations.loc[index,"Remove"] = True
                                print(pos_mutations.loc[index])
                            break
                        #If there are no sequences in the lineage continue to next parent lineage
                        else:
                            print(current_lineage_alias,nucleotides)
                            print(parent_lineage_alias,parent_nucleotides)
                    #If parent lineage does not have this mutation, stop checking parents
                    else:
                        break
    mutation_array = inner_mutations[inner_mutations.Remove ==False].reset_index(drop=True)
    print(mutation_array)
#     mutation_array = pd.concat(mutation_array,axis=0).reset_index(drop=True)
    return mutation_array


def convert_lineage_name(lineage,alias_table):
    if lineage.split(".")[0] in alias_table:
        return alias_table[lineage.split(".")[0]]+"."+".".join(lineage.split(".")[1:])
    else:
        if lineage.count(".")>3:
            for key in alias_table.keys():
                if alias_table[key] == ".".join(lineage.split(".")[:4]):
                    return key+"."+".".join(lineage.split(".")[4:])
    return lineage
            
                


def build_reference_dict(pseudo_references):
    pseudo_references = AlignIO.read(pseudo_references, "fasta")

    lineages_dictionary = {}
    for lineage_record in pseudo_references:
        lineages_dictionary[lineage_record.id] = lineage_record.seq
        split_lineage = lineage_record.id.split(".")
        if split_lineage[0] in lineage_aliases:
            split_lineage = lineage_aliases[split_lineage[0]].split(".")+split_lineage[1:]
            lineages_dictionary[".".join(split_lineage)] = lineage_record.seq
    return lineages_dictionary

def get_lineage_emergence(lineage_chunk_map):
    emergence_dict = {}
    for lineage in lineage_chunk_map.lineage.unique():
#         print(lineage)
        counts = lineage_chunk_map[lineage_chunk_map.lineage ==lineage].epi_week.value_counts().sort_index()
        #Earliest week must have at least 10 sequences, or 1% of total sequences (for lineages with very small counts)
        counts = counts[(counts > 10) | (counts>counts.sum()*0.1)]
        earliest_week = counts.index[0]
        emergence_dict[lineage] = int(earliest_week)
    return emergence_dict


def get_emergence_weeks():
    emergence_dict = {}
    unique_lineages = list(mutations.lineage.unique())
    lineage_chunks = []
    lineage_chunks_size = math.ceil(len(unique_lineages)/args.threads)
    lineage_chunks = [unique_lineages[i:i + lineage_chunks_size] for i in range(0, len(unique_lineages), lineage_chunks_size)]
    pool = multiprocessing.Pool(args.threads)
    results = pool.map(get_lineage_emergence,lineage_chunks)
    [emergence_dict.update(result) for result in results]
    for lineage in unique_lineages:
        if len(lineage.split(".")) == 1:
            emergence_dict[lineage] = 0
        else:
            #Enforce that a lineage cannot emerge before its parent and if it does, set the emergence to at earliest a week after the parent
            split_lineage = lineage.split(".")
            for j in range(len(split_lineage)-1,0,-1):
                parent_lineage = ".".join(split_lineage[:j])
                if parent_lineage in emergence_dict:
                    if emergence_dict[lineage]<=emergence_dict[parent_lineage] :
                        emergence_dict[lineage] = emergence_dict[parent_lineage]+1
                else:
                    continue
                break
    print(emergence_dict)
    return emergence_dict
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',"--input", metavar='str', type=str,help='Mutation File Path')
    parser.add_argument('-r',"--reference", metavar='str', type=str,help='Reference File Path')
    parser.add_argument('-l',"--lineage", metavar='str', type=str,help='Lineages File Path')
    parser.add_argument('-o',"--output", metavar='str', type=str,help='Output directory path')
    parser.add_argument('-t',"--type", metavar='str', type=str,help='Nucleotide or Amino Acid Mutations')
    parser.add_argument('-ln',"--level_norm", metavar='str', type=str,help='Enable level normalisation',default = "True")
    parser.add_argument('-th',"--threads", metavar='int', type=int,help='Number of threads',default=8)
    parser.add_argument('-mm',"--make_matrix", metavar='str', type=str,help='Make Signature Matrix',default="True")
    parser.add_argument('-ps',"--pseudo_reference", metavar='str', type=str,help='Path to pseudo_reference file directory',default="")
    
    
    args = parser.parse_args()

    if args.pseudo_reference != "":
        reference = SeqIO.read(args.reference, "genbank")
        lineage_aliases = cf.load_lineages(args.lineage)
        pseudo_references = build_reference_dict(f"{args.pseudo_reference}/pseudo_references.fasta")
        pseudo_reference_counts = pd.read_csv(f"{args.pseudo_reference}/pseudo_references_counts.csv")
    
    print(args)
    time = datetime.now()
    if ".csv" in args.input:
        mutations = pd.read_csv(args.input)
    else:
        mutations = pd.read_json(args.input)
    if "epi_week" not in mutations.columns:
        mutation_metadata = mutations["Sequence_Name"].str.split('|',expand=True)
#         The new headers have the hierarchical lineage in the Non-Shortcut-lineage column which is called lineage here
#         I shift the columns around so there is no need for search and replacement for aliases
        mutation_metadata.columns = ["sequence_name","GISAID_ID","Host","Clade","alias","lineage","Length","sample_date","epi_week","Country","FullCountry","Error"]
        mutation_metadata = mutation_metadata[["sequence_name","GISAID_ID","Host","Clade","lineage","alias","Length","sample_date","epi_week","Country","FullCountry","Error"]]
        mutations = pd.concat([mutations,mutation_metadata],axis=1)

    #Get rid of erroneous samples with wrong dates
    mutations = mutations.astype({'epi_week': 'int32'})
    mutations = mutations[mutations.epi_week !=-1]
    mutations = mutations[mutations.sample_date.str.contains("\?\?") == False]
    #Remove Sequences where the sample year does not match (This also removes some sequences with poorly formatted sequence names but it is a low volume)
    mutations =  mutations[mutations.sequence_name.str[-4:] == mutations.sample_date.str[:4]]
    
    mutations = mutations.sort_values(by="epi_week")
    
    #Attempt to define emergence week of lineage to filter out wrongly labelled sequences i.e B.1.1.7 did not emerge in Feb 2020 but is labelled as such
#     emergence_dict = {}
#     unique_lineages = list(mutations.lineage.unique())
#     lineage_chunks = []
#     lineage_chunks_size = math.ceil(len(unique_lineages)/(args.threads*10))
#     lineage_chunks = [unique_lineages[i:i + lineage_chunks_size] for i in range(0, len(unique_lineages), lineage_chunks_size)]
#     linage_chunk_map = [ mutations[mutations.lineage.isin(chunk)] for chunk in lineage_chunks]
#     pool = multiprocessing.Pool(args.threads)
#     results = pool.map(get_lineage_emergence,linage_chunk_map)
#     [emergence_dict.update(result) for result in results]

    emergence_dict = {}
    unique_lineages = list(mutations.lineage.unique())
    for lineage in unique_lineages:
        print(lineage)
        counts = mutations[mutations.lineage ==lineage].epi_week.value_counts().sort_index()
        counts = counts[(counts > 10) | (counts>counts.sum()*0.1)]
        earliest_week = counts.index[0]
        emergence_dict[lineage] = int(earliest_week)
        
    for lineage in unique_lineages:
        if len(lineage.split(".")) == 1:
            emergence_dict[lineage] = 0
        else:
            #Enforce that a lineage cannot emerge before its parent and if it does, set the emergence to at earliest a week after the parent
            split_lineage = lineage.split(".")
            for j in range(len(split_lineage)-1,0,-1):
                parent_lineage = ".".join(split_lineage[:j])
                if parent_lineage in emergence_dict:
                    if emergence_dict[lineage]<=emergence_dict[parent_lineage] :
                        emergence_dict[lineage] = emergence_dict[parent_lineage]+1
                else:
                    continue
                break
    mutations["emergence_week"] = [int(emergence_dict[lineage]) for lineage in mutations.lineage]
    #Keep mutations that are post-emergence
    mutations = mutations[(mutations.epi_week>=mutations.emergence_week)]
    
    
    if args.type == "Amino-Acid" or args.type == "amino-acid":
        mutations = mutations.rename(columns={'AMINO_POS': 'POS'})
    # mutations_with_duplicates = mutations
    mutations = mutations.drop_duplicates(subset=["REF","ALT","POS","lineage"])
    mutations.to_csv(f"{args.output}unique_mutations_pre_normalised.csv",index=False)
    
    print(f'Data loaded :{datetime.now()-time}')
    unique_positions = list(mutations.POS.unique())
    chunks = []
    chunk_size = math.ceil(len(unique_positions)/(args.threads*10))
    chunks = [unique_positions[i:i + chunk_size] for i in range(0, len(unique_positions), chunk_size)]
    if args.pseudo_reference != "":
        inner_mutations_map = [mutations[mutations.POS.isin(chunk)] for chunk in chunks]
        pool = multiprocessing.Pool(args.threads)
        results = pool.map(position_ratio_checker,inner_mutations_map)
        pool.close()
        pool.join()
        mutations = pd.concat(results)
        print(mutations,len(mutations))
    else:
        inner_mutations_map = [mutations[mutations.POS.isin(positions)] for chunk in chunks]
        pool = multiprocessing.Pool(args.threads)
        results = pool.map(position_checker,inner_mutations_map)
        pool.close()
        pool.join()
        mutations = pd.concat(results)
    
    #Output without level normalisation
    if args.level_norm == "True":
        print(mutations.shape)
        print("Level Normalisation")
        A_lineages = [lineage.split(".") for lineage in mutations.lineage.unique() if lineage[0] == "A" ]
        B_lineages = [lineage.split(".") for lineage in mutations.lineage.unique() if lineage[0] == "B" ]
        if len(A_lineages)!=0:
            A_lineages_depth = np.max([len(lineage) for lineage in A_lineages])
            a_levels = [[] for i in range(A_lineages_depth)]
            for lineage in A_lineages:
                depth = len(lineage)
                parent = '.'.join(lineage[:-1])
                # print(lineage,parent)
                lineage_to_level = mutations.loc[mutations.lineage == ".".join(lineage)]
                lineage_to_level.loc[:,"Parent"] = parent
                a_levels[depth-1].append(lineage_to_level)
            for i,lineage_depth in enumerate(a_levels):
                a_levels[i] = pd.concat(a_levels[i])
                a_levels[i] = a_levels[i].drop_duplicates(subset=["REF","ALT","POS","Parent"])
            a_levels = pd.concat(a_levels)
        else:
            a_levels = [pd.DataFrame()]
        if len(B_lineages)!=0:
            B_lineages_depth = np.max([len(lineage) for lineage in B_lineages])
            b_levels = [[] for i in range(B_lineages_depth)]
            for lineage in B_lineages:
                depth = len(lineage)
                parent = '.'.join(lineage[:-1])
                # print(lineage,parent)
                lineage_to_level = mutations.loc[mutations.lineage == ".".join(lineage)]
                lineage_to_level.loc[:,"Parent"] = parent
                b_levels[depth-1].append(lineage_to_level)
            for i,lineage_depth in enumerate(b_levels):
                if len(b_levels[i]) != 0:
                    b_levels[i] = pd.concat(b_levels[i])
                    b_levels[i] = b_levels[i].sort_values("epi_week").drop_duplicates(subset=["REF","ALT","POS","Parent"])
                else:
                    b_levels[i] = pd.DataFrame()
            b_levels = pd.concat(b_levels)
        else:
            b_levels = pd.DataFrame()
        mutations = pd.concat([a_levels,b_levels],axis=0)
    # mutations = mutations.drop("Keep",axis=1)
      
    #This outputs all unique mutations, but each mutation at a position will only be represented once.
    mutations.to_csv(f"{args.output}normalised_unique_mutations.csv",index=False)
    # #This bit outputs all mutations that match a unique mutation i.e does not remove duplicates
    # mutations_filtered[mutations_with_duplicates[(mutations_with_duplicates.POS == row.POS) &
    #                                              (mutations_with_duplicates.REF == row.REF) &
    #                                              (mutations_with_duplicates.ALT == row.ALT) &
    path = args.output.split(".")             
                              
    # mutations_filtered.to_csv(f'{path[0]}_non_unique.{path[1]}')
    if args.make_matrix == "True":
        matrix = cf.convert_to_sig_matrix(mutations,args.type)
        matrix.to_csv(f"{args.output}normalised_unique_mutations_signature_matrix.csv")
        
    