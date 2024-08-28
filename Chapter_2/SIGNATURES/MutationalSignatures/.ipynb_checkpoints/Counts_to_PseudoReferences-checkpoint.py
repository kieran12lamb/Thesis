import argparse
import pandas as pd
import numpy as np
import multiprocessing
from datetime import datetime
import math
import os.path
import compiled_functions as cf
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
import itertools
import json
    
def load_lineages(path):
    #Load in lineage aliases from JSON
    with open(path) as json_file:
        data = json.load(json_file)
#     lineage_dataframe = pd.concat([pd.DataFrame(data=[str(key),str(data[key])]) for key in data.keys()],axis=1).T
#     lineage_dataframe.columns = ["alias","lineage"]
#     #Remove lineages A and B and any recombinant lineages
#     lineage_dataframe["lineage"] = lineage_dataframe.lineage.astype("str")
#     lineage_dataframe = lineage_dataframe[(lineage_dataframe.lineage != "") & (lineage_dataframe.lineage.str.contains(",") == False)]
#     lineage_dataframe = lineage_dataframe.reset_index(drop=True)
    keys = list(data.keys())
    print(data)
    for key in keys:
        #Remove A and B aliases as they do not alias anything
        if data[key] == "":
            data.pop(key)
        #Remove Recombinant lineages
        elif isinstance(data[key],list):
            data.pop(key)
    return data
    

def reference_chunk(i):
    ref_chunk = list(level)[i:i+args.chunk]
    records = {}
    for lineage in ref_chunk:
        #Ignore recombinant lineages for now
        if "/" in lineages_dictionary[lineage].index.name:
            pseudo_reference = str(reference.seq)
            continue
        pseudo_reference = ""
        for col in lineages_dictionary[lineage].columns:
            #Columns have been set to tuples for some reason
            #Get index(nucleotide) with the greatest count
            max_nucleotide_count_index = lineages_dictionary[lineage][col].idxmax()
            #Get the count value of the nucleotide
            max_nucleotide_count = lineages_dictionary[lineage][col][max_nucleotide_count_index]
            #Get the total number of counts for all nucleotides
            total_nucleotide_count = lineages_dictionary[lineage][col].sum()
            #Check that the max count isn't 0 and that it is >= 75% of nucleotides counted
            if (max_nucleotide_count != 0) and ((max_nucleotide_count/total_nucleotide_count) >=0.75):
                wt = max_nucleotide_count_index
            else :
                #If the lineage is either A or B use the reference nucleotide (These are the only lineages that do not have a .)
                if "." not in lineage:
                    wt = reference.seq[int(col[0])]
                else:
                    wt = ""
                    #Get the hierarchical index name instead of the aliased one
                    lin = lineages_dictionary[lineage].index.name.split(".")
                    # Loop up through the lineage levels i.e B.1.1.7->B.1.1->B.1->B
                    for j in range(len(lin)-1,0,-1):
                        parent_lineage = '.'.join(lin[0:j])
                        if parent_lineage in reference_dict:
                            wt = reference_dict[parent_lineage].seq[int(col[0])]
                            print("Inherted new mutation")
                            break
                    if wt== "":
                        print(f'Shouldnt be happening, here is the problem lineage {".".join(lin)}')
                        triggered =True
                        wt = reference.seq[int(col[0])]
            pseudo_reference+=str(wt)
        records.update({lineage:SeqRecord(Seq(pseudo_reference),id=lineage, name= lineage)})
    return records
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',"--input", metavar='str', type=str,help='Reference Counts Path')
    parser.add_argument('-r',"--reference", metavar='str', type=str,help='Reference File Path')
    parser.add_argument('-l',"--lineage", metavar='str', type=str,help='Lineage Aliases Path')
    parser.add_argument('-c',"--chunk", metavar='int', type=int,help='Chunk Size')
    parser.add_argument('-b',"--block", metavar='int', type=int,help='Block Size')
    parser.add_argument('-th',"--threads", metavar='int', type=int,help='Thread Count')
    parser.add_argument('-o',"--output", metavar='int', type=str,help='Output File Path')  
    args = parser.parse_args()
    time = datetime.now()
    count_dataframe = pd.read_csv(args.input)
    all_lineages = count_dataframe.lineage.unique()
    lineage_aliases = load_lineages(args.lineage)
    reference = SeqIO.read(args.reference, "genbank")
    
    #A dictionary that contains a dataframe of the counts for each lineage
    #The index name of each dataframe contains the long form lineage name which is necessary for working up through the hierarchy
    lineages_dictionary = {}
    reference_dict = {}
    
    for lineage in all_lineages:
        lineage_counts = count_dataframe[count_dataframe.lineage == lineage]
        lineage_counts = lineage_counts.drop("lineage",axis=1)
        lineage_counts = lineage_counts.set_index("nucleotide",drop=True)
        split_lineage = lineage.split(".")
        if split_lineage[0] in lineage_aliases:
            split_lineage = lineage_aliases[split_lineage[0]].split(".")+split_lineage[1:]
            lineage_counts.index.name = ".".join(split_lineage)
            print(lineage_counts.index.name)
        else:
            lineage_counts.index.name = lineage
        lineages_dictionary[lineage] = lineage_counts

    max_lineage_depth = np.max([len(lineages_dictionary[lineage].index.name.split(".")) for lineage in all_lineages])
    phylogeny_levels = [[] for i in range(max_lineage_depth)]
    #Assign each lineage to its appropriate level in the tree, this ensures that lineages above are calculated before lineages below
    for lineage in lineages_dictionary.keys():
        phylogeny_levels[len(lineages_dictionary[lineage].index.name.split("."))-1].append(lineage) 
        
    for idx,level in enumerate(phylogeny_levels):
        i=0
        while(i<len(level)):
            print(f'Phylogeny Level {idx}/{len(phylogeny_levels)}, Time Taken: {datetime.now()-time}')
            pool = multiprocessing.Pool(args.threads)
            if len(level)> i+(args.chunk*args.threads):
                results = pool.map(reference_chunk,np.arange(i,i+(args.chunk*args.threads),args.chunk))
            else:
                length = len(level)
                results = pool.map(reference_chunk,np.arange(i,length,args.chunk))
            pool.close()
            pool.join()
            [reference_dict.update(dictionary) for dictionary in list(results)]
            i+=args.chunk*args.threads
        print("Level Finished")
    reference_array = reference_dict.values()
    pseudo_reference_fasta = MultipleSeqAlignment(reference_array)
    AlignIO.write(pseudo_reference_fasta,args.output, "fasta")
