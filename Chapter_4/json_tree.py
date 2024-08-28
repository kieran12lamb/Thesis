import json
import pandas as pd
from Bio import SeqIO

sequences = []
with open("Resources/Sequences/Generated_Sequences/pseudo_references.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        for node in data:
            if node["name"] == record.id:
                print()
                print(node["full_name"])
                name = record.id+"|"+node["full_name"]+"|"+record.id
                print(name)
                record.id = name
                record.description=""
        sequences.append(record)
with open("Resources/Sequences/Generated_Sequences/annotated_pseudo_references.fasta", "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")

# mutation_type= 'amino'
mutation_type= 'nucleotide'
normalised = 'Normalised'
# normalised = 'Reference'

def string_to_list_in_pairs (s):
    return ['|'.join(pair) for pair in zip(s[:-1], s[1:])]

def calculate_distances(dict1,dict2):
    pairwise_distances = {}

    for key1 in dict1.keys():
        for key2 in dict2.keys():
            
            # print(key1+"_"+key2)
            # print((dict1[key1],dict1[key1]))
            pairwise_distances[key1+"|"+key2] = abs(dict1[key1]-dict2[key2])
    return pairwise_distances

def find_node(node_name,node_list):
    for node in node_list:
        if node["name"] == node_name:
            return node

def go_to_root(node,nodes,mutations):
    local_mutation_dict = {}
    inheritence_list = []
    while node['parent'] != "Root" or node["name"] in checked_nodes: #while not at root
        node_height = node["full_name"].count(".")
        node_mutations = mutations[mutations.Sequence_Name == node["name"]] #get mutations for node
        if mutation_type == 'amino':
            node_mutations = node_mutations.dropna()
            mutation_pairs = string_to_list_in_pairs(node_mutations.Amino_Acid) #get pairs of mutations
            mutation_heights = {mut:node_height for mut in node_mutations.Amino_Acid} #get heights of mutations
        else:    
            mutation_pairs = string_to_list_in_pairs(node_mutations.Nucleotide) #get pairs of mutations
            mutation_heights = {mut:node_height for mut in node_mutations.Nucleotide} #get heights of mutations
        mutation_pair_distances = {pair:0 for pair in mutation_pairs} #get distances between pairs
        
        if len(local_mutation_dict.keys())  !=0:
            # print(len(local_mutation_dict.keys()) )
            for previous_node in local_mutation_dict.keys(): #for each previous node
                # print("Calculating pairs of mutations between lineages "+node["name"]+"and "+previous_node )
                distances = calculate_distances(mutation_heights,local_mutation_dict[previous_node]["heights"]) #calculate distances between pairs
                # print(distances)
                mutation_pair_distances = {**mutation_pair_distances, **distances}
                inheritence_list.append(previous_node)
        # print(inheritence_list)
        # print(mutation_pair_distances)
        local_mutation_dict[node["name"]] = {"pairs":mutation_pair_distances,"heights":mutation_heights} #add to dictionary
        checked_nodes.add(node["name"]) #add node to checked nodes
        # print(node["name"],node["parent"],node["alias_parent"])
        node = find_node(node["alias_parent"],nodes) #go to parent
    return local_mutation_dict

with open('Metadata/node_list.json') as f:
    data = json.load(f)

mutations = pd.read_csv("Output/Generated_Sequences/"+normalised+"/mutations.csv")
mutations["Sequence_Name"] = mutations["Sequence_Name"].str.split("|",expand=True)[0]
tree_tips = [node for node in data if len(node['children']) == 0]

checked_nodes = set()

sample_mutation_dict = {}
for node in tree_tips:
    #Ignore recombinants for now
    # print("X" in node["name"],node["full_name"])
    if "X" not in node["name"]:
        if "]" not in node["full_name"]:
            # print(node)
            node_mutations = go_to_root(node,data,mutations)
            sample_mutation_dict = {**sample_mutation_dict, **node_mutations}
            
if mutation_type == 'amino':
    with open("Output/Generated_Sequences/"+normalised+"/sample_distance_mutations_"+mutation_type+".json", "w") as outfile:
        outfile.write(json.dumps(sample_mutation_dict))
else:
    with open("Output/Generated_Sequences/"+normalised+"/sample_distance_mutations"+mutation_type+".json", "w") as outfile:
        outfile.write(json.dumps(sample_mutation_dict))

cooccurrence_distances = {}
for lineage in sample_mutation_dict.keys():
    print(sample_mutation_dict[lineage])
    for mutation in sample_mutation_dict[lineage]["pairs"].keys():
        if mutation in cooccurrence_distances.keys():
            cooccurrence_distances[mutation].append( sample_mutation_dict[lineage]["pairs"][mutation])
        else:
            cooccurrence_distances[mutation] = []
            cooccurrence_distances[mutation].append(sample_mutation_dict[lineage]["pairs"][mutation])
with open("Output/Generated_Sequences/"+normalised+"/cooccurrence_mutations_"+mutation_type+".json", "w") as outfile:
    outfile.write(json.dumps(cooccurrence_distances))
print(cooccurrence_distances)

