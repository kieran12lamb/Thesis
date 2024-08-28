import numpy as np
from datetime import datetime
from Bio.Seq import Seq
import pandas as pd
from Bio import AlignIO
from Bio.pairwise2 import identity_match
from Bio.Align import substitution_matrices
from Bio import pairwise2
import sys
import time
from sklearn.utils import resample
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from datetime import datetime
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mtick
from Bio import SeqIO
from Bio import AlignIO
import json
import multiprocessing


def createTables(mutation_bases, ref, classes, orfs, ref_seq, alignment_seq, location=False,pseudo_references=None,lineages=None):
    # Potentially change this by filtering using snpPositions to reduce dimensionality
    time = datetime.now()
    # Create empty column
    empty_column = np.array(np.tile("", len(mutation_bases)), dtype="string_")

    # Create reference sequence matrix
    if pseudo_references is None:
        ref_matrix = np.tile(np.array(ref, dtype="string_"),
                            (len(mutation_bases), 1))
    else:
        ref_matrix = np.array([pseudo_references[str(lineage)] for lineage in lineages])
    # Create matrix of previous nucleotides
    prior_nuc = np.vstack([empty_column, np.array(
        mutation_bases[:, :-1].T, dtype="string_")]).T

    # Create Matrix of next nucleotides
    next_nuc = np.vstack(
        [np.array(mutation_bases[:, 1:].T, dtype="string_"), empty_column]).T

    # Create Trimer matrix from adding prior,reference and next nucleotides
    trimer_matrix = np.char.add(np.char.add(np.array(prior_nuc, dtype="string_"), np.array(
        ref_matrix, dtype="string_")), np.array(next_nuc, dtype="string_"))

    # Create Substitution matrix
    substitution_matrix = np.char.add(np.array(
        ref_matrix, dtype="string_"), np.array(mutation_bases, dtype="string_"))

    # Create Matrix of Dashes
    dash_matrix = np.array(np.full_like(
        substitution_matrix, "-"), dtype="string_")

    # Create Sub-Trimer Matrix by adding the substitution, dash and trimer matrices
    sub_trimer_matrix = np.char.add(np.char.add(np.array(substitution_matrix, dtype="string_"), np.array(
        dash_matrix, dtype="string_")), np.array(trimer_matrix, dtype="string_"))

    # Get locations of Sub-Trimers that are real (i.e no N's, ?'s etc...)
    bool_sub_trimers = np.reshape(np.isin(np.ravel(sub_trimer_matrix), np.array(
        classes[:, 5], dtype="string_")), sub_trimer_matrix.shape)

    # Get indexes of Mutation Sites
    nucleotide_amino_acid_indexes = [np.nonzero(bool_sub_trimers)[1][np.where(
        np.nonzero(bool_sub_trimers)[0] == i)] for i in range(len(bool_sub_trimers))]

    # Get locations of Substitutions that are real (i.e no AA's, ?'s etc...)
    bool_sub_matrix = np.reshape(np.isin(np.ravel(substitution_matrix), np.array(
        classes[:, 2], dtype="string_")), substitution_matrix.shape)
    
    # Get reference peptide sequence
    if pseudo_references is None:
        reference_sequence_peptide = [np.array(list(iterative_translate("".join(orfs.iloc[i].Locations.extract(ref_seq)))), dtype="string_") for i in range(len(orfs))]
        reference_sequence_peptide = np.concatenate(reference_sequence_peptide)
        # Get amino acid changes for each ORF in each Sequence
        amino_acids = np.vstack([getAminoAcidChange(
            alignment_seq[i], reference_sequence_peptide, orfs) for i in range(len(alignment_seq))])
    else:
        pseudo_reference_sequence_peptides = np.array([np.concatenate([np.array(list(iterative_translate("".join(orfs.iloc[i].Locations.extract(Seq("".join(np.array(pseudo_references[lineage],dtype="str"))))))), dtype="string_") for i in range(len(orfs))]) for lineage in lineages],dtype="string_")
        # Get amino acid changes for each ORF in each Sequence
        amino_acids = np.vstack([getAminoAcidChange(
            alignment_seq[i], pseudo_reference_sequence_peptides[i], orfs) for i in range(len(alignment_seq))])
    # print(f'Time Taken for first cf bit: {datetime.now()-time}')
    if location == True:
        time = datetime.now()
        # Include gaps here since we are not making signatures
        gap_classes = np.array(["A-", "C-", "G-", "T-", ])
        possible_mutations = np.array(np.concatenate(
            [np.unique(classes[:, 2]), gap_classes]), dtype="string_")

        # Locate valid mutation locations
        substitution_indexes = np.reshape(np.isin(
            np.ravel(substitution_matrix), possible_mutations), sub_trimer_matrix.shape)
        substitution_indexes = np.sort(np.unique(np.concatenate([np.nonzero(
            substitution_indexes[i])for i in range(len(substitution_indexes))], axis=1)))
        # print(substitution_indexes)
        # Retrieve real substitutions
        sub_trimer_matrix = sub_trimer_matrix[:, substitution_indexes]
        if pseudo_references is None:
            reference_sequence_peptide_pair = np.char.add(
                                                    np.char.add(np.array(reference_sequence_peptide, dtype="string_"), 
                                                                np.array(np.full_like(reference_sequence_peptide, "-"), dtype="string_")),
                                                                np.array(reference_sequence_peptide, dtype="string_"))
            reference_sequence_peptide_pair = np.vstack(
                np.tile(reference_sequence_peptide_pair, (len(amino_acids), 1)))
        else:
            reference_sequence_peptide_pair = np.char.add(np.char.add(np.array(pseudo_reference_sequence_peptides, dtype="string_"), 
                                                                np.array(np.full_like(pseudo_reference_sequence_peptides, "-"), dtype="string_")),
                                                                np.array(pseudo_reference_sequence_peptides, dtype="string_"))
        # possible_amino_mutations =np.sort(np.unique(np.concatenate([np.nonzero(reference_sequence_peptide_pair!= peptide) for peptide in amino_acids],axis=1)))
        possible_amino_mutations = np.sort(np.unique(np.concatenate(
            np.nonzero(reference_sequence_peptide_pair != amino_acids))))
        amino_acids = amino_acids[:, possible_amino_mutations]
        substitution_count_matrix = pd.DataFrame(
            {"REF": [], "ALT": [], 'POS': []})
        amino_acid_count_matrix = pd.DataFrame(
            {"REF": [], "ALT": [], 'POS': []})
        # Loop for the number of Sequences
        for i in range(len(sub_trimer_matrix)):
            # Get nucleotide sequence which
            nucleotide_sequence = sub_trimer_matrix[i]
            mutations = [pd.DataFrame({"REF": [str(np.array(sub, dtype=str))[0]], "ALT":[str(np.array(sub, dtype=str))[1]], 'POS':substitution_indexes[j], "CONTEXT":str(np.array(sub, dtype=str))[3:]})
                         for j, sub in enumerate(nucleotide_sequence)
                         if (len(str(np.array(sub, dtype=str))) == 6) if (str(np.array(sub, dtype=str))[0] != str(np.array(sub, dtype=str))[1]) if (str(np.array(sub, dtype=str))[1] + str(np.array(sub, dtype=str))[3]+str(np.array(sub, dtype=str))[5] != "---")]
            # print(mutations)
            if len(mutations) != 0:
                mutations = pd.concat(mutations)
                mutations["Sequence_Name"] = alignment_seq[i].id
                valid_ref_nucleotides = mutations["REF"].isin(
                    set(["A", "C", "G", "T", "-"]))
                valid_alt_nucleotides = mutations["ALT"].isin(
                    set(["A", "C", "G", "T", "-"]))
                mutations = mutations[np.where((valid_ref_nucleotides == True) & (
                    valid_alt_nucleotides == True), True, False)]
                substitution_count_matrix = pd.concat(
                    [substitution_count_matrix, mutations])
            # Make this faster by doing bool check like for substitutions
            amino_sequence = amino_acids[i]
            amino_mutations = [pd.DataFrame({"REF": [str(np.array(sub, dtype=str))[0]], "ALT":[str(np.array(sub, dtype=str))[2]], 'POS':possible_amino_mutations[j]})
                               for j, sub in enumerate(amino_sequence)
                               if str(np.array(sub, dtype=str))[0] != str(np.array(sub, dtype=str))[2]]
            if len(amino_mutations) != 0:
                amino_mutations = pd.concat(amino_mutations)
                amino_mutations["Sequence_Name"] = alignment_seq[i].id
                invalid_alt_nucleotides = amino_mutations["ALT"].isin(["X"])
                amino_mutations = amino_mutations[invalid_alt_nucleotides == False]
                amino_acid_count_matrix = pd.concat(
                    [amino_acid_count_matrix, amino_mutations])
        # print(f'Time Taken for second cf bit: {datetime.now()-time}')
        return substitution_count_matrix, amino_acid_count_matrix
    else:
        # Retrieve real Amino Acids
        amino_acid_count_matrix = [np.array(np.unique(
            amino_acids[i], return_counts=True)) for i in range(len(amino_acids))]
        # Retrieve real Sub-Trimers
        sub_trimer_count_matrix = [np.array(np.unique(
            sub_trimer_matrix[i][bool_sub_trimers[i]], return_counts=True)) for i in range(len(sub_trimer_matrix))]
        return sub_trimer_count_matrix, amino_acid_count_matrix
    
####################### From Mutation Generator!##################################
def find_mutations(alignment,ref,pseudo_references,lineage_hierarchy_dict,Orfs,filters):
        nuc_table = []
        amino_table = []
        for record in alignment:
            if "|" in record.id and "length" in filters and int(record.id.split("|")[6]) <= filters["length"]:
                continue
            amino_acid = ''.join([iterative_translate(Seq("".join(Orfs.iloc[i].Locations.extract(record.seq)).replace("?", "N"))) for i in range(len(Orfs))])
            if "|" in record.id:
                lineage = record.id.split("|")[4]
            else:
                lineage = ""
            parent_lineage = get_parent_lineage(lineage,pseudo_references,lineage_hierarchy_dict)
            if parent_lineage in pseudo_references:
                parent_reference = pseudo_references[parent_lineage]
            else:
                parent_reference = ref.seq
            reference_amino_acid = ''.join([iterative_translate(Seq("".join(Orfs.iloc[i].Locations.extract(parent_reference)).replace("?", "N"))) for i in range(len(Orfs))])
            for idx,unit in enumerate(record.seq):
                #If there is a nucleotide mutation
                if parent_reference[idx] != unit:
                    REF = str(parent_reference[idx])
                    ALT = str(unit)
                    if REF in set(["A","G","C","T","-"]) and ALT in set(["A","G","C","T","-"]):
                        if idx == 0:
                            CONTEXT = " "+parent_reference[idx]+parent_reference[idx+1]
                        elif idx == len(record.seq)-1:
                            CONTEXT = str(parent_reference[idx-1])+str(parent_reference[idx])+" "
                        else:
                            CONTEXT = str(parent_reference[idx-1])+str(parent_reference[idx])+str(parent_reference[idx+1])
                        POS = idx
                        BIO_POS = idx+1
                        nuc_row = pd.DataFrame({"REF": REF, 
                                                "ALT": ALT, 
                                                "CONTEXT": CONTEXT, 
                                                'POS':POS, 
                                                "BIO_POS": BIO_POS,  
                                                "Sequence_Name":record.id,
                                                "ORF": assignOrf(Orfs,int(POS),"nucleotide"),
                                                "AMINO_POS":getAminoPos(Orfs,int(POS))},
                                                index=[0])
                        nuc_table.append(nuc_row)
                if idx < len(reference_amino_acid):
                    #If there is a nucleotide mutation
                    if reference_amino_acid[idx] != amino_acid[idx]:
                        REF = str(reference_amino_acid[idx])
                        ALT = str(amino_acid[idx])
                        if REF !="X" and ALT !="X":
                            if idx == 0:
                                CONTEXT = " "+reference_amino_acid[idx]+reference_amino_acid[idx+1]
                            elif idx == len(reference_amino_acid)-1:
                                CONTEXT = str(reference_amino_acid[idx-1])+str(reference_amino_acid[idx])+" "
                            else:
                                CONTEXT = str(reference_amino_acid[idx-1])+str(reference_amino_acid[idx])+str(reference_amino_acid[idx+1])
                            AMINO_POS = idx
                            AMINO_BIO_POS = idx+1
                            amino_row = pd.DataFrame({  "REF": REF, 
                                                        "ALT": ALT, 
                                                        "CONTEXT": CONTEXT, 
                                                        'AMINO_POS':AMINO_POS, 
                                                        "AMINO_BIO_POS": AMINO_BIO_POS,  
                                                        "Sequence_Name":record.id,
                                                        "ORF": assignOrf(Orfs,int(AMINO_POS),"amino-acid")},
                                                        index=[0])
                            amino_table.append(amino_row)
        nuc_table = pd.concat(nuc_table,axis=0)
        amino_table = pd.concat(amino_table,axis=0)
        return nuc_table,amino_table
    
def get_parent_lineage(lineage,pseudo_references,lineage_hierarchy_dict):
        full_lineage = ""
        if "." in lineage:
            # print(f"Has a . {lineage}")
            # Get first part of lineage name which may be an alias
            potential_alias = lineage.split(".")[0]
            # Check alias is in list of known aliases
            # print(potential_alias)
            if potential_alias in lineage_hierarchy_dict:
                # print(f"{lineage} has an alias . {potential_alias}")
                # Convert alias to hierarchical name
                full_lineage = lineage_hierarchy_dict[potential_alias]
                # print(lineage,full_lineage)
                # Insert hierarchical name into lineage
                lineage = lineage.replace(potential_alias,full_lineage)
                # Get parent lineage from newly hierarchical lineage name
                parent_lineage = ".".join(lineage.split(".")[:-1])
            else:
                parent_lineage = ".".join(lineage.split(".")[:-1])
            # Check parent is not an alias itself, as aliases are not lineages themselves i.e R == B.1.1.316 and will not be called R`
            split_parent = parent_lineage.split(".")
            i = 1
            known_parent = False
            while known_parent is False:
                # print(lineage,parent_lineage)
                if parent_lineage == "":
                    parent_lineage = "ref"
                    known_parent = True
                    # print(f'Lineage {lineage} only has the reference as a closest ancestor')
                if parent_lineage in pseudo_references:
                    known_parent = True
                    # print(parent_lineage)
                else:
                    parent_lineage = ".".join(split_parent[:-i])
                i+=1
        else:
            parent_lineage = "ref"
        return parent_lineage

def assignOrf(orf_table,position, index_type):
    if index_type == "nucleotide":
        for index in range(len(orf_table)):
            if  position >= orf_table.iloc[index].Start and position<orf_table.iloc[index].End:
                return orf_table.index[index]
        return "Non-Coding"
    elif index_type == "amino-acid":
        for index in range(len(orf_table)):
            if position >= orf_table.iloc[index].AminoStart and position<orf_table.iloc[index].AminoEnd:
                return orf_table.index[index]
        return "Non-Coding"

def getAminoPos(orf_table,position):
    # Return the amino position of a nucleotide index
    for index in range(len(orf_table)):
        # Get the ORF the nucleotide is in
        if position >= orf_table.iloc[index].Start  and position<=orf_table.iloc[index].End:
            position = int((position-orf_table.iloc[index].Offset)/3)
            return position
    # No position due to non-coding segment
    return np.nan
###################################################


# Calculate if a mutation at a given point resulted in an amino acid change
def getAminoAcidChange(alignemnt_seq, reference_sequence_peptide, orfs):
    # Align current sequences with references
    current_sequence_peptides = np.concatenate([np.array(list(iterative_translate(Seq("".join(
        orfs.iloc[i].Locations.extract(alignemnt_seq)).replace("?", "N")))), dtype="string_") for i in range(len(orfs))])
    # print(len(current_sequence_peptides),len(reference_sequence_peptide))
    amino_acid_indexes = [i for i in range(len(
        reference_sequence_peptide)) if reference_sequence_peptide[i] != current_sequence_peptides[i]]
    
    amino_acids = list(np.char.add(np.char.add(reference_sequence_peptide, np.array(
        np.full_like(reference_sequence_peptide, "-"), dtype="string_")), current_sequence_peptides))
    return amino_acids

# def iterative_translate(sequence):
#     codon = ""
#     gaps = ""
#     amino_acid = ""
#     for nuc in sequence:
#         if nuc == "-":
#             gaps += nuc
#         else:
#             codon += nuc

#         if len(codon) == 3:
#             codon = codon.replace("?", "N")
#             amino_acid += str(Seq(codon).translate())
#             codon = ""
#         if len(gaps) == 3:
#             amino_acid += "-"
#             gaps = ""
# #     # Insert final gap if there are 3 or more unused nucleotides or gaps
# #     if len(codon+gaps) >= 3:
# #         amino_acid += "-"
#     return amino_acid

def iterative_translate(sequence):
    amino_acid = ""
    for i in range(0,len(sequence)-2,3):
        codon = str(sequence[i:i+3])
        codon = codon.replace("?", "N")
        if "-" in codon:
            if codon == "---":
                amino_acid +="-"
            else:
                amino_acid+= "X"
        else:
            amino_acid += str(Seq(codon).translate())          
    return amino_acid

def extract_alignment_subset(alignment_file, metadata_file, feature, num_sequences, replacement=True):
    metadata = pd.read_csv(metadata_file)
    metadata = metadata[[feature, "sequence_name"]]
    ids = pd.Series([])
    for value in np.sort(metadata[feature].unique()):
        metadata_subset = metadata[metadata[feature] == value]
        if len(metadata_subset) < num_sequences:
            metadata_subset_ids = metadata_subset.sequence_name
        else:
            metadata_subset_ids = resample(
                metadata_subset, replace=replacement, n_samples=num_sequences).sequence_name
        ids = ids.append(metadata_subset_ids)
    align = AlignIO.read(alignment_file, "fasta")
    sequences = pd.concat([pd.Series([seq.id, str(seq.seq)])
                          for seq in align if seq.id in ids.unique()], axis=1).T
    sequences.columns = ["sequence_name", "sequence"]
    ids = pd.DataFrame({"sequence_name": ids.values})
    merge = pd.merge(left=ids, right=sequences, how="left", on="sequence_name")
    return merge

def csv_to_msa(dataframe=None, csvPath=None, outputDir=None):
    if dataframe is None and csvPath is None:
        return
    elif dataframe is None:
        align = pd.read_csv(csvPath)
    else:
        align = dataframe
    sequences = []
    for i in range(len(align)):
        sequences.append(SeqRecord(
            Seq(align.iloc[i]["sequence"]),
            id=align.iloc[i]["sequence_name"],
            description=""
        ))
    sequences = MultipleSeqAlignment(sequences)
    return sequences

def makeSubstitutionClasses():
    # Array of substitution classes
    substitutionClasses = ['CA', 'CG', 'CT', 'TA',
                           'TG', 'TC', 'AC', 'AG', 'AT', 'GA', 'GT', 'GC']
    classes = pd.DataFrame()
    for substitution in substitutionClasses:
        # Generate all possible trimers for each substitution type
        for p in itertools.product('ATGC', repeat=2):
            trimer = p[0]+substitution[0]+p[1]
            mutant_trimer = p[0]+substitution[1]+p[1]
            row = pd.Series(data=[str(Seq(trimer).translate()),
                                  trimer,
                                  substitution,
                                  str(Seq(mutant_trimer).translate()),
                                  0],
                            index=['Amino_Acid', 'Trimer', 'Substitution', 'Mutant_Amino_Acid', 'Count'])
            classes = classes.append(pd.DataFrame(row).T)
            classes = classes.reset_index(drop=True)
    classes['sub-trimer'] = classes['Substitution']+"-"+classes['Trimer']
    classes['sub-amino'] = classes['Substitution']+"-"+classes['Amino_Acid']
    classes['amino-amino'] = classes['Amino_Acid'] + \
        "-"+classes['Mutant_Amino_Acid']
    return classes

def convert_to_sig_matrix(df,df_type):
    if df_type == "nucleotide" or df_type == "Nucleotide":
        classes = makeSubstitutionClasses()
        sig_matrix = []
        if "Nucleotide_Context" not in df.columns:
            df["Sub-Trimer"] = df["REF"]+df["ALT"]+"-"+df["CONTEXT"]
        else:
            df["Sub-Trimer"] = df["Nucleotide_Context"]
        df = df[df["ALT"] !="-"]
        df = df[["Sequence_Name","Sub-Trimer"]]
        result =  df.groupby(["Sequence_Name","Sub-Trimer"]).size().unstack()
        result = pd.concat([result,pd.DataFrame(columns=classes["sub-trimer"])])
        result = result[classes["sub-trimer"]].fillna(0)
        result.index.name = 'sequence_name'
        return result
    elif df_type == "amino-acid" or df_type == "Amino-Acid" or df_type == "AminoAcid":
        df["Sub"] = df["REF"]+"-"+df["ALT"]
        df = df[["Sequence_Name","Sub"]]
        result =  df.groupby(["Sequence_Name","Sub"]).size().unstack()
        result = result.fillna(0)
        result.index.name = 'sequence_name'
        return result
        
def plot_signature(sig_type,signatures,traditional=False,plot_size=(50,25),comparableY=True):
        if sig_type == "Nucleotide":
            max_bar = np.max(np.ravel(signatures.values))
            if traditional:
                trad_subs = [ ["CA","CG","CT","TA","TC","TG"],
                        ["GT","GC","GA","AT","AG","AC"] ]
                traditional_sig = []
                for i, sub in enumerate(trad_subs[0]):
                    traditional_sig.append(np.add(signatures.iloc[signatures.index.str.startswith(sub)],signatures.iloc[signatures.index.str.startswith(trad_subs[1][i])]))
                signatures = pd.concat(traditional_sig,axis=0)
                max_bar = np.max(np.ravel(signatures.values))
                sub_groups = signatures.index.str[:2]
                palette = ['#04BBEC','black','#E42824','grey','#A0CF63','#EEC4C4']
                fig, axes = plt.subplots(len(signatures.columns),len(np.unique(sub_groups)),figsize=(plot_size[0], plot_size[1]), sharey='row')
            else:
                palette = ['#04BBEC','black','#E42824','grey','#A0CF63','#EEC4C4']
                if len(np.unique(signatures.index.str[:2]))>6:
                    palette = ["mediumorchid","orange", "brown", '#04BBEC','black','#E42824', "teal", "gold", "mediumblue",'grey','#A0CF63','#EEC4C4']
                sub_groups = signatures.index.str[:2]
                fig, axes = plt.subplots(len(signatures.columns),len(np.unique(sub_groups)),figsize=(plot_size[0], plot_size[1]), sharey='row')
            if len(signatures.columns) == 1:
                plot_data = pd.DataFrame(np.array(signatures[signatures.columns[0]]).T)
                plot_data.index = signatures.index
                plot_data.columns = ['values']
                plot_data['sub_groups'] = sub_groups
                plot_data['trimers'] = [trimer[3:] for trimer in plot_data.index]
                for j in range(len(np.unique(sub_groups))):
                    ax = axes[j]
                    ax.set_title(np.unique(sub_groups)[j])
                    sub_plot_data = plot_data[plot_data.sub_groups == np.unique(sub_groups)[j]]
                    sub_plot_data = sub_plot_data.sort_values(by="trimers")
                    sns.barplot(data=sub_plot_data, x='trimers', y="values", hue="sub_groups", palette=[palette[j]],ax=ax,edgecolor='black')
                    ax.get_legend().remove()
                    if comparableY ==False:
                        max_bar = np.max(np.ravel(np.ravel(plot_data["values"])))
                    scaling_factor = max_bar/0.02
                    ax.add_patch(plt.Rectangle((-0.75,max_bar+(0.001*scaling_factor)),25, 0.002*scaling_factor,facecolor=palette[j],clip_on=True,linewidth = 0.1))
                    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
                    ax.set_ylim(0, max_bar+(0.003*scaling_factor))
                    ax.set_xlim(left=-0.75)
                    plt.subplots_adjust(wspace=0,hspace=0.3)
                    ax.set_xticklabels(sub_plot_data.trimers, rotation=90)
                    ax.set_xlabel("")
                    ax.set_ylabel("Mutation Type Probability")
                    ax.set_title(np.unique(sub_groups)[j],fontsize=10)
                    if j !=0:
                        ax.axes.yaxis.set_visible(False)
            else:
                for i in range(len(signatures.columns)):
                    plot_data = pd.DataFrame(np.array(signatures[signatures.columns[i]]).T)
                    plot_data.index = signatures.index
                    plot_data.columns = ['values']
                    plot_data['sub_groups'] = sub_groups
                    plot_data['trimers'] = [trimer[3:] for trimer in plot_data.index]
                    # max_bar = np.max(plot_data['values'])
                    for j in range(len(np.unique(sub_groups))):
                        ax = axes[i][j]
                        ax.set_title(np.unique(sub_groups)[j])
                        sub_plot_data = plot_data[plot_data.sub_groups == np.unique(sub_groups)[j]]
                        sub_plot_data = sub_plot_data.sort_values(by="trimers")
                        sns.barplot(data=sub_plot_data, x='trimers', y="values", hue="sub_groups", palette=[palette[j]],ax=ax,edgecolor='black')
                        ax.get_legend().remove()
                        if comparableY ==False:
                            max_bar = np.max(np.ravel(np.ravel(plot_data["values"])))
                        scaling_factor = max_bar/0.02
                        ax.add_patch(plt.Rectangle((-0.75,max_bar+(0.001*scaling_factor)),25, 0.002*scaling_factor,facecolor=palette[j],clip_on=True,linewidth = 0.1))
                        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
                        ax.set_ylim(0, max_bar+(0.003*scaling_factor))
                        ax.set_xlim(left=-0.75)
                        plt.subplots_adjust(wspace=0,hspace=0.3)
                        ax.set_xticklabels(sub_plot_data.trimers, rotation=90)
                        ax.set_xlabel("")
                        ax.set_ylabel("Mutation Type Probability")
                        ax.set_title(np.unique(sub_groups)[j],fontsize=10)
                        if j !=0:
                            ax.axes.yaxis.set_visible(False)  
                plt.subplots_adjust(wspace=0,hspace=0.3)
        elif sig_type == "Amino Acid":
            cmap = ["#0433FF","#FF2603","#08F900","#000332","#FE35B6","#035301","#FFD300","#009FFF","#9B4D42","#06FBBE","#783FC1","#209698","#FFACFC","#B2CB71","#F1275C","#FE8F42","#DD3DFF","#211A00","#711354","#766C95","#00AD24"]
            max_bar = np.max(np.ravel(signatures.values))
            widths = []
            for i,acid in enumerate(np.sort(signatures.index.str[0].unique())):
                widths.append(len(signatures[signatures.index.str[0] == acid].index.str[2].unique()))
            fig, axes = plt.subplots(len(signatures.columns),len(np.unique(signatures.index.str[0])),figsize=(plot_size[0], plot_size[1]), sharey='row',gridspec_kw={'width_ratios':widths})

            for sig_index in range(len(signatures.columns)):
                plot_data = signatures[signatures.columns[sig_index]]
                summed_vals = plot_data.sum()
                for i,acid in enumerate(np.sort(plot_data.index.str[0].unique())):
                    sub_allAA = plot_data[plot_data.index.str[0] == acid]
                    sub_allAA = sub_allAA.sort_index()
                    axes[sig_index,i].bar(x=sub_allAA.index.str[2], height= sub_allAA.values,color=cmap[i])
                for i, ax in enumerate(axes[sig_index]):
                    ax.set_title(np.sort(plot_data.index.str[0].unique())[i])
                    scaling_factor = max_bar/0.02
                    ax.add_patch(plt.Rectangle((-0.75,max_bar+(0.001*scaling_factor)),
                                                widths[i]*1.005,
                                                0.002*scaling_factor,
                                                facecolor=cmap[i],
                                                linewidth = 0))
                    ax.yaxis.set_major_formatter(mtick.PercentFormatter(summed_vals))
                    ax.set_ylim(0, max_bar+(0.003*scaling_factor))
                    ax.set_xlim(left=-0.75)
                    if i!=0:
                        ax.axes.yaxis.set_visible(False)
                        ax.yaxis.set_major_formatter(mtick.PercentFormatter(summed_vals))
            plt.subplots_adjust(wspace=0,hspace=0.3)
        else:
            print("Sig Type must be Nucleotide or Amino Acid")
def load_lineages(path):
    #Load in lineage aliases from JSON
    with open(path) as json_file:
        data = json.load(json_file)
    keys = list(data.keys())
    for key in keys:
        #Remove A and B aliases as they do not alias anything
        if data[key] == "":
            data.pop(key)
        #Remove Recombinant lineages
        elif isinstance(data[key],list):
            data.pop(key)
    return data

def get_lineage_defining_mutations(pseudo_lineages_path,reference_path,lineage_path):
    reference = SeqIO.read(reference_path, "genbank")
    lineage_references = AlignIO.read(pseudo_lineages_path, "fasta")
    lineage_references_dict = {}
    [lineage_references_dict.update({ref.name:ref.seq}) for i,ref in enumerate(lineage_references)]
    lineage_aliases = load_lineages(lineage_path)
    lineage_aliases["ref"] = reference
    lineage_references_dict["ref"] = reference.seq
    lineage_mutations = []
    for i,lineage_ref in enumerate(lineage_references):
        # Skip Recombinant lineage
        if lineage_ref.id == "XA":
            continue
        split_lineage = str(lineage_ref.id).split(".")
        if len(split_lineage) == 1:
            parent_lineage = "ref"
        else:
            if split_lineage[0] in lineage_aliases:
                split_lineage = lineage_aliases[split_lineage[0]].split(".")+split_lineage[1:]
            for j in range(len(split_lineage)-1,0,-1):
                parent_lineage = '.'.join(split_lineage[0:j])
                if parent_lineage in lineage_references_dict:
                    break
        full_lineage = ".".join(split_lineage)
        [lineage_mutations.append(pd.DataFrame(index=[f'{lineage_references_dict[parent_lineage][i]}{i}{lineage_ref.seq[i]}'],data=1,columns=[full_lineage]))for i in range(len(reference)) if lineage_references_dict[parent_lineage][i] != lineage_ref[i]]
    lineage_mutations = pd.concat(lineage_mutations)
    return lineage_mutations

def calculate_convergent_mutations(pseudo_lineages_path,reference_path,lineage_path):
    reference = SeqIO.read(reference_path, "genbank")
    lineage_references = AlignIO.read(pseudo_lineages_path, "fasta")
    lineage_references_dict = {}
    [lineage_references_dict.update({ref.name:ref.seq}) for i,ref in enumerate(lineage_references)]
    lineage_aliases = load_lineages(lineage_path)
    lineage_aliases["ref"] = reference
    lineage_references_dict["ref"] = reference.seq
    lineage_mutations = []
    lineage_amino_mutations = []
    for i,lineage_ref in enumerate(lineage_references):
        # if i == 100:
        #     break
        # Skip Recombinant lineage
        if lineage_ref.id == "XA":
            continue
        split_lineage = str(lineage_ref.id).split(".")
        if len(split_lineage) == 1:
            parent_lineage = "ref"
        else:
            if split_lineage[0] in lineage_aliases:
                split_lineage = lineage_aliases[split_lineage[0]].split(".")+split_lineage[1:]
            for j in range(len(split_lineage)-1,0,-1):
                parent_lineage = '.'.join(split_lineage[0:j])
                if parent_lineage in lineage_references_dict:
                    break
        full_lineage = ".".join(split_lineage)
        parent_amino_acid = iterative_translate(lineage_references_dict[parent_lineage])
        current_lineage_amino_acid = iterative_translate(lineage_ref.seq)

        # Calculate locations of nucleotide_mutations
        [lineage_mutations.append(pd.DataFrame(index=[f'{lineage_references_dict[parent_lineage][i]}{i}{lineage_ref.seq[i]}'],data=1,columns=[full_lineage]))for i in range(len(reference)) if lineage_references_dict[parent_lineage][i] != lineage_ref[i]]

        [lineage_amino_mutations.append(pd.DataFrame(index=[f'{parent_amino_acid[i]}{i}{current_lineage_amino_acid[i]}'],data=1,columns=[full_lineage]))for i in range(len(parent_amino_acid)) if parent_amino_acid[i] != current_lineage_amino_acid[i]]
    lineage_mutations_nucleotide = pd.concat(lineage_mutations)
    lineage_amino_mutations = pd.concat(lineage_amino_mutations)


    mutation_types = ["nucleotide","amino"]
    
    for mutation_type in mutation_types:
        if mutation_type == "nucleotide":
            lineage_mutations = lineage_mutations_nucleotide
        else:
            lineage_mutations = lineage_amino_mutations
        #Remove positions where mutation has only occurred once i.e not convergent
        rotated_lineage_mutations = lineage_mutations.T
        non_unique_mutations = []
        for col in rotated_lineage_mutations.columns:
            if len(rotated_lineage_mutations[col].shape)>1:
                non_unique_mutations.append(col)
        lineage_mutations = lineage_mutations[lineage_mutations.index.isin(non_unique_mutations)]
        
        max_lineage_depth = np.max([len(lineage.split(".")) for lineage in lineage_mutations.columns])
        phylogeny_levels = [[] for i in range(max_lineage_depth)]
        #Assign each lineage to its appropriate level in the tree, this ensures that lineages above are calculated before lineages below
        for lineage in lineage_mutations.columns:
            phylogeny_levels[len(lineage.split("."))-1].append(lineage) 
        
        convergent_mutations = []
        to_be_validated_mutations = []
        for idx,phylogeny_level in enumerate(phylogeny_levels):
            possible_level_convergent = lineage_mutations[phylogeny_level].groupby(by=lineage_mutations[phylogeny_level].index).sum().sort_index()
            possible_level_convergent = possible_level_convergent.dropna(how='all')
            hidden_ancestor_mutations = possible_level_convergent[possible_level_convergent.sum(axis=1) >1]
            
            #These mutations may be convergent, but are observed multiple times at the same level so could be unobserved ancestor
            A_hidden_ancestor_mutations =hidden_ancestor_mutations[hidden_ancestor_mutations.columns[hidden_ancestor_mutations.columns.str.startswith("A")]]
            A_hidden_ancestor_mutations = A_hidden_ancestor_mutations[A_hidden_ancestor_mutations.sum(axis=1)>1].index

            #These mutations may be convergent, but are observed multiple times at the same level so could be unobserved ancestor
            B_hidden_ancestor_mutations =hidden_ancestor_mutations[hidden_ancestor_mutations.columns[hidden_ancestor_mutations.columns.str.startswith("B")]]
            B_hidden_ancestor_mutations = B_hidden_ancestor_mutations[B_hidden_ancestor_mutations.sum(axis=1)>1].index

            #These mutations are convergent as they occur at the same level but in a different branch i.e A or B of the tree
            occurs_in_A_and_B = A_hidden_ancestor_mutations.intersection(B_hidden_ancestor_mutations)
            convergent_mutations+=list(occurs_in_A_and_B)

            #These mutations need to be checked to see if they occur in other levels of the tree
            not_in_A_and_B = list(A_hidden_ancestor_mutations[A_hidden_ancestor_mutations.isin(occurs_in_A_and_B) == False]) + list(B_hidden_ancestor_mutations[B_hidden_ancestor_mutations.isin(occurs_in_A_and_B) == False] )
            to_be_validated_mutations+=not_in_A_and_B

            #These mutations are duplicates but not in the same tree level, therefore they are convergent
            real_convergent = list(possible_level_convergent[possible_level_convergent.sum(axis=1) ==1].index)
            convergent_mutations+=list(real_convergent)
        
        #Check mutations to ensure they occur at multiple levels of the tree, i.e that they are convergent across levels
        validated_convergent_mutations = list(np.unique([mutation for idx,mutation in enumerate(to_be_validated_mutations) if to_be_validated_mutations.count(mutation)>1]))
        convergent_mutations+=validated_convergent_mutations
        if mutation_type == "nucleotide":
            convergent_mutations_nucleotide = convergent_mutations
        else:
            convergent_mutations_amino = convergent_mutations

    return convergent_mutations_nucleotide,convergent_mutations_amino
    

def normalise_single_lineage(parent_record,mutations,pseudo_reference_counts):
    #Mutations should be sorted by date, with non-unique mutations removed
    lineage = mutations.index.name
    # Loop through each mutation at the position
    for pos in mutations.POS.unique():
        pos_mutations = mutations[mutations.POS == pos]
        reference_nucleotide = parent_record.seq[pos]
        for index in pos_mutations.index:
            nucleotides = pseudo_reference_counts[pseudo_reference_counts.lineage ==lineage][str(pos)]
            nucleotides = nucleotides/nucleotides.sum()
            nucleotides.index = pseudo_reference_counts[pseudo_reference_counts.lineage ==lineage].nucleotide

            parent_nucleotides = pseudo_reference_counts[pseudo_reference_counts.lineage ==parent_record.id][str(pos)]
            parent_nucleotides = parent_nucleotides/parent_nucleotides.sum()
            parent_nucleotides.index = pseudo_reference_counts[pseudo_reference_counts.lineage ==parent_record.id].nucleotide
            alt_nuc = pos_mutations.iloc[index].ALT
            #Check that lineages have sequences in the dataset i.e B.1.617 has no sequences 
            if nucleotides.size != 0 and parent_nucleotides.size != 0:
                if nucleotides.loc[alt_nuc]>=0.5 and parent_nucleotides.loc[alt_nuc] >=0.5:
                    mutations.loc[index,"Remove"] = True
            else:
                print("No Sequences")
    return mutations[mutations.Remove ==False].reset_index(drop=True)


        
def makeOrfTable(genbank_record):
    orfs=[]
    for feature in genbank_record.features:
        if feature.type =="CDS":
            orf = feature.qualifiers['gene'][0]
            for i, locations in enumerate(feature.location.parts):
                orfs.append([orf, locations.start, locations.end, i, locations])
    orfs = pd.DataFrame(orfs)
    orfs.columns = ['ORF','Start','End','Part','Locations']
    orfs = orfs.set_index("ORF")
    return orfs


def translate_with_genbank(sequence,ref):
    orfs = makeOrfTable(ref)
    translated_sequence = [np.array(list(iterative_translate("".join(orfs.iloc[i].Locations.extract(sequence))))) for i in range(len(orfs))]
    return translated_sequence


def translate_with_orfs(sequence,ref):
    orfs = makeOrfTable(ref)
    orf_order = [orf_name for orf_name in orfs.index]
    mutations = []
    for orf in range(len(orf_order)):
        for position in range(len(ref.seq[orf])):
            if ref.seq[orf][position] != sequence[orf][position]:
                mutation = f"{orf_order[orf]}:{str(ref.seq[orf][position])}{position+1}{str(sequence[orf][position])}"
                mutations.append(mutation)
    return mutations



def translate_and_get_mutations(sequence_record,ref):
    ref_protein_sequence= translate_with_genbank(ref.seq,ref)
    #Translate Current Sequence
    sequence_translation = translate_with_genbank(sequence_record.seq,ref)
    #Retrieve Mutations
    mutations = translate_with_orfs(sequence_translation,ref_protein_sequence,ref)
    return sequence_translation,mutations
    
def convert_lineage_name(lineage,alias_table):
    if lineage.split(".")[0] in alias_table:
        return alias_table[lineage.split(".")[0]]+"."+".".join(lineage.split(".")[1:])
    else:
        if lineage.count(".")>3:
            for key in alias_table.keys():
                if alias_table[key] == ".".join(lineage.split(".")[:4]):
                    return key+"."+".".join(lineage.split(".")[4:])
    return lineage
    
def get_alias(lineage,aliases):
    if lineage.split(".")[0] in aliases:
        return ".".join([aliases[lineage.split(".")[0]]] + lineage.split(".")[1:])
    elif len(lineage.split(".")) >4 and ".".join(lineage.split(".")[0:4]) in aliases:
        return ".".join([aliases[lineage.split(".")[0:4]]] + lineage.split(".")[4:])
    return lineage

def get_parents(lineage_alias):
    return [".".join(lineage_alias.split(".")[0:i+1]) for i in range(len(lineage_alias.split(".")))]

def get_closest_parent(lineage,sequences):
    parents = reversed(get_parents(lineage)[:-1])
    print(lineage)
    for parent in parents:
        print(parent)
        if parent in sequences:
            print("returned")
            return parent
    return "Ref"


def mutations_to_vcf(mutations,pseudo_ref_dict,aliases,filename):
    vcf = []
    lineage_sets = np.array_split(mutations.lineage.unique(), 8)
    star_map_input = [[mutations[mutations.lineage.isin(lineage_set)],pseudo_ref_dict,aliases] for idx,lineage_set in enumerate(lineage_sets)]
    pool = multiprocessing.Pool(8)
    results = pool.starmap(lineage_vcf_maker,star_map_input)
    pool.close()
    pool.join()
    vcf = pd.concat(list(results),axis=0)
    vcf = vcf.fillna('.')
    vcf = vcf.reindex()
    header_lines = "\n".join([
        '##fileformat=VCFv4.3',
        '##contig=<ID=0>',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">',
        '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of sequences">',
        '##FORMAT=<ID=ALT,Number=.,Type=String,Description="Genotype">',
        '##FORMAT=<ID=CONTEXT,Number=.,Type=String,Description="Nucleotide Context">',
        '\n'
        ])
    with open(filename, 'w') as vcf_file:
        vcf_file.write(header_lines)
    vcf.to_csv(filename,sep="\t",index=False,mode="a")
    
    

def lineage_vcf_maker(mutations,pseudo_ref_dict,aliases):
    lineages = mutations.lineage.unique()
    FORMAT = "ALT:CONTEXT"
    dict_list = []
    for lineage in lineages:
        print(f"Lineage {lineage}")
        lineage_muts = mutations[mutations.lineage == lineage]
        positions = lineage_muts.BIO_POS.unique()
        alias = get_alias(lineage,aliases=aliases)
        lineage_parent = get_closest_parent(alias,pseudo_ref_dict)
        parent_alias = get_alias(lineage_parent,aliases=aliases)
        for pos in positions:
            position_lineage_muts = lineage_muts[lineage_muts.BIO_POS == pos]
            data = {'#CHROM' : lineage, 
                    'POS' : pos, 
                    'ID': "SNP", 
                    'REF' : pseudo_ref_dict[parent_alias][pos-1],
                    'ALT':",".join(list(position_lineage_muts.ALT.unique())),
                    'QUAL':'.',
                    'FILTER':"Pass",
                    'INFO':";".join([f'NS={position_lineage_muts.size}',
                                     f'AF={",".join([str(allel/position_lineage_muts.ALT.value_counts().sum()) for allel in position_lineage_muts.ALT.value_counts().sort_index()])}'
                                    ]),
                    'FORMAT':FORMAT
                   }
            samples = ({name:f'{position_lineage_muts[position_lineage_muts.Sequence_Name == name].ALT.values[0]}:{position_lineage_muts[position_lineage_muts.Sequence_Name == name].CONTEXT.values[0]}'
                      for name in position_lineage_muts.Sequence_Name.unique()})
            data.update(samples)
            data = pd.DataFrame(data, index=[0])  # the `index` argument is important
            dict_list.append(data)
    lineage_list = pd.concat(dict_list,axis=0)
    return lineage_list

            
            
        
        
        
            
        
        
    
    
    