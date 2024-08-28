# import libraries
import pandas as pd
import argparse

# create hard-coded domain dfs

# make df of sars-cov-2 orf1ab protein products and their start and end aa positions within orf1ab
# positions taken from Wu et al. 2022: https://doi.org/10.1038/s41401-021-00851-w

proteins = [ "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP9", "NSP10", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16"]
proteinStarts = [1, 181, 819, 2764, 3264, 3570, 3860, 3943, 4141, 4254, 4393, 5325, 5926, 6453, 6799]
proteinEnds = [180, 818, 2763, 3263, 3569, 3859, 3942, 4140, 4253, 4392, 5324, 5925, 6452, 6798, 7096]
proteinsData = {'Protein':proteins, 'Protein Starts':proteinStarts, 'Protein Ends':proteinEnds}
proteinsdf = pd.DataFrame(data=proteinsData)

# spike
spikeDomains = [ "s-no-domain", "NTD", "s-no-domain", "RBD", "s-no-domain", "FP", "s-no-domain", "HR1", "s-no-domain", "HR2", "TMD", "CD"]
spikeDomainStarts = [1, 13, 306, 319, 542, 788, 807, 912, 985, 1163, 1213, 1237]
spikeDomainEnds = [12, 305, 318, 541, 787, 806, 911, 984, 1162, 1212, 1236, 1273]
spikeDomainData = {'Domains':spikeDomains, 'Domain Starts':spikeDomainStarts, 'Domain Ends':spikeDomainEnds}
spikeDomainDF = pd.DataFrame(data=spikeDomainData)

# envelope - https://www.nature.com/articles/s41594-020-00536-8/figures/1
eDomains = ["NTD", "TMD", "CTD"]
eDomainStarts = [1, 8, 39]
eDomainEnds = [7, 38, 75]
eDomainData = {'Domains':eDomains, 'Domain Starts':eDomainStarts, 'Domain Ends':eDomainEnds}
eDomainDF = pd.DataFrame(data=eDomainData)

# nucleocapsid - https://www.nature.com/articles/s41467-021-21953-3#Fig1
nDomains = ["NTD", "RBD", "LINK", "Dimerisation", "CTD"]
nDomainStarts = [1, 50, 174, 247, 365 ]
nDomainEnds = [49, 173, 246, 364, 419]
nDomainData = {'Domains':nDomains, 'Domain Starts':nDomainStarts, 'Domain Ends':nDomainEnds}
nDomainDF = pd.DataFrame(data=nDomainData)

# membrane - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8221690/
mDomains = ["n-terminal ecto", "TMI", "TMI-TMII Link", "TMII","TMII-TMIII Link", "TMIII", "m-no-domain", "cytoplasmic"]
mDomainStarts = [1, 20, 40, 51, 74, 78, 101, 118]
mDomainEnds = [19, 38, 50, 73, 77, 100, 117, 222]
mDomainData = {'Domains':mDomains, 'Domain Starts':mDomainStarts, 'Domain Ends':mDomainEnds}
mDomainDF = pd.DataFrame(data=mDomainData)

# orf3a: https://doi.org/10.3389/fmicb.2022.854567, https://doi.org/10.1128/mSystems.00266-20
orf3aSubDoms = [ 'orf3a-NTD', 'orf3a-TMDI-III', 'orf3a-no-domain', 'orf3a-β1–β8', 'orf3a-no-domain', 'orf3a-CTD']
orf3aSubDomsStart = [ 1, 40, 129, 145, 236, 239]
orf3aSubDomsEnd = [39, 128, 144, 235, 238, 276]
orf3aSubDomData = {'Domains':orf3aSubDoms, 'Domain Starts':orf3aSubDomsStart, 'Domain Ends':orf3aSubDomsEnd}
orf3aSubDomsDF = pd.DataFrame(data= orf3aSubDomData)

# orf7a: https://www.cell.com/iscience/pdf/S2589-0042(21)00155-3.pdf
orf7aSubDoms = ['orf7a-N-terminal', 'orf7a-Ig-like-ectodomain', 'orf7a-hydrophobic-TMD', 'orf7a-ER-retention-motif']
orf7aSubDomsStart = [1, 16, 97, 117]
orf7aSubDomsEnd = [15, 96, 116, 122]
orf7aSubDomData = {'Domains': orf7aSubDoms, 'Domain Starts': orf7aSubDomsStart, 'Domain Ends': orf7aSubDomsEnd}
orf7aSubDomsDF = pd.DataFrame(data=orf7aSubDomData)

# orf8
orf8SubDoms = ['orf8-N-terminal-sig-seq', 'orf8-Ig-like-domain']
orf8SubDomsStart = [1, 16]
orf8SubDomsEnd = [15, 122]
orf8SubDomData = {'Domains':orf8SubDoms, 'Domain Starts':orf8SubDomsStart, 'Domain Ends':orf8SubDomsEnd}
orf8SubDomsDF = pd.DataFrame(data= orf8SubDomData)

# non-structural protein domains (NSPs)
# only nsp1 will start from 1, need to add the aa position of the nsp within orf1ab to each of the others
# eg position 1 in nsp3 will be 819
# positions of NSP domains taken from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7293463/ 
nspSubDoms = [ 'NSP1', 'NSP2', 'NSP3-Ubl1', 'NSP3-HVR','NSP3-MacI', 'NSP3-no-domain', 'NSP3-MacII', 'NSP3-MacIII', 'NSP3-DPUP', 'NSP3-Ubl2', 'NSP3-PL2pro', 'NSP3-no-domain', 'NSP3-NAB', 'NSP3-Beta-SM', 'NSP3-TM', 'NSP3-no-domain', 'NSP3-Subdomain-Y', 'NSP4-no-domain', 'NSP4-TM', 'NSP4-no-domain', 'NSP5-no-domain', 'NSP5-Mpro', 'NSP5-no-domain', 'NSP6-TM', 'NSP7', 'NSP8', 'NSP9', 'NSP10', 'NSP12-no-domain', 'NSP12-NiRAN', 'NSP12-Interface', 'NSP12-RdRP', 'NSP13-ZBD', 'NSP13-no-domain', 'NSP13-Helicase', 'NSP13-no-domain', 'NSP14', 'NSP15-no-domain', 'NSP15-EndoU', 'NSP15-no-domain', 'NSP16' ]
nspSubDomsStart = [1, 181, 819, 930, 1025, 1199, 1232, 1369, 1494, 1564, 1623, 1882, 1907, 2022, 2231, 2373, 2402, 2764, 2777, 3159, 3264, 3267, 3564, 3570, 3860, 3943, 4141, 4254, 4393, 4444, 4643, 4759, 5325, 5420, 5675, 5915, 5926, 6453, 6646, 6797, 6799]
nspSubDomsEnd = [180, 818, 929, 1024, 1198, 1231, 1368, 1493, 1563, 1622, 1881, 1906, 2021, 2230, 2372, 2401, 2763, 2776, 3158, 3263, 3266, 3563, 3569, 3859, 3942, 4140, 4253, 4392, 4443, 4642, 4758, 5324, 5419, 5674, 5914, 5925, 6452, 6645, 6796, 6798, 7096]
nspSubDomnData = {'Domains':nspSubDoms, 'Domain Starts':nspSubDomsStart, 'Domain Ends':nspSubDomsEnd}
nspSubDomnDF = pd.DataFrame(data=nspSubDomnData)

# parse mutations file and process it - fix sequence name, remove Nan amino acid valies, add two columns: reference and alternative aa
def process_mutations_csv(mut_file):
    # Read in csv file containing mutations
    df = pd.read_csv(mut_file)
    
    # Get rid of duplicated sequence name information from the Sequence_Name column
    df[['Sequence_Name', 's1', 's2']] = df['Sequence_Name'].str.split("|", expand=True)
    df = df.drop(['s1', 's2'], axis=1)
    
    # Removing columns where Amino_Acid value is NaN
    df = df[~df['Amino_Acid'].isna()]
    
    # Add alt_aa position column
    alt_aa_list = df["Amino_Acid"].values.tolist()
    alternativePosList = [x[-1:] for x in alt_aa_list]

    # Add ref_aa position column
    refPosList = [x.split(":")[1][0] for x in alt_aa_list]
    df["REF_AA"] = refPosList
    df["ALT_AA"] = alternativePosList
        
    return df

# Set up argparse for command line arguments
parser = argparse.ArgumentParser(description="Process mutations file and generate SARS_smut_list.txt.")
parser.add_argument("mutations_file", metavar="mutations_file", type=str, help="Path to the mutations file in CSV format.")

# Parse command line arguments
args = parser.parse_args()

# call function with command line argument
mutations_file = args.mutations_file
df = process_mutations_csv(mutations_file)

# Function to add protein column to mutations df
def assign_protein(df, proteinsdf):
    proteinList = []

    # check which NSP the amino acid is within the protein positions for and append..
    for index, row in df.iterrows():
        if row['ORF'] == 'orf1ab':
            for _, protein_row in proteinsdf.iterrows():
                if row['AMINO_ORF_POS'] >= protein_row['Protein Starts'] and row['AMINO_ORF_POS'] <= protein_row['Protein Ends']:
                    proteinList.append(protein_row['Protein'])
                    break
        else:
            proteinList.append(row['ORF'])

    df['Protein'] = proteinList

# call function
assign_protein(df, proteinsdf)

# Function to add sliding window column to mutations df

# Need to give the name of the ORF and the dataframe that has its mutations
def SlidingWindow(df, WindowLength):
    # Store the calculated windows in a list
    windows = []
    # loop through each row of the mutation dataframe, to calculate window of each mutation
    for index, row in df.iterrows():
        window = (int)(row['AMINO_ORF_POS'] / WindowLength) + 1
        # add calculated window to list of calculated windows
        windows.append(window)

    # add calculated windows as new column 'Window' in current ORFs mutation dataframe
    df['Window'] = windows

SlidingWindow(df, 20)

# make a dataframe for each ORF containing all its mutations
# append each of the positions within an orf to a list 

# ORF1ab
orf1abdf = (df.loc[(df['ORF'] == 'orf1ab')])
orf1abdfList = list(orf1abdf.AMINO_ORF_POS)

# S
sdf = (df.loc[(df['ORF'] == 'S')])
sdfList = list(sdf.AMINO_ORF_POS)

# N
ndf = (df.loc[(df['ORF'] == 'N')])
ndfList = list(ndf.AMINO_ORF_POS)

# E
edf = (df.loc[(df['ORF'] == 'E')])
edfList = list(edf.AMINO_ORF_POS)

# M
mdf = (df.loc[(df['ORF'] == 'M')])
mdfList = list(mdf.AMINO_ORF_POS)

# ORF3a
orf3adf = (df.loc[(df['ORF'] == 'ORF3a')])
orf3aList = list(orf3adf.AMINO_ORF_POS)

# orf6
orf6df = (df.loc[(df['ORF'] == 'ORF6')])
orf6List = list(orf6df.AMINO_ORF_POS)

# orf7a
orf7adf = (df.loc[(df['ORF'] == 'ORF7a')])
orf7aList = list(orf7adf.AMINO_ORF_POS)

# orf8
orf8df = (df.loc[(df['ORF'] == 'ORF8')])
orf8List = list(orf8df.AMINO_ORF_POS)

# orf10
orf10df = (df.loc[(df['ORF'] == 'ORF10')])
orf10List = list(orf10df.AMINO_ORF_POS)

# function to append domains for each mutation
def sub_dom_function(orfList, domainDF, orfDF):
    domList = []
    for mut in orfList:
        for index, row in domainDF.iterrows():
            if mut >= row['Domain Starts'] and mut <= row['Domain Ends']:
                domList.append(row['Domains'])
                break

    orfDF['Sub-Domains'] = domList


sub_dom_function(sdfList, spikeDomainDF, sdf)
sub_dom_function(edfList, eDomainDF, edf)
sub_dom_function(ndfList, nDomainDF, ndf)
sub_dom_function(mdfList, mDomainDF, mdf)
sub_dom_function(orf3aList, orf3aSubDomsDF, orf3adf)
sub_dom_function(orf8List, orf8SubDomsDF, orf8df)
sub_dom_function(orf7aList, orf7aSubDomsDF, orf7adf)
sub_dom_function(orf1abdfList, nspSubDomnDF, orf1abdf)

# orf6 and orf10
# function to append sub-domain for ORF, where there are no sub-domains
def subdomain_list(orf, n):
    subdomainList = []
    for i in orf:
        subdomainList.append(n)
    return subdomainList

# adding sub domains for orf6 to dataframe
orf6SubDomList = subdomain_list(orf6List, 'orf6')
orf6df['Sub-Domains'] = orf6SubDomList

# adding sub domains for orf10 to dataframe
orf10SubDomList = subdomain_list(orf10List, 'orf10')
orf10df['Sub-Domains'] = orf10SubDomList

# concatenate the dfs

# final dataframe containing each mutations with annotations for sliding window, sub-domain and protein
concatenatedDF = pd.concat([orf1abdf, orf3adf, orf6df, orf7adf, orf8df, orf10df, sdf, edf, ndf, mdf], axis=0)
# convert mutation position to an integer
concatenatedDF['AMINO_ORF_POS'] = df['AMINO_ORF_POS'].astype('int')
# group the dataframe by sequence name 
concatenatedDF['pos'] = concatenatedDF.groupby('Sequence_Name').ngroup()
concatenatedDF = concatenatedDF.sort_values(by='pos').reset_index(drop=True)

# function to convert list of lists to list of strings, each element in the string separated by a pipe
# this will be used to create piped annotations 
def list_of_strings(l):
    str_l = ["|".join([str(j) for j in i]) for i in l]
    return str_l

# add column to the dataframe called annotations
# will have pipe separated annotations in the correct order to be run through WESME tool
annotation_order = ['Nucleotide','ALT_AA', 'AMINO_ORF_POS', 'REF_AA', 'Window', 'Sub-Domains', 'Protein', 'ORF']
concatenatedDF["annotations"] = list_of_strings([list(row[annotation_order]) for i,row in concatenatedDF.iterrows()])

# Create file output in right format to be run through WESME
fileContents = []

# first create header line
samples = list(concatenatedDF.Sequence_Name.unique())
headerString = 'samples' + " " + ','.join(samples)
fileContents.append(headerString)

for mutation in concatenatedDF.annotations.unique():
    filtereddf = concatenatedDF[concatenatedDF.annotations == mutation]
    mut = str(filtereddf.annotations.iloc[0]).replace(" ","-") + " "
    mutList = []
    for lineage in filtereddf.Sequence_Name.unique(): # added .unique()
        lineageindex = samples.index(lineage)
        mutList.append(str(lineageindex))
    mut = mut + (','.join(mutList))
    fileContents.append(mut)


# saving file contents list as a .txt file
with open("SARS_REF_smut_list.txt", mode="w") as file:
    file.write("\n".join(fileContents))
