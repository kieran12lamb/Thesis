import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sklearn
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
import allel
import seaborn as sns
import itertools
from datetime import datetime
import matplotlib.ticker as mtick
import MutationalSignatures.compiled_functions as cf
import multiprocessing
import math
import traceback
import os
#The Genome Signal Processor class is used to build mutational signature tables from 
#FASTA multi-sequence alignments
class SignatureMatrix:

    #Constructor for the Genome Signal Processor 
    def __init__(self, multiAlignmentPath,refPath,trinucleotide_file=None, amino_acid_file=None ,seq_count = None):
        self.amino_acid_substitution_matrix = None
        self.amino_acid_substitution_matrix = None
        if trinucleotide_file or amino_acid_file is not None:
            if trinucleotide_file is not None:
                self.trinucleotide_substitution_matrix = pd.read_csv(trinucleotide_file,index_col=False)
                self.trinucleotide_substitution_matrix = self.trinucleotide_substitution_matrix.set_index("Sequence_Name") 
            if amino_acid_file is not None:
                self.amino_acid_substitution_matrix = pd.read_csv(amino_acid_file,index_col=False)
                self.amino_acid_substitution_matrix = self.amino_acid_substitution_matrix.set_index("Sequence_Name") 
            self.ref = SeqIO.read(refPath, "genbank")
            self.orfs = self.makeOrfTable(self.ref)
            #Create the substitution classes
            self.classes = self.makeSubstitutionClasses()
            self.all_amino_acids = self.all_amino_acid_permutations()
        else:
            #Load alignment file
            self.multiAlignmentPath = multiAlignmentPath
            if seq_count is None:
                self.seq_count = sum(1 for line in open(multiAlignmentPath))/2
            else:
                self.seq_count = seq_count
            self.ref = SeqIO.read(refPath, "genbank")
            self.orfs = self.makeOrfTable(self.ref)
            #Create the substitution classes
            self.classes = self.makeSubstitutionClasses()
            self.all_amino_acids = self.all_amino_acid_permutations()

    def update(self,multiAlignmentPath):
        self.align = AlignIO.read(multiAlignmentPath, "fasta")
        self.chunk_size = 20
        self.threads=8
        align_ids = [seq.id for seq in self.align]
        align_ids_in_matrix = list(pd.Series(align_ids)[pd.Series(align_ids).isin(self.trinucleotide_substitution_matrix.index) ==False ].index)
        print(f'{len(align_ids_in_matrix)} new sequences')
        i=0
        while(i<len(align_ids_in_matrix)):
            pool = multiprocessing.Pool(self.threads)
            indexes = [align_ids_in_matrix[i:i+(self.chunk_size*self.threads)][j:j+self.chunk_size] for j in range(0, len(align_ids_in_matrix[i:]), self.chunk_size)]
            indexes = filter(None, indexes)
            results = pool.map_async(self.update_table_chunk,indexes)
            pool.close()
            pool.join()

            results = results.get()
            results = list(zip(*results))
            pre_trinucleotide_matrix = pd.concat(results[0])
            pre_trinucleotide_matrix = pre_trinucleotide_matrix.set_index("Sequence_Name",drop=True)

            pre_amino_acid_substitution_matrix = pd.concat(results[1])
            pre_amino_acid_substitution_matrix = pre_amino_acid_substitution_matrix.set_index("Sequence_Name",drop=True)

            self.trinucleotide_substitution_matrix = pd.concat([self.trinucleotide_substitution_matrix,pre_trinucleotide_matrix])
            self.amino_acid_substitution_matrix = pd.concat([self.amino_acid_substitution_matrix,pre_amino_acid_substitution_matrix])
            i=i+(self.chunk_size*self.threads)
        change = self.amino_acid_substitution_matrix.columns
        change = [val[0] != val[2] for val in change]
        change = self.amino_acid_substitution_matrix.columns[change]
        self.amino_acid_substitution_matrix = self.amino_acid_substitution_matrix[change]
        self.amino_acid_substitution_matrix = self.amino_acid_substitution_matrix.fillna(0)

    def build(self,chunk_size=20,threads=8,start_from=0,use_file = None):
        # divide alignment to free up memory
        time = datetime.now()
        unit_align = 10000
        j = start_from
        #Generate mutational signature tables
        self.trinucleotide_substitution_matrix = pd.DataFrame()
        self.amino_acid_substitution_matrix = pd.DataFrame()
        self.chunk_size = chunk_size
        total_blocks = math.ceil(self.seq_count/unit_align)
        while(j<self.seq_count):
            print(f'Block {j/unit_align}/{total_blocks}, Time Taken: {datetime.now()-time}')
            if self.seq_count > j+unit_align:
                self.align = AlignIO.read(self.multiAlignmentPath, "fasta")[j:j+unit_align]
            else:
                self.align = AlignIO.read(self.multiAlignmentPath, "fasta")[j:self.seq_count]
            i=0
            total_chunks = math.ceil(len(self.align)/(self.chunk_size*threads))
            while(i<len(self.align)):
                chunk_time = datetime.now()
                pool = multiprocessing.Pool(threads)
                if len(self.align)> i+(chunk_size*threads):
                    results = pool.map(self.table_chunk,np.arange(i,i+(chunk_size*threads),chunk_size))
                else:
                    results = pool.map(self.table_chunk,np.arange(i,len(self.align),chunk_size))
                pool.close()
                pool.join()
                # results = results.get()
                results = list(zip(*results))
                all_aa = pd.DataFrame(columns=self.all_amino_acids)
                all_aa["Sequence_Name"] = []
                trinucleotide_substitution_matrix = pd.concat(results[0])

                amino_acid_substitution_matrix = pd.concat([all_aa,pd.concat(results[1])])
                change = pd.Series(amino_acid_substitution_matrix.columns).isin(all_aa)
                # change = [val[0] != val[2] for val in change]
                change = amino_acid_substitution_matrix.columns[change]
                amino_acid_substitution_matrix = amino_acid_substitution_matrix[change]
                amino_acid_substitution_matrix = amino_acid_substitution_matrix.fillna(0)
                if use_file is not None:
                    self.file_upload(trinucleotide_substitution_matrix,amino_acid_substitution_matrix,use_file)
                else:
                    self.trinucleotide_substitution_matrix = pd.concat([self.trinucleotide_substitution_matrix,trinucleotide_substitution_matrix])
                    self.amino_acid_substitution_matrix = pd.concat([self.amino_acid_substitution_matrix,amino_acid_substitution_matrix])
                print(f'Chunk {i/(chunk_size*threads)}/{total_chunks}, Time Taken: {datetime.now()-chunk_time}')
                i=i+(chunk_size*threads)
            j+=unit_align
        
    def makeOrfTable(self, genbank_record):
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
        
    def makeSubstitutionClasses(self):
        #Array of substitution classes
        substitutionClasses = ['CA', 'CG', 'CT', 'TA','TG', 'TC', 'AC', 'AG', 'AT', 'GA', 'GT', 'GC']
        classes = pd.DataFrame()
        for substitution in substitutionClasses:
            #Generate all possible trimers for each substitution type
            for p in itertools.product('ATGC', repeat=2):
                trimer = p[0]+substitution[0]+p[1]
                mutant_trimer = p[0]+substitution[1]+p[1]
                row = pd.Series(data=[str(Seq(trimer).translate()),
                                    trimer,
                                    substitution,
                                    str(Seq(mutant_trimer).translate()),
                                    0],
                                index=['Amino_Acid', 'Trimer', 'Substitution','Mutant_Amino_Acid' ,'Count'])
                classes = classes.append(pd.DataFrame(row).T)
                classes = classes.reset_index(drop=True)
        classes['sub-trimer'] = classes['Substitution']+"-"+classes['Trimer']
        classes['sub-amino'] = classes['Substitution']+"-"+classes['Amino_Acid']
        classes['amino-amino'] = classes['Amino_Acid']+"-"+classes['Mutant_Amino_Acid']
        
        return classes
    def all_amino_acid_permutations(self):
        amino_acids = ['A', 'R', 'N', 'D','B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I','L','K','M','F','P','S','T','W','Y','V','-','*']
        all_amino_acids = []
        for pair in itertools.product(amino_acids, repeat=2):
            all_amino_acids.append(pair[0]+"-"+pair[1])
        change = [val[0] != val[2] for val in all_amino_acids]
        change = np.array(all_amino_acids)[change]
        return change
        
    #Join metadata files with the correct sequences
    def addMetadata(self,metadata_files=None):
        if metadata_files is not None:
            for metadata in metadata_files:
                if ".tsv" in metadata: 
                    metadata = pd.read_csv(metadata,sep="\t")
                else:
                    metadata = pd.read_csv(metadata)  # [['Sequence_Name', 'lineage']]
                if "strain" in metadata.columns:
                    metadata = metadata.rename(columns={"strain": "Sequence_Name"})
                if "pangolin_lineage" in metadata.columns:
                    metadata = metadata.rename(columns={"pangolin_lineage": "lineage"})
                if "date" in metadata.columns:
                    metadata = metadata.rename(columns={"date": "sample_date"})

                if "Sequence_Name" == self.trinucleotide_substitution_matrix.index.name:
                    self.trinucleotide_substitution_matrix = self.trinucleotide_substitution_matrix.merge(metadata, left_index=True,right_on='Sequence_Name', how='left')
                    self.amino_acid_substitution_matrix = self.amino_acid_substitution_matrix.merge(metadata, left_index=True,right_on='Sequence_Name', how='left')
                elif "Sequence_Name" not in self.trinucleotide_substitution_matrix.columns:
                    
                    no_Sequence_Name_metadata = metadata.drop("Sequence_Name",axis=1).columns

                    self.trinucleotide_substitution_matrix = self.trinucleotide_substitution_matrix.drop(no_Sequence_Name_metadata,axis=1)
                    self.trinucleotide_substitution_matrix = self.trinucleotide_substitution_matrix.merge(metadata, left_index=True,right_on='Sequence_Name', how='left')
                    self.trinucleotide_substitution_matrix = self.trinucleotide_substitution_matrix.set_index("Sequence_Name",drop=True)
                    
                    self.amino_acid_substitution_matrix = self.amino_acid_substitution_matrix.drop(no_Sequence_Name_metadata,axis=1)
                    self.amino_acid_substitution_matrix = self.amino_acid_substitution_matrix.merge(metadata, left_index=True,right_on='Sequence_Name', how='left')
                    self.amino_acid_substitution_matrix = self.amino_acid_substitution_matrix.set_index("Sequence_Name",drop=True)
                else:
                    self.trinucleotide_substitution_matrix = self.trinucleotide_substitution_matrix.merge(metadata, on='Sequence_Name', how='left')
                    self.amino_acid_substitution_matrix = self.amino_acid_substitution_matrix.merge(metadata, on='Sequence_Name', how='left')
                    self.amino_acid_substitution_matrix  = self.amino_acid_substitution_matrix.set_index("Sequence_Name",drop=True)
                    self.trinucleotide_substitution_matrix = self.trinucleotide_substitution_matrix.set_index("Sequence_Name",drop=True)
        else:
            if self.trinucleotide_substitution_matrix is not None:
                self.trinucleotide_substitution_matrix = self.trinucleotide_substitution_matrix.reset_index().rename(columns={'index': 'Sequence_Name'})
                metadata  = self.trinucleotide_substitution_matrix.Sequence_Name.str.split("|",expand=True,)
                metadata.columns = ["Sequence_Name","GISAID_ID","Host","Clade","lineage","Non-Shortcut-Lineage","Length","sample_date","epi_week","Country","FullCountry","Error"]
                metadata.drop("Sequence_Name",axis=1)
                self.trinucleotide_substitution_matrix = pd.concat([self.trinucleotide_substitution_matrix, metadata],axis=1)
                # self.trinucleotide_substitution_matrix.index = self.trinucleotide_substitution_matrix.Sequence_Name
                self.trinucleotide_substitution_matrix = self.trinucleotide_substitution_matrix.set_index("Sequence_Name",drop=True)
            if self.amino_acid_substitution_matrix is not None:
                self.amino_acid_substitution_matrix = self.amino_acid_substitution_matrix.reset_index().rename(columns={'index': 'Sequence_Name'})
                metadata  = self.amino_acid_substitution_matrix.Sequence_Name.str.split("|",expand=True,)
                metadata.columns = ["Sequence_Name","GISAID_ID","Host","Clade","lineage","Non-Shortcut-Lineage","Length","sample_date","epi_week","Country","FullCountry","Error"]
                self.amino_acid_substitution_matrix = pd.concat([self.amino_acid_substitution_matrix,metadata],axis=1)
                # self.amino_acid_substitution_matrix.index = self.amino_acid_substitution_matrix.Sequence_Name
                self.amino_acid_substitution_matrix = self.amino_acid_substitution_matrix.set_index("Sequence_Name",drop=True)

    def table_chunk(self,i):
        chunk_size = self.chunk_size
        if i<len(self.align) and i+chunk_size>len(self.align):
            alignment = np.char.array([list(rec) for rec in self.align[i:]])
            # alignment=np.array(alignnp[i:],dtype="string_")
            alignment_seq = self.align[i:]
        else:
            alignment = np.char.array([list(rec) for rec in self.align[i:i+chunk_size]])
            # alignment=np.array(alignnp[i:i+chunk_size],dtype="string_")
            alignment_seq = self.align[i:i+chunk_size]
        if alignment.ndim ==1:
            alignment = np.array([alignment],dtype="string_")
            alignment_seq = self.align
            
        alignment_ids = [alignment_seq[i].id for i in range(0,len(alignment_seq))]

        #Mutatation Trimer Count Table Row
        mctTrimerRow = pd.DataFrame(index=self.classes['sub-trimer'].unique()).T

        #Mutatation Amino Count Table Row
        mctAminoRow = pd.DataFrame(index=self.classes['amino-amino'].unique()).T
        
        trinucleotide_count_matrix, amino_acid_count_matrix = cf.createTables(alignment,
                                                                            np.char.array(self.ref.seq.split(maxsplit=True)),
                                                                            np.array(self.classes),
                                                                            self.orfs,
                                                                            self.ref.seq,
                                                                            alignment_seq)
        amino_acid_count_matrix = [pd.Series(index = amino_acid_count_matrix[j][0].astype(str), data = amino_acid_count_matrix[j][1].astype(int)) for j in range(len(amino_acid_count_matrix))]
        mctAminoRows = pd.concat([mctAminoRow,pd.concat(amino_acid_count_matrix,axis=1).T]).fillna(0)
        mctAminoRows['Sequence_Name'] = alignment_ids

        trinucleotide_count_matrix = [pd.Series(index = trinucleotide_count_matrix[j][0].astype(str), data = trinucleotide_count_matrix[j][1].astype(int)) for j in range(len(trinucleotide_count_matrix))]
        mctTrimerRows = pd.concat([mctTrimerRow,pd.concat(trinucleotide_count_matrix,axis=1).T]).fillna(0)
        mctTrimerRows['Sequence_Name'] = alignment_ids

        return mctTrimerRows, mctAminoRows

    def update_table_chunk(self,indexes):
        chunk_size = len(indexes)
        alignment = np.char.array([list(rec) for rec in [self.align[i] for i in indexes]])
        alignment_seq = [self.align[i] for i in indexes]
        alignment_ids = [alignment_seq[i].id for i in range(0,len(alignment_seq))]
        #Mutatation Trimer Count Table Row
        mctTrimerRow = pd.DataFrame(index=self.classes['sub-trimer'].unique()).T
        #Mutatation Amino Count Table Row
        mctAminoRow = pd.DataFrame(index=self.classes['amino-amino'].unique()).T
        trinucleotide_count_matrix, amino_acid_count_matrix = cf.createTables(alignment,
                                                                            np.char.array(self.ref.seq.split(maxsplit=True)),
                                                                            np.array(self.classes),
                                                                            self.orfs,
                                                                            self.ref.seq,
                                                                            alignment_seq)
        amino_acid_count_matrix = [pd.Series(index = amino_acid_count_matrix[j][0].astype(str), data = amino_acid_count_matrix[j][1].astype(int)) for j in range(len(amino_acid_count_matrix))]
        mctAminoRows = pd.concat([mctAminoRow,pd.concat(amino_acid_count_matrix,axis=1).T]).fillna(0)
        mctAminoRows['Sequence_Name'] = alignment_ids

        trinucleotide_count_matrix = [pd.Series(index = trinucleotide_count_matrix[j][0].astype(str), data = trinucleotide_count_matrix[j][1].astype(int)) for j in range(len(trinucleotide_count_matrix))]
        mctTrimerRows = pd.concat([mctTrimerRow,pd.concat(trinucleotide_count_matrix,axis=1).T]).fillna(0)
        mctTrimerRows['Sequence_Name'] = alignment_ids

        return mctTrimerRows, mctAminoRows

    def export(self,trinucleotide_file, amino_acid_file):
        self.trinucleotide_substitution_matrix.to_csv(path_or_buf=trinucleotide_file)
        self.amino_acid_substitution_matrix.to_csv(path_or_buf=amino_acid_file)

    def file_upload(self,substitution_count_matrix,amino_acid_count_matrix,use_file):
        if os.path.isfile(use_file[0]): 
            substitution_count_matrix.to_csv(path_or_buf=use_file[0], mode='a', header=False, index = False)
        else:
            substitution_count_matrix.to_csv(path_or_buf=use_file[0], index = False)
        if os.path.isfile(use_file[1]): 
            amino_acid_count_matrix.to_csv(path_or_buf=use_file[1], mode='a', header=False, index = False)
        else:
            amino_acid_count_matrix.to_csv(path_or_buf=use_file[1], index = False)

            
    def get_trinucleotide_matrix(self):
        return self.trinucleotide_substitution_matrix[self.trinucleotide_substitution_matrix.columns[:192]] 
    
    def get_trinucleotide_matrix_metadata(self):
        return self.trinucleotide_substitution_matrix[self.trinucleotide_substitution_matrix.columns[192:]] 
    
    def get_amino_acid_matrix(self):
        return self.amino_acid_substitution_matrix[self.amino_acid_substitution_matrix.columns[:-11]] 
    
    def get_amino_acid_matrix_metadata(self):
        return self.amino_acid_substitution_matrix[self.amino_acid_substitution_matrix.columns[-11:]] 