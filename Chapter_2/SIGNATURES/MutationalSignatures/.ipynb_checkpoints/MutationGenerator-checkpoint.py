import os
from Bio import AlignIO,SeqIO
from Bio.Seq import Seq
import numpy as np
import MutationalSignatures.compiled_functions as cf
import pandas as pd
import itertools
from functools import reduce
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import multiprocessing
from datetime import datetime
import math
import os.path
import json


class MutationGenerator:

    def __init__(self, 
                 refPath,
                 multiAlignment=None,
                 multiAlignmentPath = None,
                 seq_count=None,
                 random_sequences=False,
                 chunk_size = 20,
                 block_size = 16000,
                 threads=8,
                 use_database=None,
                 use_file=None,
                 filters=None,
                 pseudo_references=None,
                 lineage_path=None):
       #Generate mutational signature tables
        self.substitution_count_matrix = pd.DataFrame()
        self.amino_acid_count_matrix = pd.DataFrame()

        self.makeSubstitutionClasses()
        self.ref = SeqIO.read(refPath, "genbank")
        self.Orfs = self.makeOrfTable(self.ref)
        self.multiAlignment = multiAlignment
        
        self.filters = filters
        self.threads = threads
        self.chunk_size = chunk_size
        self.block_size = block_size
        
        self.use_file = use_file
        self.use_database = use_database
        if pseudo_references is not None:
            self.lineage_hierarchy_dict = self.lineage_hierarchy(lineage_path)
            self.pseudo_references = self.build_reference_dict(self.ref.seq,pseudo_references,self.lineage_hierarchy_dict)
        else:
            self.pseudo_references = {}
            self.lineage_hierarchy_dict = {}
        if multiAlignment is None:
            self.multiAlignmentPath = multiAlignmentPath
            if seq_count is None:
                self.full_align  = AlignIO.read(self.multiAlignmentPath, "fasta")
                self.seq_count = int(len(self.full_align))
                if len(self.ref.seq) != len(self.full_align[0].seq):
                    for record in self.full_align:
                        if self.ref.name in record.id:
                            self.ref.seq = record.seq.upper()
                            print("AS SEQUENCE IS FROM ALIGNMENT< TRANSLATION MAY BE INCORRECT")
        else:
            self.full_align  = self.multiAlignment
            self.align = self.multiAlignment
            self.seq_count = int(len(multiAlignment))


    def build(self):
        seq_count = self.seq_count
        chunk_size = self.chunk_size
        threads = self.threads
#         self.full_align = AlignIO.read(self.multiAlignmentPath, "fasta")
        time = datetime.now()
         # divide alignment to free up memory
        unit_align = self.block_size
        j = 0
        total_blocks = math.ceil(seq_count/unit_align)
        while(j<seq_count):
            print(f'Block {j/unit_align}/{total_blocks}, Time Taken: {datetime.now()-time}')
            if self.multiAlignment is None:
                if seq_count > j+unit_align:
                    self.align = self.full_align[j:j+unit_align]
                else:
                    self.align = self.full_align[j:seq_count]
            print(f'Block loaded :{datetime.now()-time}')
 
            sub_count_array = []
            amino_count_array = []
            i=0
            total_chunks = math.ceil(len(self.align)/(chunk_size*threads))
            while(i<len(self.align)):
                chunk_time = datetime.now()
                pool = multiprocessing.Pool(threads)
                if len(self.align)> i+(chunk_size*threads):
                    ranges = np.arange(i,i+(chunk_size*threads),chunk_size)
                    starmap_input = [self.table_chunk(i)for i in ranges]
                    results = pool.starmap(cf.find_mutations,starmap_input)
                else:
                    length = len(self.align)
                    ranges = np.arange(i,length,chunk_size)
                    starmap_input = [self.table_chunk(i)for i in ranges]
                    results = pool.starmap(cf.find_mutations,starmap_input)
                pool.close()
                pool.join()
                results = list(zip(*results))
                substitution_count_matrix = pd.concat(results[0])
                amino_acid_count_matrix = pd.concat(results[1])

                # Write to file/database every chunk
                if self.use_database is not None:
                    self.database_upload(substitution_count_matrix,amino_acid_count_matrix,self.use_database)
                elif self.use_file is not None:
                    self.file_upload(substitution_count_matrix,amino_acid_count_matrix,self.use_file)
                else:
                    self.substitution_count_matrix = pd.concat([self.substitution_count_matrix,substitution_count_matrix])
                    self.amino_acid_count_matrix = pd.concat([self.amino_acid_count_matrix,amino_acid_count_matrix])

                i=i+(chunk_size*threads)
                print(f'Chunk {i/(chunk_size*threads)}/{total_chunks}, Time Taken: {datetime.now()-chunk_time}')
            j+=unit_align
        print(f'Time Taken: {datetime.now()-time}')

    def build_reference_dict(self, ref, pseudo_references, hierarchical_lineages):
        reference_dictionary = {}
        reference_dictionary["ref"] = ref
        pseudo_references = AlignIO.read(pseudo_references, "fasta")
        for record in pseudo_references:
            reference_dictionary[record.id] = record.seq
            potential_alias = record.id.split(".")[0]
            if potential_alias in hierarchical_lineages:
                reference_dictionary[record.id.replace(potential_alias, hierarchical_lineages[potential_alias])] = record.seq
        return reference_dictionary

    def lineage_hierarchy(self,lineage_path):
        if lineage_path is None:
            return {}
        with open(lineage_path) as json_file:
            data = json.load(json_file)
        
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

    def makeSubstitutionClasses(self):
        #Array of substitution classes
        substitutionClasses = ['CA', 'CG', 'CT', 'TA','TG', 'TC', 'AC', 'AG', 'AT', 'GA', 'GT', 'GC']
        self.classes = pd.DataFrame()
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
                self.classes = self.classes.append(pd.DataFrame(row).T)
                self.classes = self.classes.reset_index(drop=True)
        self.classes['sub-trimer'] = self.classes['Substitution']+"-"+self.classes['Trimer']
        self.classes['sub-amino'] = self.classes['Substitution']+"-"+self.classes['Amino_Acid']
        self.classes['amino-amino'] = self.classes['Amino_Acid']+"-"+self.classes['Mutant_Amino_Acid']

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
        orfs = self.calculate_orf_offsets(orfs)
        starts = []
        ends = []
        for index in range(len(orfs)): 
            start = int(orfs.iloc[index]['Start'])
            end = int(orfs.iloc[index]['End'])
            starts.append(cf.getAminoPos(orfs, start))
            ends.append(cf.getAminoPos(orfs, end))
        orfs["AminoStart"] = starts
        orfs["AminoEnd"] = ends
        return orfs

    def calculate_orf_offsets(self,orf_table):
        orf_table['Offset'] = 0
        for index in range(1,len(orf_table)):
            if index == 1:
                orf_table["Offset"][index-1] = orf_table.Start[index-1]
            # if orf_table.Start[index] >orf_table.End[index-1]:
            orf_table["Offset"][index] = orf_table.Offset[index-1]+ (orf_table.Start[index] - orf_table.End[index-1])
            # else:
            #     orf_table["Offset"][index] = orf_table.Offset[index-1]
        return orf_table

    def table_chunk(self,i):
        chunk_size = self.chunk_size
        if i<len(self.align) and i+chunk_size>len(self.align):
            alignment = self.align[i:]
        else:
            alignment = self.align[i:i+chunk_size]
        if len(alignment) ==1:
            alignment = self.align
        return [alignment,self.ref,self.pseudo_references,self.lineage_hierarchy_dict,self.Orfs,self.filters]

    def database_upload(self,substitution_count_matrix,amino_acid_count_matrix,use_database):
        use_database[0].insert_many(substitution_count_matrix.to_dict('records'))
        use_database[1].insert_many(amino_acid_count_matrix.to_dict('records'))

    def file_upload(self,substitution_count_matrix,amino_acid_count_matrix,use_file):
        if os.path.isfile(use_file[0]): 
            substitution_count_matrix.to_csv(path_or_buf=use_file[0], mode='a', header=False, index = False)
        else:
            substitution_count_matrix.to_csv(path_or_buf=use_file[0], index = False)
        if os.path.isfile(use_file[1]): 
            amino_acid_count_matrix.to_csv(path_or_buf=use_file[1], mode='a', header=False, index = False)
        else:
            amino_acid_count_matrix.to_csv(path_or_buf=use_file[1], index = False)
        


