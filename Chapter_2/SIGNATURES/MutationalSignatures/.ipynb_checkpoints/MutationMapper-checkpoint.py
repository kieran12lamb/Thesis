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


class MutationMapper:

    def __init__(self, refPath,
                 multiAlignment=None,
                 multiAlignmentPath = None,
                 seq_count=None,
                 random_sequences=False,
                 chunk_size = 20,
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
        
        self.lineage_hierarchy_dict = self.lineage_hierarchy(lineage_path)

        self.use_file = use_file
        self.use_database = use_database
        if pseudo_references is not None:
            reference_dictionary = {}
            reference_dictionary["ref"] = np.array(list(self.ref.seq),dtype="string_")
            self.pseudo_references = AlignIO.read(pseudo_references, "fasta")
            for record in self.pseudo_references:
                reference_dictionary[record.id] = np.array(list(record.seq),dtype="string_")
                potential_alias = record.id.split(".")[0]
                if potential_alias in self.lineage_hierarchy_dict:
                    reference_dictionary[record.id.replace(potential_alias, self.lineage_hierarchy_dict[potential_alias])] = np.array(list(record.seq),dtype="string_")
            self.pseudo_references = reference_dictionary
        else:
            self.pseudo_references = None
            

        if multiAlignment is None:
            self.multiAlignmentPath = multiAlignmentPath
            if seq_count is None:
                align = AlignIO.read(self.multiAlignmentPath, "fasta")
                self.seq_count = int(len(align))
                if len(self.ref.seq) != len(align[0].seq):
                    for record in align:
                        if self.ref.name in record.id:
                            self.ref.seq = record.seq.upper()
                            print("AS SEQUENCE IS FROM ALIGNMENT< TRANSLATION MAY BE INCORRECT")
        else:
            if seq_count is None:
                align = self.multiAlignment
                self.align = align
                self.seq_count = int(len(align))
                if len(self.ref.seq) != len(align[0].seq):
                    print(self.ref.name)
                    for record in align:
                        if self.ref.name in record.id:
                            self.ref.seq = record.seq.upper()
                            print("AS SEQUENCE IS FROM ALIGNMENT< TRANSLATION MAY BE INCORRECT")
            self.seq_count = int(len(multiAlignment))

    def build(self):
        seq_count = self.seq_count
        chunk_size = self.chunk_size
        threads = self.threads
        time = datetime.now()
         # divide alignment to free up memory
        unit_align = 16000
        j = 0

        total_blocks = math.ceil(seq_count/unit_align)
        while(j<seq_count):
            print(f'Block {j/unit_align}/{total_blocks}, Time Taken: {datetime.now()-time}')
            if self.multiAlignment is None:
                if seq_count > j+unit_align:
                    self.align = AlignIO.read(self.multiAlignmentPath, "fasta")[j:j+unit_align]
                else:
                    self.align = AlignIO.read(self.multiAlignmentPath, "fasta")[j:seq_count]
            else:
                self.align = self.multiAlignment
            print(f'Block loaded :{datetime.now()-time}')
                
            try:
                self.alignnp = []
                self.lineages = []
                align = []
                for idx, seq in enumerate(self.align):
                    metadata = seq.id.split("|")
                    if int(metadata[6]) <= self.filters["length"]:
                        continue
                    self.alignnp.append(np.char.array(seq))
                    align.append(seq)
                    lineage =metadata[4]
                    full_lineage = ""
                    if "." in lineage:
                        print(f"Has a . {lineage}")
                        # Get first part of lineage name which may be an alias
                        potential_alias = lineage.split(".")[0]
                        # Check alias is in list of known aliases
                        print(potential_alias)
                        if potential_alias in self.lineage_hierarchy_dict:
                            print(f"{lineage} has an alias . {potential_alias}")
                            # Convert alias to hierarchical name
                            full_lineage = self.lineage_hierarchy_dict[potential_alias]
                            print(lineage,full_lineage)
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
                            print(lineage,parent_lineage)
                            if parent_lineage == "":
                                parent_lineage = "ref"
                                known_parent = True
                                print(f'Lineage {lineage} only has the reference as a closest ancestor')
                            if parent_lineage in self.pseudo_references:
                                known_parent = True
                                print(parent_lineage)
                            else:
                                parent_lineage = ".".join(split_parent[:-i])
                            i+=1
                    else:
                        parent_lineage = "ref"
                    self.lineages.append(parent_lineage)
                print(self.lineages)
                self.alignnp = np.char.array(self.alignnp)
                self.align = align
            except Exception as e:
                print("Not GISAID")
                print(e)
                self.alignnp = np.char.array([list(rec.seq.upper()) for rec in self.align])
                print(self.align)
            print(f'Block Filtered and Re-Typed :{datetime.now()-time}')
            sub_count_array = []
            amino_count_array = []
            i=0
            total_chunks = math.ceil(len(self.alignnp)/(chunk_size*threads))
            while(i<len(self.alignnp)):
                chunk_time = datetime.now()

                pool = multiprocessing.Pool(threads)
                if len(self.alignnp)> i+(chunk_size*threads):
                    table_chunk = self.table_chunk
                    results = pool.map(table_chunk,np.arange(i,i+(chunk_size*threads),chunk_size))
                else:
                    length = len(self.alignnp)
                    table_chunk = self.table_chunk
                    results = pool.map(table_chunk,np.arange(i,length,chunk_size))
                pool.close()
                pool.join()
                # results = results.get()
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


    def getdist(self,data):
        dist = [np.nan]
        for i in range(1,len(data.Pos)):
            if data.iloc[i].Pos - data.iloc[i-1].Pos == 0:
                dist.append(dist[len(dist)-1])
            else:
                dist.append(data.iloc[i].Pos - data.iloc[i-1].Pos)
        return dist

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
            # print(f'end:{end}')
            starts.append(self.getAminoPos(orfs, start))
            ends.append(self.getAminoPos(orfs, end))
            # print(self.getAminoPos(orfs, start), self.getAminoPos(orfs, end))
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

    def assignOrf(self,orf_table,position, index_type):
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
    
    def getAminoPos(self,orf_table,position):
        # Return the amino position of a nucleotide index
        for index in range(len(orf_table)):
            # Get the ORF the nucleotide is in
            if position >= orf_table.iloc[index].Start  and position<=orf_table.iloc[index].End:
                position = int((position-orf_table.iloc[index].Offset)/3)
                return position
        # No position due to non-coding segment
        return np.nan

    def table_chunk(self,i):
        chunk_size = self.chunk_size
        if i<len(self.alignnp) and i+chunk_size>len(self.alignnp):
            alignment=np.array(self.alignnp[i:],dtype="string_")
            alignment_seq = self.align[i:]
            lineage_subset = self.lineages[i:]
        else:
            alignment=np.array(self.alignnp[i:i+chunk_size],dtype="string_")
            alignment_seq = self.align[i:i+chunk_size]
            lineage_subset = self.lineages[i:i+chunk_size]
        if alignment.ndim ==1:
            alignment = np.array([alignment],dtype="string_")
            alignment_seq = self.align
            lineage_subset = self.lineages
        scm,aac  = cf.createTables( mutation_bases = alignment,
                                    ref = np.char.array(self.ref.seq.split(maxsplit=True)),
                                    classes = np.array(self.classes),
                                    orfs = self.Orfs,
                                    ref_seq = self.ref.seq,
                                    alignment_seq = alignment_seq,
                                    location = True,
                                    pseudo_references=self.pseudo_references,
                                    lineages = lineage_subset)
        scm['ORF'] = [self.assignOrf(self.Orfs,int(pos),"nucleotide") for pos in scm.POS]
        scm['AMINO_POS'] = [self.getAminoPos(self.Orfs,int(pos)) for pos in scm.POS]
        aac['ORF'] = [self.assignOrf(self.Orfs,int(pos),"amino-acid") for pos in aac.POS]
        print(scm)
        return scm, aac

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

    def lineage_hierarchy(self,lineage_path):
        lineages = pd.read_csv(lineage_path)
        
        lineages.lineage = lineages.lineage.str.replace(" ","")
        d = {}
        for row in lineages.itertuples():
            d[row.alias] = row.lineage
        return d
    