import FASTX
import BioSequences
using  GenomicAnnotations
import BioGenerics
import JSON
using ArgParse
using DataFrames
using Combinatorics
using DataStructures
using CSV
using DelimitedFiles
using Dates
# using ThreadedIterables
# Test Files to Run

# julia FastaMutationsJulia.jl --reference_file "/Users/kieranlamb/Documents/Github/Signatures_Julia/Resources/References/GISAID_REFERENCE.gb" --fasta_file "/Users/kieranlamb/Documents/Github/Signatures_Julia/Resources/PseudoReferences/pseudo_references.fasta" --output "/Users/kieranlamb/Documents/Github/Signatures_Julia/Outputs/" --alias_file "/Users/kieranlamb/Documents/Github/Signatures_Julia/Resources/LineageInformation/pango_aliases.json"

# julia FastaMutationsJulia.jl --reference_file "/home4/2191618l/Github/Signatures_Julia/Outputs/H3N2/mpox.gb" --fasta_file "/home4/2191618l/Github/Signatures_Julia/Outputs/mpox/inferred_nodes_annotated.fa" --output "/home4/2191618l/Github/Signatures_Julia/Outputs/mpox/" --alias_file "/home4/2191618l/Github/Signatures_Julia/Outputs/mpox/alias.json" --pseudoreference_folder "/home4/2191618l/Github/Signatures_Julia/Resources/Outputs/mpox/inferred_nodes_annotated.fa"

# julia FastaMutationsJulia.jl --reference_file "/home4/2191618l/Github/Signatures_Julia/Outputs/mpox/H3N2.gb" --fasta_file "/home4/2191618l/Github/Signatures_Julia/Outputs/H3N2/H3N2_Pseudoreferences.fasta" --output "/home4/2191618l/Github/Signatures_Julia/Outputs/H3N2/" --alias_file "/home4/2191618l/Github/Signatures_Julia/Outputs/H3N2/alias.json"

# julia FastaMutationsJulia.jl --reference_file "/home4/2191618l/Github/Signatures_Julia/Resources/References/GISAID_REFERENCE.gb" --fasta_file "../../../../home5/nCov/Richard/gisaid/20220128/7_aligned.fasta" --pseudoreference_folder "/home4/2191618l/Github/Signatures_Julia/Resources/PseudoReferences/" --output "/home4/2191618l/Github/Signatures_Julia/Outputs/January_2022/" --alias_file "/home4/2191618l/Github/Signatures_Julia/Resources/LineageInformation/pango_aliases.json"

# julia FastaMutationsJulia.jl --reference_file "/home4/2191618l/Github/Signatures_Julia/Outputs/mpox/mpox.gb" --fasta_file "/home4/2191618l/Github/Signatures_Julia/Outputs/mpox/mpox_alignment.fasta" --pseudoreference_folder "/home4/2191618l/Github/Signatures_Julia/Outputs/mpox/PseudoReferences/" --output "/home4/2191618l/Github/Signatures_Julia/Outputs/mpox/" --alias_file "/home4/2191618l/Github/Signatures_Julia/Outputs/mpox/alias.json"


# julia FastaMutationsJulia.jl --reference_file "/home4/2191618l/Github/Signatures_Julia/Resources/References/GISAID_REFERENCE.gb" --fasta_file "/home4/2191618l/Github/Signatures_Julia/Resources/Alignments/Deer/al_all_deer_Apr21_1.fas" --pseudoreference_folder "/home4/2191618l/Github/Signatures_Julia/Outputs/October_2022_Filtered/" --output "/home4/2191618l/Github/Signatures_Julia/Resources/Alignments/Deer/" --alias_file "/home4/2191618l/Github/Signatures_Julia/Resources/LineageInformation/pango_aliases.json" --date_file "/home4/2191618l/Github/Signatures_Julia/Outputs/October_2022_Filtered/pango_dates.csv"


function parse_commandline()
    settings = ArgParseSettings()
    @add_arg_table settings begin
        "--reference_file","-r"
            help = "The reference file in genbank format."
            arg_type = String
            default = ""
        "--fasta_file", "-f"
            help = "The alignment file."
            arg_type = String
            default = ""
        "--pseudoreference_folder", "-p"
            help = "The pseudoreference folder."
            arg_type = String
            default = ""
        "--alias_file", "-a"
            help = "The lineage alias file."
            arg_type = String
            default = ""
        "--output_file","-o"
            help = "The location for the mutations file "
            arg_type = String
            default = ""
        "--context","-c"
            help = "The context length for mutations "
            arg_type = Int
            default = 3
        "--signature_matrix","-s"
            help = "Create signature matrix output"
            arg_type = Bool
            default = true
        "--phylogenetic_normalisation","-n"
            help = "Create phylogeneticly normalised output"
            arg_type = Bool
            default = true
        "--date_file","-d"
            help = "The files containing lineage origin dates."
            arg_type = String
            default = ""
        
    end
    return parse_args(settings)
end;

function get_mutations(ref, alt,context_length,orf_list,orf_dictionary,record_sequence_name)
    mutations =  Dict([("Sequence_Name",[]),
                            ("REF",[]),
                            ("ALT",[]),
                            ("POS",[]),
                            ("BIO_POS",[]),
                            ("AMINO_ORF_POS",[]),
                            ("Nucleotide",[]),
                            ("Amino_Acid",[]),
                            ("Nucleotide_Context",[]),
                            ("Amino_Acid_Context",[]),
                            ("Codon_Position",[]),
                            ("Codon",[]),
                            ("ORF",[]),
                            ("Synonymous",[]),
                            ])
    ref = string(ref)
    alt = string(alt)
    valid_nucleotides = Set(["A","C","G","T","-"])
    for (index,nucleotide) in enumerate(ref)
        if nucleotide != alt[index] && string(nucleotide) in valid_nucleotides && string(alt[index]) in valid_nucleotides
            if string(alt[index]) != "N" && string(alt[index]) != "X"
                nucleotide_mutation = ""
                amino_acid_mutation = ""
                ref_codon = ""
                alt_codon = ""
                middle_codon_position = ""
                synonymous = true

                nucleotide_mutation = string(string(nucleotide),string(index),string(alt[index]))
                # Loop Through different orfs (Usually only 1)
                for orf in orf_list[index]
                    if (orf == "Non-Coding")
                        append!(mutations["Sequence_Name"],[record_sequence_name])
                        append!(mutations["REF"],[string(nucleotide)])
                        append!(mutations["ALT"],[string(alt[index])])
                        append!(mutations["POS"],[index-1])
                        append!(mutations["BIO_POS"],[index])
                        append!(mutations["AMINO_ORF_POS"],[""])
                        append!(mutations["Nucleotide"],[nucleotide_mutation])
                        append!(mutations["Nucleotide_Context"],[string(string(nucleotide),string(alt[index]),"-",string(get_mutation_context(index,ref,context_length,"Nucleotide"))) ])
                        append!(mutations["Amino_Acid"],[""])
                        append!(mutations["Amino_Acid_Context"],[""])
                        # append!(mutations["Codon Aware Amino Acid Context"],[""])
                        append!(mutations["Codon_Position"],[""])
                        append!(mutations["Codon"],[""])
                        append!(mutations["ORF"],[orf])
                        append!(mutations["Synonymous"],[true])
                        continue
                    else
                        #Loop through sub-orfs
                        for range in orf_dictionary[orf]
                            ref_orf_range_sequence = ref[range]
                            alt_orf_range_sequence = alt[range]

                            relative_orf_pos = 1+(index-range[1])
                            if (relative_orf_pos < 1 || relative_orf_pos > length(ref_orf_range_sequence))
                                continue
                            end
                            codon_pos = relative_orf_pos%3
                            
                            # Codon position wont matter on the last/first 3 positions
                            if (relative_orf_pos > length(ref_orf_range_sequence)-3)
                                ref_codon = ref_orf_range_sequence[length(ref_orf_range_sequence)-2:length(ref_orf_range_sequence)]
                                alt_codon = alt_orf_range_sequence[length(alt_orf_range_sequence)-2:length(alt_orf_range_sequence)]
                            elseif relative_orf_pos <=3
                                ref_codon = ref_orf_range_sequence[1:3]
                                alt_codon = alt_orf_range_sequence[1:3]
                            else
                                if codon_pos == 1
                                    ref_codon =  ref_orf_range_sequence[relative_orf_pos:relative_orf_pos+2]
                                    alt_codon = alt_orf_range_sequence[relative_orf_pos:relative_orf_pos+2]
                                elseif codon_pos == 2
                                    ref_codon =  ref_orf_range_sequence[relative_orf_pos-1:(relative_orf_pos-1)+2]
                                    alt_codon = alt_orf_range_sequence[relative_orf_pos-1:(relative_orf_pos-1)+2]
                                else
                                    ref_codon =  ref_orf_range_sequence[relative_orf_pos-2:relative_orf_pos]
                                    alt_codon = alt_orf_range_sequence[relative_orf_pos-2:relative_orf_pos]
                                    
                                end 
                            end
                            
                            if codon_pos == 1
                                middle_codon_position = relative_orf_pos+1
                                relative_orf_amino_pos =(relative_orf_pos÷3 )+1
                            elseif codon_pos == 2
                                middle_codon_position = relative_orf_pos
                                relative_orf_amino_pos =(relative_orf_pos÷3 )+1
                            else
                                middle_codon_position = relative_orf_pos-1
                                codon_pos = 3
                                # Integer divisionis ÷ in julia
                                relative_orf_amino_pos = (relative_orf_pos÷3 )
                            end


                            ref_amino = translate_with_gaps(string(ref_codon))
                            alt_amino = translate_with_gaps(string(alt_codon))
            
                            if ref_amino != alt_amino
                                amino_acid_mutation = string(string(orf)string(":")string(ref_amino),string(relative_orf_amino_pos),string(alt_amino))
                                synonymous = false
                            end
                            Nucleotide_Context = string(string(nucleotide),string(alt[index]),"-",string(get_mutation_context(relative_orf_pos,ref_orf_range_sequence,context_length,"Nucleotide"))) 
                            # Amino_Acid_Context = string(string(ref_amino),string(alt_amino),"-",string(get_mutation_context(relative_orf_pos,ref_orf_range_sequence,context_length,"Amino Acid")))
                            Codon_Aware_Amino_Acid_Context =  string(string(ref_amino),string(alt_amino),"-",string(get_mutation_context(middle_codon_position,ref_orf_range_sequence,context_length,"Amino Acid")))

                            cased_codon = String([i != codon_pos ? lowercase(c) : c for (i, c) in enumerate(ref_codon)])
                            append!(mutations["Sequence_Name"],[record_sequence_name])
                            append!(mutations["REF"],[string(nucleotide)])
                            append!(mutations["ALT"],[string(alt[index])])
                            append!(mutations["POS"],[index-1])
                            append!(mutations["BIO_POS"],[index])
                            append!(mutations["AMINO_ORF_POS"],[relative_orf_amino_pos])
                            append!(mutations["Nucleotide"],[nucleotide_mutation])
                            append!(mutations["Nucleotide_Context"],[Nucleotide_Context])
                            append!(mutations["Amino_Acid"],[amino_acid_mutation])
                            append!(mutations["Amino_Acid_Context"],[Codon_Aware_Amino_Acid_Context])
                            # append!(mutations["Codon Aware Amino Acid Context"],[Codon_Aware_Amino_Acid_Context])
                            append!(mutations["Codon_Position"],[codon_pos])
                            append!(mutations["Codon"],[string(cased_codon)])
                            append!(mutations["ORF"],[string(orf)])
                            append!(mutations["Synonymous"],[synonymous])
                        end
                    end
                end

            end
        end
    end
    return mutations
end;

function translate_with_gaps(sequence)
    # println(sequence)
    amino_acid = ""
    for i in range(1,stop=length(sequence)-2,step=3)
        codon = string(sequence[Integer(i):Integer(i+2)])
        codon = replace(codon, "?" => "N")
        if occursin("-", codon)
            if cmp(codon,"---") == 0
                amino_acid = amino_acid * "-"
            else
                amino_acid = amino_acid * "X"
            end
        elseif occursin(" ", codon)
            amino_acid = amino_acid * " "
        else
            amino_acid = amino_acid*string(BioSequences.translate(BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}(codon)))
        end
    end      
    return amino_acid
end;

function convert_lineage_name(lineage,lineage_aliases)
    if haskey(lineage_aliases, split(lineage,".")[1]) 
        alias = [lineage_aliases[split(lineage,".")[1]]]
        remainder = split(lineage,".")[2:end]
        return join( vcat(alias,remainder) ,".")
    elseif length(split(lineage,".")) >4 && haskey(lineage_aliases,join(split(lineage,".")[1:5]))
        alias = [lineage_aliases[split(lineage,".")[1:5]]]
        remainder = split(lineage,".")[5:end]
        return join(vcat(alias,remainder) ,".")
    end
    return lineage
end;

function make_pseudoreference_dict(fasta)
    pseudoreference_dict = Dict()
    for record in fasta
        pseudoreference_dict[seqheader(record)] = BioGenerics.sequence(record)
    end
    return pseudoreference_dict
end;

function make_pseudoreference_proportion_dict(csv_file)
    loaded_csv = readdlm(csv_file, ',') |> Matrix{Any}
    pseudoreference_proportion_dict = Dict()
    for (index,row) in enumerate(eachrow(loaded_csv))
        if index == 1
            continue
        else
            lineage = row[length(row)]
            nucleotide = row[1]
            counts = row[2:length(row)-1]
            if (haskey(pseudoreference_proportion_dict,lineage))
                pseudoreference_proportion_dict[lineage][nucleotide] = counts
            else
                pseudoreference_proportion_dict[lineage] = Dict()
                pseudoreference_proportion_dict[lineage][nucleotide] = counts
            end
        end
    end
    return pseudoreference_proportion_dict
end;

function get_closest_ancestor(lineage_alias,pseudoreference_dict)
    
    parents= []
    for i in range(1,stop=length(split(lineage_alias,".")),step=1)
        new_parent = join(split(lineage_alias,".")[1:i],".")
        append!(parents,[new_parent])
    end
    parents = reverse(parents)[2:end]
    for parent in parents
        println(parent)
        if haskey(pseudoreference_dict,parent)
            return parent
        end
    end
    return "Ref"
end;

function make_substitution_matrix(context)
    columns = []
    context_length = context-1
    substitutionClasses = ["CA", "CG", "CT", "TA", "TG", "TC", "AC", "AG", "AT", "GA", "GT", "GC"]
    for substitution in substitutionClasses
        # Fix this once you understand how the arguments work, cant give list must be iterators
        if context == 5
            context_combinations = vec(join.(Iterators.product("ACGT","ACGT","ACGT","ACGT")))
        elseif context == 3
            context_combinations = vec(join.(Iterators.product("ACGT","ACGT")))
        end

        for p in context_combinations
            trimer = p[1:Int(context_length/2)]*substitution[1]*p[Int((context_length/2)+1):end]
            append!(columns,[string(substitution,"-",trimer)])
        end
    end
    sort!(columns)
    substitution_matrix = DataFrame()
    empty_row = Dict{String,Any}(i => 0 for i in columns)
    empty_row["Sequence_Name"] = ""
    append!(substitution_matrix,DataFrame(empty_row))
    sort!(substitution_matrix)
    delete!(substitution_matrix, 1)
    return substitution_matrix
end;

function amino_acid_mutations(record_sequence,reference,gene,context)
    Amino_Acid_Mutations = []
    Amino_Acid_Contexts = []
    translated_record_sequence = ""
    translated_reference_sequence = ""
    if length(GenBank.locus(gene).order) > 0
        sub_orfs = GenBank.locus(gene).order
        for sub_orf in sub_orfs
            translated_reference_sequence = translated_reference_sequence*string(translate_with_gaps(reference[sub_orf]))
            translated_record_sequence = translated_record_sequence*string(translate_with_gaps(record_sequence[sub_orf]))
        end
    else
        translated_reference_sequence = translate_with_gaps(reference[GenBank.locus(gene).position])
        translated_record_sequence = translate_with_gaps(record_sequence[GenBank.locus(gene).position])
        # println(translated_record_sequence)
    end

    Amino_Acids, Amino_Acid_Context = get_mutations(translated_reference_sequence,translated_record_sequence,context,string(gene.gene))
    append!(Amino_Acid_Mutations,Amino_Acids)
    append!(Amino_Acid_Contexts,Amino_Acid_Context)
    return Amino_Acid_Mutations, Amino_Acid_Contexts
end

function map_orfs_to_position_array(main_reference)
    orf_at_position = []
    for (pos,nuc) in enumerate(main_reference.sequence)
        in_orf = false
        orf_names = Set()
        orf_list = []
        for gene in main_reference.genes
            
            if cmp(string(GenBank.feature(gene)), "CDS") == 0
                #Add product to orfs that have multiple translations
                orf_name = string(gene.gene)
                # println(orf_name)
                if in(orf_name, orf_names)
                    orf_name = string(gene.gene)*"("*string(gene.product)*")"
                end
                push!(orf_names,orf_name)
                #Check position is in orf
                if length(GenBank.locus(gene).order) > 0
                    sub_orfs = GenBank.locus(gene).order
                    for sub_orf in sub_orfs
                        # println(sub_orf)
                        if (in(pos,sub_orf))
                            # println("Position "*string(pos)*" is in "*string(orf_name))
                            push!(orf_list,orf_name)
                            in_orf = true
                        end
                    end
                else
                    if (in(pos,GenBank.locus(gene).position))
                        # println("Position "*string(pos)*" is in "*string(orf_name))
                        push!(orf_list,orf_name)
                        in_orf = true
                    end
                end
            end
        end
        if (in_orf == false)
            push!(orf_list,"Non-Coding")
        end
        push!(orf_at_position,orf_list)
    end
    return orf_at_position
end;

function make_orf_range_dict(main_reference)
    orf_dict = Dict{String,Any}()
    orf_names = Set()
    for gene in main_reference.genes
        if cmp(string(GenBank.feature(gene)), "CDS") == 0
            orf_list = []
            #Add product to orfs that have multiple translations
            orf_name = string(gene.gene)
            # println(orf_name)
            if in(orf_name, orf_names)
                orf_name = string(gene.gene)*"("*string(gene.product)*")"
            end
            push!(orf_names,orf_name)
            #Check position is in orf
            if length(GenBank.locus(gene).order) > 0
                sub_orfs = GenBank.locus(gene).order
                for sub_orf in sub_orfs
                    push!(orf_list,sub_orf)
                end
            else
                push!(orf_list,GenBank.locus(gene).position)
            end
            orf_dict[orf_name] = orf_list
        end
    end
    return orf_dict
end;

function get_mutation_context(position,ref_relative_sequence,context_length,mutation_type)
    #Calculate number of gaps to place to pad context at the start/end
    if (mutation_type == "Amino Acid")
        context_length = (context_length*3)
    end

    context_string = ""
    # println(position)
    if position < context_length
        end_distance = position-1
        padding_length = abs(Int(((context_length-1)/2)-end_distance))
        padding = " "^padding_length
        stop = Int(position+((context_length-1)/2))
        context_string = padding*string(ref_relative_sequence[position:stop])
    elseif context_length > length(ref_relative_sequence)-position 
        end_distance = length(ref_relative_sequence)-position
        padding_length = abs(Int(((context_length-1)/2)-end_distance))
        padding = " "^padding_length
        context_string =ref_relative_sequence[position:length(ref_relative_sequence)]*padding
    else
        start = Int(position-((context_length-1)/2))
        stop = Int(position+((context_length-1)/2))
        context_string = ref_relative_sequence[start:stop]
    end
    if mutation_type == "Nucleotide"
        return context_string
    elseif mutation_type == "Amino Acid"
        return translate_with_gaps(context_string)
    end
end;

function check_inheritence(row,parent_lineages,pseudo_references,lineage_aliases)
    position = int(row["POS"])
    for parent_lineage in parent_lineages
        current_parent = get_nearest_known_parent_lineage(parent_lineage,lineage_aliases)
        current_parent_base = pseudo_references[current_parent][position]
        row_base = row["REF"]
    end
    return false
end;

function check_valid_date(sample_date,sample_lineage, date_file_df, date_format,pseudoreference_dict)
    if sample_lineage == "" || sample_date == ""
        return false
    end
    
    #Filter date dataframe to retrieve lineage assignment date
    @inbounds lineage_date = DataFrame(filter(date_file_df -> date_file_df.Lineage == sample_lineage, date_file_df))
    
    if  nrow(lineage_date) == 0 || ismissing(lineage_date[:,"Earliest date"][1])
        next_best_lineage = get_closest_ancestor(sample_lineage,pseudoreference_dict)
        #Filter date dataframe to retrieve lineage assignment date
        @inbounds lineage_date = DataFrame(filter(date_file_df -> date_file_df.Lineage == next_best_lineage, date_file_df))
        if nrow(lineage_date) == 0 || ismissing(lineage_date[:,"Earliest date"][1])
            return false
        end
        lineage_date =  Date(lineage_date[:,"Earliest date"][1])
    else
        lineage_date =  Date(lineage_date[:,"Earliest date"][1])
    end
        
    valid = false
    if occursin("?", sample_date)
        valid = false
    else
        valid = Date(sample_date,date_format) >= Date(lineage_date)
    end
    # if valid == false
    #     println("Invalid Date of sequence for lineage "*sample_lineage*", origin date is "*string(lineage_date)*" but sample date is "*string(sample_date))
    # end 
    return valid
end;

function is_recombinant_lineage(lineage,lineage_aliases)
    if haskey(lineage_aliases, split(lineage,".")[1]) 
        alias = lineage_aliases[split(lineage,".")[1]]
        if alias isa String == false
            return true
        end
    end
    return false
end;

function seqheader(record)
    if (FASTX.FASTA.hasdescription(record) == false)
        record_sequence_name = BioGenerics.seqname(record)
    else
        #Sequence names sometimes have spaces 
        record_sequence_name = string(BioGenerics.seqname(record))*" "*string(FASTX.FASTA.description(record))
        record_sequence_name = replace(record_sequence_name," "=>"_")
    end
    return record_sequence_name
end

function main()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Pre-processing Stage: Parse Arguments and load files
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    parsed_args = parse_commandline()
    fasta_file = ""
    reference_file = ""
    output = ""
    pseudoreference_folder = ""
    phylogenetic_normalisation = ""
    alias_file = ""
    # header_format = ["Sequence_Name","GISAID_ID","Host","Clade","Alias","Lineage","Length","Sample_Date","Epi_Week","Country","FullCountry","Error"]
#     header_format = ["ID", "Host", "Country", "Serotype", "Sample_Date", "Strain", "GI", "aa107","aa325", "Clade", "ap2m1", "ap2mx", "sumo1","ap2m2", "casp1", "ap2m3" , "Lineage"]
#     header_format=["Sample_Name","Date","
    # header_format=["GISAID_ID","Sample_Name","Alias","Sample_Date","VOC","Lineage"]
    header_format=["GISAID_ID","Lineage","Alias"]
    context = 3
    signature_matrix = true
    last_index = 1
    date_format = "yyyy-mm-dd"
    date_file = ""

    for (arg,val) in parsed_args
        if arg == "fasta_file"
            fasta_file = val
        elseif arg == "reference_file"
            reference_file = val
        elseif arg == "pseudoreference_folder"
            pseudoreference_folder = val
        elseif arg == "alias_file"
            alias_file = val
        elseif arg == "header_format"
            header_format = val
        elseif arg == "context"
            context = val
        elseif arg == "signature_matrix"
            signature_matrix = val
        elseif arg == "phylogenetic_normalisation"
            phylogenetic_normalisation = val
        elseif arg == "date_file"
            date_file = val
        else
            output = val
        end
    end
     
    substitution_matrix = make_substitution_matrix(context)
    empty_row = Dict{String,Any}(i => 0 for i in names(substitution_matrix))
    main_reference = GenomicAnnotations.readgbk(reference_file)[1]
    fasta = open(FASTX.FASTA.Reader,fasta_file)
    mutation_set = ""
    lineage_aliases = Dict()
    pseudoreference_dict = Dict()
    record_list = []

    #Read lineage aliases
    if alias_file != ""
        lineage_aliases = JSON.parse(open(alias_file))
        filter!(p->(p[2] !=""),lineage_aliases)
    end

    if pseudoreference_folder != ""
        pseudoreference_dict = make_pseudoreference_dict(open(FASTX.FASTA.Reader,pseudoreference_folder*"pseudo_references.fasta"))
        pseudoreference_dict["Ref"] = main_reference.sequence
        #Not used anymore
        #pseudoreference_proportion_dict = make_pseudoreference_proportion_dict(pseudoreference_folder*"pseudo_references_counts.csv")
    end

    if date_file != ""
        #Read in the date file that lets us filter for mutations that are wrongly annotated
        date_file_df = DataFrame(CSV.File(date_file,dateformat=date_format))
    end
    
    orf_at_position = map_orfs_to_position_array(main_reference)
    orf_range_dict = make_orf_range_dict(main_reference)


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Stage 1: Main Loop where mutations for each sequence are calculated, processed and written to files
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    for (index,record) in enumerate(fasta)
        
        #Get Sequence and Sequence Name
        record_sequence_name = seqheader(record)
        record_sequence = BioGenerics.sequence(record)
        
        #Get Metadata from header 
        header_metadata = split(record_sequence_name,"|")
        
        if (length(header_metadata) == length(header_format))
            if "Lineage" in header_format
                println("Lineage here")
                lineage = header_metadata[findfirst( x -> x == "Lineage", header_format)]
                #Check if date is valid for the lineage (If it isnt, skip to the next sequence)
                if "Sample_Date" in header_format && date_file != ""
                    sample_date = header_metadata[findfirst( x -> x == "Sample_Date", header_format)]
                    # println(sample_date)
                    if "Alias" in header_format
                        alias_lineage_name = header_metadata[findfirst( x -> x == "Alias", header_format)]
                        if alias_lineage_name == lineage
                            lineage = convert_lineage_name(alias_lineage_name,lineage_aliases)
                        end
    
                        if check_valid_date(sample_date,alias_lineage_name,date_file_df,date_format,pseudoreference_dict) == false
                            continue
                        end
                    else
                        if check_valid_date(sample_date,lineage,date_file_df,date_format,pseudoreference_dict) == false
                            continue
                        end
                    end
                end
            end
        else
            println(header_metadata)
            continue
        end

       #Use single reference
        if cmp(pseudoreference_folder, "") == 0
            reference = string(main_reference.sequence)
            #Get the mutations for the sequence 
            mutations = get_mutations(reference,record_sequence,context,orf_at_position,orf_range_dict,record_sequence_name)
            mutations_df = DataFrame(Sequence_Name = mutations["Sequence_Name"],
                                    REF = mutations["REF"],
                                    ALT = mutations["ALT"],
                                    POS = Nucleotide = mutations["POS"],
                                    BIO_POS = Nucleotide = mutations["BIO_POS"],
                                    AMINO_ORF_POS = Nucleotide = mutations["AMINO_ORF_POS"],
                                    Nucleotide = mutations["Nucleotide"],
                                    Nucleotide_Context = mutations["Nucleotide_Context"],
                                    Amino_Acid = mutations["Amino_Acid"],
                                    Amino_Acid_Context = mutations["Amino_Acid_Context"],
                                    Codon_Position = mutations["Codon_Position"],
                                    Codon = mutations["Codon"],
                                    ORF = mutations["ORF"],
                                    Synonymous = mutations["Synonymous"])
            #No Ambiguous mutations since reference known
            mutations_df[!, :Ambiguous] .= false
        #Pseudo-referencing is being used (Recommended)
        else
            if is_recombinant_lineage(lineage,lineage_aliases)
                parent_lineages = lineage_aliases[lineage]
                mutations_df = ""
                for parent_lineage in parent_lineages
                    reference = string(pseudoreference_dict[parent_lineage]) 
                    #Get the mutations for the sequence 
                    mutations = get_mutations(reference,record_sequence,context,orf_at_position,orf_range_dict,record_sequence_name)
                    inner_mutations_df = DataFrame(Sequence_Name = mutations["Sequence_Name"],
                                        REF = mutations["REF"],
                                        ALT = mutations["ALT"],
                                        POS = Nucleotide = mutations["POS"],
                                        BIO_POS = Nucleotide = mutations["BIO_POS"],
                                        AMINO_ORF_POS = Nucleotide = mutations["AMINO_ORF_POS"],
                                        Nucleotide = mutations["Nucleotide"],
                                        Nucleotide_Context = mutations["Nucleotide_Context"],
                                        Amino_Acid = mutations["Amino_Acid"],
                                        Amino_Acid_Context = mutations["Amino_Acid_Context"],
                                        Codon_Position = mutations["Codon_Position"],
                                        Codon = mutations["Codon"],
                                        ORF = mutations["ORF"],
                                        Synonymous = mutations["Synonymous"])
                    #Concatenate dataframes only theres been at least one mutations_df dataframe, otherwise initialise it
                    if mutations_df == "" 
                        mutations_df = inner_mutations_df
                    else
                        append!(mutations_df,inner_mutations_df)
                    end
                end
                sort!(mutations_df,[:POS])
                # Remove duplicate mutations that are present regardless of parent (these are mutatations we know to be real)
                mutations_df[!,"Duplicate"] = nonunique(mutations_df,["REF","ALT","POS"])
                #Non-Ambiguous mutations can be identifiedas the first duplicated mutation (cant be identified as unique, since once duplicates are removed we can no longer tell ambiguity)
                mutations_df[!,"Ambiguous"] = vcat([false],[mutations_df[i,"Duplicate"] == false  && mutations_df[i+1,"Duplicate"] == true ? false : true  for i in 1:nrows(mutations_df)-1])
                filter!([:Duplicate] => x -> x == false, mutations_df)
                #Remove mutations that are present only in one parent (These mutations are unlikely to be real i.e more parsimonious to assume inheritance from other parental lineage)
                #These mutations will have different REF nucleotides, but the same ALT, meaning they can be uniquely identified since only one parent will have the mutation (the other wont as it has been passed on)
                mutations_df[!,"Duplicate"] = nonunique(mutations_df,["REF","POS",])
                filter!([:Duplicate] => x -> x == true, mutations_df)
                #Remove Duplicate Filter Column
                select!(mutations_df, Not(:Duplicate))
            else
                println(get_closest_ancestor(lineage,pseudoreference_dict))
                reference = string(pseudoreference_dict[get_closest_ancestor(lineage,pseudoreference_dict)]) 
                println(get_closest_ancestor(lineage,pseudoreference_dict),lineage)
                 if (length(record_sequence) != length(reference) )
                    println("Sequences are not aligned!")
                    println(lineage)
                    println(length(record_sequence))
                    println(length(reference))
                    continue
                end
                #Get the mutations for the sequence 
                mutations = get_mutations(reference,record_sequence,context,orf_at_position,orf_range_dict,record_sequence_name)
                mutations_df = DataFrame(Sequence_Name = mutations["Sequence_Name"],
                                        REF = mutations["REF"],
                                        ALT = mutations["ALT"],
                                        POS = Nucleotide = mutations["POS"],
                                        BIO_POS = Nucleotide = mutations["BIO_POS"],
                                        AMINO_ORF_POS = Nucleotide = mutations["AMINO_ORF_POS"],
                                        Nucleotide = mutations["Nucleotide"],
                                        Nucleotide_Context = mutations["Nucleotide_Context"],
                                        Amino_Acid = mutations["Amino_Acid"],
                                        Amino_Acid_Context = mutations["Amino_Acid_Context"],
                                        Codon_Position = mutations["Codon_Position"],
                                        Codon = mutations["Codon"],
                                        ORF = mutations["ORF"],
                                        Synonymous = mutations["Synonymous"])
                #No Ambiguous mutations since reference known
                mutations_df[!, :Ambiguous] .= false
            end
        end
        
        if mutation_set == ""
            mutation_set = mutations_df
        else
            append!(mutation_set,mutations_df)
        end

        if signature_matrix == true
            new_row = copy(empty_row)
            for substitution_category in mutations["Nucleotide_Context"]
                if haskey(new_row,substitution_category)
                    new_row[substitution_category] = new_row[substitution_category]+1 
                end
            end
            new_row["Sequence_Name"] = record_sequence_name
            append!(substitution_matrix,DataFrame(new_row))
        end

        if index%50000 == 0
            println("Sequence "*string(index))
            if index == 50000
                write_type = "w"
                append_type = false
            else
                write_type = "a"
                append_type = true
            end
            println(mutation_set)

            CSV.write(output*"mutations.csv", mutation_set,append=append_type)

            if signature_matrix == true
                substitution_matrix = select!(substitution_matrix,:Sequence_Name, Not([:Sequence_Name]))
                CSV.write(output*"signature_matrix.csv", substitution_matrix,append=append_type)
            end
            mutation_set = ""
            substitution_matrix =make_substitution_matrix(context)
            last_index=index
        end 
    end

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Stage 2: Write to file any remainder mutations
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if last_index == 1
        write_type = "w"
        append_type = false
    else
        write_type = "a"
        append_type = true
    end
    CSV.write(output*"mutations.csv", mutation_set,append=append_type)
    
    if signature_matrix == true
        substitution_matrix = select!(substitution_matrix,:Sequence_Name, Not([:Sequence_Name]))
        CSV.write(output*"signature_matrix.csv", substitution_matrix,append=append_type)
    end

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#    Stage 3: Remove Duplicates and Filter Mutations 
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#     De-duplicate mutations and keep oldest known de-novo mutations in each lineage
#     Read in file containing all the mutations we have saved
#     df = DataFrame(CSV.File(output*"mutations.csv"))
    df = CSV.read(output*"mutations.csv", DataFrame)
    df[!,:Sequence_Header] = df[!,:Sequence_Name]
    
    if ("Sample_Name" in header_format)

        #Split Sequence Name and Expand Columns
        transform!(df, :Sequence_Name => ByRow(x -> split(x, '|')) => header_format)

        #Remove Unique (sort dates in order before removing duplicate mutations)
        sort!(df,[:Sample_Date])
        df[!,"Duplicate"] = nonunique(df,["REF","ALT","POS","Lineage"])
        filter!([:Duplicate] => x -> x == false, df)

        #Remove Duplicate Filter Column
        select!(df, Not(:Duplicate))

        #Write Normalised mutations
        CSV.write(output*"normalised_mutations.csv", df)


        normalised_substitution_dict = Dict{String,Dict}()
        normalised_substitution_matrix = make_substitution_matrix(context)

        #Write normalised substitution matrix
        for (index,row) in enumerate(eachrow(df))
            if haskey(normalised_substitution_dict,row.Sequence_Header) == false
                dict_row = Dict{String,Any}(i => 0 for i in names(normalised_substitution_matrix))
                dict_row["Sequence_Name"] = row.Sequence_Header
                normalised_substitution_dict[row.Sequence_Header] = dict_row
            end
            sub = row.Nucleotide_Context
            if sub in(names(normalised_substitution_matrix))
                normalised_substitution_dict[row.Sequence_Header][sub] = normalised_substitution_dict[row.Sequence_Header][sub]+=1
            end
        end

        for key in collect(keys(normalised_substitution_dict))
            row = DataFrame(normalised_substitution_dict[key])
            append!(normalised_substitution_matrix,row)
        end

        #Write Normalised mutations
        CSV.write(output*"normalised_signature_matrix.csv", normalised_substitution_matrix)
    end
end


main()

