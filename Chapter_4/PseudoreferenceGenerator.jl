import FASTX
import BioSequences
using  GenomicAnnotations
import BioGenerics
import JSON
using ArgParse
using DataFrames
using CSV
using DelimitedFiles

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
        "--alias_file", "-a"
            help = "The lineage alias file."
            arg_type = String
            default = ""
        "--output_folder","-o"
            help = "The location for the mutations file "
            arg_type = String
            default = ""
    end
    return parse_args(settings)
end;


function get_mutations(ref_seq,record_seq)
    frequency_table = DataFrame(A = zeros(length(ref_seq)), C = zeros(length(ref_seq)), G = zeros(length(ref_seq)), T = zeros(length(ref_seq)), - = zeros(length(ref_seq)))
    valid_bases = Set(["A","C","G","T","-"])
    
    ref_seq = string(ref_seq)
    record_seq = string(record_seq)
    
    for (index,ref_base) in enumerate(ref_seq)
        record_base = string(record_seq[index])
        if(record_base in valid_bases)
            frequency_table[index,record_base] = 1
        end
    end
    return frequency_table
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

function is_recombinant_lineage(lineage,lineage_aliases)
    if haskey(lineage_aliases, split(lineage,".")[1]) 
        alias = lineage_aliases[split(lineage,".")[1]]
        if alias isa String == false
            return true
        end
    end
    return false
end;

function get_wild_type_base(row,index,Frequency_Dataframe,lineage_aliases,main_reference,pseudoreference_dict,summed_rows,Lineage)
    
    #Recombinant lineages will have a list of lineages, so must detect this before converting lineage names
    is_recombinant = is_recombinant_lineage(Lineage,lineage_aliases)
    #Unknown wild type base
    wild_type_base = ""
    sum_of_bases = summed_rows[index]
    base_frequency = 0 
    base = ""
    for base_name in names(Frequency_Dataframe)
        if row[base_name]> base_frequency
            base_frequency = row[base_name]
            base = base_name
        end
    end

    #If base frequencey is greater than 0.75, assign it as the wild type
    if sum_of_bases!= 0 && base_frequency >=0.75
        wild_type_base = base
    #If not, we need to guess the base using parental lineages
    else
        #Assume sequence parent is the reference sequence
        parental_lineages = []
        #Check if sequence is a recombinant lineage (i.e has multiple parents)
        if is_recombinant
            parental_lineages = lineage_aliases[Lineage]
            for (parent_index,parent) in enumerate(parental_lineages)
                parental_lineages[parent_index] = convert_lineage_name(parent,lineage_aliases)
            end
        #Check if sequence lineage is a direct descendant from reference sequence
        elseif (occursin(".",convert_lineage_name(Lineage,lineage_aliases)))
            wild_type_base = string(main_reference.sequence[index])
        #If sequence is not recombinant or directly descendant from the reference, get its full name
        else
            parental_lineages = [convert_lineage_name(Lineage,lineage_aliases)]
        end

        #Use parent lineages to find appropriate wild type
        if length(parental_lineages) != 0
            #Use the closest availible wild type base(If a parent lineage has to go all the way back to the reference,
            #while another only has to go back one step, assume the closest traversal is better)
            wild_types = []
            wild_type_heights = []
            for (parent_index,parent) in enumerate(parental_lineages)
                max_height = count(i->(i=='.'),parent)
                for parent_height in 1:max_height
                    current_parent = join(split(parent,".")[1:end-parent_height],".")
                    if haskey(pseudoreference_dict,current_parent)
                        append!(wild_types,pseudoreference_dict[current_parent][index])
                        append!(wild_type_heights,parent_height)
                    end
                end
            end
            if length(wild_types) == 0
                wild_type_base = main_reference.sequence[index]
            else
                wild_type_height_ordering = sort!(DataFrame(Parent_Base = wild_types, Height = wild_type_heights),[:Height],rev = false)
                wild_type_base = string(wild_type_height_ordering[1,Parent_Base])
            end
        end 
    end
    return wild_type_base
end

function find_lineage_height(lineage_name,lineage_aliases)
    height = 0
    #If the lineage is recombinant, find the longest known parent
    if is_recombinant_lineage(lineage_name,lineage_aliases)
        #Set max parental lineage name and height (Starts at 0)
        longest_parental_lineage = ""
        longest_lineage_height = count(i->(i=='.'), longest_parental_lineage)

        #Get recombinant parental lineages
        parental_lineages = lineage_aliases[lineage_name]

        #Find max height of recombinant parents
        for (parent_index,parent) in enumerate(parental_lineages)
            parental_lineages[parent_index] = convert_lineage_name(parent,lineage_aliases)
            lineage_height = count(i->(i=='.'), parental_lineages[parent_index]) + count(i->(i=='.'),lineage_name)
            if lineage_height > longest_lineage_height
                longest_parental_lineage = parental_lineages[parent_index]
            end
        end 
        height = longest_lineage_height 
    #If not recombinant, 
    else
        height = count(i->(i=='.'), convert_lineage_name(lineage_name,lineage_aliases))
    end
    return height
end



function main()
    
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #Pre-processing Stage: Parse Arguments and load files
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    parsed_args = parse_commandline()
    fasta_file = ""
    reference_file = ""
    output_folder = ""
    alias_file = ""
    
    for (arg,val) in parsed_args
        if arg == "fasta_file"
            fasta_file = val
        elseif arg == "reference_file"
            reference_file = val
        elseif arg == "alias_file"
            alias_file = val
        else
            output_folder = val
        end
    end
    
#     header_format = ["Sequence_Name","GISAID_ID","Host","Clade","Alias","Lineage","Length","Sample_Date","Epi_Week","Country","FullCountry","Error"]
#     ["Sequence_Name","GISAID_ID","Host","Clade","Alias","Lineage","Length","Sample_Date","Epi_Week","Country","FullCountry","Error"]
#     USA/CO-CDC-MMB08755899/2021|EPI_ISL_2646617|HUMAN|GK|AY.44|B.1.617.2.44|29787|2021-05-11|72|USA|NORTH_AMERICA/USA/COLORADO| 
    # header_format=["GISAID_ID","Sample_Name","Lineage","Sample_Date","VOC","Alias"]
    header_format=["Sample_Name","Lineage","Alias"]
    
    main_reference = GenomicAnnotations.readgbk(reference_file)[1]
    fasta = open(FASTX.FASTA.Reader,fasta_file)
    
    #Read lineage aliases
    if alias_file != ""
        lineage_aliases = JSON.parse(open(alias_file))
        filter!(p->(p[2] !=""),lineage_aliases)
    end
    
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #Stage 1: Get Base frequencies for all lineages
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Lineage_Base_Frequency_Dict = Dict{String,DataFrame}()
    for (index,record) in enumerate(fasta)
        if (FASTX.FASTA.hasdescription(record) == false)
            record_sequence_name = BioGenerics.seqname(record)
        else
            #Sequence names sometimes have spaces 
            record_sequence_name = string(BioGenerics.seqname(record))*" "*string(FASTX.FASTA.description(record))
            record_sequence_name = replace(record_sequence_name," "=>"_")
        end
        

        #Current Fasta Sequence
        record_sequence = BioGenerics.sequence(record)
        println("Sequence name is "*record_sequence_name)
        if "Alias" in header_format
            record_lineage = split(record_sequence_name, "|")[findall( x -> x == "Alias", header_format)[1]]
        elseif "Lineage" in header_format
            record_lineage = split(record_sequence_name, "|")[findall( x -> x == "Lineage", header_format)[1]]
        end
        
        if record_lineage == ""
            continue
        end
 
        #Check if lineage has been discovered already
        if (length(findall( x -> x == record_lineage, collect(keys(Lineage_Base_Frequency_Dict)))) == 0)
            println(record_lineage)
            Lineage_Base_Frequency_Dict[record_lineage] = DataFrame(A = zeros(length(main_reference.sequence)), C = zeros(length(main_reference.sequence)), G = zeros(length(main_reference.sequence)), T = zeros(length(main_reference.sequence)), - = zeros(length(main_reference.sequence))) 
        end
        
        #Update base frequencies for lineage
        Lineage_Base_Frequency_Dict[record_lineage] = Lineage_Base_Frequency_Dict[record_lineage] .+ get_mutations(main_reference.sequence,record_sequence)
    end
    println("Calculated Lineage Base Frequencies")
    
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #Stage 2: Build Hierarchy List of lineages (We need this so that lineages earlier in the tree are built before later ones)
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    all_lineages_list = collect(keys(Lineage_Base_Frequency_Dict))
    lineage_heights = []
    for lineage_name in all_lineages_list
        append!(lineage_heights,find_lineage_height(lineage_name,lineage_aliases))
    end
    height_ordering = sort!(DataFrame(Lineage = all_lineages_list, Height = lineage_heights),[:Height],rev = false) 
    println("Ordered Lineages by Height")
    
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #Stage 3: Build Pseudoreference Fasta file from frequencies
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    #Empty dictionary of pseudoreferences
    pseudoreference_dict = Dict{String,String}()
    for Lineage in height_ordering[!,:Lineage]
        
        #Fill NaNs with 0
        Lineage_Base_Frequency_Dict[Lineage] .= ifelse.(isnan.(Lineage_Base_Frequency_Dict[Lineage]), 0, Lineage_Base_Frequency_Dict[Lineage])
        
        #Get Frequency dataframe
        Frequency_Dataframe = Lineage_Base_Frequency_Dict[Lineage]
        
        #Calculate Base proportions
        summed_rows = sum.(eachrow(Frequency_Dataframe))
        
        Frequency_Dataframe = Frequency_Dataframe./summed_rows
        
        #Pseudoreference string for sequence
        pseudoreference_sequence = ""
        
        #Build reference loop
        for (index,row) in enumerate(eachrow(Frequency_Dataframe))
            pseudoreference_sequence = pseudoreference_sequence*string(get_wild_type_base(row,index,Frequency_Dataframe,lineage_aliases,main_reference,pseudoreference_dict,summed_rows,Lineage))
        end
#         println("$Lineage - "*string(length(pseudoreference_sequence)))
        pseudoreference_dict[Lineage] = pseudoreference_sequence
    end
    println("Constructed Pseudoreferences")
    
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #Stage 4: Write results to Fasta File and CSV for frequencies
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    FASTX.FASTA.Writer(open("$output_folder/pseudo_references.fasta", "w"))  do fasta_file_writer
        for (Lineage, Pseudoreference) in pseudoreference_dict
            Lineage_Base_Frequency_Dict[Lineage][!,:POS] = 1:length(main_reference.sequence)
            Lineage_Base_Frequency_Dict[Lineage][!,:POS ] = string.(Lineage_Base_Frequency_Dict[Lineage][:,:POS ])
            Lineage_Base_Frequency_Dict[Lineage] = permutedims(Lineage_Base_Frequency_Dict[Lineage],"POS")
            Lineage_Base_Frequency_Dict[Lineage][!, :Lineage] .= Lineage
            write(fasta_file_writer, FASTX.FASTA.Record(Lineage, Pseudoreference))
        end
    end
    CSV.write("$output_folder/pseudo_references_counts.csv", reduce(vcat, values(Lineage_Base_Frequency_Dict)))
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

end

main()