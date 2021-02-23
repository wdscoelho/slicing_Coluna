############################################################################################

# Thesis title: Modeling and Optimization for 5G Network Design
# Ph.D. student: Wesley DA SILVA COELHO
# Thesis Director : Stefano SECCI (Conservatoire National des Arts et MÃ©tiers - CEDRIC)
# Thesis Supervisors : Amal BENHAMICHE and Nancy PERROT (Orange Labs)
############################################################################################

#                   INCLUDE PACKAGES"
using DataStructures
using NBInclude
using LightGraphs 
using MetaGraphs

using JuMP,MathOptInterface 

struct VNF
           id::Int8
           typee::Symbol
           treatment_capacity::Float16
           dataVolPerUE::Float16
           cpu_request::Float16
           ram_request::Float16
           storage_request::Float16
           compression::Float16 
           amount_of_traffic_sent_to_g::Vector{Float16}
           UE_dependency::Int64 
    
    
           VNF(id::Int64,typee::Symbol,treatment_capacity::Float64,dataVolPerUE::Float64,cpu_request::Float64,ram_request::Float64,storage_request::Float64, compression::Float16, amount_of_traffic_sent_to_g::Vector{Float16},UE_dependency::Int64) = new(id,typee,treatment_capacity,dataVolPerUE,cpu_request,ram_request,storage_request,compression,amount_of_traffic_sent_to_g,UE_dependency)

end

struct SliceRequest
        id::Int8 
        TotalAmountUE::Int16
        minBandwidth::Float16
        setVNFtoInstall::String
        VNFs_To_Connect::String
        VNF_sharing_file::String
        VNF_sharing::Dict
        set_VNFs_to_install::Vector{Int16}
        set_commodities::Vector{Dict}
        maxLatencyDataLayer::Float16
        nodeSharing::Vector{Int8}
        
        
        SliceRequest(id::Int8,TotalAmountUE::Int16,minBandwidth::Float16,setVNFtoInstall::String, VNFs_To_Connect::String, VNF_sharing_file::String,VNF_sharing::Dict,set_VNFs_to_install::Vector{Int16},set_commodities::Vector{Dict},maxLatencyDataLayer::Float16,nodeSharing::Vector{Int8}) = 
                 new(id,TotalAmountUE,minBandwidth,setVNFtoInstall,VNFs_To_Connect, VNF_sharing_file,VNF_sharing,set_VNFs_to_install,set_commodities,maxLatencyDataLayer,nodeSharing)
end

mutable struct Instance
        set_VNFs::Vector{VNF}
        physical_network::MetaDiGraph
        setSlices::Vector{SliceRequest}
        VNF_connection::Dict
        number_of_NFs::Int16
        number_of_AN_based_NFs::Int16
        number_of_CN_based_NFs::Int16
        maxLatencyBetweenFunctions::Array{Any,1}
        slice_cont::Array{Any,1}
        band_link_total::Float16
        maxLatencyBetweenDU_NFS1::Float16
        Instance(set_VNFs::Vector{VNF}, physical_network::MetaDiGraph, setSlices::Vector{SliceRequest}) =  new(set_VNFs,physical_network,setSlices)
end

mutable struct Parameters
    solver::String
    boost::Bool
    instance_name::String
    relaxation::Bool
    warming_up::Bool
    number_of_lazy_cuts::Int64
    objf::String
    formulation::String
    class::String
    best_sol::Float64

    number_of_phy_nodes::Int64
    graph_density::Float16
    link_capacity_ratio::String
    node_capacity_ratio::String
    number_of_slice_requests::Int64
    number_of_commodities_per_NS::Int64
    number_of_CP_NFSs::Int64
    number_of_DP_NFSs::Int64
    prob_sharing_NFS::Float16
    prob_sharing_node::Float16

    tree_size::Int64
    number_of_int_solutions::Int64
    final_lower_bound::Float16
    first_lower_bound::Float16
    final_upper_bound::Float16
    first_upper_bound::Float16
    node_last_int_sol::Int64
    final_status::MathOptInterface.TerminationStatusCode
    LP_solution::Float16
    MILP_solution::Float16
    root_gap::Float16
    
    number_of_variables_gamma::Int64
    number_of_frac_variables_gamma::Int64
    number_of_base_relax_variables_gamma::Int64
    number_of_base_final_variables_gamma::Int64
    
    number_of_variables_x::Int64
    number_of_frac_variables_x::Int64
    number_of_base_relax_variables_x::Int64
    number_of_base_final_variables_x::Int64
    number_of_variables_y::Int64
    number_of_frac_variables_y::Int64
    number_of_base_relax_variables_y::Int64
    number_of_base_final_variables_y::Int64
    number_of_variables_w::Int64
    number_of_frac_variables_w::Int64
    number_of_base_relax_variables_w::Int64
    number_of_base_final_variables_w::Int64
    number_of_variables_z::Int64
    number_of_frac_variables_z::Int64
    number_of_base_relax_variables_z::Int64
    number_of_base_final_variables_z::Int64
    
    total_MILP_time::Float16
    total_LP_time::Float16
    root_node_processing_time::Float16
    MILP_relaxation_time::Float16
    branch_and_bound_time::Float16
    vi_creation_time::Float16
    number_constraints_before_VI::Int64
    number_constraints_after_VI::Int64
    valid_inequalities::String
    simplex_iterations::Int64
    instance_type::String
    instance_number::Int64
    start_branching_upper_bound::Float16
    
    cut_class::String
    cut_time::Float16
    number_of_cuts::Int64
    lazy_time::Float64
    best_solution_value::Float64
    
    Parameters() =  new()
end     



function max_Int(a::Int64,b::Int64)
    if a< b
        return b
    else
        return a
    end
end

function get_Instance(input_folder::String,instance_number::Int64)
    
   #Here, we get all proprities needed to represent the instance
        set,dpNFS,cpNFS = get_VNFs(joinpath(input_folder,"instance_$(instance_number)_NFS_types.dat"))
        instance = Instance(set,
                            get_physical_network(joinpath(input_folder,"instance_$(instance_number)_physical_network.dat")),
                            get_NS_requests(joinpath(input_folder,"instance_$(instance_number).dat")))

        instance.VNF_connection = get_VNFs_To_Connect(instance)
        instance.number_of_AN_based_NFs = dpNFS
        instance.number_of_CN_based_NFs = cpNFS
        instance.maxLatencyBetweenFunctions = get_latency_between_functions(joinpath(input_folder,"instance_$(instance_number)_latency_btw_functions.dat"))
        instance.number_of_NFs = 0
        for s in 1:length(instance.setSlices)
            instance.number_of_NFs+=length(instance.setSlices[s].set_commodities)
        end
        for u in 1: props(instance.physical_network)[:number_nodes] 
            if props(instance.physical_network, u)[:node_type] == "non_access"
                instance.number_of_NFs += 1
            end 
        end
        
        return instance
    
end

function get_NS_requests(file::String) 
    #Auxiliar variables
    setSlices = Vector{SliceRequest}()
   # We open a text file found on the recieved path "file::String"
    open(file) do file
        
        #we read line by line
        for ln in eachline(file)

            #Here, we storage each word / number in a position of aux_string vector
            aux_string = split(ln)
            
            #error treatment: if it is a line with nothing written, we return to the first isntruction of this loop
            if length(aux_string) == 0
                continue #going dic=rectly to the first line of the loop
            end
            
            #getting number of slice_requests 
            if aux_string[1] == "number_of_requests" 
                number_of_requests = Meta.parse(aux_string[2]) # Meta.parse tries to convert strings into a most suitable type
            
            #getting slice requests and their attributes
            elseif aux_string[1]  == "slice_request" 
                VNFs_to_install_instance_name, aux_set_VNFs_to_install =  get_VNFs_to_install(Meta.parse(aux_string[10]),parse(Int8,aux_string[2]))    

                # adding its attributes
                push!(setSlices, SliceRequest(parse(Int8,aux_string[4]),parse(Int16,aux_string[6]),
                                         parse(Float16,aux_string[8]),Meta.parse(aux_string[10]),Meta.parse(aux_string[12]),
                                         Meta.parse(aux_string[14]), get_VNF_sharing(Meta.parse(aux_string[14])),
                                        aux_set_VNFs_to_install,get_commodities(Meta.parse(aux_string[16])), parse(Float16,aux_string[18]),get_node_sharing(Meta.parse(aux_string[20]),parse(Int64,aux_string[2]))))
           
                end# end of conditions   
               
        end#end for - we have read all lines from file
     
    end#closing file
    return setSlices 

end # end of function get_NS_requests

function get_node_sharing(file::String,slice_id::Int64)
    #auxiliar variables
    my_vector = Vector{Int8}()
    line_count = 0
    # We open a text file found on the recieved path "file::String"
    open(file) do file
        #we read line by line
        for ln in eachline(file)
            line_count = line_count+1
            if line_count==slice_id
                aux_string = split(ln)
                for n in 1:length(aux_string)
                    push!(my_vector, parse(Int8, aux_string[n]))
                end
            end
        end
    end
    return my_vector
end

function get_VNF_sharing(file::String) 
    
    #auxiliar variables
    number_of_slices = 0
    number_VNFs = 0
    aux_Matrix1 = Vector{ Vector{Int8} } 
    aux_Matrix = aux_Matrix1(undef,0)
    aux_dic = Dict()    
   
    # We open a text file found on the recieved path "file::String"
    open(file) do file
  
        #we read line by line
        for ln in eachline(file)

            #Here, we storage each word / number in a position of aux_string vector
            aux_string = split(ln)

            #error treatment: if it is a line with nothing written or with ";", we return to the first isntruction of this loop
            if length(aux_string) == 0
                continue #going directly to the first line of the loop
            end
            
            if aux_string[1]  == ";" 
                    aux_dic["$(length(aux_dic)+1)"] = aux_Matrix
                    aux_Matrix = aux_Matrix1(undef,0)
                    continue #going directly to the first line of the loop
            end

            #getting number of VNFs and Slices 
            if aux_string[1] == "number_of_VNFs" 
                number_VNFs = parse(Int64,aux_string[2]) # Meta.parse tries to convert strings into a most suitable type

            elseif aux_string[1] == "number_of_slices" 
                number_of_slices = parse(Int8,aux_string[2]) # Meta.parse tries to convert strings into a most suitable type

            else
                #creating Sharing Matrix
                aux_vector = Array{Bool,1}(undef, number_VNFs)          

                for (n, f) in enumerate(aux_string)
                    aux_vector[n] = parse(Int8,f)
                end  
                 push!(aux_Matrix,aux_vector)       

            end# end of conditions   

        end#end for - we have read all lines from file   

    end#closing file
    
    return aux_dic 

end # end of function get_VNF_sharing

function get_physical_network(file::String) 

    physical_network = MetaDiGraph() # It is this graph that we return in the end of this function. We use "Meta" to have several attributes on each node and each edge.
    set_prop!(physical_network, :type, "all")
   # We open a text file found on the recieved path "file::String"
    open(file) do file
        
        #we read line by line
        for ln in eachline(file)
            
            #Here, we storage each word / number in a position of aux_string vector
            aux_string = split(ln)
            
            #error treatment: if it is a line with nothing written, we return to the first instruction of this loop
            if length(aux_string) == 0
                continue #going directly to the first line of the loop
            end
            
            #getting number of nodes, set it as a attributes of physical_network and add this number of nodes to the graph
            if aux_string[1] == "number_of_nodes" 
                MetaGraphs.set_prop!(physical_network, :number_nodes, Meta.parse(aux_string[2])) # Meta.parse tries to convert strings into a most suitable type
                MetaGraphs.add_vertices!(physical_network, Meta.parse(aux_string[2]))
                
            #getting number of edges set it as a attributes of physical_network graph    
            elseif aux_string[1]  == "number_of_arcs" 
                MetaGraphs.set_prop!(physical_network, :number_of_arcs, Meta.parse(aux_string[2]))
            
            #getting nodes and their attributes
            elseif aux_string[1]  == "node" 
                
                #we add its attributes
                MetaGraphs.set_props!(physical_network, Meta.parse(aux_string[2]), Dict(:node_id =>Meta.parse(aux_string[4]), :node_type =>aux_string[6] , :longitude =>Meta.parse(aux_string[8]), :latitude =>Meta.parse(aux_string[10]), :coverage_radius =>Meta.parse(aux_string[12]), :ram_capacity =>Meta.parse(aux_string[14]), :cpu_capacity =>Meta.parse(aux_string[16]), :storage_capacity =>Meta.parse(aux_string[18]), :networking_capacity =>Meta.parse(aux_string[20]), :ram_cost =>Meta.parse(aux_string[22]), :cpu_cost =>Meta.parse(aux_string[24]), :storage_cost =>Meta.parse(aux_string[26]), :networking_cost =>Meta.parse(aux_string[28]), :ue_capacity =>Meta.parse(aux_string[30]), :availability =>Meta.parse(aux_string[32]), :internal_delay =>Meta.parse(aux_string[34]), :cost_node =>Meta.parse(aux_string[36]) ))
           
            
            #getting edges and their attributes. Here, we can add an edge and its attributes with the same function
            elseif aux_string[1]  == "arc"            
                MetaGraphs.add_edge!(physical_network, Meta.parse(aux_string[6]), Meta.parse(aux_string[8]), Dict(:edge_id=>Meta.parse(aux_string[4]), :source =>Meta.parse(aux_string[6]), :target =>Meta.parse(aux_string[8]), :type =>aux_string[10], :delay =>Meta.parse(aux_string[12]), :max_bandwidth =>Meta.parse(aux_string[14]), :availability=>Meta.parse(aux_string[16]), :cost =>Meta.parse(aux_string[18])))
           
            elseif aux_string[1]  == "node_capacity_types"
                aux_vec = Vector{Symbol}()                
                MetaGraphs.set_prop!(physical_network, :number_node_capacity_types, parse(Int8,aux_string[2]))              
                for i in 1:parse(Int8,aux_string[2])
                    push!(aux_vec, Meta.parse(aux_string[i+2]))
                end
                MetaGraphs.set_prop!(physical_network, :node_capacity_types, aux_vec)              
                
         end# end of conditions   
                
        end#end for - we have read all lines from file
     
        end#closing file
    
    return physical_network 

end # end of function get_physical_network


function get_VNFs(file::String) 

    #Auxiliar variables
    number_VNFs = 0
    number_of_AN_based_NFs = 0
    number_of_CN_based_NFs = 0
    set_of_VNFs = Vector{VNF}()
    
   # We open a text file found on the recieved path "file::String"
    open(file) do file
        
        #we read line by line
        for ln in eachline(file)

            #Here, we storage each word / number in a position of aux_string vector
            aux_string = split(ln)
            
            #error treatment: if it is a line with nothing written, we return to the first isntruction of this loop
            if length(aux_string) == 0
                continue #going dic=rectly to the first line of the loop
            end
            
            #getting number of VNFs 
            if aux_string[1] == "number_of_VNFs" 
                number_VNFs = Meta.parse(aux_string[2]) # Meta.parse tries to convert strings into a most suitable type
            
            elseif aux_string[1]  == "number_of_DP_NFs" 
                number_of_AN_based_NFs =  parse(Int64,aux_string[2])
            
            elseif aux_string[1]  == "number_of_CP_NFs" 
                number_of_CN_based_NFs = parse(Int64,aux_string[2])
                        
            #getting VNFs and their attributes
            elseif aux_string[1]  == "NFS" 
                
                #we add its attributes
                t = Meta.parse(aux_string[6])
                aux_vector = Array{Float16,1}(undef, number_VNFs)    
                for i in 1:number_VNFs 
                    aux_vector[i] = Meta.parse(aux_string[21+i])
                end
                push!(set_of_VNFs, VNF(Meta.parse(aux_string[4]),Meta.parse(aux_string[6]),Meta.parse(aux_string[8]),Meta.parse(aux_string[10]),
                                         Meta.parse(aux_string[12]),Meta.parse(aux_string[14]),Meta.parse(aux_string[16]), parse(Float16,aux_string[18]),aux_vector,parse(Int64,aux_string[20])))
            end# end of conditions   
                
        end#end for - we have read all lines from file
     
    end#closing file
    
    return set_of_VNFs,number_of_AN_based_NFs,number_of_CN_based_NFs 

end # end of function get_VNFs

function get_latency_between_functions(file::String)
    
    number_of_functions = 0
    aux_Matrix1 = Vector{ Vector{Float16} } 
    #aux_Matrix_return =  aux_Matrix01(undef,0)
    aux_Matrix_return = []
    # We open a text file found on the recieved path "file::String"
   
    open(file) do file
  
        #we read line by line
        for ln in eachline(file)

            #Here, we storage each word / number in a position of aux_string vector
            aux_string = split(ln)

            #error treatment: if it is a line with nothing written or with ";", we return to the first isntruction of this loop
            if length(aux_string) == 0
                continue #going directly to the first line of the loop
            end
            
            #getting number of VNFs 
            if aux_string[1] == "number_of_NFSs" 
                number_of_functions = parse(Int64,aux_string[2]) # Meta.parse tries to convert strings into a most suitable type

            else
                #we add its attributes
                aux_vector = Array{Float16,1}(undef, number_of_functions)          

                for (n, f) in enumerate(aux_string)
                    aux_vector[n] = parse(Float16,f)
                end  
                 push!(aux_Matrix_return,aux_vector)       

            end# end of conditions   

        end#end for - we have read all lines from file   

    end#closing file
    
    return aux_Matrix_return 

end

function get_VNFs_to_install(file::String,slice_id::Int8) 
    set_VNF_to_install_instance = ""
    number_of_VNFs_to_install = 0
    line_count = 0
    # We open a text file found on the recieved path "file::String"
    open(file) do file
        #we read line by line
        for ln in eachline(file)
                line_count = line_count+1
                #Here, we storage each word / number in a position of aux_string vector
                aux_string = split(ln)

                #error treatment: if it is a line with nothing written or with ";", we return to the first isntruction of this loop
                if length(aux_string) == 0
                    continue #going directly to the first line of the loop
                end

                #getting micro areas
                if aux_string[1] == "set_VNF_to_install" 
                    set_VNF_to_install_instance = Meta.parse(aux_string[2]) # Meta.parse tries to convert strings into a most suitable type

                elseif aux_string[1] == "number_of_VNFs_to_install" 
                    number_of_VNFs_to_install = parse(Int16,aux_string[2]) # Meta.parse tries to convert strings into a most suitable type
                         
                elseif line_count - 3 == slice_id
                    #we add its attributes
                    aux_vector = Array{Int16,1}(undef, number_of_VNFs_to_install)          

                    for (n, f) in enumerate(aux_string)
                        aux_vector[n] = parse(Int16,f)
                    end  

                    return set_VNF_to_install_instance , aux_vector 

                end# end of conditions   

        end#end for - we have read all lines from file   

    end#closing file
    

end # end of function get_VNF_sharing

function get_VNFs_To_Connect(instance::Instance)
    
    #auxiliar variables
    number_VNFs = 0
    aux_Matrix1 = Vector{ Vector{Bool} } 
    aux_Matrix = aux_Matrix1(undef,0)
    aux_dic = Dict()    
    
    # We open a text file found on the recieved path "file::String"
    for each_slice in instance.setSlices
       file =  each_slice.VNFs_To_Connect

        open(file) do file
  
            #we read line by line
            for ln in eachline(file)

                #Here, we storage each word / number in a position of aux_string vector
                aux_string = split(ln)

                #error treatment: if it is a line with nothing written or with ";", we return to the first isntruction of this loop
                if length(aux_string) == 0
                    continue #going directly to the first line of the loop
                end

                #getting number of VNFs 
                if aux_string[1] == "number_of_NFSs" 
                    number_VNFs = parse(Int64,aux_string[2]) # Meta.parse tries to convert strings into a most suitable type

                else
                    #creating Conection Matrix
                    aux_vector = Array{Bool,1}(undef, number_VNFs)          

                    for (n, f) in enumerate(aux_string)
                        aux_vector[n] = parse(Bool,f)
                    end  
                     push!(aux_Matrix,aux_vector)  

                end# end of conditions   

            end#end for - we have read all lines from file   
             
            aux_dic["$(each_slice.id)"] = aux_Matrix
            aux_Matrix = aux_Matrix1(undef,0) 
       
        end#closing file
    end#end for each slice
   
    return aux_dic 
end#end of function


function get_commodities(file::String) 
    
    #auxiliar variables
    aux_vector = Vector{Dict}()
    number_of_commodities = Int16
   
   # We open a text file found on the recieved path "file::String"
    open(file) do file
        
        #we read line by line
        for ln in eachline(file)

            #Here, we storage each word / number in a position of aux_string vector
            aux_string = split(ln)
            
            #error treatment: if it is a line with nothing written, we return to the first isntruction of this loop
            if length(aux_string) == 0
                continue #going dic=rectly to the first line of the loop
            end
            
            #getting number of commodities 
            if aux_string[1] == "number_of_commodities" 
                number_of_commodities = Meta.parse(aux_string[2]) # Meta.parse tries to convert strings into a most suitable type
           
            #getting commodities and their attributes
            elseif aux_string[1]  == "commodity_id" 
                
                #we add its attributes
                aux_dict = Dict()
                aux_dict["commodity_id"] = parse(Int64, aux_string[2])
                aux_dict["origin_node"] = parse(Int64, aux_string[4])
                aux_dict["target_node"] = parse(Int64, aux_string[6])
                aux_dict["volume_of_data"] = parse(Float16, aux_string[8])
                
                push!(aux_vector, aux_dict)
            end# end of conditions   
        end#end for - we have read all lines from file
     
    end#closing file
    
    return aux_vector 

end # end of function get_VNF_sharing



function create_my_model(instance::Instance)
    my_model =Model(CPLEX.Optimizer)
    S = 1:length(instance.setSlices)
    N= 1:instance.number_of_NFs
    SS =  1:length(instance.setSlices)
 #-------------------------Variables--------------------
    @variable(my_model, gamma[s in S, k in 1:length(instance.setSlices[s].set_commodities), a in 1:props(instance.physical_network)[:number_of_arcs],f in 1:length(instance.set_VNFs)+1, g in 1:length(instance.set_VNFs)+1], Bin)
    @variable(my_model, 0<=z[s in S,f in 1:instance.number_of_AN_based_NFs]<=1)
    @variable(my_model, x[s in S,f in 1:length(instance.set_VNFs), m in N, u in 1: props(instance.physical_network)[:number_nodes]], Bin)
    @variable(my_model, w[s in S,f in 1:length(instance.set_VNFs), m in N, u in 1: props(instance.physical_network)[:number_nodes]] >= 0)
    @variable(my_model, y[f in 1:length(instance.set_VNFs), m in N, u in 1: props(instance.physical_network)[:number_nodes]] >=0, Int)

    
    
    #---------------------------------------------
    # Lower Bound
    #---------------------------------------------)
    dim_f = []
     for f in 1:length(instance.set_VNFs)
        dim = 0.0 
        for  s in S
            if f!=1 && instance.set_VNFs[f].typee != :cp 
                dim += sum(instance.setSlices[s].set_VNFs_to_install[f]*instance.set_VNFs[f-1].compression*k["volume_of_data"]  for k in instance.setSlices[s].set_commodities)/instance.set_VNFs[f].treatment_capacity                
            elseif f== 1  && instance.set_VNFs[f].typee != :cp 
                dim += sum(instance.setSlices[s].set_VNFs_to_install[f]*k["volume_of_data"] for  k in instance.setSlices[s].set_commodities)/instance.set_VNFs[f].treatment_capacity                
            else
                dim += instance.setSlices[s].set_VNFs_to_install[f]*instance.set_VNFs[f].dataVolPerUE*instance.setSlices[s].TotalAmountUE/instance.set_VNFs[f].treatment_capacity
            end
        end
        push!(dim_f,ceil(dim))    
    end
                         
     LB =  sum(dim_f[i] for i in 1:length(dim_f))
     my_Ys = zero(AffExpr)
     add_to_expression!(my_Ys, sum(y[f,m,u] for m in N ,f in 1:length(instance.set_VNFs),   u in 1:props(instance.physical_network)[:number_nodes]))

     @constraint(my_model, lb, my_Ys >= LB)          
    #---------------------------------------------
    # auxiliar vectors
    #---------------------------------------------
    
   non_access = Vector{Int64}()
   for  u in 1: props(instance.physical_network)[:number_nodes]
       if props(instance.physical_network, u)[:node_type] == "non_access"
          push!(non_access, u)
        end
    end
    O_s = Vector{Vector{Int64}}() 
    non_O_s =Vector{Vector{Int64}}() 
    for s in 1:length(instance.setSlices)
        O = Vector{Int64}() 
        non_O = Vector{Int64}()
        for k in instance.setSlices[s].set_commodities
             push!(O,k["origin_node"])                            
        end
        push!(O_s,O)
       for u in 1: props(instance.physical_network)[:number_nodes]
            if  props(instance.physical_network, u)[:node_type] == "access" && (u in O_s[s]) == false
                push!(non_O,u)
            end
        end
        push!(non_O_s,non_O)
    end
    #---------------------------------------------
    # SPLIT SELECTION
    #---------------------------------------------)
    @constraint(my_model, split[s in S, f in 1:instance.number_of_AN_based_NFs-1], z[s,f]<=z[s,f+1])
    #---------------------------------------------
    # DIMENSIONING - CP NFSs
    #---------------------------------------------
    @constraint(my_model, dim_cp[s in S, f in (instance.number_of_AN_based_NFs+1):(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),n in  N,u in non_access], w[s,f,n,u]  == instance.set_VNFs[f].dataVolPerUE*instance.setSlices[s].TotalAmountUE*x[s,f,n,u] /instance.set_VNFs[f].treatment_capacity)
 
    #---------------------------------------------
    # DIMENSIONING - DISTRIBUTED DP NFSs
    #---------------------------------------------
    @constraint(my_model,dim_dis_dp[n in N,s in S,k in 1:length(instance.setSlices[s].set_commodities), f in 2:instance.number_of_AN_based_NFs], w[s,f,n,instance.setSlices[s].set_commodities[k]["origin_node"]]  == instance.set_VNFs[f-1].compression*instance.setSlices[s].set_commodities[k]["volume_of_data"]*x[s,f,n,instance.setSlices[s].set_commodities[k]["origin_node"]]/instance.set_VNFs[f].treatment_capacity)
    @constraint(my_model,dimen_dis_dp[n in N,s in S,k in 1:length(instance.setSlices[s].set_commodities)], w[s,1,n,instance.setSlices[s].set_commodities[k]["origin_node"]] == instance.setSlices[s].set_commodities[k]["volume_of_data"]*x[s,1,n,instance.setSlices[s].set_commodities[k]["origin_node"]]/instance.set_VNFs[1].treatment_capacity)

    #---------------------------------------------
    # DIMENSIONING - CENTRALIZED DP NFSs
    #---------------------------------------------
    
    @constraint(my_model, dim_cen_dp[n in N,s in S,  f in 2:instance.number_of_AN_based_NFs, u in non_access],w[s,f,n,u] == instance.setSlices[s].set_VNFs_to_install[f]*(sum(instance.set_VNFs[f-1].compression*instance.setSlices[s].set_commodities[k]["volume_of_data"] for k in 1:length(instance.setSlices[s].set_commodities))*x[s,f,n,u])/instance.set_VNFs[f].treatment_capacity)                
    @constraint(my_model, dim_cenn_dp[n in N,s in S, u in non_access], w[s,1,n,u] ==  instance.setSlices[s].set_VNFs_to_install[1]*(sum(instance.setSlices[s].set_commodities[k]["volume_of_data"] for k in 1:length(instance.setSlices[s].set_commodities))*x[s,1,n,u])/instance.set_VNFs[1].treatment_capacity)                
 
                                           
    #---------------------------------------------
    # PLACEMENT - DISTRIBUTED DP NFSs 
    #---------------------------------------------

        
    @constraint(my_model, plac_dp_in[s in S,f in 1:instance.number_of_AN_based_NFs, u in O_s[s]],sum(x[s,f,m,u] for m in N)== 1 - z[s,f])
    @constraint(my_model, plac_dp_out[s in S,f in 1:instance.number_of_AN_based_NFs],sum(x[s,f,m,u] for m in N,u in non_O_s[s])== 0)
    O_s = nothing
    non_O_s = nothing
            

    #---------------------------------------------
    # PLACEMENT - CENTRALIZED DP NFSs     
    #---------------------------------------------
    Matrix = Vector{Vector{AffExpr}}()     
    for s in S 
        my_s = Vector{AffExpr}()
        push!(Matrix,my_s)
        for f in 1:instance.number_of_AN_based_NFs
            expression = zero(AffExpr)
            push!(Matrix[s], expression)
            for u in non_access
                add_to_expression!(Matrix[s][f], sum(x[s,f,m,u] for m in N)) 
            end
        end
    end
    @constraint(my_model, plac_cent_dp[s in S,  f in 1:instance.number_of_AN_based_NFs], Matrix[s][f] == z[s,f])  
    Matrix = nothing                           
                
    #---------------------------------------------
    # PLACEMENT - CP NFSs
    #---------------------------------------------

        
    Matrix = Vector{Vector{AffExpr}}() 
    for s in S
        my_s = Vector{AffExpr}()
        push!(Matrix,my_s)
        my_s = nothing
        for f in 1:(instance.number_of_AN_based_NFs + instance.number_of_CN_based_NFs)
            expression = zero(AffExpr)
            push!(Matrix[s], expression)
            if f > instance.number_of_AN_based_NFs
                for u in non_access
                    add_to_expression!(Matrix[s][f], sum(x[s,f,m,u] for m in N)) 
                end
            end
 
        end
    end
        
    @constraint(my_model, plac_cp_1[s in S, f in (1+instance.number_of_AN_based_NFs):(instance.number_of_AN_based_NFs + instance.number_of_CN_based_NFs)],Matrix[s][f] == instance.setSlices[s].set_VNFs_to_install[f]) 
    
    Matrix = nothing
    #---------------------------------------------
    # VIRTUAL ISOLATION ================== not in Axis
    #---------------------------------------------
    no_S = Vector()
    for s in SS, t in SS
        if t != s
            push!(no_S, t)
        end
    end
    
     @constraint( my_model, virtual_iso[m in N,s in SS, t in no_S[s], u in 1:props(instance.physical_network)[:number_nodes],   f in 1:length(instance.set_VNFs),g in 1:length(instance.set_VNFs)], x[s,f,m,u] + x[t,g,m,u]  <= 1 + instance.setSlices[s].VNF_sharing["$(t)"][f][g]*instance.setSlices[t].VNF_sharing["$(s)"][g][f])  
 
    #---------------------------------------------
    # Physical ISOLATION ================== not in Axis
    #---------------------------------------------         
    @constraint(my_model, phy_iso[s in SS, t in no_S[s], f in 1:length(instance.set_VNFs), g in 1:length(instance.set_VNFs),u in non_access],sum(x[s,f,n,u] for n in N) + sum(x[t,g,m,u] for m in N)<=instance.setSlices[s].nodeSharing[t]*instance.setSlices[t].nodeSharing[s] +1 )

    #---------------------------------------------
    # PACKING  ================== not in Axis
    #--------------------------------------------
     @constraint(my_model, packing[m in N,f in 1:length(instance.set_VNFs),  u in 1: props(instance.physical_network)[:number_nodes]],sum(w[s,f,m,u] for s in S)  <= y[f,m,u])
     @constraint(my_model, one_NF[n in N, s in SS ,t in no_S[s], f in 1:length(instance.set_VNFs), g in 1:length(instance.set_VNFs), u in 1:props(instance.physical_network)[:number_nodes]-1,v in u+1: props(instance.physical_network)[:number_nodes]], x[s,f,n,u] + x[t,g,n,v] <=1)
     no_S = nothing
    #---------------------------------------------
    # NODE CAPACITY ================== not in Axis
    #---------------------------------------------
    @constraint(my_model, cap_ram[u in 1: props(instance.physical_network)[:number_nodes]], sum(instance.set_VNFs[f].ram_request*y[f,m,u] for  m in N,f in 1:length(instance.set_VNFs)) <= props(instance.physical_network,u)[:ram_capacity])
    @constraint(my_model, cap_storage[u in 1: props(instance.physical_network)[:number_nodes]], sum(instance.set_VNFs[f].storage_request*y[f,m,u] for  m in N, f in 1:length(instance.set_VNFs)) <= props(instance.physical_network,u)[:storage_capacity])
    @constraint(my_model, cap_cpu[u in 1: props(instance.physical_network)[:number_nodes]], sum(instance.set_VNFs[f].cpu_request*y[f,m,u] for f in 1:length(instance.set_VNFs), m in N) <= props(instance.physical_network,u)[:cpu_capacity])

     
    #---------------------------------------------
    #  LATENCY
    #---------------------------------------------
    
    end_to_end_delay = Vector{Vector{AffExpr}}()
    pair_delay = Vector{Vector{Vector{Vector{AffExpr}}}}()
    for s in S
        aux_vector = Vector{AffExpr}()
        push!(end_to_end_delay,aux_vector)
        aux_vector1 = Vector{Vector{Vector{AffExpr}}}()
        push!(pair_delay,aux_vector1)
        for k in 1:length(instance.setSlices[s].set_commodities)
            ex = zero(AffExpr)
           aux_vector2 = Vector{Vector{AffExpr}}()
            push!(pair_delay[s],aux_vector2)
            push!(end_to_end_delay[s],ex)
            for f in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
                aux_vector3 = Vector{AffExpr}()
                push!(pair_delay[s][k],aux_vector3)
                for g in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
                    ex1 = zero(AffExpr)
                    push!(pair_delay[s][k][f],ex1)
                    for a in edges(instance.physical_network)
                        delay_a = get_prop(instance.physical_network,a,:delay)
                        if f < instance.number_of_AN_based_NFs && g == f+1
                            add_to_expression!(end_to_end_delay[s][k], delay_a*gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g])
                        end 
                        if f == 1 && g == 1
                            add_to_expression!(end_to_end_delay[s][k], delay_a * (gamma[s,k,get_prop(instance.physical_network,a,:edge_id),length(instance.set_VNFs)+1,1] +gamma[s,k,get_prop(instance.physical_network,a,:edge_id),instance.number_of_AN_based_NFs,length(instance.set_VNFs)+1]))
                         end 
                        if instance.VNF_connection["$(s)"][f][g] == true
                            add_to_expression!(pair_delay[s][k][f][g], delay_a * gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g]) 
                        end
                    end
                add_to_expression!(pair_delay[s][k][f][g],- instance.maxLatencyBetweenFunctions[f][g])  
                end
            end
           add_to_expression!(end_to_end_delay[s][k], -  instance.setSlices[s].maxLatencyDataLayer) 
        end
    end

    @constraint(my_model, delay_pair_NFS[s in S, k in 1:length(instance.setSlices[s].set_commodities),f in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),g in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs) ], pair_delay[s][k][f][g] <= 0)  
    @constraint(my_model,end_2_end_delay[s in S, k in 1:length(instance.setSlices[s].set_commodities) ], end_to_end_delay[s][k] <= 0) 
    end_to_end_delay = nothing
    pair_delay = nothing

    
     #---------------------------------------------
    # FLOW - DP
    #---------------------------------------------
    flow_DP = Vector{Vector{Vector{Vector{AffExpr}}}}()
    flow_DP_ori = Vector{Vector{Vector{Vector{AffExpr}}}}()
    flow_DP_targ = Vector{Vector{Vector{Vector{AffExpr}}}}()
    for s in S
        a = Vector{Vector{Vector{AffExpr}}}()
        push!(flow_DP,a)
        aa = Vector{Vector{Vector{AffExpr}}}()
        push!(flow_DP_ori,a)
        aaa = Vector{Vector{Vector{AffExpr}}}()
        push!(flow_DP_targ,a)
       for k in 1:length(instance.setSlices[s].set_commodities)
            b = Vector{Vector{AffExpr}}()
            push!(flow_DP[s],b)   
            bb = Vector{Vector{AffExpr}}()
            push!(flow_DP_ori[s],b)
            bbb = Vector{Vector{AffExpr}}()
            push!(flow_DP_targ[s],b)
            for f in 1:(instance.number_of_AN_based_NFs-1)
                c = Vector{AffExpr}()
                push!(flow_DP[s][k],c)
                cc = Vector{AffExpr}()
                push!(flow_DP_ori[s][k],c)  
                ccc = Vector{AffExpr}()
                push!(flow_DP_targ[s][k],c)  
                for u in 1: props(instance.physical_network)[:number_nodes]
                    ex0 = zero(AffExpr)
                    push!(flow_DP[s][k][f],ex0) 
                    exx0 = zero(AffExpr)
                    push!(flow_DP_ori[s][k][f],exx0) 
                    exxx0 = zero(AffExpr)
                    push!(flow_DP_targ[s][k][f],exxx0) 
                    for a in edges(instance.physical_network)
                        if src(a) == u
                           add_to_expression!(flow_DP_targ[s][k][f][u], gamma[s,k,get_prop(instance.physical_network,a,:edge_id),instance.number_of_AN_based_NFs,length(instance.set_VNFs)+1])
                           add_to_expression!(flow_DP_ori[s][k][f][u], gamma[s,k,get_prop(instance.physical_network,a,:edge_id),length(instance.set_VNFs)+1,1])
                           add_to_expression!(flow_DP[s][k][f][u], gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,f+1])
                        elseif dst(a) == u
                            add_to_expression!(flow_DP_targ[s][k][f][u],- gamma[s,k,get_prop(instance.physical_network,a,:edge_id),instance.number_of_AN_based_NFs,length(instance.set_VNFs)+1])
                            add_to_expression!(flow_DP_ori[s][k][f][u],- gamma[s,k,get_prop(instance.physical_network,a,:edge_id),length(instance.set_VNFs)+1,1] )
                            add_to_expression!(flow_DP[s][k][f][u], -gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,f+1])
                        end
                    end
                    if props(instance.physical_network, u)[:node_type] != "access"
                        add_to_expression!(flow_DP[s][k][f][u],- sum(x[s,f,m,u] for m in N) +sum(x[s,f+1,m,u] for m in N))
                    elseif u!=instance.setSlices[s].set_commodities[k]["origin_node"]
                        add_to_expression!(flow_DP[s][k][f][u], 0)
                    else               
                         add_to_expression!(flow_DP[s][k][f][u],- z[s,f+1] + z[s,f])
                    end  
                    if  u==instance.setSlices[s].set_commodities[k]["target_node"]
                        add_to_expression!(flow_DP_targ[s][k][f][u], 1)
                        add_to_expression!(flow_DP_ori[s][k][f][u], 0)        

                    elseif  u==instance.setSlices[s].set_commodities[k]["origin_node"]
                        add_to_expression!(flow_DP_targ[s][k][f][u], - 1+z[s,instance.number_of_AN_based_NFs])
                        add_to_expression!(flow_DP_ori[s][k][f][u],- z[s,1])


                    elseif  props(instance.physical_network, u)[:node_type] == "non_access" 
                        add_to_expression!(flow_DP_targ[s][k][f][u],- sum(x[s,instance.number_of_AN_based_NFs,m,u] for m in N))
                        add_to_expression!(flow_DP_ori[s][k][f][u], sum(x[s,1,m,u] for m in N))

                    else 
                        add_to_expression!(flow_DP_targ[s][k][f][u], 0)
                        add_to_expression!(flow_DP_ori[s][k][f][u], 0)
                    end
                end
            end
        end
    end
    
    @constraint(my_model, dp_flow[s in S, k in 1:length(instance.setSlices[s].set_commodities), f in 1:instance.number_of_AN_based_NFs,u in 1: props(instance.physical_network)[:number_nodes]],flow_DP[s][k][f][u]==0 )
    @constraint(my_model, dp_flow_origin[s in S, k in 1:length(instance.setSlices[s].set_commodities), f in 1:instance.number_of_AN_based_NFs,u in 1: props(instance.physical_network)[:number_nodes]],flow_DP_ori[s][k][f][u]==0 )
    @constraint(my_model, dp_flow_target[s in S, k in 1:length(instance.setSlices[s].set_commodities), f in 1:instance.number_of_AN_based_NFs,u in 1: props(instance.physical_network)[:number_nodes]],flow_DP_targ[s][k][f][u]==0 )

    flow_DP = nothing
    flow_DP_ori = nothing
    flow_DP_targ = nothing
    
    #---------------------------------------------
    # FLOW - BETWEEN CP NFSs
    #---------------------------------------------
    flow_CP = Vector{Vector{Vector{Vector{AffExpr}}}}()
    for s in S
        a = Vector{Vector{Vector{AffExpr}}}()
        push!(flow_CP,a)
        for f in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
            b = Vector{Vector{AffExpr}}()
            push!(flow_CP[s],b)
            for g in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
                c = Vector{AffExpr}()
                push!(flow_CP[s][f],c)
                for u in 1: props(instance.physical_network)[:number_nodes]
                    ex = zero(AffExpr)
                    push!(flow_CP[s][f][g],ex)
                    if f> instance.number_of_AN_based_NFs&& g> instance.number_of_AN_based_NFs&& f!=g && instance.VNF_connection["$(s)"][f][g] == true && instance.setSlices[s].set_VNFs_to_install[g] == true && instance.setSlices[s].set_VNFs_to_install[f] == true
                        for a in edges(instance.physical_network)
                          if src(a) == u
                               add_to_expression!(flow_CP[s][f][g][u], gamma[s,1,get_prop(instance.physical_network,a,:edge_id),f,g])         
                            end
                            if dst(a) == u
                               add_to_expression!(flow_CP[s][f][g][u], - gamma[s,1,get_prop(instance.physical_network,a,:edge_id),f,g])          
                            end
                        end
                        add_to_expression!(flow_CP[s][f][g][u], - sum(x[s,f,m,u] for m in N) + sum(x[s,g,m,u] for m in N))
                            end #if
                        end
                    end
                end
            end# s
    @constraint(my_model, cp_flow[s in S, f in (instance.number_of_AN_based_NFs+1):(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),g in (instance.number_of_AN_based_NFs+1):(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),u in 1: props(instance.physical_network)[:number_nodes]],flow_CP[s][f][g][u]==0 )
    flow_CP = nothing    
    #---------------------------------------------
    # FLOW - BETWEEN DP AND CP NFSs
    #---------------------------------------------               
    flow_CP_DP = Vector{Vector{Vector{Vector{Vector{AffExpr}}}}}()
    for s in S
        a = Vector{Vector{Vector{Vector{AffExpr}}}}()
        push!(flow_CP_DP,a)
        for  k in 1:length(instance.setSlices[s].set_commodities)
             b = Vector{Vector{Vector{AffExpr}}}()
            push!(flow_CP_DP[s],b)
            for  f in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
                c = Vector{Vector{AffExpr}}()
                push!(flow_CP_DP[s][k],c)
                for g in 1:instance.number_of_AN_based_NFs
                    d = Vector{AffExpr}()
                    push!(flow_CP_DP[s][k][f],d)
                    for u in 1: props(instance.physical_network)[:number_nodes]
                        ex1 = zero(AffExpr)
                        push!(flow_CP_DP[s][k][f][g],ex1)
                        if f>instance.number_of_AN_based_NFs && instance.VNF_connection["$(s)"][f][g] == true && instance.setSlices[s].set_VNFs_to_install[g] == true && instance.setSlices[s].set_VNFs_to_install[f] == true
                            for a in edges(instance.physical_network)
                                if src(a) == u
                                   add_to_expression!(flow_CP_DP[s][k][f][g][u], gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g])          
                                end
                                if dst(a) == u
                                   add_to_expression!(flow_CP_DP[s][k][f][g][u], - gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g])         
                                end
                            end
                            if props(instance.physical_network, u)[:node_type] == "non_access"           
                                add_to_expression!(flow_CP_DP[s][k][f][g][u],- sum(x[s,f,m,u] for m in N) + sum(x[s,g,m,u] for m in N))

                            elseif props(instance.physical_network, u)[:node_type] == "access" && props(instance.physical_network, u)[:node_id] ==  instance.setSlices[s].set_commodities[k]["origin_node"]           
                                add_to_expression!(flow_CP_DP[s][k][f][g][u],- z[s,g]+1)
                            else
                              add_to_expression!(flow_CP_DP[s][k][f][g][u], - 0)
                            end
                        end# end if
                    end
                end
            end
        end    
    end
    @constraint(my_model, cp_dp_flow[s in S,k in 1:length(instance.setSlices[s].set_commodities), f in (instance.number_of_AN_based_NFs+1):(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),g in 1:instance.number_of_AN_based_NFs,u in 1: props(instance.physical_network)[:number_nodes]],flow_CP_DP[s][k][f][g][u]==0 )
    
    flow_CP_DP = nothing
    flow_CP_DP1 = Vector{Vector{Vector{Vector{Vector{AffExpr}}}}}()
    for s in S
        a = Vector{Vector{Vector{Vector{AffExpr}}}}()
        push!(flow_CP_DP1,a)
        for  k in 1:length(instance.setSlices[s].set_commodities)
             b = Vector{Vector{Vector{AffExpr}}}()
            push!(flow_CP_DP1[s],b)
            for  f in 1:instance.number_of_AN_based_NFs
                c = Vector{Vector{AffExpr}}()
                push!(flow_CP_DP1[s][k],c)
                for g in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
                    d = Vector{AffExpr}()
                    push!(flow_CP_DP1[s][k][f],d)
                    for u in 1: props(instance.physical_network)[:number_nodes]
                        ex1 = zero(AffExpr)
                        push!(flow_CP_DP1[s][k][f][g],ex1)
                        if g>instance.number_of_AN_based_NFs &&instance.VNF_connection["$(s)"][f][g] == true && instance.setSlices[s].set_VNFs_to_install[g] == true && instance.setSlices[s].set_VNFs_to_install[f] == true
                            for a in edges(instance.physical_network)
                                if src(a) == u
                                   add_to_expression!(flow_CP_DP1[s][k][f][g][u], gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g])          
                                end
                                if dst(a) == u
                                   add_to_expression!(flow_CP_DP1[s][k][f][g][u], - gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g])         
                                end
                            end 
                            
                            if props(instance.physical_network, u)[:node_type] == "non_access"           
                                add_to_expression!(flow_CP_DP1[s][k][f][g][u],- sum(x[s,f,m,u] for m in N) + sum(x[s,g,m,u] for m in N))

                            elseif props(instance.physical_network, u)[:node_type] == "access" && props(instance.physical_network, u)[:node_id] ==  instance.setSlices[s].set_commodities[k]["origin_node"]           
                                add_to_expression!(flow_CP_DP1[s][k][f][g][u], z[s,f]-1)
                            else
                              add_to_expression!(flow_CP_DP1[s][k][f][g][u], - 0)
                            end
                        end# end if
                    end
                end
            end
        end    
    end
    @constraint(my_model, dp_cp_flow[s in S,k in 1:length(instance.setSlices[s].set_commodities), f in 1:instance.number_of_AN_based_NFs,g in (instance.number_of_AN_based_NFs+1):(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),u in 1: props(instance.physical_network)[:number_nodes]],flow_CP_DP1[s][k][f][g][u]==0 )
    flow_CP_DP1 = nothing

    

    #---------------------------------------------
    # LINK CAPACITY ================== not in Axis                             
    #---------------------------------------------
    cap_a = Vector{AffExpr}()
    for a in edges(instance.physical_network)           
        sum_gamma = zero(AffExpr)            
        for s in 1:length(instance.setSlices), k in 1:length(instance.setSlices[s].set_commodities)
            add_to_expression!(sum_gamma,gamma[s,k,get_prop(instance.physical_network,a,:edge_id),length(instance.set_VNFs)+1,1]*instance.setSlices[s].set_commodities[k]["volume_of_data"])
            add_to_expression!(sum_gamma,gamma[s,k,get_prop(instance.physical_network,a,:edge_id),instance.number_of_AN_based_NFs,length(instance.set_VNFs)+1]*instance.setSlices[s].set_commodities[k]["volume_of_data"]*instance.set_VNFs[instance.number_of_AN_based_NFs].compression)
            for f in 1:(instance.number_of_AN_based_NFs+instance.number_of_CN_based_NFs), g in 1:(instance.number_of_AN_based_NFs+instance.number_of_CN_based_NFs)
                if f!=g 
                    if f <= instance.number_of_AN_based_NFs && g <= instance.number_of_AN_based_NFs && f==g-1                      
                        add_to_expression!(sum_gamma,gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g]*instance.setSlices[s].set_commodities[k]["volume_of_data"]*instance.set_VNFs[f].compression)
                    end                            
                    if f>instance.number_of_AN_based_NFs && g>instance.number_of_AN_based_NFs
                        add_to_expression!(sum_gamma,gamma[s,1,get_prop(instance.physical_network,a,:edge_id),f,g]*instance.set_VNFs[f].amount_of_traffic_sent_to_g[g]*(instance.setSlices[s].TotalAmountUE/length(instance.setSlices[s].set_commodities)))
                    elseif f>instance.number_of_AN_based_NFs && g<=instance.number_of_AN_based_NFs
                        add_to_expression!(sum_gamma,gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g]*instance.set_VNFs[f].amount_of_traffic_sent_to_g[g]*(instance.setSlices[s].TotalAmountUE/length(instance.setSlices[s].set_commodities)))
                    elseif f<=instance.number_of_AN_based_NFs && g>instance.number_of_AN_based_NFs
                        add_to_expression!(sum_gamma,gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g]*instance.set_VNFs[f].amount_of_traffic_sent_to_g[g]*(instance.setSlices[s].TotalAmountUE/length(instance.setSlices[s].set_commodities)))

                    end
                end

            end
        end
        add_to_expression!(sum_gamma, - get_prop(instance.physical_network,a,:max_bandwidth))
        push!(cap_a, sum_gamma)
    end
    @constraint(my_model,cap_link[a in 1:length(cap_a)], cap_a[a] <= 0)

    cap_a = nothing
    
   
    

    @objective(my_model, Min , sum((props(instance.physical_network,u)[:ram_cost]*instance.set_VNFs[f].ram_request+props(instance.physical_network,u)[:storage_cost]*instance.set_VNFs[f].storage_request+props(instance.physical_network,u)[:cpu_cost]*instance.set_VNFs[f].cpu_request)*y[f,m,u] for  m in N, f in 1:length(instance.set_VNFs), u in 1: props(instance.physical_network)[:number_nodes]) )                    

    optimize!(my_model)

    return objective_value(my_model)



end

function create_COLUNA(instance::Instance)
    coluna = nothing
    #my_model = Model(CPLEX.Optimizer
    coluna = optimizer_with_attributes(
    Coluna.Optimizer,
    "params" => Coluna.Params(
        solver = Coluna.Algorithm.TreeSearchAlgorithm(maxnumnodes =20) 
    ),
    "default_optimizer" => CPLEX.Optimizer 
)


    @axis(S , 1:length(instance.setSlices))
    my_model = BlockModel(coluna)
    N= 1:instance.number_of_NFs
    SS =  1:length(instance.setSlices)

    #-------------------------Variables--------------------
    @variable(my_model, gamma[s in S, k in 1:length(instance.setSlices[s].set_commodities), a in 1:props(instance.physical_network)[:number_of_arcs],f in 1:length(instance.set_VNFs)+1, g in 1:length(instance.set_VNFs)+1], Bin)
    @variable(my_model, 0<=z[s in S,f in 1:instance.number_of_AN_based_NFs]<=1)
    @variable(my_model, x[s in S,f in 1:length(instance.set_VNFs), m in N, u in 1: props(instance.physical_network)[:number_nodes]], Bin)
    @variable(my_model, w[s in S,f in 1:length(instance.set_VNFs), m in N, u in 1: props(instance.physical_network)[:number_nodes]] >= 0)
    @variable(my_model, y[f in 1:length(instance.set_VNFs), m in N, u in 1: props(instance.physical_network)[:number_nodes]] >=0, Int)

    
    
    #---------------------------------------------
    # Lower Bound
    #---------------------------------------------)
    dim_f = []
     for f in 1:length(instance.set_VNFs)
        dim = 0.0 
        for  s in S
            if f!=1 && instance.set_VNFs[f].typee != :cp 
                dim += sum(instance.setSlices[s].set_VNFs_to_install[f]*instance.set_VNFs[f-1].compression*k["volume_of_data"]  for k in instance.setSlices[s].set_commodities)/instance.set_VNFs[f].treatment_capacity                
            elseif f== 1  && instance.set_VNFs[f].typee != :cp 
                dim += sum(instance.setSlices[s].set_VNFs_to_install[f]*k["volume_of_data"] for  k in instance.setSlices[s].set_commodities)/instance.set_VNFs[f].treatment_capacity                
            else
                dim += instance.setSlices[s].set_VNFs_to_install[f]*instance.set_VNFs[f].dataVolPerUE*instance.setSlices[s].TotalAmountUE/instance.set_VNFs[f].treatment_capacity
            end
        end
        push!(dim_f,ceil(dim))    
    end
                         
     LB =  sum(dim_f[i] for i in 1:length(dim_f))
     my_Ys = zero(AffExpr)
     add_to_expression!(my_Ys, sum(y[f,m,u] for m in N ,f in 1:length(instance.set_VNFs),   u in 1:props(instance.physical_network)[:number_nodes]))

     @constraint(my_model, lb, my_Ys >= LB)          
    #---------------------------------------------
    # auxiliar vectors
    #---------------------------------------------
    
   non_access = Vector{Int64}()
   for  u in 1: props(instance.physical_network)[:number_nodes]
       if props(instance.physical_network, u)[:node_type] == "non_access"
          push!(non_access, u)
        end
    end
    O_s = Vector{Vector{Int64}}() 
    non_O_s =Vector{Vector{Int64}}() 
    for s in 1:length(instance.setSlices)
        O = Vector{Int64}() 
        non_O = Vector{Int64}()
        for k in instance.setSlices[s].set_commodities
             push!(O,k["origin_node"])                            
        end
        push!(O_s,O)
       for u in 1: props(instance.physical_network)[:number_nodes]
            if  props(instance.physical_network, u)[:node_type] == "access" && (u in O_s[s]) == false
                push!(non_O,u)
            end
        end
        push!(non_O_s,non_O)
    end
    #---------------------------------------------
    # SPLIT SELECTION
    #---------------------------------------------)
    @constraint(my_model, split[s in S, f in 1:instance.number_of_AN_based_NFs-1], z[s,f]<=z[s,f+1])
    #---------------------------------------------
    # DIMENSIONING - CP NFSs
    #---------------------------------------------
    @constraint(my_model, dim_cp[s in S, f in (instance.number_of_AN_based_NFs+1):(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),n in  N,u in non_access], w[s,f,n,u]  == instance.set_VNFs[f].dataVolPerUE*instance.setSlices[s].TotalAmountUE*x[s,f,n,u] /instance.set_VNFs[f].treatment_capacity)
 
    #---------------------------------------------
    # DIMENSIONING - DISTRIBUTED DP NFSs
    #---------------------------------------------
    @constraint(my_model,dim_dis_dp[n in N,s in S,k in 1:length(instance.setSlices[s].set_commodities), f in 2:instance.number_of_AN_based_NFs], w[s,f,n,instance.setSlices[s].set_commodities[k]["origin_node"]]  == instance.set_VNFs[f-1].compression*instance.setSlices[s].set_commodities[k]["volume_of_data"]*x[s,f,n,instance.setSlices[s].set_commodities[k]["origin_node"]]/instance.set_VNFs[f].treatment_capacity)
    @constraint(my_model,dimen_dis_dp[n in N,s in S,k in 1:length(instance.setSlices[s].set_commodities)], w[s,1,n,instance.setSlices[s].set_commodities[k]["origin_node"]] == instance.setSlices[s].set_commodities[k]["volume_of_data"]*x[s,1,n,instance.setSlices[s].set_commodities[k]["origin_node"]]/instance.set_VNFs[1].treatment_capacity)

    #---------------------------------------------
    # DIMENSIONING - CENTRALIZED DP NFSs
    #---------------------------------------------
    
    @constraint(my_model, dim_cen_dp[n in N,s in S,  f in 2:instance.number_of_AN_based_NFs, u in non_access],w[s,f,n,u] == instance.setSlices[s].set_VNFs_to_install[f]*(sum(instance.set_VNFs[f-1].compression*instance.setSlices[s].set_commodities[k]["volume_of_data"] for k in 1:length(instance.setSlices[s].set_commodities))*x[s,f,n,u])/instance.set_VNFs[f].treatment_capacity)                
    @constraint(my_model, dim_cenn_dp[n in N,s in S, u in non_access], w[s,1,n,u] ==  instance.setSlices[s].set_VNFs_to_install[1]*(sum(instance.setSlices[s].set_commodities[k]["volume_of_data"] for k in 1:length(instance.setSlices[s].set_commodities))*x[s,1,n,u])/instance.set_VNFs[1].treatment_capacity)                
 
                                           
    #---------------------------------------------
    # PLACEMENT - DISTRIBUTED DP NFSs 
    #---------------------------------------------

        
    @constraint(my_model, plac_dp_in[s in S,f in 1:instance.number_of_AN_based_NFs, u in O_s[s]],sum(x[s,f,m,u] for m in N)== 1 - z[s,f])
    @constraint(my_model, plac_dp_out[s in S,f in 1:instance.number_of_AN_based_NFs],sum(x[s,f,m,u] for m in N,u in non_O_s[s])== 0)
    O_s = nothing
    non_O_s = nothing
            

    #---------------------------------------------
    # PLACEMENT - CENTRALIZED DP NFSs     
    #---------------------------------------------
    Matrix = Vector{Vector{AffExpr}}()     
    for s in S 
        my_s = Vector{AffExpr}()
        push!(Matrix,my_s)
        for f in 1:instance.number_of_AN_based_NFs
            expression = zero(AffExpr)
            push!(Matrix[s], expression)
            for u in non_access
                add_to_expression!(Matrix[s][f], sum(x[s,f,m,u] for m in N)) 
            end
        end
    end
    @constraint(my_model, plac_cent_dp[s in S,  f in 1:instance.number_of_AN_based_NFs], Matrix[s][f] == z[s,f])  
    Matrix = nothing                           
                
    #---------------------------------------------
    # PLACEMENT - CP NFSs
    #---------------------------------------------

        
    Matrix = Vector{Vector{AffExpr}}() 
    for s in S
        my_s = Vector{AffExpr}()
        push!(Matrix,my_s)
        my_s = nothing
        for f in 1:(instance.number_of_AN_based_NFs + instance.number_of_CN_based_NFs)
            expression = zero(AffExpr)
            push!(Matrix[s], expression)
            if f > instance.number_of_AN_based_NFs
                for u in non_access
                    add_to_expression!(Matrix[s][f], sum(x[s,f,m,u] for m in N)) 
                end
            end
 
        end
    end
        
    @constraint(my_model, plac_cp_1[s in S, f in (1+instance.number_of_AN_based_NFs):(instance.number_of_AN_based_NFs + instance.number_of_CN_based_NFs)],Matrix[s][f] == instance.setSlices[s].set_VNFs_to_install[f]) 
    
    Matrix = nothing
    #---------------------------------------------
    # VIRTUAL ISOLATION ================== not in Axis
    #---------------------------------------------
    no_S = Vector()
    for s in SS, t in SS
        if t != s
            push!(no_S, t)
        end
    end
    
     @constraint( my_model, virtual_iso[m in N,s in S, t in no_S[s], u in 1:props(instance.physical_network)[:number_nodes],   f in 1:length(instance.set_VNFs),g in 1:length(instance.set_VNFs)], x[s,f,m,u] + x[t,g,m,u]  <= 1 + instance.setSlices[s].VNF_sharing["$(t)"][f][g]*instance.setSlices[t].VNF_sharing["$(s)"][g][f])  
 
    #---------------------------------------------
    # Physical ISOLATION ================== not in Axis
    #---------------------------------------------         
    @constraint(my_model, phy_iso[s in SS, t in no_S[s], f in 1:length(instance.set_VNFs), g in 1:length(instance.set_VNFs),u in non_access],sum(x[s,f,n,u] for n in N) + sum(x[t,g,m,u] for m in N)<=instance.setSlices[s].nodeSharing[t]*instance.setSlices[t].nodeSharing[s] +1 )

    #---------------------------------------------
    # PACKING  ================== not in Axis
    #--------------------------------------------
     @constraint(my_model, packing[m in N,f in 1:length(instance.set_VNFs),  u in 1: props(instance.physical_network)[:number_nodes]],sum(w[s,f,m,u] for s in S)  <= y[f,m,u])
     @constraint(my_model, one_NF[n in N, s in SS ,t in no_S[s], f in 1:length(instance.set_VNFs), g in 1:length(instance.set_VNFs), u in 1:props(instance.physical_network)[:number_nodes]-1,v in u+1: props(instance.physical_network)[:number_nodes]], x[s,f,n,u] + x[t,g,n,v] <=1)
     no_S = nothing
    #---------------------------------------------
    # NODE CAPACITY ================== not in Axis
    #---------------------------------------------
    @constraint(my_model, cap_ram[u in 1: props(instance.physical_network)[:number_nodes]], sum(instance.set_VNFs[f].ram_request*y[f,m,u] for  m in N,f in 1:length(instance.set_VNFs)) <= props(instance.physical_network,u)[:ram_capacity])
    @constraint(my_model, cap_storage[u in 1: props(instance.physical_network)[:number_nodes]], sum(instance.set_VNFs[f].storage_request*y[f,m,u] for  m in N, f in 1:length(instance.set_VNFs)) <= props(instance.physical_network,u)[:storage_capacity])
    @constraint(my_model, cap_cpu[u in 1: props(instance.physical_network)[:number_nodes]], sum(instance.set_VNFs[f].cpu_request*y[f,m,u] for f in 1:length(instance.set_VNFs), m in N) <= props(instance.physical_network,u)[:cpu_capacity])

     
    #---------------------------------------------
    #  LATENCY
    #---------------------------------------------
    
    end_to_end_delay = Vector{Vector{AffExpr}}()
    pair_delay = Vector{Vector{Vector{Vector{AffExpr}}}}()
    for s in S
        aux_vector = Vector{AffExpr}()
        push!(end_to_end_delay,aux_vector)
        aux_vector1 = Vector{Vector{Vector{AffExpr}}}()
        push!(pair_delay,aux_vector1)
        for k in 1:length(instance.setSlices[s].set_commodities)
            ex = zero(AffExpr)
           aux_vector2 = Vector{Vector{AffExpr}}()
            push!(pair_delay[s],aux_vector2)
            push!(end_to_end_delay[s],ex)
            for f in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
                aux_vector3 = Vector{AffExpr}()
                push!(pair_delay[s][k],aux_vector3)
                for g in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
                    ex1 = zero(AffExpr)
                    push!(pair_delay[s][k][f],ex1)
                    for a in edges(instance.physical_network)
                        delay_a = get_prop(instance.physical_network,a,:delay)
                        if f < instance.number_of_AN_based_NFs && g == f+1
                            add_to_expression!(end_to_end_delay[s][k], delay_a*gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g])
                        end 
                        if f == 1 && g == 1
                            add_to_expression!(end_to_end_delay[s][k], delay_a * (gamma[s,k,get_prop(instance.physical_network,a,:edge_id),length(instance.set_VNFs)+1,1] +gamma[s,k,get_prop(instance.physical_network,a,:edge_id),instance.number_of_AN_based_NFs,length(instance.set_VNFs)+1]))
                         end 
                        if instance.VNF_connection["$(s)"][f][g] == true
                            add_to_expression!(pair_delay[s][k][f][g], delay_a * gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g]) 
                        end
                    end
                add_to_expression!(pair_delay[s][k][f][g],- instance.maxLatencyBetweenFunctions[f][g])  
                end
            end
           add_to_expression!(end_to_end_delay[s][k], -  instance.setSlices[s].maxLatencyDataLayer) 
        end
    end

    @constraint(my_model, delay_pair_NFS[s in S, k in 1:length(instance.setSlices[s].set_commodities),f in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),g in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs) ], pair_delay[s][k][f][g] <= 0)  
    @constraint(my_model,end_2_end_delay[s in S, k in 1:length(instance.setSlices[s].set_commodities) ], end_to_end_delay[s][k] <= 0) 
    end_to_end_delay = nothing
    pair_delay = nothing

    
     #---------------------------------------------
    # FLOW - DP
    #---------------------------------------------
    flow_DP = Vector{Vector{Vector{Vector{AffExpr}}}}()
    flow_DP_ori = Vector{Vector{Vector{Vector{AffExpr}}}}()
    flow_DP_targ = Vector{Vector{Vector{Vector{AffExpr}}}}()
    for s in S
        a = Vector{Vector{Vector{AffExpr}}}()
        push!(flow_DP,a)
        aa = Vector{Vector{Vector{AffExpr}}}()
        push!(flow_DP_ori,a)
        aaa = Vector{Vector{Vector{AffExpr}}}()
        push!(flow_DP_targ,a)
       for k in 1:length(instance.setSlices[s].set_commodities)
            b = Vector{Vector{AffExpr}}()
            push!(flow_DP[s],b)   
            bb = Vector{Vector{AffExpr}}()
            push!(flow_DP_ori[s],b)
            bbb = Vector{Vector{AffExpr}}()
            push!(flow_DP_targ[s],b)
            for f in 1:(instance.number_of_AN_based_NFs-1)
                c = Vector{AffExpr}()
                push!(flow_DP[s][k],c)
                cc = Vector{AffExpr}()
                push!(flow_DP_ori[s][k],c)  
                ccc = Vector{AffExpr}()
                push!(flow_DP_targ[s][k],c)  
                for u in 1: props(instance.physical_network)[:number_nodes]
                    ex0 = zero(AffExpr)
                    push!(flow_DP[s][k][f],ex0) 
                    exx0 = zero(AffExpr)
                    push!(flow_DP_ori[s][k][f],exx0) 
                    exxx0 = zero(AffExpr)
                    push!(flow_DP_targ[s][k][f],exxx0) 
                    for a in edges(instance.physical_network)
                        if src(a) == u
                           add_to_expression!(flow_DP_targ[s][k][f][u], gamma[s,k,get_prop(instance.physical_network,a,:edge_id),instance.number_of_AN_based_NFs,length(instance.set_VNFs)+1])
                           add_to_expression!(flow_DP_ori[s][k][f][u], gamma[s,k,get_prop(instance.physical_network,a,:edge_id),length(instance.set_VNFs)+1,1])
                           add_to_expression!(flow_DP[s][k][f][u], gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,f+1])
                        elseif dst(a) == u
                            add_to_expression!(flow_DP_targ[s][k][f][u],- gamma[s,k,get_prop(instance.physical_network,a,:edge_id),instance.number_of_AN_based_NFs,length(instance.set_VNFs)+1])
                            add_to_expression!(flow_DP_ori[s][k][f][u],- gamma[s,k,get_prop(instance.physical_network,a,:edge_id),length(instance.set_VNFs)+1,1] )
                            add_to_expression!(flow_DP[s][k][f][u], -gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,f+1])
                        end
                    end
                    if props(instance.physical_network, u)[:node_type] != "access"
                        add_to_expression!(flow_DP[s][k][f][u],- sum(x[s,f,m,u] for m in N) +sum(x[s,f+1,m,u] for m in N))
                    elseif u!=instance.setSlices[s].set_commodities[k]["origin_node"]
                        add_to_expression!(flow_DP[s][k][f][u], 0)
                    else               
                         add_to_expression!(flow_DP[s][k][f][u],- z[s,f+1] + z[s,f])
                    end  
                    if  u==instance.setSlices[s].set_commodities[k]["target_node"]
                        add_to_expression!(flow_DP_targ[s][k][f][u], 1)
                        add_to_expression!(flow_DP_ori[s][k][f][u], 0)        

                    elseif  u==instance.setSlices[s].set_commodities[k]["origin_node"]
                        add_to_expression!(flow_DP_targ[s][k][f][u], - 1+z[s,instance.number_of_AN_based_NFs])
                        add_to_expression!(flow_DP_ori[s][k][f][u],- z[s,1])


                    elseif  props(instance.physical_network, u)[:node_type] == "non_access" 
                        add_to_expression!(flow_DP_targ[s][k][f][u],- sum(x[s,instance.number_of_AN_based_NFs,m,u] for m in N))
                        add_to_expression!(flow_DP_ori[s][k][f][u], sum(x[s,1,m,u] for m in N))

                    else 
                        add_to_expression!(flow_DP_targ[s][k][f][u], 0)
                        add_to_expression!(flow_DP_ori[s][k][f][u], 0)
                    end
                end
            end
        end
    end
    
    @constraint(my_model, dp_flow[s in S, k in 1:length(instance.setSlices[s].set_commodities), f in 1:instance.number_of_AN_based_NFs,u in 1: props(instance.physical_network)[:number_nodes]],flow_DP[s][k][f][u]==0 )
    @constraint(my_model, dp_flow_origin[s in S, k in 1:length(instance.setSlices[s].set_commodities), f in 1:instance.number_of_AN_based_NFs,u in 1: props(instance.physical_network)[:number_nodes]],flow_DP_ori[s][k][f][u]==0 )
    @constraint(my_model, dp_flow_target[s in S, k in 1:length(instance.setSlices[s].set_commodities), f in 1:instance.number_of_AN_based_NFs,u in 1: props(instance.physical_network)[:number_nodes]],flow_DP_targ[s][k][f][u]==0 )

    flow_DP = nothing
    flow_DP_ori = nothing
    flow_DP_targ = nothing
    
    #---------------------------------------------
    # FLOW - BETWEEN CP NFSs
    #---------------------------------------------
    flow_CP = Vector{Vector{Vector{Vector{AffExpr}}}}()
    for s in S
        a = Vector{Vector{Vector{AffExpr}}}()
        push!(flow_CP,a)
        for f in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
            b = Vector{Vector{AffExpr}}()
            push!(flow_CP[s],b)
            for g in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
                c = Vector{AffExpr}()
                push!(flow_CP[s][f],c)
                for u in 1: props(instance.physical_network)[:number_nodes]
                    ex = zero(AffExpr)
                    push!(flow_CP[s][f][g],ex)
                    if f> instance.number_of_AN_based_NFs&& g> instance.number_of_AN_based_NFs&& f!=g && instance.VNF_connection["$(s)"][f][g] == true && instance.setSlices[s].set_VNFs_to_install[g] == true && instance.setSlices[s].set_VNFs_to_install[f] == true
                        for a in edges(instance.physical_network)
                          if src(a) == u
                               add_to_expression!(flow_CP[s][f][g][u], gamma[s,1,get_prop(instance.physical_network,a,:edge_id),f,g])         
                            end
                            if dst(a) == u
                               add_to_expression!(flow_CP[s][f][g][u], - gamma[s,1,get_prop(instance.physical_network,a,:edge_id),f,g])          
                            end
                        end
                        add_to_expression!(flow_CP[s][f][g][u], - sum(x[s,f,m,u] for m in N) + sum(x[s,g,m,u] for m in N))
                            end #if
                        end
                    end
                end
            end# s
    @constraint(my_model, cp_flow[s in S, f in (instance.number_of_AN_based_NFs+1):(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),g in (instance.number_of_AN_based_NFs+1):(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),u in 1: props(instance.physical_network)[:number_nodes]],flow_CP[s][f][g][u]==0 )
    flow_CP = nothing    
    #---------------------------------------------
    # FLOW - BETWEEN DP AND CP NFSs
    #---------------------------------------------               
    flow_CP_DP = Vector{Vector{Vector{Vector{Vector{AffExpr}}}}}()
    for s in S
        a = Vector{Vector{Vector{Vector{AffExpr}}}}()
        push!(flow_CP_DP,a)
        for  k in 1:length(instance.setSlices[s].set_commodities)
             b = Vector{Vector{Vector{AffExpr}}}()
            push!(flow_CP_DP[s],b)
            for  f in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
                c = Vector{Vector{AffExpr}}()
                push!(flow_CP_DP[s][k],c)
                for g in 1:instance.number_of_AN_based_NFs
                    d = Vector{AffExpr}()
                    push!(flow_CP_DP[s][k][f],d)
                    for u in 1: props(instance.physical_network)[:number_nodes]
                        ex1 = zero(AffExpr)
                        push!(flow_CP_DP[s][k][f][g],ex1)
                        if f>instance.number_of_AN_based_NFs && instance.VNF_connection["$(s)"][f][g] == true && instance.setSlices[s].set_VNFs_to_install[g] == true && instance.setSlices[s].set_VNFs_to_install[f] == true
                            for a in edges(instance.physical_network)
                                if src(a) == u
                                   add_to_expression!(flow_CP_DP[s][k][f][g][u], gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g])          
                                end
                                if dst(a) == u
                                   add_to_expression!(flow_CP_DP[s][k][f][g][u], - gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g])         
                                end
                            end
                            if props(instance.physical_network, u)[:node_type] == "non_access"           
                                add_to_expression!(flow_CP_DP[s][k][f][g][u],- sum(x[s,f,m,u] for m in N) + sum(x[s,g,m,u] for m in N))

                            elseif props(instance.physical_network, u)[:node_type] == "access" && props(instance.physical_network, u)[:node_id] ==  instance.setSlices[s].set_commodities[k]["origin_node"]           
                                add_to_expression!(flow_CP_DP[s][k][f][g][u],- z[s,g]+1)
                            else
                              add_to_expression!(flow_CP_DP[s][k][f][g][u], - 0)
                            end
                        end# end if
                    end
                end
            end
        end    
    end
    @constraint(my_model, cp_dp_flow[s in S,k in 1:length(instance.setSlices[s].set_commodities), f in (instance.number_of_AN_based_NFs+1):(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),g in 1:instance.number_of_AN_based_NFs,u in 1: props(instance.physical_network)[:number_nodes]],flow_CP_DP[s][k][f][g][u]==0 )
    
    flow_CP_DP = nothing
    flow_CP_DP1 = Vector{Vector{Vector{Vector{Vector{AffExpr}}}}}()
    for s in S
        a = Vector{Vector{Vector{Vector{AffExpr}}}}()
        push!(flow_CP_DP1,a)
        for  k in 1:length(instance.setSlices[s].set_commodities)
             b = Vector{Vector{Vector{AffExpr}}}()
            push!(flow_CP_DP1[s],b)
            for  f in 1:instance.number_of_AN_based_NFs
                c = Vector{Vector{AffExpr}}()
                push!(flow_CP_DP1[s][k],c)
                for g in 1:(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs)
                    d = Vector{AffExpr}()
                    push!(flow_CP_DP1[s][k][f],d)
                    for u in 1: props(instance.physical_network)[:number_nodes]
                        ex1 = zero(AffExpr)
                        push!(flow_CP_DP1[s][k][f][g],ex1)
                        if g>instance.number_of_AN_based_NFs &&instance.VNF_connection["$(s)"][f][g] == true && instance.setSlices[s].set_VNFs_to_install[g] == true && instance.setSlices[s].set_VNFs_to_install[f] == true
                            for a in edges(instance.physical_network)
                                if src(a) == u
                                   add_to_expression!(flow_CP_DP1[s][k][f][g][u], gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g])          
                                end
                                if dst(a) == u
                                   add_to_expression!(flow_CP_DP1[s][k][f][g][u], - gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g])         
                                end
                            end 
                            
                            if props(instance.physical_network, u)[:node_type] == "non_access"           
                                add_to_expression!(flow_CP_DP1[s][k][f][g][u],- sum(x[s,f,m,u] for m in N) + sum(x[s,g,m,u] for m in N))

                            elseif props(instance.physical_network, u)[:node_type] == "access" && props(instance.physical_network, u)[:node_id] ==  instance.setSlices[s].set_commodities[k]["origin_node"]           
                                add_to_expression!(flow_CP_DP1[s][k][f][g][u], z[s,f]-1)
                            else
                              add_to_expression!(flow_CP_DP1[s][k][f][g][u], - 0)
                            end
                        end# end if
                    end
                end
            end
        end    
    end
    @constraint(my_model, dp_cp_flow[s in S,k in 1:length(instance.setSlices[s].set_commodities), f in 1:instance.number_of_AN_based_NFs,g in (instance.number_of_AN_based_NFs+1):(instance.number_of_CN_based_NFs+instance.number_of_AN_based_NFs),u in 1: props(instance.physical_network)[:number_nodes]],flow_CP_DP1[s][k][f][g][u]==0 )
    flow_CP_DP1 = nothing

    

    #---------------------------------------------
    # LINK CAPACITY ================== not in Axis                             
    #---------------------------------------------
    cap_a = Vector{AffExpr}()
    for a in edges(instance.physical_network)           
        sum_gamma = zero(AffExpr)            
        for s in 1:length(instance.setSlices), k in 1:length(instance.setSlices[s].set_commodities)
            add_to_expression!(sum_gamma,gamma[s,k,get_prop(instance.physical_network,a,:edge_id),length(instance.set_VNFs)+1,1]*instance.setSlices[s].set_commodities[k]["volume_of_data"])
            add_to_expression!(sum_gamma,gamma[s,k,get_prop(instance.physical_network,a,:edge_id),instance.number_of_AN_based_NFs,length(instance.set_VNFs)+1]*instance.setSlices[s].set_commodities[k]["volume_of_data"]*instance.set_VNFs[instance.number_of_AN_based_NFs].compression)
            for f in 1:(instance.number_of_AN_based_NFs+instance.number_of_CN_based_NFs), g in 1:(instance.number_of_AN_based_NFs+instance.number_of_CN_based_NFs)
                if f!=g 
                    if f <= instance.number_of_AN_based_NFs && g <= instance.number_of_AN_based_NFs && f==g-1                      
                        add_to_expression!(sum_gamma,gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g]*instance.setSlices[s].set_commodities[k]["volume_of_data"]*instance.set_VNFs[f].compression)
                    end                            
                    if f>instance.number_of_AN_based_NFs && g>instance.number_of_AN_based_NFs
                        add_to_expression!(sum_gamma,gamma[s,1,get_prop(instance.physical_network,a,:edge_id),f,g]*instance.set_VNFs[f].amount_of_traffic_sent_to_g[g]*(instance.setSlices[s].TotalAmountUE/length(instance.setSlices[s].set_commodities)))
                    elseif f>instance.number_of_AN_based_NFs && g<=instance.number_of_AN_based_NFs
                        add_to_expression!(sum_gamma,gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g]*instance.set_VNFs[f].amount_of_traffic_sent_to_g[g]*(instance.setSlices[s].TotalAmountUE/length(instance.setSlices[s].set_commodities)))
                    elseif f<=instance.number_of_AN_based_NFs && g>instance.number_of_AN_based_NFs
                        add_to_expression!(sum_gamma,gamma[s,k,get_prop(instance.physical_network,a,:edge_id),f,g]*instance.set_VNFs[f].amount_of_traffic_sent_to_g[g]*(instance.setSlices[s].TotalAmountUE/length(instance.setSlices[s].set_commodities)))

                    end
                end

            end
        end
        add_to_expression!(sum_gamma, - get_prop(instance.physical_network,a,:max_bandwidth))
        push!(cap_a, sum_gamma)
    end
    @constraint(my_model,cap_link[a in 1:length(cap_a)], cap_a[a] <= 0)

    cap_a = nothing
    function my_lazy_constraint_callback(cb_data)
        TOL = 1e-6                                                                                        
        Threads.@threads for u in 1: props(instance.physical_network)[:number_nodes]
            if props(instance.physical_network, u)[:node_type] != "app"
                for  s in 1:length(instance.setSlices)-1, t in s+1:length(instance.setSlices), f in 1:length(instance.set_VNFs), g in 1:length(instance.set_VNFs), n in 1:instance.number_of_NFs
                    #---------------------------------------------
                    # VIRTUAL ISOLATION
                    #---------------------------------------------
                    if instance.setSlices[s].VNF_sharing["$(t)"][f][g]*instance.setSlices[t].VNF_sharing["$(s)"][g][f]<=0+TOL && callback_value(cb_data, x[s,f,n,u])*callback_value(cb_data, x[t,g,n,u])>0 + TOL 
                        con = @build_constraint(x[s,f,n,u] + x[t,g,n,u]  <= 1)
                        MOI.submit(my_model, MOI.LazyConstraint(cb_data), con)
                        parameters.number_of_lazy_cuts += 1                                                                                            
                    end 
                    #---------------------------------------------
                    # Physical ISOLATION
                    #---------------------------------------------
                    if props(instance.physical_network, u)[:node_type] == "non_access"
                        for m in 1:instance.number_of_NFs
                            if instance.setSlices[s].nodeSharing[t]*instance.setSlices[t].nodeSharing[s]<=0+TOL && callback_value(cb_data, x[s,f,m,u])*callback_value(cb_data, x[t,g,m,u])>0 + TOL                                                                                             
                                con = @build_constraint(x[s,f,m,u] + x[t,g,n,u]  <= 1)
                                MOI.submit(my_model, MOI.LazyConstraint(cb_data), con) 
                                parameters.number_of_lazy_cuts += 1  

                            end
                        end
                    end
                #----------------------------------------------------------------------------                                                                                           
                #packing   
                #-------------------------------------                                                                                           
                    for v in 1:props(instance.physical_network)[:number_nodes] 
                        if props(instance.physical_network, v)[:node_type] != "app"&&  v == u+1 && callback_value(cb_data, x[s,f,n,u])*callback_value(cb_data, x[t,g,n,v])>0 + TOL
                            con = @build_constraint(x[s,f,n,u] + x[t,g,n,v]  <= 1)
                            MOI.submit(my_model, MOI.LazyConstraint(cb_data), con) 
                            parameters.number_of_lazy_cuts += 1 

                        end
                    end                                                                                                
                end 
            end # if
        end # end of all u  

    end # end of my_lazy_constraint_callback function
                                                                 
   #MOI.set(my_model, MOI.LazyConstraintCallback(), my_lazy_constraint_callback)       

    
    @objective(my_model, Min , sum((props(instance.physical_network,u)[:ram_cost]*instance.set_VNFs[f].ram_request+props(instance.physical_network,u)[:storage_cost]*instance.set_VNFs[f].storage_request+props(instance.physical_network,u)[:cpu_cost]*instance.set_VNFs[f].cpu_request)*y[f,m,u] for  m in N, f in 1:length(instance.set_VNFs), u in 1: props(instance.physical_network)[:number_nodes]) )                    
    @dantzig_wolfe_decomposition(my_model, decomposition, S)
    master = getmaster(decomposition)
    subproblems = getsubproblems(decomposition)

    optimize!(my_model)

    return objective_value(my_model)


end
