import pandas as pd
from Impact_score_calculation import unique_list_of_list_constructor
from Impact_score_calculation import network_table_reader
from Impact_score_calculation import Node


def subnetwork_inclusion_initiator(input_node:Node,network_node_objects_dict: dict,directionality_reaction:str,directionality_other:str,distance_level_limit:int):
    input_node.inclusion_in_subnetwork = True
    if distance_level_limit > 0:
        # initiating the inclusion_in_subnetwork parameter for neighbours connected via non-reaction type network interactions
        if directionality_other == "bidirectional":
            next_nodes = input_node.next_nodes
            previous_nodes = input_node.previous_nodes
            non_reaction_neighbouring_nodes = unique_list_of_list_constructor(next_nodes + previous_nodes)  # remove duplicates
        else:  # if directionality_other == "unidirectional"
            non_reaction_neighbouring_nodes = input_node.next_nodes

        for element in non_reaction_neighbouring_nodes:
            neighbour_id = element[0]
            neighbour = network_node_objects_dict[neighbour_id]
            neighbour.inclusion_in_subnetwork = True #levels should be zero, so there is no need to reset them again

        # initiating the inclusion_in_subnetwork parameter for neighbours connected via reaction type network interactions
        if directionality_reaction == "bidirectional":
            next_nodes = input_node.next_nodes
            previous_nodes = input_node.previous_nodes
            reaction_neighbouring_nodes = unique_list_of_list_constructor(next_nodes + previous_nodes)  # remove duplicates
        else:  # if directionality_reaction == "unidirectional"
            reaction_neighbouring_nodes = input_node.next_nodes

        for element in reaction_neighbouring_nodes:
            neighbour_id = element[0]
            neighbour = network_node_objects_dict[neighbour_id]
            neighbour.inclusion_in_subnetwork = True #levels should be zero, so there is no need to reset them again

def subnetwork_level_and_inclusion_elongator(input_node:Node,previous_level_assigned_node: Node,network_node_objects_dict: dict,directionality_reaction:str,directionality_other:str,distance_level_limit:int):
    # validation inputs ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if type(input_node) != Node or type(previous_level_assigned_node) != Node:
        raise Exception("The first two inputs should be Node objects.")

    if type(network_node_objects_dict) != dict:
        raise Exception("The network should be contained in a dictionary: {node_id:node_object}.")
    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception("Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")

    possible_directionalities = ["bidirectional", "unidirectional"]
    if (directionality_reaction not in possible_directionalities) or (directionality_other not in possible_directionalities):
        raise Exception("The two directionionality inputs can only be bidirectional or unidirectional.")

    if (type(distance_level_limit) != int or distance_level_limit < 0):
        raise Exception("The distance_level_limit should be an intiger equal to or above zero.")

    # network extender  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if input_node.level < distance_level_limit:  # keep extending till level limit is reached

        # extending through non-reaction type network interactions
        if directionality_other == "bidirectional":
            next_nodes = input_node.next_nodes
            previous_nodes = input_node.previous_nodes
            neighbouring_nodes = unique_list_of_list_constructor(next_nodes + previous_nodes)  # remove duplicates
        else:  # if directionality_other == "unidirectional"
            neighbouring_nodes = input_node.next_nodes

        for element in neighbouring_nodes:
            neighbour_id = element[0]
            interaction_type_with_neighbour = element[1]
            neighbour = network_node_objects_dict[neighbour_id]

            if neighbour_id == previous_level_assigned_node.id:  # pass over the node that this function was previously called on. Otherwise you could get stuck in a loop.
                continue

            if interaction_type_with_neighbour != "reaction":
                if interaction_type_with_neighbour == "identical_id_connection":
                    # the same level should be assigned to a group of 'identical nodes' (meaning group of nodes that have different variants of the same id).
                    neighbour.level = input_node.level
                    neighbour.inclusion_in_subnetwork = True
                    # the next node will be the start point for assigning weights to the next layer of nodes
                    subnetwork_level_and_inclusion_elongator(neighbour, input_node, network_node_objects_dict,directionality_reaction, directionality_other,distance_level_limit)

                else:
                    potential_new_level = input_node.level + 1
                    # if the next node can be reached via several ways from the gene of interest, we should assign the lowest level to this next node (if the level of the next node has been changed before. Remember that the default lvl == 0)
                    if (potential_new_level >= neighbour.level) and neighbour.inclusion_in_subnetwork:
                        continue
                    neighbour.level = potential_new_level
                    neighbour.inclusion_in_subnetwork = True
                    # if better level was found --> continue exploring network in that direction
                    # the next node will be the start point for assigning levels to the next layer of nodes
                    subnetwork_level_and_inclusion_elongator(neighbour, input_node, network_node_objects_dict,directionality_reaction, directionality_other, distance_level_limit)

        # extending through reaction type network interactions
        if directionality_reaction == "bidirectional":
            next_nodes = input_node.next_nodes
            previous_nodes = input_node.previous_nodes
            neighbouring_nodes = unique_list_of_list_constructor(next_nodes + previous_nodes)  # remove duplicates
        else:  # if directionality_reaction == "unidirectional"
            neighbouring_nodes = input_node.next_nodes

        for element in neighbouring_nodes:
            neighbour_id = element[0]
            interaction_type_with_neighbour = element[1]
            neighbour = network_node_objects_dict[neighbour_id]

            if neighbour_id == previous_level_assigned_node.id:  # pass over the node that this function was previously called on. Otherwise you could get stuck in a loop.
                continue

            if interaction_type_with_neighbour == "reaction":  # interaction type == "reaction": These interaction types have to be handled differently, because you want to assign weights to the metabolites, not the genes.
                if (input_node.type == "gene" or input_node.type == "group_gene") and neighbour.type == "compound":  # should only be passed for the first gene when you enter the reaction network from another type of network
                    # calculate new potential weight and the accompanying level
                    potential_new_level = input_node.level + 1
                    # if the next node can be reached via several ways from the gene of interest, we should assign the lowest level to this next node (if the level of the next node has been changed before. Remember that the default lvl == 0)
                    if (potential_new_level >= neighbour.level) and neighbour.inclusion_in_subnetwork:
                        continue

                    neighbour.level = potential_new_level
                    neighbour.inclusion_in_subnetwork = True
                    # if better level was found --> continue exploring network in that direction
                    # the next node will be the start point for assigning weights to the next layer of nodes
                    subnetwork_level_and_inclusion_elongator(neighbour, input_node, network_node_objects_dict,directionality_reaction, directionality_other, distance_level_limit)

                elif input_node.type == "compound" and (neighbour.type == "gene" or neighbour.type == "group_gene"):
                    # calculate new potential weight and the accompanying level
                    potential_new_level = input_node.level
                    # if the next node can be reached via several ways from the gene of interest, we should assign the lowest level to this next node (if the level of the next node has been changed before. Remember that the default lvl == 0)
                    if (potential_new_level >= neighbour.level) and neighbour.inclusion_in_subnetwork:
                        continue

                    neighbour.level = potential_new_level
                    neighbour.inclusion_in_subnetwork = True
                    # if better level was found --> continue exploring network in that direction
                    # the next node will be the start point for assigning weights to the next layer of nodes
                    subnetwork_level_and_inclusion_elongator(neighbour, input_node, network_node_objects_dict,directionality_reaction, directionality_other, distance_level_limit)


                elif input_node.type == "compound" and neighbour.type == "compound":
                    raise Exception("Error: compound-compound edge defined as reaction.")
                elif (input_node.type == "gene" or input_node.type == "group_gene") and (
                        neighbour.type == "gene" or neighbour.type == "group_gene"):
                    raise Exception("Error: gene-gene edge defined as reaction.")
                else:  # gene & group_gene
                    raise Exception("Error: odd edge defined as reaction. Probably, group - group.")

def subnetwork_minimal_node_list_constructor(network_node_objects_dict: dict) -> list:
    minimal_node_list = []
    for node_id in network_node_objects_dict:
        node = network_node_objects_dict[node_id]
        if node.inclusion_in_subnetwork:
            minimal_node_list.append(node_id)
    return minimal_node_list

def subnetwork_table_constructor_from_node_list(minimal_node_list: list,reverse_interaction_doubled:bool,path_inputfile_results_impact_analysis:str,impact_value_type:str,impact_threshold:float,max_step_distance:int,output_directory_and_filename:str):


def subnetwork_table_constructor_from_network_file(path_inputfile_network_1:str,reverse_interaction_doubled:bool,path_inputfile_results_impact_analysis:str,impact_value_type:str,impact_threshold:float,max_step_distance:int,output_directory_and_filename: str):
    """
    method: build subnetwork of a provided
    """
    network_node_dict = network_table_reader(path_inputfile_network_1,reverse_interaction_doubled)
    results_impact_analysis_pd = pd.read_csv(path_inputfile_results_impact_analysis,sep="\t")

    #INPUT VALIDATION ------------------------------------------------------------------------------------------------------------------------------------------------------------------
    topological_independent_total_impact_score
    total_impact_score

    #CREATING THE SUBNETWORKS -----------------------------------------------------------------------------------------------------------------------------------------------------------
    #create subnetwork for each node of interest that scored above the impact threshold
    rows = results_impact_analysis_pd.shape[0]
    for row in range(rows):
        NOI_id = results_impact_analysis_pd.iloc[rows,0]
        NOI_network_id = results_impact_analysis_pd.iloc[rows,1]
        total_impact_score = results_impact_analysis_pd.iloc[rows,2]
        topological_independent_total_impact_score = results_impact_analysis_pd.iloc[rows,4]

        if ((impact_value_type == 'total_impact_score') and (total_impact_score>=impact_threshold)) or ((impact_value_type == 'topological_independent_total_impact_score') and (topological_independent_total_impact_score>=impact_threshold)):


            """
        Node object:
        self.id = id
        self.type = type
        self.previous_nodes = previous_nodes [id,type,info]
        self.next_nodes = next_nodes [id,type,info]
        self.weight = weight
        self.level = 0
        self.changed = False
        self.changed_omics_type = []
            
        This method reads in a tsv file of the following format:
        Column 0: "source_id"
        Column 1: "source_type"
        Column 2: "target_id"
        Column 3: "target_type"
        Column 4: "interaction_type"
        Column 5: "interaction_id"
        Column 6: "interaction_info"
        """

def cytoscape_node_table_nodetype_extension_constructor(path_inputfile_network_1:str,output_directory_and_filename: str):
    """
    This method make a table containing column that lists all the node ids and column that contains all the node types (gene or compound). It can be used to extend the node table in cytoscape (this enables colouring based on type).
    RISK: seperation of types is based on compound ids containing "cpd"
    """
    # INPUT VALIDATON --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #read tsv files containing the two networks
    network_1 = pd.read_csv(path_inputfile_network_1,sep="\t")

    #check format of first seven columns
    required_column_names = ["source_id","source_type","target_id","target_type","interaction_type","interaction_id","interaction_info"]

    network_1_column_names = list(network_1.columns)
    for index in range(7):
        if network_1_column_names[index] != required_column_names[index]:
            raise Exception("Please check the format of the input network tsv files. The first seven columns should be named: source_id,source_type,target_id,target_type,interaction_type,interaction_id,interaction_info.")

    #MAKING THE NODE TABLE EXTENSION (TYPE) ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    sources_set = set(list(network_1["source_id"]))
    target_set = set(list(network_1["target_id"]))
    total_set = sources_set + target_set
    total_ids_list = list(total_set)
    total_type_list = []
    for id in total_ids_list:
        if "cpd" in id:
            total_type_list.append("compound")
        else:
            total_type_list.append("gene")

    node_table_dict = {"node_id":total_ids_list,"node_type":total_type_list}
    node_table_pd = pd.DataFrame(node_table_dict)
    node_table_pd.to_csv(output_directory_and_filename,sep="\t",index=False)

def cytoscape_node_table_general_extension_constructor(path_inputfile_node_info:str,column_index_unmatched_ids:int,output_directory_and_filename: str):
    print("nothing")
