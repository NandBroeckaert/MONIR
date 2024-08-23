import pandas as pd
from Impact_score_calculation import network_table_reader
from Impact_score_calculation import Node
from Impact_score_calculation import node_level_resetter
from Impact_score_calculation import node_inclusion_in_subnetwork_resetter
from Impact_score_calculation import direct_neighbours_id_interactiontype_interactioninfo_list_of_lists_constructor
from Impact_score_calculation import find_similar_network_node_id

def subnetwork_inclusion_initiator(input_node:Node,network_node_objects_dict: dict,directionality_reaction:str,directionality_other:str):
    """
    This method sets the inclusion_in_subnetwork attribute of the requested direct neighbours of the input_node to True.
    :param input_node:
    :param network_node_objects_dict: The dictionary that contains all the network node objects (format: {node_id:node_object}).
    :param directionality_reaction: Determines for which neighbours that are connected via reaction type interactions the inclusion_in_subnetwork attributes will be changed. Two possible options are available: "bidirectional" (up- and downstream neighbours) and "unidirectional" (only downstream neighbours)".
    :param directionality_other: Determines for which neighbours that are connected via non-reaction type interactions the inclusion_in_subnetwork attributes will be changed. Two possible options are available: "bidirectional" (up- and downstream neighbours) and "unidirectional" (only downstream neighbours)".
    """
    # the  input node inclusion_in_subnetwork should also be set to true, given that it must be included in the subnetwork later
    input_node.inclusion_in_subnetwork = True
    # selecting the requested neighbours
    neighbouring_nodes = direct_neighbours_id_interactiontype_interactioninfo_list_of_lists_constructor(input_node,directionality_reaction,directionality_other)
    # initiating the inclusion_in_subnetwork parameter for neighbours connected via non-reaction type network interactions
    for element in neighbouring_nodes:
        neighbour_id = element[0]
        neighbour = network_node_objects_dict[neighbour_id]
        neighbour.inclusion_in_subnetwork = True #levels should be zero, so there is no need to reset them again


def subnetwork_inclusion_elongator(input_node:Node,previous_level_assigned_node: Node,network_node_objects_dict: dict,directionality_reaction:str,directionality_other:str,distance_level_limit:int):
    """
    This method sets the inclusion_in_subnetwork attribute of the requested direct neighbours of the input_node to True and also updates the level attribute of these nodes (excluding the previous_level_assigned_node).
    :param input_node:
    :param previous_level_assigned_node:
    :param network_node_objects_dict: The dictionary that contains all the network node objects (format: {node_id:node_object}).
    :param directionality_reaction: Determines for which neighbours that are connected via reaction type interactions the inclusion_in_subnetwork attributes will be changed. Two possible options are available: "bidirectional" (up- and downstream neighbours) and "unidirectional" (only downstream neighbours)".
    :param directionality_other: Determines for which neighbours that are connected via non-reaction type interactions the inclusion_in_subnetwork attributes will be changed. Two possible options are available: "bidirectional" (up- and downstream neighbours) and "unidirectional" (only downstream neighbours)".
    :param distance_level_limit: The method will continue exploring the network until the level of the input node is equal to or larger than this distance_level_limit.
    """
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

        neighbouring_nodes = direct_neighbours_id_interactiontype_interactioninfo_list_of_lists_constructor(input_node,directionality_reaction,directionality_other)
        for element in neighbouring_nodes:
            neighbour_id = element[0]
            interaction_type_with_neighbour = element[1]
            neighbour = network_node_objects_dict[neighbour_id]

            # pass over the node that this function was previously called on. Otherwise you could get stuck in a loop.
            if neighbour_id == previous_level_assigned_node.id:
                continue

            # extending through non-reaction type network interactions
            if interaction_type_with_neighbour != "reaction":
                if interaction_type_with_neighbour == "identical_id_connection":
                    # the same level should be assigned to a group of 'identical nodes' (meaning group of nodes that have different variants of the same id).
                    neighbour.level = input_node.level
                    neighbour.inclusion_in_subnetwork = True
                    # the next node will be the start point for assigning weights to the next layer of nodes
                    subnetwork_inclusion_elongator(neighbour, input_node, network_node_objects_dict,directionality_reaction, directionality_other,distance_level_limit)

                else:
                    potential_new_level = input_node.level + 1
                    # if the next node can be reached via several ways from the gene of interest, we should assign the lowest level to this next node (if the level of the next node has been changed before. Remember that the default lvl == 0)
                    if (potential_new_level >= neighbour.level) and neighbour.inclusion_in_subnetwork:
                        continue
                    neighbour.level = potential_new_level
                    neighbour.inclusion_in_subnetwork = True
                    # if better level was found --> continue exploring network in that direction
                    # the next node will be the start point for assigning levels to the next layer of nodes
                    subnetwork_inclusion_elongator(neighbour, input_node, network_node_objects_dict,directionality_reaction, directionality_other, distance_level_limit)

            # extending through reaction type network interactions
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
                    subnetwork_inclusion_elongator(neighbour, input_node, network_node_objects_dict,directionality_reaction, directionality_other, distance_level_limit)

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
                    subnetwork_inclusion_elongator(neighbour, input_node, network_node_objects_dict,directionality_reaction, directionality_other, distance_level_limit)


                elif input_node.type == "compound" and neighbour.type == "compound":
                    raise Exception("Error: compound-compound edge defined as reaction.")
                elif (input_node.type == "gene" or input_node.type == "group_gene") and (
                        neighbour.type == "gene" or neighbour.type == "group_gene"):
                    raise Exception("Error: gene-gene edge defined as reaction.")
                else:  # gene & group_gene
                    raise Exception("Error: odd edge defined as reaction. Probably, group - group.")


def subnetwork_node_list_constructor(network_node_objects_dict: dict,distance_level_limit:int) -> list:
    """
    This method will compile a list of all the nodes for which the inclusion_in_subnetwork attribute is set to True and the level attribute is smaller than the distance_level_limit.
    :param network_node_objects_dict: The dictionary that contains all the network node objects (format: {node_id:node_object}).
    :param distance_level_limit:
    :return: A list of node ids (str).
    """
    if distance_level_limit <= 0:
        raise Exception("The distance_level_limit should be a positive integer.")

    node_id_list = []
    for node_id in network_node_objects_dict:
        node = network_node_objects_dict[node_id]
        if node.inclusion_in_subnetwork and node.level < distance_level_limit:
            node_id_list.append(node_id)
    return node_id_list


def single_subnetwork_table_constructor(full_network_pd:pd.DataFrame,node_id_list: list,incl_neighbours:bool)-> pd.DataFrame:
    """
    This method makes a subnetwork table from the full network table based on the node_id_list.
    :param full_network_pd: A pandas dataframe containing the full network table.
        Format dataframe:
            Column 0: "source_id"
            Column 1: "source_type"
            Column 2: "target_id"
            Column 3: "target_type"
            Column 4: "interaction_type"
            Column 5: "interaction_id"
            Column 6: "interaction_info"
    :param node_id_list: The minimal list of node ids that should be included in the subnetwork table.
    :param incl_neighbours: If set to true, the direct neighbours of nodes in the node_id list will be included in the subnetwork table.
    :return: A pandas dataframe containing the subnetwork table.
    """
    #CREATE SUBNETWORK TABLE
    if incl_neighbours:
        subnetwork_pd = full_network_pd[full_network_pd['source_id'].isin(node_id_list) or full_network_pd['target_id'].isin(node_id_list)]
    else:
        subnetwork_pd = full_network_pd[full_network_pd['source_id'].isin(node_id_list) and full_network_pd['target_id'].isin(node_id_list)]
    return subnetwork_pd


def subnetwork_table_constructor_from_network_file_and_impact_dataframe(path_inputfile_network:str,reverse_interaction_doubled:bool,directionality_reaction:str,directionality_other:str,results_impact_analysis_pd:pd.DataFrame,impact_value_type:str,impact_threshold:float,distance_level_limit:int,distance_level_limit_based_on_impact_results:bool,incl_neighbours:bool,output_directory_and_filename: str):
    """
    This method makes a file, containing a subnetwork, for each node of interest (NOI) that exceeds the specified impact threshold.
    :param path_inputfile_network: path to the input network file.
        Format table in the input network file:
            Column 0: "source_id"
            Column 1: "source_type"
            Column 2: "target_id"
            Column 3: "target_type"
            Column 4: "interaction_type"
            Column 5: "interaction_id"
            Column 6: "interaction_info"
    :param reverse_interaction_doubled: If set to true, reversible reaction edges are contained in the network table in two directions (A to B, B to A).
    :param directionality_reaction: Determines for which neighbours that are connected via reaction type interactions the inclusion_in_subnetwork attributes will be changed. Two possible options are available: "bidirectional" (up- and downstream neighbours) and "unidirectional" (only downstream neighbours)".
    :param directionality_other: Determines for which neighbours that are connected via non-reaction type interactions the inclusion_in_subnetwork attributes will be changed. Two possible options are available: "bidirectional" (up- and downstream neighbours) and "unidirectional" (only downstream neighbours)".
    :param results_impact_analysis_pd: This dataframe contains the results of the impact analysis.
        Format impact analysis dataframe:
            Column 0: "NOI_id"
            Column 1: "NOI_network_id"
            Column 2: "total_impact_score"
            Column 3: "hypothetical_maximum_total_impact_score"
            Column 4: "topological_independent_total_impact_score"
            Column 5: "contributing_nodes"
            Column 6: "contributing_nodes_sub_scores"
            Column 7: "contributing_nodes_max_level"
            Column 8: "contributing_nodes_level"
            Column 9: "top_impacting_node_type"
            Column 10: "node_types"
            Column 11: "node_types_sub_scores"
    :param impact_value_type: Determines which type of impact value will be used for the impact threshold. Two options are available: "total_impact_score" and "topological_independent_total_impact_score".
    :param impact_threshold: Only subnetworks will be made for nodes of interest with an impact value above this threshold (float).
    :param distance_level_limit: Only nodes with a level under this threshold will be included in the subnetwork table. The minimal value is one. In practice, the distance_level_limit is the number of reactions (metabolite -> gene -> metabolite) or gene interactions (gene -> gene) that will be contained in the subnetwork, starting from the direct neighbours of the NOI.
    :param distance_level_limit_based_on_impact_results:If set to True, the given distance_level_limit will be disregarded and the 'contributing_nodes_max_level + 1' will be used instead.
    :param incl_neighbours: If set to true, the direct neighbours of nodes in the node_id list will be included in the subnetwork table.
    :param output_directory_and_filename: The path to the output directory and filename for the output file.
    """
    network_pd = pd.read_csv(path_inputfile_network, sep="\t")
    network_node_dict = network_table_reader(path_inputfile_network,reverse_interaction_doubled)

    #create subnetwork for each node of interest that scored above the impact threshold
    rows = results_impact_analysis_pd.shape[0]
    for row in range(rows):
        NOI_network_id = results_impact_analysis_pd.iloc[rows,1]
        total_impact_score = results_impact_analysis_pd.iloc[rows,2]
        topological_independent_total_impact_score = results_impact_analysis_pd.iloc[rows,4]
        contributing_nodes_max_level = results_impact_analysis_pd.iloc[rows,7]

        if distance_level_limit_based_on_impact_results:
            distance_level_limit = contributing_nodes_max_level

        if ((impact_value_type == 'total_impact_score') and (total_impact_score>=impact_threshold)) or ((impact_value_type == 'topological_independent_total_impact_score') and (topological_independent_total_impact_score>=impact_threshold)):
            NOI = network_node_dict[NOI_network_id]

            #reset inclusion and level attributes of the node objects
            node_level_resetter(network_node_dict)
            node_inclusion_in_subnetwork_resetter(network_node_dict)

            #initiate the requested neighbouring nodes
            subnetwork_inclusion_initiator(NOI,network_node_dict,directionality_reaction,directionality_other,distance_level_limit)

            #determine requested neighbours of starting node
            neighbouring_nodes = direct_neighbours_id_interactiontype_interactioninfo_list_of_lists_constructor(NOI,directionality_reaction,directionality_other)

            #elongate inclusion attributes in the direction of the requested neighbours
            for element in neighbouring_nodes:
                start_node_id = element[0]
                start_node = network_node_dict[start_node_id]
                subnetwork_inclusion_elongator(start_node,NOI,network_node_dict,directionality_reaction,directionality_other,distance_level_limit)

            #make list of nodes that should be included in the subnetwork
            subnetwork_node_id_list = subnetwork_node_list_constructor(network_node_dict)

            #make subnetwork pandas dataframe
            subnetwork_pd = single_subnetwork_table_constructor(network_pd,subnetwork_node_id_list,incl_neighbours)

            #write to output directory
            subnetwork_pd.to_csv(output_directory_and_filename+"_"+NOI_network_id, sep="\t", index=False)


def subnetwork_table_constructor(path_inputfile_network:str,reverse_interaction_doubled:bool,directionality_reaction:str,directionality_other:str,path_inputfile_results_impact_analysis:str,impact_value_type:str,impact_threshold:float,distance_level_limit:int,distance_level_limit_based_on_impact_results:bool,incl_neighbours:bool,output_directory_and_filename: str):
    """
    This method makes a file, containing a subnetwork, for each node of interest (NOI) that exceeds the specified impact threshold.
    :param path_inputfile_network: path to the input network file.
        Format table in the input network file:
            Column 0: "source_id"
            Column 1: "source_type"
            Column 2: "target_id"
            Column 3: "target_type"
            Column 4: "interaction_type"
            Column 5: "interaction_id"
            Column 6: "interaction_info"
    :param reverse_interaction_doubled: If set to true, reversible reaction edges are contained in the network table in two directions (A to B, B to A).
    :param directionality_reaction: A string that determines which nodes that are connected via reaction type interactions will be included in the subnetwork.
            If unidirectional is selected, only downstream nodes of the speciefied interaction type will be included."
            The default setting is bidirectional.
    :param directionality_other: A string that determines which nodes that are connected via non-reaction type interactions will be included in the subnetwork.
            If unidirectional is selected, only downstream nodes of the speciefied interaction type will be included."
            The default setting is bidirectional.
    :param path_inputfile_results_impact_analysis: path to the file that contains the results of the impact analysis.
        Format impact analysis table:
            Column 0: "NOI_id"
            Column 1: "NOI_network_id"
            Column 2: "total_impact_score"
            Column 3: "hypothetical_maximum_total_impact_score"
            Column 4: "topological_independent_total_impact_score"
            Column 5: "contributing_nodes"
            Column 6: "contributing_nodes_sub_scores"
            Column 7: "contributing_nodes_max_level"
            Column 8: "contributing_nodes_level"
            Column 9: "top_impacting_node_type"
            Column 10: "node_types"
            Column 11: "node_types_sub_scores"
    :param impact_value_type: Determines which type of impact value will be used for the impact threshold. Two options are available: "total_impact_score" and "topological_independent_total_impact_score".
    :param impact_threshold: Only subnetworks will be made for nodes of interest with an impact value above this threshold (float).
    :param distance_level_limit: Only nodes with a level under this threshold will be included in the subnetwork table. The minimal value is one. In practice, the distance_level_limit is the number of reactions (metabolite -> gene -> metabolite) or gene interactions (gene -> gene) that will be contained in the subnetwork, starting from the direct neighbours of the NOI.
    :param distance_level_limit_based_on_impact_results:If set to True, the given distance_level_limit will be disregarded and the 'contributing_nodes_max_level + 1' will be used instead.
    :param incl_neighbours: If set to true, the direct neighbours of nodes in the node_id list will be included in the subnetwork table.
    :param output_directory_and_filename: The path to the output directory and filename for the output file.
    """
    results_impact_analysis_pd = pd.read_csv(path_inputfile_results_impact_analysis,sep="\t")

    #create subnetwork for each node of interest that scored above the impact threshold
    subnetwork_table_constructor_from_network_file_and_impact_dataframe(path_inputfile_network,reverse_interaction_doubled,directionality_reaction,directionality_other,results_impact_analysis_pd,impact_value_type,impact_threshold,distance_level_limit,distance_level_limit_based_on_impact_results,incl_neighbours,output_directory_and_filename)


def cytoscape_node_table_nodeshortid_nodetype_extension_constructor(path_inputfile_network:str,output_directory_and_filename: str):
    """
    This method makes a table containing column that lists all the (network) node ids and column that contains the short node ids (part before '_') and a column that contains all the node types (gene or compound). It can be used to extend the node table in cytoscape (this enables colouring based on type).
    RISK: separation of types is based on compound ids containing "cpd:"
    :param path_inputfile_network: The path to the input network file.
        Format table in the input network file:
            Column 0: "source_id"
            Column 1: "source_type"
            Column 2: "target_id"
            Column 3: "target_type"
            Column 4: "interaction_type"
            Column 5: "interaction_id"
            Column 6: "interaction_info"
    :param output_directory_and_filename: The path to the output directory and filename for the output file.
    """
    # INPUT VALIDATON --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #read tsv files containing the network
    network = pd.read_csv(path_inputfile_network,sep="\t")

    #check format of first seven columns
    required_column_names = ["source_id","source_type","target_id","target_type","interaction_type","interaction_id","interaction_info"]

    network_1_column_names = list(network.columns)
    for index in range(7):
        if network_1_column_names[index] != required_column_names[index]:
            raise Exception("Please check the format of the input network tsv files. The first seven columns should be named: source_id,source_type,target_id,target_type,interaction_type,interaction_id,interaction_info.")

    #MAKING THE NODE TABLE EXTENSION (TYPE) ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    sources_set = set(list(network["source_id"]))
    target_set = set(list(network["target_id"]))
    total_set = sources_set + target_set
    total_ids_list = list(total_set)
    total_type_list = []
    total_shortid_list = []

    for id in total_ids_list:
        short_id = id.split('_')[0]
        total_shortid_list.append(short_id)

        if "cpd:" in id:
            total_type_list.append("compound")
        else:
            total_type_list.append("gene")

    node_table_dict = {"node_id":total_ids_list,"short_node_id":total_shortid_list,"node_type":total_type_list}
    node_table_pd = pd.DataFrame(node_table_dict)
    node_table_pd.to_csv(output_directory_and_filename,sep="\t",index=False)


def cytoscape_node_table_general_extension_constructor(path_inputfile_network:str,reverse_interaction_doubled:bool,path_inputfile_node_annotations:str,column_index_ids_annotation_inputfile:int,output_directory_and_filename: str):
    """
    This method will add a column containing network node ids to your annotation table. This table can then be used to extend the node table in cytoscape with the info in the annotation file.
    :param path_inputfile_network: The path to the input network file.
    :param reverse_interaction_doubled: If set to true, reversible reaction edges are contained in the network table in two directions (A to B, B to A).
    :param path_inputfile_node_annotations: The path to the input node annotations file.
    :param column_index_ids_annotation_inputfile: The index of the column containing short or same node ids that are in the network. (e.g. network_id = P06675_pathway1, annotation_id = P06675)
    :param output_directory_and_filename: The path to the output directory and filename for the output file.
    """
    network_dict = network_table_reader(path_inputfile_network,reverse_interaction_doubled)
    node_annotations_pd = pd.read_csv(path_inputfile_node_annotations,sep="\t")

    network_node_id_list = []
    rows = node_annotations_pd.shape[0]
    for row in range(rows):
        annotation_node_id = node_annotations_pd.iloc[row,column_index_ids_annotation_inputfile]
        network_node_id = find_similar_network_node_id(network_dict,annotation_node_id)

        if network_node_id == None:
            network_node_id = "No similar node_id in network"

        network_node_id_list.append(network_node_id)
    node_annotations_pd.insert(0, "network_node_id",network_node_id_list, True)
    node_annotations_pd.to_csv(output_directory_and_filename,sep="\t",index=False)