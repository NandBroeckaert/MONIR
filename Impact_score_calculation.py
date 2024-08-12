"""These methods aim to:
1) Read the tsv file containing all the required information concerning the biological
network (e.g. a transcriptional network, metabolic network, or a combination of the two) and construct a network out of node objects.
2) Calculate an impact scores of genes in the network. Here, the impact score represents a sort of probability that the gene in question plays a role in the omics changes observed in the network.
"""

import pandas as pd
from Subnetwork_selector_for_visualisation import subnetwork_table_constructor_from_network_file_and_impact_dataframe

class Node:
    """
    This class will be used to represent the nodes/network inside python.
    Important to note is that the previous and next nodes variables are list of lists.
    Each element in the list consists of three elements: the id of the next/previous node and the interaction type and info with that node
    """

    def __init__(self, id, type, previous_nodes, next_nodes, weight, level, changed, changed_omics_type,inclusion_in_subnetwork):
        """
        :param id: a string that represents the identifier of the node
            Note: compounds should have a name that starts with 'cpd:' (e.g. cpd:C06730).
            Note: genes that are part of a KEGG pathway will have the following format: genename_pathwayid (e.g. PA2507_rn:R05299)
            Note: a gene that is part of several KEGG pathways will be stored in multiple nodes. Their ids will be very similar, but the second part of the id will be unique due to the unique KEGG pathway id (e.g. gene1_reaction1, gene1_reaction2, gene1_reaction3,...)
        :param type: a string that represents the type of node
            Note: three options are available: gene, group_gene, compound
        :param previous_nodes: a list of lists that contains information on the previous nodes (meaning the sources of directed edges towards this node).
        These lists in this list consist of three parts: the id of the previous node, the interaction type and the interaction info (see KGG_network_construction module for more info on interaction info).
            Format: [[previous_node1_id, previous_node1_interaction_type, previous_node1_interaction_info],[previous_node2_id, previous_node2_interaction_type, previous_node2_interaction_info],...]
            Note: there are seven possible interaction types:chemical,reaction,ECrel,PPrel,GErel,PCrel,identical_id_connection
            Note: possible interaction info: e.g. reversible,irreversible,..   see KEGG KGML webpage for more info about futher details about the possible interaction types
        :param next_nodes: a list of list that contains information on the next nodes (meaning the targets of directed edges starting from this node).
            Format: same as for the previous_nodes parameter.
        :param weight: this float value represent the contribution of this node to the impact score of a certain node in the network
        :param level: this integer value represent the number of steps this node is removed from a certain node in the network
            Note: default value = 0
            Note: direct neighbours of a node have a level of zero (so level assignment starts at zero)
            Note: in a chemical reaction, the level is only raised at each next node that represents a compound.
            Note: you need a second variable if you want to select all nodes within a certain number of steps from a certain node, if you didn't calculate the level value for all nodes in the network (because you will have a pattern like this: start node - 0,1,2,3,...,X,0,0,0,...)
        :param changed: set to true if omics changes were measured for this node.
        :param changed_omics_type: a string which represents the type of omics data for which a change was measured for this node.
            Note: types that can be used for a more specific version of the impact score calculation: proteomics,transcriptomics,methylation,ubiquitination,glycosylation,phosphorylation.
        :param inclusion_in_subnetwork: a boolean parameter that is used to make subnetworks of the large network for downstream visualisation (see subnetwork_selector_for_visualisation.py).
        """

        self.id = id
        self.type = type
        self.previous_nodes = previous_nodes
        self.next_nodes = next_nodes
        self.weight = weight
        self.level = 0
        self.changed = False
        self.changed_omics_type = []
        self.inclusion_in_subnetwork = False


def unique_list_of_list_constructor(non_unique_list_of_lists: list) -> list:
    """
    This function will take a list of lists and return a unique list of lists.
    """
    unique_list_of_lists = [list(x) for x in set(tuple(x) for x in non_unique_list_of_lists)]
    return unique_list_of_lists


def network_table_reader(path_inputfile_network: str, reverse_interaction_doubled: bool) -> dict:
    """
    This method reads in a tsv file containing network information and returns a dictionary of node objects that represent the network {node_id:node_object}.

    :param path_inputfile_network: path to the input file. A tsv file containing all the network information.
        Format network tsv file:
            Column 0: source_id (str)
            Column 1: source_type (str)
            Column 2: target_id (str)
            Column 3: target_type (str)
            Column 4: interaction_type (str) (Note: there are six possible interaction types:chemical,reaction,ECrel,PPrel,GErel,PCrel)
            Column 5: interaction_id (str) (Note: only if present)
            Column 6: interaction_info (str) (e.g. reversible,irreversible,..   see KEGG KGML webpage for more info about futher details about the possible interaction types)
    :param reverse_interaction_doubled: If set to true, reversible reaction edges are contained in the network table in two directions (A to B, B to A).
    :return: a dictionary of node objects, representing a network.
        Format: {node_id:node_object}.
        Note: a node's target previous and next nodes lists should not contain duplicates. A node can only refer once to another node!
        Note: a node that is part of a reversible reaction will be contained in the previous and next nodes lists of nodes that are part of this interaction.
    """
    # validation input ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if (type(path_inputfile_network) != str) or (type(reverse_interaction_doubled) != bool):
        raise Exception("Invalid input types detected. Your first input should be a string that describes the path to the network file. The second input is a boolean type.")
    try:
        network_table = pd.read_csv(path_inputfile_network, sep="\t")
    except:
        raise Exception("No file was found. Please check your path to the file.")

    # check names of first seven columns
    required_column_names = ["source_id", "source_type", "target_id", "target_type", "interaction_type",
                             "interaction_id", "interaction_info"]
    network_column_names = list(network_table.columns)
    for index in range(7):
        if network_column_names[index] != required_column_names[index]:
            raise Exception(
                "Please check the format of the input network tsv files. The first seven columns should be named: source_id,source_type,target_id,target_type,interaction_type,interaction_id,interaction_info.")

    # make first seven columns the correct type -----------------------------------------------------------------------------------------------------------------------------------------------------------------
    network_table["source_id"].astype(str)
    network_table["source_type"].astype(str)
    network_table["target_id"].astype(str)
    network_table["target_type"].astype(str)
    network_table["interaction_type"].astype(str)
    network_table["interaction_id"].astype(str)
    network_table["interaction_info"].astype(str)

    # make dictionary of all nodes and their types
    node_id_sources = network_table["source_id"].tolist()
    node_id_targets = network_table["target_id"].tolist()
    node_ids = node_id_sources + node_id_targets

    node_type_sources = network_table["source_type"].tolist()
    node_type_targets = network_table["target_type"].tolist()
    node_types = node_type_sources + node_type_targets

    dict_node_id_to_type = dict(zip(node_ids,node_types))  # keys is the unique list of nodes. The value (/type) for each key the last value in the node_types list for a particular key (/id)

    # Construct network (construct all node objects)
    network_node_objects_dict = {}
    for node_id in dict_node_id_to_type:
        node_type = dict_node_id_to_type[node_id]

        next_nodes = []
        previous_nodes = []

        if node_id in node_id_sources:

            next_nodes_id = network_table.loc[network_table["source_id"] == node_id, "target_id"].tolist()
            next_nodes_interaction_types = network_table.loc[network_table["source_id"] == node_id, "interaction_type"].tolist()
            next_nodes_interaction_info = network_table.loc[network_table["source_id"] == node_id, "interaction_info"].tolist()

            for i in range(len(next_nodes_id)):
                next_nodes.append([next_nodes_id[i], next_nodes_interaction_types[i], next_nodes_interaction_info[i]])

                if next_nodes_interaction_types[i] == "identical_id_connection":
                    previous_nodes.append([next_nodes_id[i], next_nodes_interaction_types[i], next_nodes_interaction_info[i]])

            # A node that is part of a reversible reaction be contained in the previous and next nodes lists of nodes that are part of this interaction.
            if reverse_interaction_doubled == False:
                next_nodes_id_reversible_KEGG_reactions = network_table.loc[(network_table["source_id"] == node_id) & (network_table["interaction_info"] == "reversible"), "target_id"].tolist()
                next_nodes_interaction_types_reversible_KEGG_reactions = network_table.loc[(network_table["source_id"] == node_id) & (network_table["interaction_info"] == "reversible"), "interaction_type"].tolist()
                next_nodes_interaction_info_reversible_KEGG_reactions = network_table.loc[(network_table["source_id"] == node_id) & (network_table["interaction_info"] == "reversible"), "interaction_info"].tolist()

                extension_previous_nodes = []
                for i in range(len(next_nodes_id_reversible_KEGG_reactions)):
                    extension_previous_nodes.append([next_nodes_id_reversible_KEGG_reactions[i],
                                                     next_nodes_interaction_types_reversible_KEGG_reactions[i],
                                                     next_nodes_interaction_info_reversible_KEGG_reactions[i]])
                previous_nodes += extension_previous_nodes

        if node_id in node_id_targets:

            previous_nodes_id = network_table.loc[network_table["target_id"] == node_id, "source_id"].tolist()
            previous_nodes_interaction_types = network_table.loc[network_table["target_id"] == node_id, "interaction_type"].tolist()
            previous_nodes_interaction_info = network_table.loc[network_table["target_id"] == node_id, "interaction_info"].tolist()

            for i in range(len(previous_nodes_id)):
                previous_nodes.append(
                    [previous_nodes_id[i], previous_nodes_interaction_types[i], previous_nodes_interaction_info[i]])

                if previous_nodes_interaction_types[i] == "identical_id_connection":
                    next_nodes.append([previous_nodes_id[i], previous_nodes_interaction_types[i], previous_nodes_interaction_info[i]])

            # A node that is part of a reversible reaction be contained in the previous and next nodes lists of nodes that are part of this interaction.
            if reverse_interaction_doubled == False:
                previous_nodes_id_reversible_KEGG_reactions = network_table.loc[(network_table["target_id"] == node_id) & (network_table["interaction_info"] == "reversible"), "source_id"].tolist()
                previous_nodes_interaction_types_reversible_KEGG_reactions = network_table.loc[(network_table["target_id"] == node_id) & (network_table["interaction_info"] == "reversible"), "interaction_type"].tolist()
                previous_nodes_interaction_info_reversible_KEGG_reactions = network_table.loc[(network_table["target_id"] == node_id) & (network_table["interaction_info"] == "reversible"), "interaction_info"].tolist()

                extension_next_nodes = []
                for i in range(len(previous_nodes_id_reversible_KEGG_reactions)):
                    extension_next_nodes.append([previous_nodes_id_reversible_KEGG_reactions[i],previous_nodes_interaction_types_reversible_KEGG_reactions[i],previous_nodes_interaction_info_reversible_KEGG_reactions[i]])
                next_nodes += extension_next_nodes

        # Remove duplicates. A node's target previous and next nodes lists should be not contain duplicates. A node can only refer once to another node!
        previous_nodes = unique_list_of_list_constructor(previous_nodes)
        next_nodes = unique_list_of_list_constructor(next_nodes)

        # make network object
        network_node_objects_dict[node_id] = Node(node_id, node_type, previous_nodes, next_nodes, 0, 0, False, "")

    return network_node_objects_dict


def identical_node_list_constructor(input_node: Node, identical_nodes: list, network_node_objects_dict: dict):
    """
    This method will grow an empty list to a list of Node objects that are connected to the original node via identical_id_connection type interactions.

    :param input_node: the starting node of the search
    :param identical_nodes: if provided with an empty list, the method will turn it into a list that contains all Node objects that are connected to the original node via identical_id_connection type interactions.
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    """
    # validation input -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if type(input_node) != Node:
        raise Exception("Ivalid first input. It should be a node object.")
    # the validity of the dictionary is checked only in the special_node_degree_calculator method. If you use this method directly, you have to make sure that the node dictionary is valid {node_id:node_object}.

    # construct a list of all identical id nodes to the first node ---------------------------------------------------------------------------------------------------------------------------------------------
    if input_node not in identical_nodes:
        identical_nodes.append(input_node)

    neighbours = unique_list_of_list_constructor(input_node.previous_nodes + input_node.next_nodes)
    for id_interaction_type_pair in neighbours:
        neighbour_node_id = id_interaction_type_pair[0]
        neighbour_node_interaction_type = id_interaction_type_pair[1]
        neighbour_node = network_node_objects_dict[neighbour_node_id]

        if (neighbour_node_interaction_type == "identical_id_connection") and (neighbour_node not in identical_nodes):
            identical_nodes.append(neighbour_node)
            identical_node_list_constructor(neighbour_node, identical_nodes, network_node_objects_dict)


def identical_node_id_list_constructor(input_node: Node, network_node_objects_dict: dict):
    """
    This method will return a list of Node objects that are connected to the original node via identical_id_connection type interactions.
    Note: this method does not require you to provide an empty list. The identical_node_list_constructor method does (see identical_nodes paramter).

    :param input_node: the starting node of the search
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    :return: a list of Node objects that are connected to the original node via identical_id_connection type interactions.
    """
    identical_node_list = []
    identical_node_list_constructor(input_node, identical_node_list, network_node_objects_dict)

    identical_node_id_list = []
    for identical_node in identical_node_list:
        identical_node_id_list.append(identical_node.id)
    return identical_node_id_list


def id_dict_and_list_of_identical_node_groups_constructor(network_node_objects_dict: dict) -> list:
    """
    This method will create and return a list that contains two things:
        element 0) a list of all the ids of nodes that have an identical_id_connection interaction type with another node.
        element 1) a dictionary in which the keys are the shortest id variant inside a group of nodes that are connected by identical_id_connections, and the values are a list of ids of that group

    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    """
    # validation input -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if type(network_node_objects_dict) != dict:
        raise Exception("The network should be contained in a dictionary: {node_id:node_object}.")
    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception("Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")

    # creating the list and dictionary----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    list_identical_variant_ids = []
    dict_identical_variant_groups = {}

    for node in network_node_objects_dict.values():
        if node.id not in list_identical_variant_ids:
            id_list_of_group_of_identical_variants_around_this_node = identical_node_id_list_constructor(node,network_node_objects_dict)
            if len(id_list_of_group_of_identical_variants_around_this_node) > 1:  # this list will always include the node on which the 'identical_node_id_list_constructor' method was called. As such, you only have identical variants if len > 1.
                sorted_id_list_of_group_of_identical_variants_around_this_node = sorted(id_list_of_group_of_identical_variants_around_this_node,key=len)  # sort so that the shortest id is the first element (this is the true gene/node id)
                dict_identical_variant_groups[sorted_id_list_of_group_of_identical_variants_around_this_node[0]] = sorted_id_list_of_group_of_identical_variants_around_this_node[1::]  # the first element is the true id, the rest are longer variants.
                list_identical_variant_ids += sorted_id_list_of_group_of_identical_variants_around_this_node

    return list_identical_variant_ids, dict_identical_variant_groups


def direct_neighbours_id_interactiontype_interactioninfo_list_of_lists_constructor(input_node:Node,directionality_reaction:str,directionality_other:str) -> list:
    """
    This method will make a list of lists containing info (id, interaction type, interaction type) about the direct neighbours of the input node.

    :param input_node: The starting node of the search.
    :param directionality_reaction: A string that determines which neighbours that are connected to the input node via an interaction of the reaction interaction type will be included in the list of lists.
        Note: two options are available: unidirectional or bidirectional. If unidirectional is selected, only neighbours that are downstream of the input node are included.
    :param directionality_other: A string that determines which neighbours that are connected to the input node via an interaction of the non-reaction interaction type will be included in the list of lists.
        Note: two options are available: unidirectional or bidirectional. If unidirectional is selected, only neighbours that are downstream of the input node are included.
    :return: a list of lists containing info (id, interaction type, interaction type) about the direct neighbours of the input node.
        Format: [[neighbour1_id, neighbour1_inputnode_interactiontype, neighbour1_inputnode_interactioninfo],...]
    """
    previous_nodes = input_node.previous_nodes
    neighbouring_nodes = input_node.next_nodes  # next/downstream nodes are always included. only the previous/upstream nodes are optional.
    for element in previous_nodes:
        previous_node_interaction_type = element[1]
        if (directionality_reaction == "bidirectional") and (previous_node_interaction_type == "reaction"):
            neighbouring_nodes.append(element)
        if (directionality_other == "bidirectional") and (previous_node_interaction_type != "reaction"):
            neighbouring_nodes.append(element)
    neighbouring_nodes = unique_list_of_list_constructor(neighbouring_nodes)
    return neighbouring_nodes


def id_list_of_neighbours_constructor(input_node: Node, network_node_objects_dict: dict, direction_of_neighbours: str,
                                      interaction_type: str) -> list:
    """
    This method will return a list of all the ids of the direct neighbours of a group of nodes that are connected via identical_id_connection type interactions.

    :param input_node: The direct neighbour of the group of nodes that are connected via identical_id_connection type interactions and contains this input node will be returned
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    :param direction_of_neighbours: a string that determines which neighbours will be included in the output based on how the nodes are positioned and connected to the input node. Three options are available: "previous" (aka upstream nodes), "next" (aka downstream nodes), "all".
    :param interaction_type: a string that determines which neighbours will be included in the output based on the type of interaction between the neighbours and the input node. Three options are available: "reaction", "other" (all the non-reaction interaction types), "all".
    :return: a list of node identifiers
    """
    # validation input -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if type(input_node) != Node:
        raise Exception("Ivalid first input. It should be a node object.")

    if type(network_node_objects_dict) != dict:
        raise Exception("The network should be contained in a dictionary: {node_id:node_object}.")
    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception("Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")

    possible_direction_of_neighbours = ["previous", "next", "all"]
    if direction_of_neighbours not in possible_direction_of_neighbours:
        raise Exception("Invalid direction_of_neighbours was specified. Only use one of the following: previous, next, all.")

    possible_interaction_types = ["reaction", "other", "all"]
    if interaction_type not in possible_interaction_types:
        raise Exception("Invalid interaction type was specified. Only use one of the following: reaction, other, all")

    # create a list of neighbour ids------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    identical_node_list = []
    identical_node_list_constructor(input_node, identical_node_list, network_node_objects_dict)

    previous_node_id_interaction_type_interaction_info_for_all_identical_nodes = []
    next_node_id_interaction_type_interaction_info_for_all_identical_nodes = []
    for identical_node in identical_node_list:
        previous_node_id_interaction_type_interaction_info_for_all_identical_nodes += identical_node.previous_nodes
        next_node_id_interaction_type_interaction_info_for_all_identical_nodes += identical_node.next_nodes

    previous_node_id_interaction_type_interaction_info_for_all_identical_nodes = unique_list_of_list_constructor(
        previous_node_id_interaction_type_interaction_info_for_all_identical_nodes)
    next_node_id_interaction_type_interaction_info_for_all_identical_nodes = unique_list_of_list_constructor(
        next_node_id_interaction_type_interaction_info_for_all_identical_nodes)

    # filter out the required ids
    neighbours_set = set()

    if direction_of_neighbours == "previous" or direction_of_neighbours == "all":
        # retrieve set of ids of previous nodes that aren't part of identical_id_connection (normally there should not be multiple connections in the same direction between two nodes)
        previous_nodes_ids_set = set()
        for id_interaction_type_pair in previous_node_id_interaction_type_interaction_info_for_all_identical_nodes:
            previous_node_id = id_interaction_type_pair[0]
            previous_node_interaction_type = id_interaction_type_pair[1]
            if previous_node_interaction_type != "identical_id_connection":
                if previous_node_interaction_type == "reaction" and interaction_type == "reaction":
                    previous_nodes_ids_set.add(previous_node_id)
                elif previous_node_interaction_type != "reaction" and interaction_type == "other":
                    previous_nodes_ids_set.add(previous_node_id)
                elif interaction_type == "all":
                    previous_nodes_ids_set.add(previous_node_id)

        neighbours_set = neighbours_set.union(previous_nodes_ids_set)

    if direction_of_neighbours == "next" or direction_of_neighbours == "all":
        # retrieve set of ids of next nodes that aren't part of identical_id_connection (normally there should not be multiple connections in the same direction between two nodes)
        next_nodes_ids_set = set()
        for id_interaction_type_pair in next_node_id_interaction_type_interaction_info_for_all_identical_nodes:
            next_node_id = id_interaction_type_pair[0]
            next_node_interaction_type = id_interaction_type_pair[1]
            if next_node_interaction_type != "identical_id_connection":
                if next_node_interaction_type == "reaction" and interaction_type == "reaction":
                    next_nodes_ids_set.add(next_node_id)
                elif next_node_interaction_type != "reaction" and interaction_type == "other":
                    next_nodes_ids_set.add(next_node_id)
                elif interaction_type == "all":
                    next_nodes_ids_set.add(next_node_id)

        neighbours_set = neighbours_set.union(next_nodes_ids_set)

    return list(neighbours_set)


def special_node_degree_calculator(input_node: Node, network_node_objects_dict: dict, degree_interaction_type: str,degree_type) -> int:
    """
    This method calculates the in, out or total degree of a note for specific degree interaction types (reaction or other) of a node.

    :param input_node: The input node
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    :param degree_interaction_type: only connections of the specified interaction type will be considered in the degree calculation. Three options are available: "reaction", "other" (aka non-reaction types), "all" (reaction + non-reaction).
    :param degree_type: The type of degree. Three options are available: "indegree", "outdegree", "total".
    :return: The in, out or total degree of a node (based on the specified interaction type)
    """
    # validation input -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if type(input_node) != Node:
        raise Exception("Ivalid first input. It should be a node object.")

    if type(network_node_objects_dict) != dict:
        raise Exception("The network should be contained in a dictionary: {node_id:node_object}.")
    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception(
                "Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")

    possible_interaction_type = ["reaction", "other", "all"]
    if degree_interaction_type not in possible_interaction_type:
        raise Exception(
            "Invalid degree_interaction_type was specified. Only use one of the following: reaction, other, all. ")

    possible_degree_types = ["indegree", "outdegree", "total"]
    if degree_type not in possible_degree_types:
        raise Exception("Invalid degree_type was specified. Only use one of the following: indegree, outdegree, total.")

    # calculate the degree ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    identical_node_list = []
    identical_node_list_constructor(input_node, identical_node_list, network_node_objects_dict)

    # retrieve set of ids of previous nodes that aren't part of identical_id_connection (normally there should not be multiple connections in the same direction between two nodes)
    previous_node_id_interaction_type_interaction_info_for_all_identical_nodes = []
    next_node_id_interaction_type_interaction_info_for_all_identical_nodes = []
    for identical_node in identical_node_list:
        previous_node_id_interaction_type_interaction_info_for_all_identical_nodes += identical_node.previous_nodes
        next_node_id_interaction_type_interaction_info_for_all_identical_nodes += identical_node.next_nodes

    previous_node_id_interaction_type_interaction_info_for_all_identical_nodes = unique_list_of_list_constructor(
        previous_node_id_interaction_type_interaction_info_for_all_identical_nodes)
    next_node_id_interaction_type_interaction_info_for_all_identical_nodes = unique_list_of_list_constructor(
        next_node_id_interaction_type_interaction_info_for_all_identical_nodes)

    previous_nodes_ids_set = set()
    for id_interaction_type_pair in previous_node_id_interaction_type_interaction_info_for_all_identical_nodes:
        previous_node_id = id_interaction_type_pair[0]
        previous_node_interaction_type = id_interaction_type_pair[1]

        if previous_node_interaction_type != "identical_id_connection":
            if degree_interaction_type == "reaction" and previous_node_interaction_type == "reaction":
                previous_nodes_ids_set.add(previous_node_id)
            elif degree_interaction_type == "other" and previous_node_interaction_type != "reaction":
                previous_nodes_ids_set.add(previous_node_id)
            elif degree_interaction_type == "all":
                previous_nodes_ids_set.add(previous_node_id)

    # retrieve set of ids of next nodes that aren't part of identical_id_connection (normally there should not be multiple connections in the same direction between two nodes)
    next_nodes_ids_set = set()
    for id_interaction_type_pair in next_node_id_interaction_type_interaction_info_for_all_identical_nodes:
        next_node_id = id_interaction_type_pair[0]
        next_node_interaction_type = id_interaction_type_pair[1]

        if next_node_interaction_type != "identical_id_connection":
            if degree_interaction_type == "reaction" and next_node_interaction_type == "reaction":
                next_nodes_ids_set.add(next_node_id)
            elif degree_interaction_type == "other" and next_node_interaction_type != "reaction":
                next_nodes_ids_set.add(next_node_id)
            elif degree_interaction_type == "all":
                next_nodes_ids_set.add(next_node_id)

    # calculate degree
    if degree_type == "indegree":
        degree = len(previous_nodes_ids_set)
    elif degree_type == "outdegree":
        degree = len(next_nodes_ids_set)
    elif degree_type == "total":
        all_directly_connected_nodes = previous_nodes_ids_set.union(next_nodes_ids_set)
        degree = len(all_directly_connected_nodes)

    return degree


def node_centrality_modifier_calculator(input_node: Node, network_node_objects_dict: dict, degree_type_other: str, degree_type_reaction: str) -> float:
    """
    This method calculates the centrality modifier of a node.

    :param input_node: node for which the modifier is calculated
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    :param degree_type_other: The type of degree calculated for the non-reaction type interactions. Three options are available: "indegree", "outdegree", "total".
    :param degree_type_reaction: The type of degree calculated for the reaction type interactions. Three options are available: "indegree", "outdegree", "total".
    :return: the centrality modifier of the node (float)
        Note: The reasoning behind this modifier is the following: if a node has more connections, the observed omics change for that node is less likely to be caused by the node for which we are calculating the impact score. Hence, the weight that is contributed to the next node should be lower.
        Note: The modifier is equal to the sum of the degree for reaction and non-reaction type interactions minus one. -1 because there is already a centrality penalty included in the weight inherited by the previous node. As such, we don't want to penalize twice.
        Note: The modifier has a minimum value of 1.
    """

    # validation input -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if type(input_node) != Node:
        raise Exception("This method should be given a node object as its first input.")
    # the validity of the other inputs are checked in the special_node_degree_calculator method.

    # calculate degree ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    degree_other = special_node_degree_calculator(input_node, network_node_objects_dict,"other" ,degree_type_other)
    degree_reaction = special_node_degree_calculator(input_node, network_node_objects_dict,"reaction" ,degree_type_reaction)

    degree = degree_other + degree_reaction

    if degree > 1:  # -1 because there is already a centrality penalty included in the weight inherited by the previous node. As such, we don't want to penalize twice.
        modifier = degree - 1
    else:
        modifier = 1
    return modifier


def node_weight_distribution_initiator(input_node: Node, network_node_objects_dict: dict, start_weight_value: float,
                                       directionality_reaction: str, directionality_other: str):
    """
    This method assigns a starting weight to certain direct neighbour nodes of the group of nodes that are connected via identical_id_connection type interactions and contain this input node

    :param input_node: a starting weight will be assigned to certain direct neighbour nodes of the group of nodes that are connected via identical_id_connection type interactions and contain this input node.
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    :param start_weight_value: weight (float value) that will be assigned to the neighbouring nodes.
    :param directionality_reaction: A string that determines which neighbours that are connected to the input node via an interaction of the reaction interaction type will be assigned the starting weight.
        Note: two options are available: unidirectional or bidirectional. If unidirectional is selected, only neighbours that are downstream of the input node are included.
    :param directionality_other: A string that determines which neighbours that are connected to the input node via an interaction of the non-reaction interaction type will be assigned the starting weight.
        Note: two options are available: unidirectional or bidirectional. If unidirectional is selected, only neighbours that are downstream of the input node are included.
    """
    # check validity of inputs ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    possible_directionalities = ["unidirectional", "bidirectional"]
    if (directionality_reaction not in possible_directionalities) or (
            directionality_other not in possible_directionalities):
        raise Exception("Invalid directionality was specified. Only use one of the following: downstream, bidirectional.")
    if type(start_weight_value) != float or start_weight_value < 0:
        raise Exception("The starting weight should be a float equal to or above zero.")
    if type(input_node) != Node:
        raise Exception("Invalid object type for the first input. Please input a Node object as the first input for this method.")
    if type(network_node_objects_dict) != dict:
        raise Exception("The network should be contained in a dictionary: {node_id:node_object}.")
    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception("Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")

    # assign start weights to surrounding nodes, depending on directionality: ----------------------------------------------------------------------------------------------------------------------------------
    # find ids of the neighbouring nodes for reaction interaction types
    if directionality_reaction == "bidirectional":
        neighbour_id_list_reaction = id_list_of_neighbours_constructor(input_node, network_node_objects_dict, "all",
                                                                       "reaction")
    else:
        neighbour_id_list_reaction = id_list_of_neighbours_constructor(input_node, network_node_objects_dict, "next",
                                                                       "reaction")

    # find ids of the neighbouring nodes for other interactions types
    if directionality_other == "bidirectional":
        neighbour_id_list_other = id_list_of_neighbours_constructor(input_node, network_node_objects_dict, "all", "other")
    else:
        neighbour_id_list_other = id_list_of_neighbours_constructor(input_node, network_node_objects_dict, "next", "other")

    # merge two neighbour lists
    neighbour_id_list = list(set(neighbour_id_list_reaction + neighbour_id_list_other))

    # assign starting weight to neighbours
    for id in neighbour_id_list:
        neighbour_node = network_node_objects_dict[id]
        neighbour_node.weight = start_weight_value


def node_weight_distribution_elongator(
        input_node: Node,
        previous_weight_assigned_node: Node,
        network_node_objects_dict: dict,
        directionality_reaction: str,
        directionality_other: str,
        centrality_modification: bool,
        missingness_modification: bool,
        missingness_modification_step_penalty: float,
        distance_modification: bool,
        distance_modification_step_penalty: float,
        distance_level_limit: int):
    """
    After assigning starting weights, this method will be used to update the weight of the other nodes in the network (following the options that were selected in the input variables).

    :param input_node:
    :param previous_weight_assigned_node:
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    :param directionality_reaction:
    :param directionality_other:
    :param centrality_modification:
    :param missingness_modification:
    :param missingness_modification_step_penalty:
    :param distance_modification:
    :param distance_modification_step_penalty:
    :param distance_level_limit:
    """
    """
    Updates the weights of the nodes in the network.
    For each node, weights are first distributed to the neighbouring non-reaction linked nodes.
    Afterward, weights are assigned to neighbouring nodes connected to the original node via reaction interactions.
    This allows for more flexibility in terms of the directionality that weights are assigned (e.g. unidirectional for transcriptional interactions and bidirectional for reaction interactions)
    """
    # validation inputs ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if type(input_node) != Node or type(previous_weight_assigned_node) != Node:
        raise Exception("The first two inputs should be Node objects.")

    if type(network_node_objects_dict) != dict:
        raise Exception("The network should be contained in a dictionary: {node_id:node_object}.")
    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception(
                "Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")
    possible_directionalities = ["bidirectional", "unidirectional"]

    if (directionality_reaction not in possible_directionalities) or (directionality_other not in possible_directionalities):
        raise Exception("The two directionionality inputs can only be bidirectional or unidirectional.")

    if (type(centrality_modification) != bool) or (type(missingness_modification) != bool) or (type(distance_modification) != bool):
        raise Exception("The centrality_modification, missingness_modification and distance_modification inputs should be of the boolean type.")

    if ((type(missingness_modification_step_penalty) != float) or (type(distance_modification_step_penalty) != float)):
        raise Exception("The penalty inputs should be of the float type.")

    if (type(distance_level_limit) != int or distance_level_limit < 0):
        raise Exception("The distance_level_limit should be an intiger equal to or above zero.")


    # weight distribution --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if input_node.weight != 0 and input_node.level < distance_level_limit:  # check if we need to calculate weights of neighbouring nodes

        neighbouring_nodes = direct_neighbours_id_interactiontype_interactioninfo_list_of_lists_constructor(input_node,directionality_reaction,directionality_other)
        for element in neighbouring_nodes:
            neighbour_id = element[0]
            interaction_type_with_neighbour = element[1]
            neighbour = network_node_objects_dict[neighbour_id]

            # pass over the node that this function was previously called on. Otherwise you could get stuck in a loop.
            if neighbour_id == previous_weight_assigned_node.id:
                continue

            # distribute weight through non-reaction type network interactions
            if interaction_type_with_neighbour != "reaction":
                if interaction_type_with_neighbour == "identical_id_connection":
                    # the same weight and level should be assigned to a group of 'identical nodes' (meaning group of nodes that have different variants of the same id).
                    neighbour.weight = input_node.weight
                    neighbour.level = input_node.level

                    # the next node will be the start point for assigning weights to the next layer of nodes
                    node_weight_distribution_elongator(neighbour, input_node, network_node_objects_dict,directionality_reaction,directionality_other,centrality_modification,missingness_modification,missingness_modification_step_penalty,distance_modification,distance_modification_step_penalty,distance_level_limit)

                else:
                    # calculate new potential weight and the accompanying level
                    neighbour_centrality_modifier = node_centrality_modifier_calculator(neighbour,network_node_objects_dict ,"indegree", "total")
                    potential_new_weight = (input_node.weight / (centrality_modification * neighbour_centrality_modifier+(not centrality_modification))) - (distance_modification * distance_modification_step_penalty) - (missingness_modification * (not input_node.changed) * missingness_modification_step_penalty)
                    potential_new_level = input_node.level + 1

                    # if the next node can be reached via several ways from the gene of interest, we should assign the highest weight and lowest level to this next node
                    if potential_new_weight < 0:
                        potential_new_weight = 0
                    if potential_new_weight > neighbour.weight:
                        neighbour.weight = potential_new_weight
                        neighbour.level = potential_new_level

                        # if better weight was found --> continue exploring network in that direction
                        # the next node will be the start point for assigning weights to the next layer of nodes
                        node_weight_distribution_elongator(neighbour, input_node, network_node_objects_dict,directionality_reaction,directionality_other,centrality_modification,missingness_modification,missingness_modification_step_penalty,distance_modification,distance_modification_step_penalty,distance_level_limit)

            # distribute weights through reaction type network interactions
            if interaction_type_with_neighbour == "reaction":  # interaction type == "reaction": These interaction types have to be handled differently, because you want to assign weights to the metabolites, not the genes.
                if (input_node.type == "gene" or input_node.type == "group_gene") and neighbour.type == "compound":  # should only be passed for the first gene when you enter the reaction network from another type of network
                    # calculate new potential weight and the accompanying level
                    neighbour_centrality_modifier = node_centrality_modifier_calculator(neighbour,network_node_objects_dict,"indegree", "total")
                    potential_new_weight = (input_node.weight / (centrality_modification * neighbour_centrality_modifier+(not centrality_modification))) - (distance_modification * distance_modification_step_penalty) - (missingness_modification * (not input_node.changed) * missingness_modification_step_penalty)
                    potential_new_level = input_node.level + 1

                    # if the next node can be reached via several ways from the gene of interest, we should assign the highest weight and lowest level to this next node
                    if potential_new_weight < 0:
                        potential_new_weight = 0
                    if potential_new_weight > neighbour.weight:
                        neighbour.weight = potential_new_weight
                        neighbour.level = potential_new_level

                        # if better weight was found --> continue exploring network in that direction
                        # the next node will be the start point for assigning weights to the next layer of nodes
                        node_weight_distribution_elongator(neighbour, input_node, network_node_objects_dict,directionality_reaction,directionality_other,centrality_modification,missingness_modification,missingness_modification_step_penalty,distance_modification,distance_modification_step_penalty,distance_level_limit)

                elif input_node.type == "compound" and (neighbour.type == "gene" or neighbour.type == "group_gene"):

                    next_nodes_of_neighbour_node_id_list = [x[0] for x in neighbour.next_nodes]  # get next nodes ids
                    previous_nodes_of_neighbour_node_id_list = [x[0] for x in neighbour.previous_nodes]  # get previous nodes ids
                    unique_previous_nodes_of_neighbour_node_id_list = list(set(previous_nodes_of_neighbour_node_id_list) - set(next_nodes_of_neighbour_node_id_list))  # get unique previous nodes ids (not in next nodes)
                    neighbour_nodes_of_neighbour_node_id_list = next_nodes_of_neighbour_node_id_list + unique_previous_nodes_of_neighbour_node_id_list  # make list of next nodes and unique previous nodes ids

                    for second_order_neighbour_node_id in neighbour_nodes_of_neighbour_node_id_list:
                        second_order_neighbour_node = network_node_objects_dict[second_order_neighbour_node_id]

                        interaction_type_between_neighbour_and_second_order_neighbour_node = ""
                        next_and_previous_nodes_of_neighbour_list = neighbour.next_nodes + neighbour.previous_nodes #list could contain duplicates (not important here)
                        for id_interaction_type_pair in next_and_previous_nodes_of_neighbour_list:
                            id_node = id_interaction_type_pair[0]
                            interaction_type = id_interaction_type_pair[1]
                            if id_node == second_order_neighbour_node_id:
                                interaction_type_between_neighbour_and_second_order_neighbour_node = interaction_type

                        if (second_order_neighbour_node.type == "compound") and (interaction_type_between_neighbour_and_second_order_neighbour_node == "reaction"): #make sure that the second order neighbour is a compound and part of the reaction
                            # calculate new potential weight and the accompanying level second_order_neighbour
                            second_order_neighbour_centrality_modifier = node_centrality_modifier_calculator(second_order_neighbour_node,network_node_objects_dict,"indegree", "total")
                            potential_new_weight = (input_node.weight / (centrality_modification * second_order_neighbour_centrality_modifier+(not centrality_modification))) - (distance_modification * distance_modification_step_penalty) - (missingness_modification * (not input_node.changed) * missingness_modification_step_penalty)
                            potential_new_level = input_node.level + 1

                            if potential_new_weight < 0:
                                potential_new_weight = 0
                            if potential_new_weight > second_order_neighbour_node.weight:
                                second_order_neighbour_node.weight = potential_new_weight
                                second_order_neighbour_node.level = potential_new_level

                                # if you don't want to explore upstream reactions, you should not call the function on upstream nodes
                                if (second_order_neighbour_node_id in unique_previous_nodes_of_neighbour_node_id_list) and (directionality_reaction != "bidirectional"):
                                    continue

                                # if better weight was found --> continue exploring network in that direction
                                # the next node will be the start point for assigning weights to the next layer of nodes
                                node_weight_distribution_elongator(second_order_neighbour_node, neighbour,network_node_objects_dict,directionality_reaction,directionality_other,centrality_modification,missingness_modification,missingness_modification_step_penalty,distance_modification,distance_modification_step_penalty,distance_level_limit)

                elif input_node.type == "compound" and neighbour.type == "compound":
                    raise Exception("Error: compound-compound edge defined as reaction.")
                elif (input_node.type == "gene" or input_node.type == "group_gene") and (
                        neighbour.type == "gene" or neighbour.type == "group_gene"):
                    raise Exception("Error: gene-gene edge defined as reaction.")
                else:  # gene & group_gene
                    raise Exception("Error: odd edge defined as reaction. Probably, group - group.")


def node_weight_resetter(network_node_objects_dict: dict):
    """
    Sets the weight of all the nodes in the network to zero.
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    """
    if type(network_node_objects_dict) != dict:
        raise Exception("The network should be contained in a dictionary: {node_id:node_object}.")

    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception(
                "Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")
        node.weight = 0


def node_level_resetter(network_node_objects_dict: dict):
    """
    sets the level parameter for all the nodes in the network to zero.
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    """

    if type(network_node_objects_dict) != dict:
        raise Exception("The network should be contained in a dictionary: {node_id:node_object}.")

    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception(
                "Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")
        node.level = 0

def node_inclusion_in_subnetwork_resetter(network_node_objects_dict: dict):
    """
    Sets the inclusion_in_subnetwork parameter for all the nodes in the network to False.
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    """

    if type(network_node_objects_dict) != dict:
        raise Exception("The network should be contained in a dictionary: {node_id:node_object}.")

    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception("Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")
        node.inclusion_in_subnetwork = False

def node_changed_and_omics_type_resetter(network_node_objects_dict: dict):
    """
    Sets the changes and changed_omics_type parameter for all the nodes in the network to False and an empty list, respectively.
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    """
    for node_id in network_node_objects_dict.keys():
        # Input validation
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception(
                "Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")

        # reset changed and omics_type attributes
        node.changed = False
        node.changed_omics_type = []


def find_similar_network_node_id(input_node_id: str, network_node_objects_dict: dict):
    """
    This method will attempt to find a node in the network with the same or a very similar id to the input_node_id.
    :param input_node_id: a string that represents a gene of compound or group_gene that is possibly in the network.
        Note: a compound should have 'cpd:' at the start of the id.
        Note: this method assumes that your id has the following format 'name_something' (e.g. PA2507_rn:R05299).
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    :return: a string that represents the id of the node in the network.
        Note: if no similar id is found, it will return None.
    """
    # input validation
    if type(input_node_id) != str:
        raise Exception("The first input should be a string. It should be the id of a node.")

    if type(network_node_objects_dict) != dict:
        raise Exception("The network should be contained in a dictionary: {node_id:node_object}.")
    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception("Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")

    # find node in network with a similar id to the input node id -------------------------------------------------------------------------------------------------------------------------------------
    # check for an exact match --> if exists: return the input_node_id
    if input_node_id in network_node_objects_dict:
        return input_node_id
    # check whether there is a node_id in the network dictionary with the format: input_node_id+'_'+something.
    else:
        for node_id in network_node_objects_dict:
            node_id_short = node_id.split('_')[0]
            if input_node_id == node_id_short:
                return node_id
        # return none if no similar id was found in the network dict
        return None


def node_changed_and_omics_type_updater(path_inputfile_node_omics_info: str, network_node_objects_dict: dict):
    """
    This method updates the 'changed' and 'changed_omics_type' attributes of the nodes and their identical variants (nodes connected via identical_id_connection type interactions) in the network based on info from a tsv file.

    :param path_inputfile_node_omics_info: path to the tsv file that contains the multi-omics info about nodes in the network.
        Format of tsv file:
            column 0: node_id (string)
            column 1: changed (bool)
            column 2: changed_omics_type (string) (seperated by '$')
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    """
    # input validation ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if (type(path_inputfile_node_omics_info) != str) or (type(network_node_objects_dict) != dict):
        raise Exception("Invalid input types detected. Your first input should be a string that describes the path to the node omics info file. The second input should be a dictionary that contains the network {node_id:node_object}.")
    try:  # read in table with multi-omics information about the compounds and genes in the network contained in the dictionary
        node_omics_table = pd.read_csv(path_inputfile_node_omics_info, sep="\t")
    except:
        raise Exception("No file was found. Please check your path to the file.")

    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception(
                "Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")

    # check names of three seven columns
    required_column_names = ["node_id", "changed", "changed_omics_type"]
    node_omics_table_column_names = list(node_omics_table.columns)
    for index in range(3):
        if node_omics_table_column_names[index] != required_column_names[index]:
            raise Exception("Please check the format of the input network tsv files. The first three columns should be named: node_id,changed,changed_omics_type.")

    # make first three columns the correct type
    node_omics_table["node_id"].astype(str)
    node_omics_table["changed"].astype(bool)
    node_omics_table["changed_omics_type"].astype(str)

    # update the 'changed' and 'changed_omics_type' attributes of the nodes in the network dictionary -----------------------------------------------------------------------------------------------------
    rows = node_omics_table.shape[0]

    for row in range(rows):
        node_id = node_omics_table.iloc[row, 0]
        changed = node_omics_table.iloc[row, 1]
        changed_omics_type_str = node_omics_table.iloc[row, 2]
        changed_omics_type_list = changed_omics_type_str.split("$")

        similar_network_node_id = find_similar_network_node_id(node_id, network_node_objects_dict)
        if similar_network_node_id != None:  # check if there is a similar node_id in the network dictionary
            similar_node = network_node_objects_dict[similar_network_node_id]
            similar_node.changed = changed
            similar_node.changed_omics_type = changed_omics_type_list

            # ensure that identical variants of the true node id have a False changed value. This will ensure that only the weight of a genes in a variant group are only counted once.
            list_of_node_ids_identical_to_similar_network_node_id = identical_node_id_list_constructor(similar_node,network_node_objects_dict)
            list_of_node_ids_identical_to_similar_network_node_id.remove(similar_network_node_id)  # remove the similar_network_node_id. Otherwise, you won't count the weight of this node.

            for identical_variant_id in list_of_node_ids_identical_to_similar_network_node_id:
                identical_variant = network_node_objects_dict[identical_variant_id]
                identical_variant.changed = False
                identical_variant.changed_omics_type = changed_omics_type_list


def node_changed_and_omics_type_maximum_background_updater(network_node_objects_dict: dict):
    """
    This method sets the changed and changed_omics_types attributes of all the nodes in the network to true and all possible omics layers that can influence the impact score, respectively.
        Note: This method is used to calculate the maximum impact score for each node.
    :param network_node_objects_dict: the dictionary of all Node objects in the network. Format: {node_id:node_object}.
    """
    for node in network_node_objects_dict.values():
        node.changed = True
        node.changed_omics_type = ["metabolomics","proteomics","transcriptomics","methylation","ubiquitination","glycosylation","phosphorylation","dephosphorylation"]


def impact_calculator(network_node_objects_dict: dict, interaction_specific: bool) -> pd.DataFrame:
    """
    Output: a dataframe containing
    - total impact_score
    - nodes that contribute to impact score
    - sub impact scores of nodes that contribute to impact score
    - levels of contibuting nodes
    - max level of a contributing node
    - main contributing node type
    - node types
    - subscores per node type

    :param network_node_objects_dict:
    :param interaction_specific:
    :return:

    """
    # Don't count indentical variants (ensure via the node_changed_and_omics_type_updater method)
    # Only count the weight if the right omics info is available.

    # calculate all the output values
    total_impact_score = 0
    contributing_nodes = []
    contributing_nodes_sub_scores = []
    contributing_nodes_level = []
    contributing_nodes_max_level = 0
    node_types = []
    node_types_sub_scores = []

    for node_id in network_node_objects_dict.keys():
        node = network_node_objects_dict[node_id]
        if type(node_id) != str or type(node) != Node:
            raise Exception("Wrong object type as key or value in dictionary. Please input a dictionary with {node_id:Node object}.")

        add_node = False

        if node.weight > 0 and node.changed:
            if interaction_specific: #check if the right omics layer is changing
                if (node.type == "gene" or node.type == "group_gene"):
                    changed_omics_type = node.changed_omics_type

                    for i in len(node.previous_nodes):
                        interaction_type = node.previous_node[i][1]
                        interaction_info = node.previous_node[i][2]

                        if (interaction_type == "GErel") and (("proteomics" in changed_omics_type) or ("transcriptomics" in changed_omics_type)):
                            add_node = True
                            break
                        elif ("methylation" in interaction_info) and ("methylation" in changed_omics_type):
                            add_node = True
                            break
                        elif ("ubiquitination" in interaction_info) and ("ubiquitination" in changed_omics_type):
                            add_node = True
                            break
                        elif ("glycosylation" in interaction_info) and ("glycosylation" in changed_omics_type):
                            add_node = True
                            break
                        elif (("phosphorylation" in interaction_info) or ("dephosphorylation" in interaction_info)) and ("phosphorylation" in changed_omics_type):
                            add_node = True
                            break
                        else:
                            continue
                            # The following interaction types,info is not assigned to a certain omics type.
                            # As such, weight of these nodes won't be added to the impact score.
                            # "missing interaction"
                            # "dissociation"
                            # "binding/association"
                            # "state change"
                            # "indirect effect"
                            # "inhibition"
                            # "activation"


                    for i in len(node.next_nodes):
                        interaction_type = node.next_nodes[i][1]
                        interaction_info = node.next_nodes[i][2]

                        if (interaction_type == "GErel") and (("proteomics" in changed_omics_type) or ("transcriptomics" in changed_omics_type)):
                            add_node = True
                            break
                        elif ("methylation" in interaction_info) and (("proteomics" in changed_omics_type) or ("transcriptomics" in changed_omics_type)):
                            add_node = True
                            break
                        elif ("ubiquitination" in interaction_info) and (("proteomics" in changed_omics_type) or ("transcriptomics" in changed_omics_type)):
                            add_node = True
                            break
                        elif ("glycosylation" in interaction_info) and (("proteomics" in changed_omics_type) or ("transcriptomics" in changed_omics_type)):
                            add_node = True
                            break
                        elif (("phosphorylation" in interaction_info) or ("dephosphorylation" in interaction_info)) and (("proteomics" in changed_omics_type) or ("transcriptomics" in changed_omics_type)):
                            add_node = True
                            break
                        else:
                            continue
                            # The following interaction types,info is not assigned to a certain omics type.
                            # As such, weight of these nodes won't be added to the impact score.
                            # "missing interaction"
                            # "dissociation"
                            # "binding/association"
                            # "state change"
                            # "indirect effect"
                            # "inhibition"
                            # "activation"
                            # "repression"
                            # "expression"

                else:
                    if "metabolomics" in node.changed_omics_type:
                        add_node = True
            else:
                add_node = True

        if add_node:
            total_impact_score += node.weight
            contributing_nodes.append(node.id)
            contributing_nodes_sub_scores.append(str(node.weight))
            contributing_nodes_level.append(node.level)

            if node.type in node_types:
                index = node_types.index(node.type)
                node_types_sub_scores[index] += node.weight
            else:
                node_types.append(node.type)
                node_types_sub_scores.append(node.weight)

    # find type of top impacting node
    if len(node_types_sub_scores) == 0:
        top_impacting_node_type = ""
    else:
        top_impacting_node_type = node_types[node_types_sub_scores.index(max(node_types_sub_scores))]

    # find max level of a contributing node:
    if len(node_types_sub_scores) != 0:
        contributing_nodes_max_level = max(contributing_nodes_level)

    # convert lists into strings with $ as separator
    separator = "$"
    contributing_nodes = separator.join(contributing_nodes)
    contributing_nodes_sub_scores = separator.join(contributing_nodes_sub_scores)
    contributing_nodes_level = separator.join(contributing_nodes_level)
    node_types = separator.join(node_types)
    node_types_sub_scores = [str(x) for x in node_types_sub_scores]
    node_types_sub_scores = separator.join(node_types_sub_scores)

    print([str(total_impact_score)])
    print([contributing_nodes])
    print([contributing_nodes_sub_scores])
    print([top_impacting_node_type])
    print([node_types])
    print([node_types_sub_scores])

    # make the dataframe
    node_impact_dict = {"total_impact_score": [str(total_impact_score)], "contributing_nodes": [contributing_nodes],
                        "contributing_nodes_sub_scores": [contributing_nodes_sub_scores],
                        "contributing_nodes_max_level":[str(contributing_nodes_max_level)],"contributing_nodes_level":contributing_nodes_level,
                        "top_impacting_node_type": [top_impacting_node_type], "node_types": [node_types],
                        "node_types_sub_scores": [node_types_sub_scores]}
    node_impact_df = pd.DataFrame.from_dict(node_impact_dict)
    return node_impact_df


def single_node_impact_assessor(
        NOI_id: str,
        network_node_objects_dict: dict,
        directionality_reaction: str,
        directionality_other: str,
        interaction_specific: bool,
        start_weight_value: float,
        centrality_modification: bool,
        missingness_modification: bool,
        missingness_modification_step_penalty: float,
        distance_modification: bool,
        distance_modification_step_penalty: float,
        distance_level_limit: int
) -> pd.DataFrame:
    """
        Output: a dataframe containing
        - total impact_score
        - nodes that contribute to impact score
        - sub impact scores of nodes that contribute to impact score
        - main contributing node type
        - node types
        - subscores per node type
        """
    similar_NOI_id = find_similar_network_node_id(NOI_id,network_node_objects_dict)
    if similar_NOI_id == None:
        similar_NOI_id = "No similar node_id in network"
        result_NOI_dict = {"total_impact_score": [""], "contributing_nodes": [""],
                            "contributing_nodes_sub_scores": [""],
                            "top_impacting_node_type": [""], "node_types": [""],
                            "node_types_sub_scores": [""]}
        result_NOI_pd = pd.DataFrame.from_dict(result_NOI_dict)
    else:
        similar_NOI_node = network_node_objects_dict[similar_NOI_id]

        # reset the level and weight for the impact analysis
        node_weight_resetter(network_node_objects_dict)
        node_level_resetter(network_node_objects_dict)

        # assign starting weights to the nodes flanking the node of interest (NOI)
        node_weight_distribution_initiator(similar_NOI_node, network_node_objects_dict, start_weight_value, directionality_reaction,directionality_other)

        # continue assigning weights to nodes away from the NOI in the specified directions
        # the node_weight_distribution_elongator method requires the id of a previous node. In this case, we need to provide the ids of the NOI variants that are flanking the nodes that were assigned weights earlier.
        # so, we construct a list (of lists) contains elements which consist of the id of a starting node (flank NOI variants) and the closest variant of the NOI.
        if directionality_reaction == 'bidirectional':
            weight_distribution_start_node_id_reaction_list = id_list_of_neighbours_constructor(similar_NOI_node,network_node_objects_dict,'all', "reaction")
        else:
            weight_distribution_start_node_id_reaction_list = id_list_of_neighbours_constructor(similar_NOI_node,network_node_objects_dict,'next', "reaction")

        if directionality_other == 'bidirectional':
            weight_distribution_start_node_id_other_list = id_list_of_neighbours_constructor(similar_NOI_node,network_node_objects_dict,'all', "other")
        else:
            weight_distribution_start_node_id_other_list = id_list_of_neighbours_constructor(similar_NOI_node,network_node_objects_dict,'next', "other")


        weight_distribution_start_node_id_list = list(set(weight_distribution_start_node_id_reaction_list + weight_distribution_start_node_id_other_list))
        weight_distribution_start_node_id_paired_with_NOI_variant_idslist = []  # this list (of lists) contains elements which consist of the id of a starting node (flank NOI variants) and the closest variant of the NOI.

        for id in weight_distribution_start_node_id_list:
            node = network_node_objects_dict[id]
            neighbour_info = unique_list_of_list_constructor(node.previous_nodes + node.next_nodes)
            for element in neighbour_info:
                neighbour_id = element[0]
                if similar_NOI_id in neighbour_id:
                    weight_distribution_start_node_id_paired_with_NOI_variant_idslist.append([id, neighbour_id])

        # Distribute weights throughout the rest of the network
        for id_pair in weight_distribution_start_node_id_paired_with_NOI_variant_idslist:
            starting_node = network_node_objects_dict[id_pair[0]]
            flanking_NOI_variant_node = network_node_objects_dict[id_pair[1]]

            node_weight_distribution_elongator(
                starting_node,
                flanking_NOI_variant_node,
                network_node_objects_dict,
                directionality_reaction,
                directionality_other,
                centrality_modification,
                missingness_modification,
                missingness_modification_step_penalty,
                distance_modification,
                distance_modification_step_penalty,
                distance_level_limit)

        result_NOI_pd = impact_calculator(network_node_objects_dict, interaction_specific)

    result_NOI_pd.insert(0,"NOI_id",[NOI_id],True)
    result_NOI_pd.insert(1,"NOI_network_id",[similar_NOI_id],True)

    return result_NOI_pd


def general_node_impact_assessor(
        path_inputfile_network: str,
        reverse_interaction_doubled: bool,
        path_inputfile_node_omics_info: str,
        path_inputfile_nodes_of_interest: str,
        directionality_reaction: str,
        directionality_other: str,
        interaction_specific: bool,
        start_weight_value: float,
        centrality_modification: bool,
        missingness_modification: bool,
        missingness_modification_step_penalty: float,
        distance_modification: bool,
        distance_modification_step_penalty: float,
        path_output_directory_and_filename: str,
        distance_level_limit=int(100000000),
        include_subnetwork_for_visualisation=True):
    """
    This method is used to perform an impact analysis for a given network (tsv file), list of nodes (tsv file) and multi-omics information (tsv file).

    :param path_inputfile_network: path to the input file. A tsv file containing all the network information.
        Format network tsv file:
            Column 0: source_id (str)
            Column 1: source_type (str)
            Column 2: target_id (str)
            Column 3: target_type (str)
            Column 4: interaction_type (str) (Note: there are six possible interaction types:chemical,reaction,ECrel,PPrel,GErel,PCrel)
            Column 5: interaction_id (str) (Note: only if present)
            Column 6: interaction_info (str) (e.g. reversible,irreversible,..   see KEGG KGML webpage for more info about futher details about the possible interaction types)
    :param reverse_interaction_doubled: If set to true, reversible reaction edges are contained in the network table in two directions (A to B, B to A).
    :param path_inputfile_node_omics_info: path to the tsv file that contains the multi-omics info about nodes in the network.
        Format of tsv file:
            column 0: node_id (string)
            column 1: changed (bool)
            column 2: changed_omics_type (string) (seperated by '$')
    :param path_inputfile_nodes_of_interest:
    :param directionality_reaction:
    :param directionality_other:
    :param interaction_specific:
    :param start_weight_value:
    :param centrality_modification:
    :param missingness_modification:
    :param missingness_modification_step_penalty:
    :param distance_modification:
    :param distance_modification_step_penalty:
    :param path_output_directory_and_filename:
    :param distance_level_limit:
    :param include_subnetwork_for_visualisation:
    :return:
    """
    """
        Output: a dataframe containing
        - total impact_score
        - nodes that contribute to impact score
        - sub impact scores of nodes that contribute to impact score
        - main contributing node type
        - node types
        - subscores per node type
        """
    # read in a list of nodes of interest
    if type(path_inputfile_nodes_of_interest) != str:
        raise Exception(
            "Invalid input type detected. Your first input should be a string that describes the path to the network file")
    try:
        nodes_of_interest_ids_table = pd.read_csv(path_inputfile_nodes_of_interest, sep="\t")
        nodes_of_interest_ids_table['node_id'].astype(str)
        nodes_of_interest_ids_list = nodes_of_interest_ids_table['node_id'].tolist()
    except:
        raise Exception("No file was found. Please check your path to the file.")

    # construct the network
    network_node_objects_dict = network_table_reader(path_inputfile_network, reverse_interaction_doubled)

    print('STANDARD --------------------------------------------------------------------------')
    #STANDARD IMPACT ANALYSIS
    # add omics information to nodes in the network
    node_changed_and_omics_type_updater(path_inputfile_node_omics_info, network_node_objects_dict)
    # perform the standard impact analysis for all nodes of interest using the available omics data
    total_standard_impact_analysis_table = pd.DataFrame(columns=['NOI_id','NOI_network_id','total_impact_score', 'contributing_nodes', 'contributing_nodes_sub_scores', 'top_impacting_node_type','node_types', 'node_types_sub_scores'])
    for node_id in nodes_of_interest_ids_list:
        standard_node_impact_table = single_node_impact_assessor(
            node_id,
            network_node_objects_dict,
            directionality_reaction,
            directionality_other,
            interaction_specific,
            start_weight_value,
            centrality_modification,
            missingness_modification,
            missingness_modification_step_penalty,
            distance_modification,
            distance_modification_step_penalty,
            distance_level_limit)
        total_standard_impact_analysis_table = pd.concat([total_standard_impact_analysis_table, standard_node_impact_table])

    print('MAXIMUM BACKGROUND --------------------------------------------------------------------------')
    #MAXIMUM BACKGROUND IMPACT ANALYSIS
    # change the omics data to the omics background for the calculation of the maximal node impact values
    node_changed_and_omics_type_maximum_background_updater(network_node_objects_dict)
    # perform the impact analysis to calculate the maximum impact score for all nodes of interest

    total_maximum_impact_analysis_table = pd.DataFrame(
        columns=['NOI_id', 'NOI_network_id', 'total_impact_score', 'contributing_nodes',
                 'contributing_nodes_sub_scores','contributing_nodes_max_level','contributing_nodes_level','top_impacting_node_type', 'node_types', 'node_types_sub_scores'])
    for node_id in nodes_of_interest_ids_list:
        maximum_node_impact_table = single_node_impact_assessor(
            node_id,
            network_node_objects_dict,
            directionality_reaction,
            directionality_other,
            interaction_specific,
            start_weight_value,
            centrality_modification,
            missingness_modification,
            missingness_modification_step_penalty,
            distance_modification,
            distance_modification_step_penalty,
            distance_level_limit)
        total_maximum_impact_analysis_table = pd.concat([total_maximum_impact_analysis_table, maximum_node_impact_table])

    print('TOPOLOGICAL INDEPENDENT IMPACT SCORES --------------------------------------------------------------------------')
    #CALCULATE TOPOLOGICAL INDEPENDENT IMPACT SCORES

    maximum_total_impact_score_list = total_maximum_impact_analysis_table["total_impact_score"].tolist()
    standard_total_impact_score_list = total_standard_impact_analysis_table["total_impact_score"].tolist()
    topological_independent_impact_score_list = []
    for i in range(len(standard_total_impact_score_list)):
        standard_total_impact_score = standard_total_impact_score_list[i]
        maximum_total_impact_score = maximum_total_impact_score_list[i]

        if (len(standard_total_impact_score) >= 1) and (len(maximum_total_impact_score) >=1):
            standard_total_impact_score = float(standard_total_impact_score)
            maximum_total_impact_score = float(maximum_total_impact_score)

            topological_independent_impact_score = standard_total_impact_score/maximum_total_impact_score
            topological_independent_impact_score_list.append(str(topological_independent_impact_score))
        else:
            topological_independent_impact_score_list.append('')

    #WRITING OUTPUT OF IMPACT ANALYSIS
    # add columns to output table
    total_extended_impact_analysis_table = total_standard_impact_analysis_table
    total_extended_impact_analysis_table.insert(loc=3,column="hypothetical_maximum_total_impact_score",value=maximum_total_impact_score_list)
    total_extended_impact_analysis_table.insert(loc=4, column="topological_independent_total_impact_score",value=topological_independent_impact_score_list)

    # write to total impact table to tsv file
    total_extended_impact_analysis_table.to_csv(path_output_directory_and_filename, sep="\t", index=False)

    #GENERATE SUBNETWORKS FOR VISUALIZING RESULTS OF IMPACT ANALYSIS
    if include_subnetwork_for_visualisation:
        subnetwork_table_constructor_from_network_file_and_impact_dataframe(path_inputfile_network,reverse_interaction_doubled,directionality_reaction,directionality_other,total_extended_impact_analysis_table,"total_impact_score",start_weight_value,distance_level_limit,True,True)


# analyzing acetylomics wih metabolomics  -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#14-1
general_node_impact_assessor("C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/KEGG_network_all_metabolic_reactions_pae_17052024.txt",False,'C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/14_1_MONIT_omics_format.txt','C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/14_1_MONIT_NOI_KEGG_format.txt','bidirectional','bidirectional',True,10.0,True,True,0.5,True,0.5,'C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Output/14_1_impact_NOI_acetylomcs_OMICS_metabolomics_v2.txt',10000)

#LUZ19
general_node_impact_assessor("C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/KEGG_network_all_metabolic_reactions_pae_17052024.txt",False,'C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/LUZ19_MONIT_omics_format.txt','C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/LUZ19_MONIT_NOI_KEGG_format.txt','bidirectional','bidirectional',True,10.0,True,True,0.5,True,0.5,'C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Output/LUZ19_impact_NOI_acetylomcs_OMICS_metabolomics_v2.txt',10000)

#PEV2
general_node_impact_assessor("C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/KEGG_network_all_metabolic_reactions_pae_17052024.txt",False,'C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/PEV2_MONIT_omics_format.txt','C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/PEV2_MONIT_NOI_KEGG_format.txt','bidirectional','bidirectional',True,10.0,True,True,0.5,True,0.5,'C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Output/PEV2_impact_NOI_acetylomcs_OMICS_metabolomics_v2.txt',10000)

#PhiKZ
general_node_impact_assessor("C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/KEGG_network_all_metabolic_reactions_pae_17052024.txt",False,'C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/PhiKZ_MONIT_omics_format.txt','C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/PhiKZ_MONIT_NOI_KEGG_format.txt','bidirectional','bidirectional',True,10.0,True,True,0.5,True,0.5,'C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Output/PhiKZ_impact_NOI_acetylomcs_OMICS_metabolomics_v2.txt',10000)

#YuA
general_node_impact_assessor("C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/KEGG_network_all_metabolic_reactions_pae_17052024.txt",False,'C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/YUA_MONIT_omics_format.txt','C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Input/formatted/YuA_MONIT_NOI_KEGG_format.txt','bidirectional','bidirectional',True,10.0,True,True,0.5,True,0.5,'C:/Nand_phd/data_en_analyse/WP3/Integration_acetylomics_metabolomics/Output/YuA_impact_NOI_acetylomcs_OMICS_metabolomics_v2.txt',10000)


