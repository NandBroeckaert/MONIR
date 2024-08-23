import argparse

def subparser_KEGG_network_constructor(subparsers):
    """
    This method creates the subparser for the KEGG_network_constructor_from_string method.
    :param subparsers: object to which the subparser should be added.
    """
    parser_KEGG_network_constructor = subparsers.add_parser(
        'KEGG_network_constructor',
        help = "This method will create tsv file containing network information of the selected KEGG pathways.")
    parser_KEGG_network_constructor.add_argument(
        "-p",
        "--pathways",
        required=True,
        type=str,
        help='A string containing the KEGG pathway ids of a specific organism in the KEGG database separated by the dollar sign, which will be included in the network. '
             'If you want to include all KEGG metabolic reactions, write all_metabolic_reactions.'
    )
    parser_KEGG_network_constructor.add_argument(
        "-c",
        "--organism_code",
        required=True,
        type=str,
        help='The code that KEGG uses for a specific organism (e.g. pae for Pseudomonas aeruginosa PAO1).'
    )
    parser_KEGG_network_constructor.add_argument(
        "-t",
        "--selected_interaction_types",
        required=True,
        type=str,
        help= 'A string of all the interaction types that need to be included in the network table separated by the dollar sign. '
              'This network may contain the following types of interactions: chemical, reaction, ECrel, PPrel, GErel, PCrel. '
              'The following types cannot be selected together: (chemical & reaction),(reaction & ECrel).'
    )
    parser_KEGG_network_constructor.add_argument(
        "-o",
        "--path_outputfile",
        required=True,
        type=str,
        help='The path to the output file (includes file name) (e.g. C:/test/all_metabolic_reactions_pae.tsv).'
    )
    parser_KEGG_network_constructor.add_argument(
        "-r",
        "--reverse_interactions",
        type=bool,
        action='store_true',
        help='If this option is specified, reversible reaction edges will be contained in the network table in two directions (A to B, B to A). False by default.'
    )

def subparser_network_merger(subparsers):
    """
    This method creates the subparser for the network_merger method.
    :param subparsers: object to which the subparser should be added.
    """

def subparser_node_impact_assessor(subparsers):
    """
    This method creates the subparser for the general_node_impact_assessor method.
    :param subparsers: object to which the subparser should be added.
    """
    parser_node_impact_assessor = subparsers.add_parser(
        'node_impact_assessor',
        help="This method is used to perform an impact analysis for a given network (tsv file), list of nodes (tsv file) and multi-omics information (tsv file)."
             "Additionally, subnetworks for visualsation of high scoring nodes of interest can automatically be made.")
    parser_node_impact_assessor.add_argument(
        "-n",
        "--path_inputfile_network",
        required=True,
        type=str,
        help="The path to a tsv file containing all the network information."
    )
    parser_node_impact_assessor.add_argument(
        "-m",
        "--path_inputfile_node_omics_info",
        required=True,
        type=str,
        help="The path to the tsv file that contains the multi-omics info about nodes in the network."
    )
    parser_node_impact_assessor.add_argument(
        "-t",
        "--path_inputfile_nodes_of_interest",
        required=True,
        type=str,
        help="The path to the tsv file that contains the list of identifiers for which you want to do the impact analysis (nodes of interest)."
    )
    parser_node_impact_assessor.add_argument(
        "-r",
        "--reverse_interactions",
        type=bool,
        action='store_true',
        help="If this option is specified, it is assumed that reversible reaction edges are contained in the provided network tsv file in two directions (A to B, B to A). False by default."
    )
    parser_node_impact_assessor.add_argument(
        "-a",
        "--directionality_reaction",
        type=str,
        choices=["unidirectional","bidirectional"],
        default="bidirectional",
        help="A string that determines which nodes that are connected (directly and indirectly) to the node of interest (NOI) with reaction type interactions will be taken into account for the impact analysis."
             "If unidirectional is selected, only downstream nodes of the speciefied interaction type will possibly contribute to the impact of a NOI. The default setting is bidirectional."

    )
    parser_node_impact_assessor.add_argument(
        "-c",
        "--directionality_other",
        type=str,
        choices=["unidirectional","bidirectional"],
        default="unidirectional",
        help="A string that determines which nodes that are connected (directly and indirectly) to the node of interest (NOI) with non-reaction type interactions will be taken into account for the impact analysis."
             "If unidirectional is selected, only downstream nodes of the speciefied interaction type will possibly contribute to the impact of a NOI. The default setting is unidirectional."
    )
    parser_node_impact_assessor.add_argument(
        "-i",
        "--interaction_specific",
        type=bool,
        action='store_true',
        help="If specified, the type of available omics data will be taken into account when calculating the impact score of the node of interest."
             "False by default."
    )
    parser_node_impact_assessor.add_argument(
        "-s",
        "--start_weight_value",
        type=float,
        default=float(10),
        help="The weight (float value) that will be assigned to the direct neighbour nodes of the group of nodes that are connected via identical_id_connection type interactions and contain the node of interest."
             "10 by default."
    )
    parser_node_impact_assessor.add_argument(
        "-c",
        "--centrality_modification",
        type=bool,
        action='store_false',
        help="If specified, no centrality modifier will used for calculating impact scores."
    )
    parser_node_impact_assessor.add_argument(
        "-g",
        "--missingness_modification",
        type=bool,
        action='store_false',
        help="If specified, no missingness modifier will used for calculating impact scores."
    )
    parser_node_impact_assessor.add_argument(
        "-e",
        "--missingness_modification_step_penalty",
        type=float,
        default=float(0.3),
        help="A float value that represents the penalty for missing data."
             "0.3 by default"
    )
    parser_node_impact_assessor.add_argument(
        "-d",
        "--distance_modification",
        type=bool,
        action='store_false',
        help="If specified, no distance modifier will used for calculating impact scores."
    )
    parser_node_impact_assessor.add_argument(
        "-f",
        "--distance_modification_step_penalty",
        type=float,
        default=float(0.3),
        help="A float value that represents the penalty for distance to the node for which you want to calculate the impact score."
             "0.3 by default."
    )
    parser_node_impact_assessor.add_argument(
        "-l",
        "--distance_level_limit",
        type=int,
        default=int(100000000),
        help="The weight distribution process starts at zero and will be halted when at this node level. Only nodes till this level can contribute to the impact scores."
             "100000000 by default."
    )
    parser_node_impact_assessor.add_argument(
        "-o",
        "--path_output_directory_and_filename",
        required=True,
        type=str,
        help="The path to the output tsv file that will be created that contains the results of the impact analysis."
    )
    parser_node_impact_assessor.add_argument(
        "-j",
        "--include_subnetwork_for_visualisation",
        type=bool,
        action='store_false',
        help="If this option is specified, subnetworks that are useful for visualization will not be generated. By defualt subnetworks for all NOIs hat have an total impact score above the start_weight_value."
    )


def subparser_subnetwork_table_constructor(subparsers):
    """
    This method creates the subparser for the subnetwork_table_constructor method.
    :param subparsers: object to which the subparser should be added.
    """
    parser_subnetwork_table_constructor = subparsers.add_parser(
        'subnetwork_constructor',
        help="This method makes a tsv file, containing a subnetwork, for each node of interest (NOI) in an impact analysis that exceeds the specified impact threshold.")
    parser_subnetwork_table_constructor.add_argument(
        "-n",
        "--path_inputfile_network",
        required=True,
        type=str,
        help="The path to a tsv file containing all the network information."
    )
    parser_subnetwork_table_constructor.add_argument(
        "-i",
        "--path_inputfile_results_impact_analysis",
        required=True,
        type=str,
        help="The path to the tsv file that contains the results of the impact analysis."
    )
    parser_subnetwork_table_constructor.add_argument(
        "-r",
        "--reverse_interactions",
        type=bool,
        action='store_true',
        help="If this option is specified, it is assumed that reversible reaction edges are contained in the provided network tsv file in two directions (A to B, B to A). False by default."
    )
    parser_subnetwork_table_constructor.add_argument(
        "-a",
        "--directionality_reaction",
        type=str,
        choices=["unidirectional","bidirectional"],
        default="bidirectional",
        help="A string that determines which nodes that are connected via reaction type interactions will be included in the subnetwork."
             "If unidirectional is selected, only downstream nodes of the speciefied interaction type will be included."
             "The default setting is bidirectional."
    )
    parser_subnetwork_table_constructor.add_argument(
        "-c",
        "--directionality_other",
        type=str,
        choices=["unidirectional","bidirectional"],
        default="unidirectional",
        help="A string that determines which nodes that are connected via non-reaction type interactions will be included in the subnetwork."
             "If unidirectional is selected, only downstream nodes of the speciefied interaction type will be included."
             "The default setting is bidirectional."
    )
    parser_subnetwork_table_constructor.add_argument(
        "-t",
        "--threshold_impact_value_type",
        type=str,
        choices=["total_impact_score","topological_independent_total_impact_score"],
        default="total_impact_score",
        help="Determines which type of impact value will be used for the impact threshold. Two options are available: total_impact_score and topological_independent_total_impact_score."
             "'total_impact_score' by default."
    )
    parser_subnetwork_table_constructor.add_argument(
        "-b",
        "--impact_threshold_value",
        type=float,
        required=True,
        help="Only subnetworks will be made for nodes of interest with an impact value above this threshold (float)"
    )
    parser_subnetwork_table_constructor.add_argument(
        "-l",
        "--distance_level_limit",
        type=int,
        default =int(5),
        help="Only nodes with a level under this threshold will be included in the subnetwork table. The minimal value is one."
             "In practice, the distance_level_limit is the number of reactions (metabolite -> gene -> metabolite) or gene interactions (gene -> gene) that will be contained in the subnetwork, starting from the direct neighbours of the NOI."
             "The default setting is 5."
    )
    parser_subnetwork_table_constructor.add_argument(
        "-d",
        "--distance_level_limit_based_on_impact_results",
        type=bool,
        action='store_false',
        help= "If this option is not specified, the given distance_level_limit will be disregarded and set to the 'contributing_nodes_max_level + 1'."
    )
    parser_subnetwork_table_constructor.add_argument(
        "-e",
        "--incl_neighbours",
        type=bool,
        action='store_false',
        help="If this option is not specified, the direct neighbours of nodes in the node_id list will be included in the subnetwork table."
    )
    parser_subnetwork_table_constructor.add_argument(
        'o',
        '--output_directory_and_filename',
        type=str,
        help="The path to the output tsv file that will be created (includes the name of the file)."
    )

def subparser_cytoscape_node_table_nodeshortid_nodetype_extension_constructor(subparsers):
    """
    This method creates the subparser for the cytoscape_node_table_nodeshortid_nodetype_extension_constructor method.
    :param subparsers: object to which the subparser should be added.
    """
    parser_node_table_id_and_type_extender = subparsers.add_parser(
        'node_table_id_and_type_extender',
        help="This method makes a tsv file, containing a subnetwork, for each node of interest (NOI) in an impact analysis that exceeds the specified impact threshold.")
    parser_node_table_id_and_type_extender.add_argument(
        "-n",
        "--path_inputfile_network",
        required=True,
        type=str,
        help="The path to a tsv file containing all the network information."
    )
    parser_node_table_id_and_type_extender.add_argument(
        'o',
        '--output_directory_and_filename',
        type=str,
        help="The path to the output tsv file that will be created (includes the name of the file)."
    )

def subparser_cytoscape_node_table_general_extension_constructor(subparsers):
    """
    This method creates the subparser for the cytoscape_node_table_general_extension_constructor method.
    :param subparsers: object to which the subparser should be added.
    """
    parser_annotation_table_id_extender = subparsers.add_parser(
        'annotation_table_id_extender',
        help="This method makes a tsv file, containing a subnetwork, for each node of interest (NOI) in an impact analysis that exceeds the specified impact threshold."
             "This table can then be used to extend the node table in cytoscape with the info in the annotation file.")
    parser_annotation_table_id_extender.add_argument(
        "-n",
        "--path_inputfile_network",
        required=True,
        type=str,
        help="The path to a tsv file containing all the network information."
    )
    parser_annotation_table_id_extender.add_argument(
        "-r",
        "--reverse_interactions",
        type=bool,
        action='store_true',
        help="If this option is specified, it is assumed that reversible reaction edges are contained in the provided network tsv file in two directions (A to B, B to A). False by default."
    )
    parser_annotation_table_id_extender.add_argument(
        "-a",
        "--path_inputfile_node_annotations",
        required=True,
        type=str,
        help="The path to a tsv file containing node annotations."
    )
    parser_annotation_table_id_extender.add_argument(
        "-i",
        "--column_index_ids_annotation_inputfile",
        required=True,
        type = int,
        help = "The index of the column containing short or same node ids that are in the network. (e.g. network_id = P06675_pathway1, annotation_id = P06675)"
    )
    parser_annotation_table_id_extender.add_argument(
        'o',
        '--output_directory_and_filename',
        type=str,
        help="The path to the output tsv file that will be created (includes the name of the file)."
    )

def get_parser():
    """
    This method creates the full parser.
    :return: the parser
    """
    # create the top-level parser
    parser = argparse.ArgumentParser(
        prog="MONIT",
        description="Program which helps (i) select genes that potentially drive the observed multi-omics changes and (ii) visualize the results."
    )
    parser.add_argument("-v","--version", action="version", version="MONIT 1.0")

    # create the subparsers
    subparsers = parser.add_subparsers(
        dest="subcommand",
        required=True)

        #list of subparsers
    subparser_KEGG_network_constructor(subparsers)
    subparser_network_merger(subparsers)
    subparser_node_impact_assessor(subparsers)
    subparser_subnetwork_table_constructor(subparsers)
    subparser_cytoscape_node_table_nodeshortid_nodetype_extension_constructor(subparsers)
    subparser_cytoscape_node_table_general_extension_constructor(subparsers)

    return parser

