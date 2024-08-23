import sys
import argparse
from parsers import get_parser
from KEGG_network_construction import (
    KEGG_organism_all_metabolic_pathway_retriever,
    KEGG_network_constructor_from_string,
    network_merger)
from Impact_score_calculation import general_node_impact_assessor
from Subnetwork_selector_for_visualisation import(
    subnetwork_table_constructor,
    cytoscape_node_table_nodeshortid_nodetype_extension_constructor,
    cytoscape_node_table_general_extension_constructor)

def main():
    """
    Entry point of the program.
    The commandline input will be given to the parser and the requested subcommand will be executed.
    """
    input = sys.argv[1:]
    parser = get_parser()
    args = parser.parse_args(input)
    execute_requested_subcommand(args)

def execute_requested_subcommand(args):
    """
    Based on the provided argparse arguments, a specific subcommand will be executed.
    :param args: argparse arguments
    """
    if args.subcommand == "KEGG_network_constructor":
        print("Building the requested KEGG network.")

        if args.pathways == "all_metabolic_reactions":
            all_metabolic_pathways = KEGG_organism_all_metabolic_pathway_retriever(args.organism_code)

            KEGG_network_constructor_from_string(
                KEGG_pathways_string = all_metabolic_pathways,
                organism_code = args.organism_code,
                selected_KEGG_interaction_types_string = args.selected_interaction_types,
                path_outputfile = args.path_outputfile,
                reverse_interaction_doubler = args.reverse_interactions)

        else:
            KEGG_network_constructor_from_string(
                KEGG_pathways_string = args.pathways,
                organism_code = args.organism_code,
                selected_KEGG_interaction_types_string = args.selected_interaction_types,
                path_outputfile = args.path_outputfile,
                reverse_interaction_doubler = args.reverse_interactions)


    elif args.subcommand == "network_merger":
        print("Merging the two specified networks.")
        network_merger()

    elif args.subcommand == "node_impact_assessor":
        print("Performing the impact analysis for the requested nodes.")
        general_node_impact_assessor(
            path_inputfile_network=args.path_inputfile_network,
            reverse_interaction_doubled=args.reverse_interactions,
            path_inputfile_node_omics_info=args.ath_inputfile_node_omics_info,
            path_inputfile_nodes_of_interest=args.path_inputfile_nodes_of_interest,
            directionality_reaction=args.directionality_reaction,
            directionality_other=args.directionality_other,
            interaction_specific=args.interaction_specific,
            start_weight_value=args.start_weight_value,
            centrality_modification=args.centrality_modification,
            missingness_modification=args.missingness_modification,
            missingness_modification_step_penalty=args.missingness_modification_step_penalty,
            distance_modification=args.distance_modification,
            distance_modification_step_penalty=args.distance_modification_step_penalty,
            path_output_directory_and_filename=args.path_output_directory_and_filename,
            distance_level_limit=args.distance_level_limit,
            include_subnetwork_for_visualisation=args.include_subnetwork_for_visualisation)

    elif args.subcommand == "subnetwork_constructor":
        print("Generating files containing subnetworks for downstream visualisation.")
        subnetwork_table_constructor(
            path_inputfile_network=args.path_inputfile_network,
            reverse_interaction_doubled=args.reverse_interactions,
            directionality_reaction=args.directionality_reaction,
            directionality_other=args.directionality_other,
            path_inputfile_results_impact_analysis=args.path_inputfile_results_impact_analysis,
            impact_value_type=args.threshold_impact_value_type,
            impact_threshold=args.impact_threshold_value,
            distance_level_limit=args.distance_level_limit,
            distance_level_limit_based_on_impact_results=args.distance_level_limit_based_on_impact_results,
            incl_neighbours=args.incl_neighbours,
            output_directory_and_filename=args.output_directory_and_filename,)

    elif args.subcommand == "node_table_id_and_type_extender":
        print("Generating a new file containing with all the node ids and types in the network.")
        cytoscape_node_table_nodeshortid_nodetype_extension_constructor(
            path_inputfile_network=args.path_inputfile_network,
            output_directory_and_filename=args.output_directory_and_filename)

    elif args.subcommand == "annotation_table_id_extender":
        print("Generating a new file with an additional node id column.")
        cytoscape_node_table_general_extension_constructor(
            path_inputfile_network=args.path_inputfile_network,
            reverse_interaction_doubled=args.reverse_interactions,
            path_inputfile_node_annotations=args.path_inputfile_node_annotations,
            column_index_ids_annotation_inputfile=args.column_index_ids_annotation_inputfile,
            output_directory_and_filename=args.output_directory_and_filename)

if __name__ == '__main__':
    main()