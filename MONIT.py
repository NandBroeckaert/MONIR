import argparse
import KEGG_network_construction
import Impact_score_calculation
import Subnetwork_selector_for_visualisation


"""
to be added:
KEGG_organism_all_metabolic_pathway_retriever
network_merger
general_node_impact_assessor
subnetwork_table_constructor
cytoscape_node_table_extension_constructor
"""
#sub-command functions



#create the top-level parser
parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(required = True)

#create the parser for the KEGG_organism_all_metabolic_pathway_retriever method
parser_KEGG_organism_all_metabolic_pathway_retriever = subparsers.add_parser('KEGG_organism_all_metabolic_pathway_retriever',help='This method creates a string ')

#create the parser for the KEGG_network_costructor_from_string method
parser_KEGG_network_constructor_from_string = subparsers.add_parser('KEGG_network_constructor_from_string',help="This method creates a '.tsv' file containing a network of the specified KEGG pathways.")
parser_KEGG_network_constructor_from_string.add_argument('KEGG_pathways',help='A string of KEGG pathway ids seperated by a $ sign. (example: pae00010)',type = str)
parser_KEGG_network_constructor_from_string.add_argument('organism_code',help='The organism code from the KEGG database. (example: pae)',type = str)
parser_KEGG_network_constructor_from_string.add_argument('KEGG_interaction_types',help='A string of possible KEGG interaction types that you want to include. Each type is seperated from the next by a $ sign. Possible interaction type: chemical,reaction,ECrel,PPrel,GErel,PCrel. The chemical and reaction types cannot be selected together. (example: reaction$GErel',type = str)
parser_KEGG_network_constructor_from_string.add_argument('-r','--reverse_interaction_doubler',help='Set as True if you want reversible reactions to be specified in both directions in the network table.',nargs='?',default=False,choices=[True,False],type=bool)
    #connect it up to function

