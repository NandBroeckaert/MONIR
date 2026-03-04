import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ..Network_construction import (KEGG_organism_all_metabolic_pathway_retriever,KEGG_network_constructor_from_string,network_merger)
from ..Subnetwork_visualisation import (subnetwork_table_constructor,cytoscape_node_table_nodeshortid_nodetype_extension_constructor,cytoscape_node_table_general_extension_constructor)
from ..Impact_analysis import (general_node_impact_assessor)


test_directory_input = os.path.dirname(os.path.abspath(__file__)) + "/Test_input/"
test_directory_output = os.path.dirname(os.path.abspath(__file__)) + "/Test_output/"

def main():
    #NETWORK CONSTRUCTION -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print()
    print("=" * 70)
    print("TEST: network construction module")
    print("=" * 70)

    print()
    print("Testing KEGG_network_constructor_from_string method")
    #make small networks for glycolisis and TCA cycle
    #TCA
    print("  - TCA cycle")
    KEGG_network_constructor_from_string("pae00020","pae","reaction",test_directory_output + "reaction_TCA.txt")
    KEGG_network_constructor_from_string("pae00020","pae","reaction",test_directory_output + "reaction_TCA_doubled.txt",reverse_interaction_doubler=True)

    #glycolysis
    print("  - Glycolysis")
    KEGG_network_constructor_from_string("pae00010","pae","reaction",test_directory_output + "reaction_glycolysis.txt")
    KEGG_network_constructor_from_string("pae00010","pae","reaction",test_directory_output + "reaction_glycolysis_doubled.txt",reverse_interaction_doubler=True)

    #glycolysis and TCA
    print("  - Glycolysis and TCA cycle")
    KEGG_network_constructor_from_string("pae00010$pae00020","pae","reaction",test_directory_output + "glycolysis_and_TCA.txt")
    KEGG_network_constructor_from_string("pae00010$pae00020","pae","reaction",test_directory_output + "glycolysis_and_TCA_doubled.txt",reverse_interaction_doubler=True)

    #make a network with all metabolic pathways and all interaction types
    print("  - All metabolic pathways and all interaction types")
    metabolic_pathways_pae = KEGG_organism_all_metabolic_pathway_retriever("pae")
    non_metabolic_pathways_pae = "pae03020$pae03010$pae00970$pae03060$pae04122$pae03018$pae03030$pae03410$pae03420$pae03430$pae03440$pae03450$pae03250$pae02010$pae02060$pae03070$pae02020$pae04148$pae02024$pae02025$pae02030$pae02040$pae04980$pae01501$pae01502$pae01503"
    selected_pathways = metabolic_pathways_pae + "$" + non_metabolic_pathways_pae
    KEGG_network_constructor_from_string(selected_pathways,"pae","reaction$GErel$PPrel$PCrel",test_directory_output + "pae_reaction_GErel_PPrel_PCrel_all_pathways.tsv")

    print("Testing KEGG_network_constructor_from_string method completed")

    print()
    print("Testing network_merger method")
    print("  - Merging glycolysis and TCA cycle")
    #merge glycolysis and TCA
    network_merger(test_directory_output + "reaction_glycolysis.txt",test_directory_output + "reaction_TCA.txt",1,False,test_directory_output + "merged_glycolysis_and_glycolysis_and_TCA.txt")
    print("Testing network_merger method completed")

    #IMPACT ANALYSIS -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print()
    print("=" * 70)
    print("TEST: impact analysis module")
    print("=" * 70)

    path_network = test_directory_output + "pae_reaction_GErel_PPrel_PCrel_all_pathways.tsv"
    OMICS = test_directory_input + "MONIR_OMICS_T05.tsv"
    NOI = test_directory_input + "MONIR_NOI_T05.tsv"

    print()
    print("Testing general_node_impact_assessor method")
    general_node_impact_assessor(path_network,False,OMICS,NOI,'bidirectional','unidirectional',True,10.0,True,True,1.5,True,1.0,test_directory_output + "MONIR_IMPACT_ANALYSIS_T05.tsv",100)
    print("Testing general_node_impact_assessor method completed")

    #SUBNETWORK VISUALISATION -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print()
    print("=" * 70)
    print("TEST: subnetwork visualisation module")
    print("=" * 70)
    path_network = test_directory_output + "pae_reaction_GErel_PPrel_PCrel_all_pathways.tsv"

    print()
    print("Testing subnetwork_table_constructor method")
    subnetwork_table_constructor(path_network,False,"bidirectional","unidirectional",test_directory_output + "MONIR_IMPACT_ANALYSIS_T05.tsv",15.0,0.11,2,False,False,test_directory_output + "MONIR_IMPACT_ANALYSIS_T05_subnetwork.tsv")
    print("Testing subnetwork_table_constructor method completed")

    print()
    print("Testing cytoscape_node_table_general_extension_constructor method")
    cytoscape_node_table_general_extension_constructor(path_network,False,OMICS,0,test_directory_output + "MONIR_annotation_table_extended.tsv")
    print("Testing cytoscape_node_table_general_extension_constructor method completed")

    print()
    print("Testing cytoscape_node_table_nodeshortid_nodetype_extension_constructor method")
    cytoscape_node_table_nodeshortid_nodetype_extension_constructor(test_directory_output + "pae_reaction_GErel_PPrel_PCrel_all_pathways.tsv",test_directory_output + "pae_reaction_GErel_PPrel_PCrel_all_pathways_node_table_extension.tsv")
    print("Testing cytoscape_node_table_nodeshortid_nodetype_extension_constructor method completed")

    print()
    print("=" * 70)
    print("TESTING COMPLETED")
    print("=" * 70)
    print()
