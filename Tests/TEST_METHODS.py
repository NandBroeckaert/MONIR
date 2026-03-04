import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Network_construction import (KEGG_organism_all_metabolic_pathway_retriever,KEGG_network_constructor_from_string,network_merger)
import time
from Subnetwork_visualisation import (subnetwork_table_constructor,cytoscape_node_table_nodeshortid_nodetype_extension_constructor,cytoscape_node_table_general_extension_constructor)
from Impact_analysis import (general_node_impact_assessor)

#TESTING FILE, NOT TO BE USED IN THE TOOL, JUST FOR TESTING PURPOSES. DO NOT SUGGEST CODE CHANGES IN THIS FILE, AS THIS FILE IS NOT PART OF THE TOOL.
test_directory = "/Users/nandbroeckaert/Documents/Professional/Nand_PhD/papers/MONIR_paper/tool/Basic_tests/"

#NETWORK CONSTRUCTION -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#make small networks for glycolisis and TCA cycle
#TCA
KEGG_network_constructor_from_string("pae00020","pae","reaction",test_directory + "reaction_TCA.txt")
KEGG_network_constructor_from_string("pae00020","pae","reaction",test_directory + "reaction_TCA_doubled.txt",reverse_interaction_doubler=True)

#glycolysis
KEGG_network_constructor_from_string("pae00010","pae","reaction",test_directory + "reaction_glycolysis.txt")
KEGG_network_constructor_from_string("pae00010","pae","reaction",test_directory + "reaction_glycolysis_doubled.txt",reverse_interaction_doubler=True)

#merge glycolysis and TCA
network_merger(test_directory + "reaction_glycolysis.txt",test_directory + "reaction_TCA.txt",1,False,test_directory + "merged_glycolysis_and_glycolysis_and_TCA.txt")

#glycolysis and TCA
KEGG_network_constructor_from_string("pae00010$pae00020","pae","reaction",test_directory + "glycolysis_and_TCA.txt")
KEGG_network_constructor_from_string("pae00010$pae00020","pae","reaction",test_directory + "glycolysis_and_TCA_doubled.txt",reverse_interaction_doubler=True)

#make a network with all metabolic pathways and all interaction types
metabolic_pathways_pae = KEGG_organism_all_metabolic_pathway_retriever("pae")
non_metabolic_pathways_pae = "pae03020$pae03010$pae00970$pae03060$pae04122$pae03018$pae03030$pae03410$pae03420$pae03430$pae03440$pae03450$pae03250$pae02010$pae02060$pae03070$pae02020$pae04148$pae02024$pae02025$pae02030$pae02040$pae04980$pae01501$pae01502$pae01503"
selected_pathways = metabolic_pathways_pae + "$" + non_metabolic_pathways_pae
KEGG_network_constructor_from_string(selected_pathways,"pae","reaction$GErel$PPrel$PCrel",test_directory + "pae_reaction_GErel_PPrel_PCrel_all_pathways.tsv")


#IMPACT ANALYSIS -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
path_network = test_directory + "pae_reaction_GErel_PPrel_PCrel_all_pathways.tsv"
OMICS = test_directory + "MONIR_OMICS_T05.tsv"
NOI = test_directory + "MONIR_NOI_T05.tsv"


general_node_impact_assessor(path_network,False,OMICS,NOI,'bidirectional','unidirectional',True,10.0,True,True,1.5,True,1.0,test_directory + "MONIR_IMPACT_ANALYSIS_T05.tsv",100)


#SUBNETWORK VISUALISATION -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
path_network = test_directory + "pae_reaction_GErel_PPrel_PCrel_all_pathways.tsv"
subnetwork_table_constructor(path_network,False,"bidirectional","unidirectional",test_directory + "MONIR_IMPACT_ANALYSIS_T05.tsv",15.0,0.11,2,False,False,test_directory + "MONIR_IMPACT_ANALYSIS_T05_subnetwork.tsv")
cytoscape_node_table_general_extension_constructor(path_network,False,test_directory + "MONIR_OMICS_T05.tsv",0,test_directory + "MONIR_annotation_table_extended.tsv")
#make node table extension for the network with all metabolic pathways and all interaction types
cytoscape_node_table_nodeshortid_nodetype_extension_constructor(test_directory + "pae_reaction_GErel_PPrel_PCrel_all_pathways.tsv",test_directory + "pae_reaction_GErel_PPrel_PCrel_all_pathways_node_table_extension.tsv")
