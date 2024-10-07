from KEGG_network_construction import (KEGG_organism_all_metabolic_pathway_retriever,KEGG_network_constructor_from_string,network_merger)
from Subnetwork_selector_for_visualisation import (subnetwork_table_constructor,cytoscape_node_table_nodeshortid_nodetype_extension_constructor,cytoscape_node_table_general_extension_constructor)
from Impact_score_calculation import (general_node_impact_assessor)

#testing -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#metabolic_pathways_pae = KEGG_organism_all_metabolic_pathway_retriever("pae")
#non_metabolic_pathways_pae = "pae03020$pae03010$pae00970$pae03060$pae04122$pae03018$pae03030$pae03410$pae03420$pae03430$pae03440$pae03450$pae03250$pae02010$pae02060$pae03070$pae02020$pae04148$pae02024$pae02025$pae02030$pae02040$pae04980$pae01501$pae01502$pae01503"
#selected_pathways = metabolic_pathways_pae + "$" + non_metabolic_pathways_pae
#KEGG_network_constructor_from_string(selected_pathways,"pae","reaction$GErel$PPrel$PCrel","C:/test/pae_reaction_GErel_PPrel_PCrel_all_pathways_02102024.tsv")
#cytoscape_node_table_nodeshortid_nodetype_extension_constructor("C:/test/pae_reaction_GErel_PPrel_PCrel_all_pathways_02102024.tsv","C:/test/pae_reaction_GErel_PPrel_PCrel_all_pathways_node_table_extension_02102024.tsv")
#network_merger("C:/test/pae_reaction_GErel_PPrel_PCrel_all_pathways_02102024.tsv","C:/test/RegulomePAdb_allinteractions_15012022_reformatted.tsv",1,False,"C:/test/merged_pae_KEGG_RegulomePA.tsv")
#cytoscape_node_table_nodeshortid_nodetype_extension_constructor("C:/test/merged_pae_KEGG_RegulomePA.tsv","C:/test/merged_pae_KEGG_RegulomePA_node_table_extension_02102024.tsv")

#KEGG_network_constructor_from_string("pae00020","pae","reaction","C:/test/reaction_TCA.txt")
#KEGG_network_constructor_from_string("pae00020","pae","reaction","C:/test/reaction_TCA_doubled.txt",reverse_interaction_doubler=True)
#KEGG_network_constructor_from_string("pae00010$pae00020","pae","reaction","C:/test/glycolysis_and_TCA.txt")
#KEGG_network_constructor_from_string("pae00010$pae00020","pae","reaction","C:/test/glycolysis_and_TCA_doubled.txt",reverse_interaction_doubler=True)
#KEGG_network_constructor_from_string("pae00010","pae","reaction","C:/test/reaction_glycolysis.txt")
#KEGG_network_constructor_from_string("pae00010","pae","reaction","C:/test/reaction_glycolysis_doubled.txt",reverse_interaction_doubler=True)

#KEGG_network_constructor_from_string(metabolic_pathways_pae,"pae","reaction","C:/test/reaction_all_metabolic_pae_doubled.txt",reverse_interaction_doubler=True)


#KEGG_network_constructor_from_string("pae00010","pae","chemical","C:/test/chemical_glycolysis.txt")
#KEGG_network_constructor_from_string("pae00010","pae","chemical","C:/test/chemical_glycolysis_doubled.txt",reverse_interaction_doubler=True)
#KEGG_network_constructor_from_string(metabolic_pathways_pae,"pae","chemical","C:/test/chemical_all_metabolic_pae.txt")
#KEGG_network_constructor_from_string(metabolic_pathways_pae,"pae","chemical","C:/test/chemical_all_metabolic_pae_doubled.txt",reverse_interaction_doubler=True)

#KEGG_network_constructor_from_string("pae00010","pae","ECrel","C:/test/ECrel_glycolysis_doubled.txt",reverse_interaction_doubler=True)

#merging networks
#glycolysis and TCA
#network_merger("C:/test/reaction_glycolysis_FALSE.txt","C:/test/reaction_glycolysis.txt",1,False,"C:/test/merged_glycolysis_reaction_02102024.txt")

#glycolysis and (glycolysis and TCA)
#network_merger("C:/test/glycolysis.txt","C:/test/TCA.txt",False,1,False,"C:/test/merged_glycolysis_and_glycolysis_and_TCA.txt")


# analyzing acetylomics wih metabolomics and proteomics  -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#overview paths
network = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/1_Network_creation/02102024_run1/merged_pae_KEGG_RegulomePA.tsv"

ph141_NOI = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/0_Input_data_preprocessing/NOI/14_1_MONIT_NOI_KEGG_format_07102024.txt"
ph141_omics = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/0_Input_data_preprocessing/OMICS/14_1_MONIT_omics_format_07102024.txt"

LUZ19_NOI = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/0_Input_data_preprocessing/NOI/LUZ19_MONIT_NOI_KEGG_format_07102024.txt"
LUZ19_omics = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/0_Input_data_preprocessing/OMICS/LUZ19_MONIT_omics_format_07102024.txt"

PEV2_NOI ="C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/0_Input_data_preprocessing/NOI/PEV2_MONIT_NOI_KEGG_format_07102024.txt"
PEV2_omics = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/0_Input_data_preprocessing/OMICS/PEV2_MONIT_omics_format_07102024.txt"

PhiKZ_NOI = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/0_Input_data_preprocessing/NOI/PhiKZ_MONIT_NOI_KEGG_format_07102024.txt"
PhiKZ_omics = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/0_Input_data_preprocessing/OMICS/PhiKZ_MONIT_omics_format_07102024.txt"

YuA_NOI = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/0_Input_data_preprocessing/NOI/YuA_MONIT_NOI_KEGG_format_07102024.txt"
YuA_omics = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/0_Input_data_preprocessing/OMICS/YUA_MONIT_omics_format_07102024.txt"

output_directory = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/2_Results_impact_analysis/"



#14-1
general_node_impact_assessor(network,False,ph141_omics,ph141_NOI,'bidirectional','unidirectional',True,10.0,True,True,0.5,True,0.5,output_directory+'14_1_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',10000)
#LUZ19
general_node_impact_assessor(network,False,LUZ19_omics,LUZ19_NOI,'bidirectional','unidirectional',True,10.0,True,True,0.5,True,0.5,output_directory+'LUZ19_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',10000)
#PEV2
general_node_impact_assessor(network,False,PEV2_omics,PEV2_NOI,'bidirectional','unidirectional',True,10.0,True,True,0.5,True,0.5,output_directory+'PEV2_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',10000)
#PhiKZ
general_node_impact_assessor(network,False,PhiKZ_omics,PhiKZ_NOI,'bidirectional','unidirectional',True,10.0,True,True,0.5,True,0.5,output_directory+'PhiKZ_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',10000)
#YuA
general_node_impact_assessor(network,False,YuA_omics,YuA_NOI,'bidirectional','unidirectional',True,10.0,True,True,0.5,True,0.5,output_directory+'YuA_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',10000)