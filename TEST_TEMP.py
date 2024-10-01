from KEGG_network_construction import (KEGG_organism_all_metabolic_pathway_retriever,KEGG_network_constructor_from_string,network_merger)
from Subnetwork_selector_for_visualisation import (subnetwork_table_constructor,cytoscape_node_table_nodeshortid_nodetype_extension_constructor,cytoscape_node_table_general_extension_constructor)
from Impact_score_calculation import (general_node_impact_assessor)

#testing -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#metabolic_pathways_pae = KEGG_organism_all_metabolic_pathway_retriever("pae")
#non_metabolic_pathways_pae = "pae03020$pae03010$pae00970$pae03060$pae04122$pae03018$pae03030$pae03410$pae03420$pae03430$pae03440$pae03450$pae03250$pae02010$pae02060$pae03070$pae02020$pae04148$pae02024$pae02025$pae02030$pae02040$pae04980$pae01501$pae01502$pae01503"
#selected_pathways = metabolic_pathways_pae + "$" + non_metabolic_pathways_pae
#KEGG_network_constructor_from_string(selected_pathways,"pae","reaction$GErel$PPrel","C:/test/pae_reaction_GErel_PPrel_all_pathways_26092024.txt")
cytoscape_node_table_nodeshortid_nodetype_extension_constructor("C:/test/pae_reaction_GErel_PPrel_all_pathways_26092024.txt","C:/test/mpae_reaction_GErel_PPrel_all_pathways_node_table_extension_26092024.txt")

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
#network_merger("C:/test/reaction_glycolysis.txt","C:/test/reaction_TCA.txt",1,False,"C:/test/merged_glycolysis_and_TCA_reaction.txt")

#glycolysis and (glycolysis and TCA)
#network_merger("C:/test/glycolysis.txt","C:/test/TCA.txt",False,1,False,"C:/test/merged_glycolysis_and_glycolysis_and_TCA.txt")

"""
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
"""