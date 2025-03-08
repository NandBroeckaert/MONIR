from KEGG_network_construction import (KEGG_organism_all_metabolic_pathway_retriever,KEGG_network_constructor_from_string,network_merger)
import time
from Subnetwork_selector_for_visualisation import (subnetwork_table_constructor,cytoscape_node_table_nodeshortid_nodetype_extension_constructor,cytoscape_node_table_general_extension_constructor)
from Impact_score_calculation import (general_node_impact_assessor)

#MONIR paper -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
selected_pathways_combinations = [
    "pae00010",
    "pae00010$pae00020",
    "pae00010$pae00020$pae00030",
    "pae00010$pae00020$pae00030$pae00040$pae00051$pae00052$pae00053$pae00500$pae00520$pae00541$pae00620$pae00630$pae00640$pae00650$pae00660$pae00562",
    "pae00010$pae00020$pae00030$pae00040$pae00051$pae00052$pae00053$pae00500$pae00520$pae00541$pae00620$pae00630$pae00640$pae00650$pae00660$pae00562$pae00190$pae00710$pae00720$pae00680$pae00910$pae00920",
    "pae00010$pae00020$pae00030$pae00040$pae00051$pae00052$pae00053$pae00500$pae00520$pae00541$pae00620$pae00630$pae00640$pae00650$pae00660$pae00562$pae00190$pae00710$pae00720$pae00680$pae00910$pae00920$pae00061$pae00071$pae00074$pae00561$pae00564$pae00565$pae00600$pae00592$pae01040",
    "pae00010$pae00020$pae00030$pae00040$pae00051$pae00052$pae00053$pae00500$pae00520$pae00541$pae00620$pae00630$pae00640$pae00650$pae00660$pae00562$pae00190$pae00710$pae00720$pae00680$pae00910$pae00920$pae00061$pae00071$pae00074$pae00561$pae00564$pae00565$pae00600$pae00592$pae01040$pae00230$pae00240",
    "pae00010$pae00020$pae00030$pae00040$pae00051$pae00052$pae00053$pae00500$pae00520$pae00541$pae00620$pae00630$pae00640$pae00650$pae00660$pae00562$pae00190$pae00710$pae00720$pae00680$pae00910$pae00920$pae00061$pae00071$pae00074$pae00561$pae00564$pae00565$pae00600$pae00592$pae01040$pae00230$pae00240$pae00250$pae00260$pae00270$pae00280$pae00290$pae00300$pae00310$pae00220$pae00330$pae00340$pae00350$pae00360$pae00380$pae00400",
    "pae00010$pae00020$pae00030$pae00040$pae00051$pae00052$pae00053$pae00500$pae00520$pae00541$pae00620$pae00630$pae00640$pae00650$pae00660$pae00562$pae00190$pae00710$pae00720$pae00680$pae00910$pae00920$pae00061$pae00071$pae00074$pae00561$pae00564$pae00565$pae00600$pae00592$pae01040$pae00230$pae00240$pae00250$pae00260$pae00270$pae00280$pae00290$pae00300$pae00310$pae00220$pae00330$pae00340$pae00350$pae00360$pae00380$pae00400$pae00410$pae00430$pae00440$pae00450$pae00460$pae00470$pae00480",
    "pae00010$pae00020$pae00030$pae00040$pae00051$pae00052$pae00053$pae00500$pae00520$pae00541$pae00620$pae00630$pae00640$pae00650$pae00660$pae00562$pae00190$pae00710$pae00720$pae00680$pae00910$pae00920$pae00061$pae00071$pae00074$pae00561$pae00564$pae00565$pae00600$pae00592$pae01040$pae00230$pae00240$pae00250$pae00260$pae00270$pae00280$pae00290$pae00300$pae00310$pae00220$pae00330$pae00340$pae00350$pae00360$pae00380$pae00400$pae00410$pae00430$pae00440$pae00450$pae00460$pae00470$pae00480$pae00540$pae00550$pae00552$pae00543$pae00730$pae00740$pae00750$pae00760$pae00770$pae00780$pae00785$pae00790$pae00670$pae00860$pae00130",
    "pae00010$pae00020$pae00030$pae00040$pae00051$pae00052$pae00053$pae00500$pae00520$pae00541$pae00620$pae00630$pae00640$pae00650$pae00660$pae00562$pae00190$pae00710$pae00720$pae00680$pae00910$pae00920$pae00061$pae00071$pae00074$pae00561$pae00564$pae00565$pae00600$pae00592$pae01040$pae00230$pae00240$pae00250$pae00260$pae00270$pae00280$pae00290$pae00300$pae00310$pae00220$pae00330$pae00340$pae00350$pae00360$pae00380$pae00400$pae00410$pae00430$pae00440$pae00450$pae00460$pae00470$pae00480$pae00540$pae00550$pae00552$pae00543$pae00730$pae00740$pae00750$pae00760$pae00770$pae00780$pae00785$pae00790$pae00670$pae00860$pae00130$pae00900$pae00903$pae00907$pae00523$pae01053$pae00946$pae00332$pae00261$pae00521$pae00525$pae00401$pae00405$pae00999$pae00997$pae00362$pae00627$pae00364$pae00625$pae00361$pae00623$pae00622$pae00633$pae00643$pae00930$pae00626",
    "pae00010$pae00020$pae00030$pae00040$pae00051$pae00052$pae00053$pae00500$pae00520$pae00541$pae00620$pae00630$pae00640$pae00650$pae00660$pae00562$pae00190$pae00710$pae00720$pae00680$pae00910$pae00920$pae00061$pae00071$pae00074$pae00561$pae00564$pae00565$pae00600$pae00592$pae01040$pae00230$pae00240$pae00250$pae00260$pae00270$pae00280$pae00290$pae00300$pae00310$pae00220$pae00330$pae00340$pae00350$pae00360$pae00380$pae00400$pae00410$pae00430$pae00440$pae00450$pae00460$pae00470$pae00480$pae00540$pae00550$pae00552$pae00543$pae00730$pae00740$pae00750$pae00760$pae00770$pae00780$pae00785$pae00790$pae00670$pae00860$pae00130$pae00900$pae00903$pae00907$pae00523$pae01053$pae00946$pae00332$pae00261$pae00521$pae00525$pae00401$pae00405$pae00999$pae00997$pae00362$pae00627$pae00364$pae00625$pae00361$pae00623$pae00622$pae00633$pae00643$pae00930$pae00626$pae03020$pae03010$pae00970$pae03008$pae03060$pae04122$pae03018$pae03030$pae03410$pae03420$pae03430$pae03440$pae03450$pae03250$pae02010$pae02060$pae03070$pae02020$pae04146$pae04148$pae02024$pae02025$pae02030$pae02040$pae04981$pae04980$pae01501$pae01502$pae01503",
]
i = 1
for pathway_combination in selected_pathways_combinations:
    start_time = time.time()

    KEGG_network_constructor_from_string(pathway_combination,"pae","reaction$GErel$PPrel$PCrel","C:/test/runtime_network_"+str(i)+".tsv")

    stop_time = time.time()
    elapsed_time = stop_time - start_time
    print("combination "+str(i)+":"+str(elapsed_time))

    i += 1
"""


#MONIR paper -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#visualsatie volledig netwerk
network = "C:/Nand_phd/data_en_analyse/MONIR_testcase_05032025/MONIR_NETWORK_merged_pae_KEGG_RegulomePA_30012025.tsv"
cytoscape_node_table_nodeshortid_nodetype_extension_constructor(network,"C:/Nand_phd/data_en_analyse/MONIR_testcase_05032025/network_node_table_extension.tsv")

#annotation_table = "C:/Nand_phd/data_en_analyse/MONIR_testcase_05032025/MONIR_annotation_table.tsv"
#cytoscape_node_table_general_extension_constructor(network,False,annotation_table,0,"C:/Nand_phd/data_en_analyse/MONIR_testcase_05032025/MONIR_annotation_table_extended.tsv")

#analyse

OMICS_T05 = "C:/Nand_phd/data_en_analyse/MONIR_testcase_05032025/MONIR_OMICS_T05.tsv"
OMICS_T10 = "C:/Nand_phd/data_en_analyse/MONIR_testcase_05032025/MONIR_OMICS_T10.tsv"
OMICS_T15 = "C:/Nand_phd/data_en_analyse/MONIR_testcase_05032025/MONIR_OMICS_T15.tsv"

OMICS_list = [OMICS_T05,OMICS_T10,OMICS_T15]

NOI_T05 = "C:/Nand_phd/data_en_analyse/MONIR_testcase_05032025/MONIR_NOI_T05.tsv"
NOI_T10 = "C:/Nand_phd/data_en_analyse/MONIR_testcase_05032025/MONIR_NOI_T10.tsv"
NOI_T15 = "C:/Nand_phd/data_en_analyse/MONIR_testcase_05032025/MONIR_NOI_T15.tsv"
NOI_list = [NOI_T05,NOI_T10,NOI_T15]

timepoints = ["T05","T10","T15"]
for i in range(len(timepoints)):
    OMICS = OMICS_list[i]
    NOI = NOI_list[i]
    timepoint = timepoints[i]
    print(OMICS)
    print(NOI)

    general_node_impact_assessor(network,False,OMICS,NOI,'bidirectional','unidirectional',True,10.0,True,True,1.5,True,1.0,"C:/Nand_phd/data_en_analyse/MONIR_testcase_05032025/MONIR_IMPACT_ANALYSIS_LUZ19transcription_"+timepoint+'_05032025.tsv',10000)

#visualisate subnetwerken


"""
#gp13 paper -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#network_merger("C:/test/runtime_network_12.tsv","C:/test/RegulomePAdb_allinteractions_15012022_reformatted.tsv",1,False,"C:/test/merged_pae_KEGG_RegulomePA_30012025.tsv")
network = "C:/Nand_phd/data_en_analyse/gp13_paper_27012025/MONIR_NETWORK_merged_pae_KEGG_RegulomePA_30012025.tsv"
OMICS_T000 = "C:/Nand_phd/data_en_analyse/gp13_paper_27012025/MONIR_OMICS_metabolomics_CTRLvsGP13_TP000_30012025.txt"
OMICS_T015 = "C:/Nand_phd/data_en_analyse/gp13_paper_27012025/MONIR_OMICS_metabolomics_CTRLvsGP13_TP015_30012025.txt"
OMICS_T030 = "C:/Nand_phd/data_en_analyse/gp13_paper_27012025/MONIR_OMICS_metabolomics_CTRLvsGP13_TP030_30012025.txt"
OMICS_T045 = "C:/Nand_phd/data_en_analyse/gp13_paper_27012025/MONIR_OMICS_metabolomics_CTRLvsGP13_TP045_30012025.txt"
OMICS_T060 = "C:/Nand_phd/data_en_analyse/gp13_paper_27012025/MONIR_OMICS_metabolomics_CTRLvsGP13_TP060_30012025.txt"
OMICS_T090 = "C:/Nand_phd/data_en_analyse/gp13_paper_27012025/MONIR_OMICS_metabolomics_CTRLvsGP13_TP090_30012025.txt"
OMICS_T120 = "C:/Nand_phd/data_en_analyse/gp13_paper_27012025/MONIR_OMICS_metabolomics_CTRLvsGP13_TP120_30012025.txt"
OMICS_Tall = "C:/Nand_phd/data_en_analyse/gp13_paper_27012025/MONIR_OMICS_metabolomics_CTRLvsGP13_allTP_30012025.txt"
OMICS_list = [OMICS_T000,OMICS_T015,OMICS_T030,OMICS_T045,OMICS_T060,OMICS_T090,OMICS_T120,OMICS_Tall]
NOI = "C:/Nand_phd/data_en_analyse/gp13_paper_27012025/MONIR_NOI_acetylomics_gp13vsCTRL_SigProt_LP75_PEP05_SIG05_28012025.txt"
timepoints = ["TP000","TP015","TP030","TP045","TP060","TP090","T120","Tall"]
i = 0
for OMICS in OMICS_list:
    general_node_impact_assessor(network,False,OMICS,NOI,'bidirectional','unidirectional',True,10.0,True,True,1.0,True,0.5,"C:/Nand_phd/data_en_analyse/gp13_paper_27012025/MONIR_IMPACT_ANALYSIS_gp13vCTRL_"+timepoints[i]+'_31012025.txt',10000)
    i +=1

"""






#testing random--------------------------------------------------------------------------------------------------------------------------------------------------------------------
#metabolic_pathways_pae = KEGG_organism_all_metabolic_pathway_retriever("pae")
#non_metabolic_pathways_pae = "pae03020$pae03010$pae00970$pae03060$pae04122$pae03018$pae03030$pae03410$pae03420$pae03430$pae03440$pae03450$pae03250$pae02010$pae02060$pae03070$pae02020$pae04148$pae02024$pae02025$pae02030$pae02040$pae04980$pae01501$pae01502$pae01503"
#selected_pathways = metabolic_pathways_pae + "$" + non_metabolic_pathways_pae
#KEGG_network_constructor_from_string(selected_pathways,"pae","reaction$GErel$PPrel$PCrel","C:/test/pae_reaction_GErel_PPrel_PCrel_all_pathways_02122024.tsv")
#cytoscape_node_table_nodeshortid_nodetype_extension_constructor("C:/test/pae_reaction_GErel_PPrel_PCrel_all_pathways_02102024.tsv","C:/test/pae_reaction_GErel_PPrel_PCrel_all_pathways_node_table_extension_02102024.tsv")
#network_merger("C:/test/pae_reaction_GErel_PPrel_PCrel_all_pathways_02122024.tsv","C:/test/RegulomePAdb_allinteractions_15012022_reformatted.tsv",1,False,"C:/test/merged_pae_KEGG_RegulomePA_02122024.tsv")
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

"""
# analyzing acetylomics wih metabolomics and proteomics  -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#overview paths

network = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/07102024/1_Network_creation/merged_pae_KEGG_RegulomePA.tsv"

PEV2_NOI ="C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/07102024/0_Input_data_preprocessing/NOI/PEV2_MONIT_NOI_KEGG_format_07102024.txt"
PEV2_omics = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/07102024/0_Input_data_preprocessing/OMICS/PEV2_MONIT_omics_format_07102024.txt"

output_directory = "C:/test/"


#14-1
general_node_impact_assessor(network,False,ph141_omics,ph141_NOI,'bidirectional','unidirectional',True,10.0,True,True,0.5,True,0.5,output_directory+'14_1_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',10000)
#LUZ19
general_node_impact_assessor(network,False,LUZ19_omics,LUZ19_NOI,'bidirectional','unidirectional',True,10.0,True,True,0.5,True,0.5,output_directory+'LUZ19_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',10000)

#PEV2
#general_node_impact_assessor(network,False,PEV2_omics,PEV2_NOI,'bidirectional','unidirectional',True,10.0,True,True,0.5,True,0.5,output_directory+'PEV2_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',10000)

#PhiKZ
general_node_impact_assessor(network,False,PhiKZ_omics,PhiKZ_NOI,'bidirectional','unidirectional',True,10.0,True,True,0.5,True,0.5,output_directory+'PhiKZ_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',10000)
#YuA
general_node_impact_assessor(network,False,YuA_omics,YuA_NOI,'bidirectional','unidirectional',True,10.0,True,True,0.5,True,0.5,output_directory+'YuA_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',10000)

#generating subnetworks
output_directory_subnetworks = "C:/Nand_phd/data_en_analyse/Phage_acetylomics_multi_omics_impact_analysis_2024/3_Subnetwork_analysis/"


#14_1
subnetwork_table_constructor(network,False,"bidirectional","bidirectional",output_directory+'14_1_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',"total_impact_score",10,5,False,False,output_directory_subnetworks+'141/141_5_08102024')
#LUZ19
subnetwork_table_constructor(network,False,"bidirectional","bidirectional",output_directory+'LUZ19_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',"total_impact_score",10,5,False,False,output_directory_subnetworks+'LUZ19/LUZ19_5_08102024')
#PEV2
subnetwork_table_constructor(network,False,"bidirectional","bidirectional",output_directory+'PEV2_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',"total_impact_score",10,5,False,False,output_directory_subnetworks+'PEV2/PEV2_5_08102024')
#phiKZ
subnetwork_table_constructor(network,False,"bidirectional","bidirectional",output_directory+'PhiKZ_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',"total_impact_score",10,5,False,False,output_directory_subnetworks+'PhiKZ/PhiKZ_5_08102024')
#YuA
subnetwork_table_constructor(network,False,"bidirectional","bidirectional",output_directory+'YuA_impact_NOI_acetylomcs_OMICS_metabolomics_07102024.txt',"total_impact_score",10,5,False,False,output_directory_subnetworks+'YuA/YuA_5_08102024')
"""

#merge proteomics and metabolomics
#import pandas as pd

#omics_141 = pd.read_csv("C:/Users/nandb/Downloads/14_1_MONIT_omics_format_07102024.txt",sep="\t")
#omics_LUZ19 = pd.read_csv("C:/Users/nandb/Downloads/LUZ19_MONIT_omics_format_07102024.txt",sep="\t")
#omics_PEV2 = pd.read_csv("C:/Users/nandb/Downloads/PEV2_MONIT_omics_format_07102024.txt",sep="\t")
#omics_PhiKZ = pd.read_csv("C:/Users/nandb/Downloads/PhiKZ_MONIT_omics_format_07102024.txt",sep="\t")
#omics_YUA = pd.read_csv("C:/Users/nandb/Downloads/YUA_MONIT_omics_format_07102024.txt",sep="\t")

#merged_omics = pd.merge(omics_141,omics_LUZ19,on = 'node_id',how = 'outer')
#merged_omics = pd.merge(omics_PEV2,merged_omics,on = 'node_id',how = 'outer')
#merged_omics = pd.merge(omics_PhiKZ,merged_omics,on = 'node_id',how = 'outer')
#merged_omics = pd.merge(omics_YUA,merged_omics,on = 'node_id',how = 'outer')


#merged_omics.to_csv("C:/Users/nandb/Downloads/merged_omics_node_table",sep="\t",index=False)