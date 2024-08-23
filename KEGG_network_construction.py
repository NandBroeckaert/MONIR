import pandas as pd
import requests
import time
import io
import xml.etree.ElementTree as ET


def KEGG_organism_pathway_retriever(KEGG_organism_code:str) -> list:
    """
    This method returns a list containing all the pathway identifiers of the specfied organism.
    :param KEGG_organism_code: the code that KEGG uses for a specific organism (e.g. 'pae' for Pseudomonas aeruginosa PAO1)
    :return: a list containing the pathway identifiers of the specified organism
    """
    KEGG_API_url_pathway_list = "https://rest.kegg.jp/list/pathway/{organism_code}".format(organism_code=KEGG_organism_code)
    response = requests.get(KEGG_API_url_pathway_list, timeout=5) #retrieve data from KEGG database

    if response.status_code != 200:
        raise Exception("There was a problem connecting to the KEGG database. Please check the organism code in the KEGG database.")

    response_info = str(response.content)               #turn content into string
    response_info = response_info.split("\\")           #split string into separate elements
    response_info = response_info[0::2]                 #only retain uneven elements in list. These are the pathway ids
    response_info[0] = response_info[0].split("'")[1]   #remove 'b' from the first element
    response_info = response_info[:-1]                  #remove last element. This is also not a pathway id
    for i in range(1,len(response_info)):               #remove first letter 'n' from all elements (except first one, since it doesn't have it). This 'n' is not part of the organism code.
        response_info[i]=response_info[i][1:]

    if len(response_info)<1:
        raise Exception(
            "No data for this organism code could be retrieved. Please check whether the organism is in the KEGG database.")
    return response_info


def KEGG_organism_all_metabolic_pathway_retriever(KEGG_organism_code:str) -> str:
    """
    This method returns a string containing all the metabolic pathway identifiers of the specified organism separated by the delimiter '$'.
    :param KEGG_organism_code: the code that KEGG uses for a specific organism (e.g. 'pae' for Pseudomonas aeruginosa PAO1)
    :return: a string containing all the metabolic pathway identifiers of the specified organism separated by the delimiter '$'
    """
    #all hypothetically possible metabolic pathways
    Carbohydrate_metabolism_pathway_numbers = ['00010','00020','00030','00040','00051','00052','00053','00500','00520','00620','00630','00640','00650','00660','00562']
    Energy_metabolism_pathway_numbers = ['00190','00195','00196','00710','00720','00680','00910','00920']
    Lipid_metabolism_pathway_numbers = ['00061','00062','00071','00073','00100','00120','00121','00140','00561','00564','00565','00600','00590','00591','00592','01040']
    Nucleotide_metabolism_pathway_numbers = ['00230','00240']
    Amino_acid_metabolism_pathway_numbers = ['00250','00260','00270','00280','00290','00300','00310','00220','00330','00340','00350','00360','00380','00400']
    Metabolism_of_other_amino_acids_pathway_numbers = ['00410','00430','00440','00450','00460','00470','00480']
    Glycan_biosynthesis_and_metabolism_pathway_numbers = ['00510','00513','00512','00515','00514','00532','00534','00533','00531','00563','00601','00603','00604','00511','00540','00542','00541','00550','00552','00571','00572','00543']
    Metabolism_of_cofactors_and_vitamins_pathway_numbers = ['00730','00740','00750','00760','00770','00780','00785','00790','00670','00830','00860','00130']
    Metabolism_of_terpenoids_and_polyketides_pathway_numbers = ['00900','00902','00909','00904','00906','00905','00981','00908','00903','00907','01052','00522','01051','01059','01056','01057','00253','00523','01054','01053','01055']
    Biosynthesis_of_other_secondary_metabolites_pathway_numbers = ['00940','00945','00941','00944','00942','00943','00946','00901','00403','00950','00960','00996','00232','00965','00966','00402','0311','00332','00261','00331','00521','00524','00525','00401','00404','00333','00254','00998','00999','00997']
    Xenobiotics_biodegradation_and_metabolism_pathway_numbers = ['00362','00627','00364','00625','00361','00623','00622','00633','00642','00643','00791','00930','00363','00621','00626','00624','00365','00984','00980','00982','00983']

    all_possible_metabolic_pathway_numbers = Carbohydrate_metabolism_pathway_numbers+Energy_metabolism_pathway_numbers+Lipid_metabolism_pathway_numbers+Nucleotide_metabolism_pathway_numbers+Amino_acid_metabolism_pathway_numbers+Metabolism_of_other_amino_acids_pathway_numbers+Glycan_biosynthesis_and_metabolism_pathway_numbers+Metabolism_of_cofactors_and_vitamins_pathway_numbers+Metabolism_of_terpenoids_and_polyketides_pathway_numbers+Biosynthesis_of_other_secondary_metabolites_pathway_numbers+Xenobiotics_biodegradation_and_metabolism_pathway_numbers
    all_possible_metabolic_pathways_numbers = list(set(all_possible_metabolic_pathway_numbers))
    all_possible_metabolic_pathway_id_list = []
    for number in all_possible_metabolic_pathways_numbers:
        pathway_id = KEGG_organism_code+number
        all_possible_metabolic_pathway_id_list.append(pathway_id)

    #all pathways for that organism in KEGG database
    all_pathways_for_specified_organism = KEGG_organism_pathway_retriever(KEGG_organism_code)

    #filter out pathway numbers that aren't present for this organism
    all_metabolic_pathway_ids_for_this_organism = set(all_pathways_for_specified_organism).intersection(set(all_possible_metabolic_pathway_id_list))

    #join elements intro string seperated by '$'
    separator = "$"
    all_metabolic_pathway_ids_for_this_organism_string = separator.join(all_metabolic_pathway_ids_for_this_organism)
    return all_metabolic_pathway_ids_for_this_organism_string


def KEGG_pathway_list_decoder(KEGG_pathways: str,organism_code: str) -> list:
    """
    This method makes a list of KEGG pathway ids for a specific organism based on the KEGG_pathways string and checks whether the list contains valid pathway ids.
    :param KEGG_pathways: a string containing KEGG pathway ids of the specified organism separated by the delimiter '$'
    :param organism_code: the code that KEGG uses for a specific organism (e.g. 'pae' for Pseudomonas aeruginosa PAO1)
    :return: a list containing the KEGG pathway ids that are in the KEGG_pathways string
    """

    #remove spaces and split the string in a list of pathway ids
    KEGG_pathways.replace(" ","")
    list_pathways = KEGG_pathways.split("$")

    #remove duplicates
    list_pathways = list(set(list_pathways))

    #check whether the pathway ids contain "path:", are from the right organism & whether they are in the list of pathways that belong to that organism
    all_pathways = KEGG_organism_pathway_retriever(organism_code)
    for pathway in list_pathways:
        if ":" in pathway:
            raise Exception("Input contains pathway containing 'path:'. Please remove these ex")
        if pathway not in all_pathways:
            raise Exception("Input contains pathway that does not belong to the specified organism. Check for wrong pathway ids or pathways from other organisms.")

    return list_pathways


def KEGG_interacion_type_decoder(selected_KEGG_interaction_types_string:str) -> list:
    """
    This method makes a list of interaction types based on the selected_KEGG_interaction_types_string string & checks whether the specified types are valid.
    :param selected_KEGG_interaction_types_string: a string containing KEGG interaction types separated by the delimiter '$'.
        There are six possible interaction types:chemical,reaction,ECrel,PPrel,GErel,PCrel.
        The chemical type and reaction type can not be selected together.
        The ECrel type and reaction type can not be selected together.
    :return: a list containing the interaction types specified in the selected_KEGG_interaction_types_string
    """

    #remove spaces and split the string in a list of interaction types
    selected_KEGG_interaction_types_string.replace(" ","")
    interaction_types_list = selected_KEGG_interaction_types_string.split("$")

    # remove duplicates
    interaction_types_list = list(set((interaction_types_list)))

    # the KEGG database also has a maplink type, but this is not usefull for network construction, hence it is not included.
    possible_interaction_types = ["chemical","reaction","ECrel","PPrel","GErel","PCrel"]

    #check for invalid interaction types
    for interaction_type in interaction_types_list:
        if interaction_type not in possible_interaction_types:
            raise Exception("A non valid interaction type was selected. Please make a selection out of the following types: chemical, reaction, ECrel, PPrel, GErel, PCrel")

    #the chemical type and reaction type can not be selected together.
    if "chemical" in interaction_types_list and "reaction" in interaction_types_list:
        raise Exception("The chemical type and reaction type can not be selected together.")

    #the reaction type and ECrel type can not be selected together.
    if "ECrel" in interaction_types_list and "reaction" in interaction_types_list:
        raise Exception("The ECrel type and reaction type can not be selected together.")

    return interaction_types_list


def kgml_reader(path_to_pathway_kgml_file) -> list:
    """
    This method will construct a list of dictionaries containing all the information on KEGG entries, relations and reactions contained in a KGML file.

    :param path_to_pathway_kgml_file: the path to the kgml file that contains all the info about a certain pathway
    :return: A list containing the following dictionaries:
        0) entry_dict
            # format of entry_dict per entry type:
                # compound - id:[type (str),name (list)]
                # gene - id:[type (str),name (list),reaction (str)]
                # group - id:[type (str), name(list), components (list)]          the name of a group == ["undefined"]
        1) genes_reaction_dict
            # format:
                # entry_reaction_KEGG_id:entry_names (list of gene names)
        2) relation_dict
            # format:
                # id:[entry_id_1 (str), entry_id_2 (str), type (str), interaction_info (str)]
        3) reaction_dict
            # format:
                # reaction_KEGG_id:[reaction_type,reaction_substrate_ids,reaction_products_ids]
    """
    tree = ET.parse(path_to_pathway_kgml_file)
    root = tree.getroot()

    entry_dict = {}
    genes_reaction_dict = {}
    relation_dict = {}
    reaction_dict = {}

    relation_count = 0  # used to make unique id for relations  (kglm does not provide one)

    for child in root:

        if child.tag == "entry":
            entry_id = child.attrib["id"]

            entry_type = child.attrib["type"]
            entry_name = child.attrib["name"].split(" ")
            dict_value = [entry_type, entry_name]
            if entry_type == "gene" and "reaction" in child.attrib:
                entry_reaction_KEGG_id = child.attrib["reaction"] #note: don't load map pathways. There one gene entry may refer to several reactions.
                dict_value.append(entry_reaction_KEGG_id)
                genes_reaction_dict[entry_reaction_KEGG_id] = entry_name

            if entry_type == "group":
                for grandchild in child:
                    components = []
                    if grandchild.tag == "component":
                        components.append(grandchild.attrib["id"])

                dict_value.append(components)

            # format of entry_dict per entry type:
            # compound - id:[type (str),name (list)]
            # gene - id:[type (str),name (list),reaction (str)]
            # group - id:[type (str), name(list), components (list)]          the name of a group == ["undefined"]
            entry_dict[entry_id] = dict_value

        elif child.tag == "relation":
            relation_id = "relation_" + str(relation_count)  # make unique id for relations  (kglm does not provide one)

            relation_entry_id_1 = child.attrib["entry1"]
            relation_entry_id_2 = child.attrib["entry2"]
            relation_type = child.attrib["type"]

            relation_types_with_additional_interaction_information_list = ["ECrel", "PPrel", "GErel"]
            relation_interaction_info_list = []
            if relation_type in relation_types_with_additional_interaction_information_list:
                for grandchild in child:
                    if grandchild.tag == "subtype":
                        relation_info = grandchild.attrib["name"]
                        relation_info.replace(" ", "_")
                        relation_interaction_info_list.append(relation_info)
                separator = "$"
                relation_interaction_info_string = separator.join(relation_interaction_info_list)
            else:
                relation_interaction_info_string = "NA"

            relation_list = [relation_entry_id_1, relation_entry_id_2, relation_type, relation_interaction_info_string]
            # format of relation_dict:
            # id:[entry_id_1 (str), entry_id_2 (str), type (str), interaction_info (str)]
            relation_dict[relation_id] = relation_list
            relation_count += 1

        elif child.tag == "reaction":

            reaction_KEGG_id = child.attrib['name']
            reaction_info = child.attrib['type']  # reversible or irreversible
            reaction_substrate_entry_ids = []
            reaction_products_entry_ids = []
            for grandchild in child:
                if grandchild.tag == "substrate":
                    reaction_substrate_entry_ids.append(grandchild.attrib["id"])
                elif grandchild.tag == "product":
                    reaction_products_entry_ids.append(grandchild.attrib["id"])

            # format of reaction_dict:
            # reaction_KEGG_id:[reaction_type,reaction_substrate_ids,reaction_products_ids]
            reaction_dict[reaction_KEGG_id] = [reaction_info, reaction_substrate_entry_ids, reaction_products_entry_ids]

    # return a list containing the different dictionaries
    return [entry_dict, genes_reaction_dict, relation_dict, reaction_dict]


def network_builder(kgml_information_list:list,organism_code: str,selected_interaction_types: list, reverse_interaction_doubler: bool) -> pd.DataFrame:
    """
    This method uses the dictionaries made by kgml_reader() to build a pandas dataframe that contains the information of a pathway in a network table format.

    :param kgml_information_list: A list containing info about a specific pathway/kgml files. This info is stored in the following dictionaries:
        Element 0: entry_dict
            # format of entry_dict per entry type:
                # compound - id:[type (str),name (list)]
                # gene - id:[type (str),name (list),reaction (str)]
                # group - id:[type (str), name(list), components (list)]          the name of a group == ["undefined"]
        Element 1: genes_reaction_dict
            # format:
                # entry_reaction_KEGG_id:entry_names (list of gene names)
        Element 2: relation_dict
            # format:
                # id:[entry_id_1 (str), entry_id_2 (str), type (str), interaction_info (str)]
        Element 3: reaction_dict
            # format:
                # reaction_KEGG_id:[reaction_type,reaction_substrate_ids,reaction_products_ids]
    :param organism_code: the code that KEGG uses for a specific organism (e.g. 'pae' for Pseudomonas aeruginosa PAO1)
    :param selected_interaction_types: a list of all the interaction types that need to be included in the network pandas dataframe.
        This network may contain the following types of interactions: chemical, reaction, ECrel, PPrel, GErel, PCrel
        note: chemical and reaction can not be selected together
        note: reaction and ECrel can not be selected together
    :param reverse_interaction_doubler: if set to true, reversible reaction edges will be contained in the network table in two directions (A to B, B to A).
    :return: a pandas dataframe that contains all the network information.
        Format network dataframe:
            Column 0: source_id (str)
            Column 1: source_type (str)
            Column 2: target_id (str)
            Column 3: target_type (str)
            Column 4: interaction_type (str) (Note: there are six possible interaction types:chemical,reaction,ECrel,PPrel,GErel,PCrel)
            Column 5: interaction_id (str) (Note: only if present)
            Column 6: interaction_info (str) (e.g. reversible,irreversible,..   see KEGG KGML webpage for more info about futher details about the possible interaction types)
    """

    # lists from kgml_reader
    entry_dict = kgml_information_list[0]
    genes_reaction_dict = kgml_information_list[1]
    relation_dict = kgml_information_list[2]
    reaction_dict = kgml_information_list[3]

    # incorporate reactions
    network_source_ids = []
    network_source_types = []
    network_target_ids = []
    network_target_types = []
    network_interaction_types = []
    network_interaction_ids = []
    network_interaction_infos = []

    if "reaction" in selected_interaction_types or "chemical" in selected_interaction_types:
        for reaction in reaction_dict:

            reaction_id = reaction  # interaction id
            reaction_info = reaction_dict[reaction][0]  # interaction info: reversible or irreversible

            reaction_substrate_entry_ids = reaction_dict[reaction][1]  # list of all substrate compounds (entry ids)
            reaction_product_entry_ids = reaction_dict[reaction][2]  # list of all product compounds (entry ids)

            # previous lists of substrate, product and enzyme entry ids must be converted to KEGG ids
            reaction_substrate_KEGG_ids = []  # list of all substrate compounds (KEGG ids)
            for substrate_entry_id in reaction_substrate_entry_ids:
                substrate_names = entry_dict[substrate_entry_id][1]  # one entry id may represent several compounds
                substrate_type = entry_dict[substrate_entry_id][0]

                if substrate_type == "compound":
                    for compound in substrate_names:  # one entry id may represent several compounds
                        reaction_substrate_KEGG_ids.append(str(compound))
                else:
                    raise Exception("Error in network_builder. Non-compound detected as a substrate of a reaction")

            reaction_product_KEGG_ids = []  # list of all product compounds (KEGG ids)
            for product_entry_id in reaction_product_entry_ids:  # list of all product compounds (KEGG ids)
                product_names = entry_dict[product_entry_id][1]  # one entry id may represent several compounds
                product_type = entry_dict[product_entry_id][0]

                if product_type == "compound":
                    for compound in product_names:  # one entry id may represent several compounds
                        reaction_product_KEGG_ids.append(str(compound))
                else:
                    raise Exception("Error in network_builder. Non-compound detected as a product of a reaction")

            # list of all genes that can catalyze this reaction (KEGG ids)
            reaction_enzyme_KEGG_ids = genes_reaction_dict[reaction]

            # NETWORK TYPE1: chemical network (only metabolites)
            if "chemical" in selected_interaction_types:
                reaction_type = "chemical"
                for substrate in reaction_substrate_KEGG_ids:
                    for product in reaction_product_KEGG_ids:
                        network_source_ids.append(substrate)
                        network_source_types.append("compound")
                        network_target_ids.append(product)
                        network_target_types.append("compound")

                        network_interaction_types.append(reaction_type)
                        network_interaction_ids.append(reaction_id)
                        network_interaction_infos.append(reaction_info)

                        if reaction_info == "reversible" and reverse_interaction_doubler:
                            network_source_ids.append(product)
                            network_source_types.append("compound")
                            network_target_ids.append(substrate)
                            network_target_types.append("compound")

                            network_interaction_types.append(reaction_type)
                            network_interaction_ids.append(reaction_id)
                            network_interaction_infos.append(reaction_info)

            # NETWORK TYPE2: metabolic network (metabolites and enzymes:unique enzyme code per reaction)
            if "reaction" in selected_interaction_types:
                reaction_type = "reaction"  # interaction type
                for substrate in reaction_substrate_KEGG_ids:
                    for gene in reaction_enzyme_KEGG_ids:
                        network_source_ids.append(substrate)
                        network_source_types.append("compound")
                        network_target_ids.append((gene + "_" + reaction).replace(organism_code+":",""))  # Some genes can catalyse several reactions. We give those genes a unique source/target id per reaction. Otherwise, there would be a bias in our impact score analysis.
                        network_target_types.append("gene")
                        network_interaction_types.append(reaction_type)
                        network_interaction_ids.append(reaction_id)
                        network_interaction_infos.append(reaction_info)

                        if reaction_info == "reversible" and reverse_interaction_doubler:
                            network_source_ids.append((gene + "_" + reaction).replace(organism_code+":", ""))
                            network_source_types.append("gene")
                            network_target_ids.append(substrate)
                            network_target_types.append("compound")
                            network_interaction_types.append(reaction_type)
                            network_interaction_ids.append(reaction_id)
                            network_interaction_infos.append(reaction_info)

                for product in reaction_product_KEGG_ids:
                    for gene in reaction_enzyme_KEGG_ids:
                        network_source_ids.append((gene + "_" + reaction).replace(organism_code+":", ""))
                        network_source_types.append("gene")
                        network_target_ids.append(product)
                        network_target_types.append("compound")
                        network_interaction_types.append(reaction_type)
                        network_interaction_ids.append(reaction_id)
                        network_interaction_infos.append(reaction_info)

                        if reaction_info == "reversible" and reverse_interaction_doubler:
                            network_source_ids.append(product)
                            network_source_types.append("compound")
                            network_target_ids.append((gene + "_" + reaction).replace(organism_code+":", ""))
                            network_target_types.append("gene")
                            network_interaction_types.append(reaction_type)
                            network_interaction_ids.append(reaction_id)
                            network_interaction_infos.append(reaction_info)

    if "ECrel" in selected_interaction_types or "PPrel" in selected_interaction_types or "GErel" in selected_interaction_types or "PCrel" in selected_interaction_types:
        # incorporate relations
        for relation in relation_dict:
            relation_id = relation
            relation_entry_1_entry_id = relation_dict[relation][0]
            relation_entry_1_type = entry_dict[relation_entry_1_entry_id][0]
            relation_entry_2_entry_id = relation_dict[relation][1]
            relation_entry_2_type = entry_dict[relation_entry_2_entry_id][0]
            relation_type = relation_dict[relation][2]  # interaction type
            relation_interaction_info = relation_dict[relation][3]  # interaction info

            # NETWORK TYPE 3,4,5,6
            # only add interaction types that were selected to the network
            if relation_type in selected_interaction_types:

                # first we need to obtain the names of the entry ids in the relation. Some of these nodes could be groups, as such we need retrieve the names of the genes inside group entries.
                if relation_entry_1_type == "compound" or relation_entry_1_type == "gene":
                    relation_entry_1_names = entry_dict[relation_entry_1_entry_id][1]
                    relation_entry_1_names = [x.replace(organism_code+":", "") for x in relation_entry_1_names]

                if relation_entry_2_type == "compound" or relation_entry_2_type == "gene":
                    relation_entry_2_names = entry_dict[relation_entry_2_entry_id][1]
                    relation_entry_2_names = [x.replace(organism_code+":", "") for x in relation_entry_2_names]

                if relation_entry_1_type == "group":
                    relation_entry_1_names = []
                    relation_entry_1_type = "group_gene"
                    component_entry_ids = entry_dict[relation_entry_1_entry_id][2]
                    for component_entry_id in component_entry_ids:
                        relation_entry_1_names = relation_entry_1_names + entry_dict[component_entry_id][1]
                        if entry_dict[component_entry_id][0] != "gene":
                            raise Exception("group contains things other than genes")
                    relation_entry_1_names = [x.replace(organism_code+":", "") for x in relation_entry_1_names]

                if relation_entry_2_type == "group":
                    relation_entry_2_names = []
                    relation_entry_2_type = "group_gene"
                    component_entry_ids = entry_dict[relation_entry_2_entry_id][2]
                    for component_entry_id in component_entry_ids:
                        relation_entry_2_names = relation_entry_2_names + entry_dict[component_entry_id][1]
                        if entry_dict[component_entry_id][0] != "gene":
                            raise Exception("group contains things other than genes")
                    relation_entry_2_names = [x.replace(organism_code+":", "") for x in relation_entry_2_names]

                # add interactions to network table
                for start_entry_name in relation_entry_1_names:
                    for end_entry_name in relation_entry_2_names:
                        network_source_ids.append(start_entry_name)
                        network_source_types.append(relation_entry_1_type)
                        network_target_ids.append(end_entry_name)
                        network_target_types.append(relation_entry_2_type)
                        network_interaction_types.append(relation_type)
                        network_interaction_ids.append(relation_id)
                        network_interaction_infos.append(relation_interaction_info)


    # NETWORK TABLE HAS BEEN CONSTRUCTED
    # WRITE TABLE TO TSV FILE
    pathway_network_dict = {"source_id": network_source_ids, "source_type": network_source_types,
                            "target_id": network_target_ids, "target_type": network_target_types,
                            "interaction_type": network_interaction_types,
                            "interaction_id": network_interaction_ids,
                            "interaction_info": network_interaction_infos}
    pathway_network_table = pd.DataFrame.from_dict(pathway_network_dict)

    return pathway_network_table


def KEGG_network_constructor_from_list(KEGG_pathways_list: list,organism_code: str,KEGG_interaction_types:list,path_outputfile: str, reverse_interaction_doubler: bool):
    """
    This method will create tsv file containing network information of the KEGG pathways described in the KEGG pathways list.
        Note: this method assumes that the pathway list is totally correct (see method: KEGG_pathway_list_decoder).
        Note: this method assumes that the list of selected interaction types is totally correct (see method: KEGG_interacion_type_decoder).

    :param KEGG_pathways_list: a list of KEGG pathway ids of a specific organism in the KEGG database
    :param organism_code: the code that KEGG uses for a specific organism (e.g. 'pae' for Pseudomonas aeruginosa PAO1)
    :param KEGG_interaction_types: a list of all the interaction types that need to be included in the network table.
        This network may contain the following types of interactions: chemical, reaction, ECrel, PPrel, GErel, PCrel
        note: chemical and reaction can not be selected together
        note: reaction and ECrel can not be selected together
    :param path_outputfile: path to the output file (includes file name) (e.g. C:/test/all_metabolic_reactions_pae.txt)
    :param reverse_interaction_doubler: if set to true, reversible reaction edges will be contained in the network table in two directions (A to B, B to A).
    :return: a tsv file containing info on KEGG pathways in a network format.
        Format network tsv file:
            Column 0: source_id (str)
            Column 1: source_type (str)
            Column 2: target_id (str)
            Column 3: target_type (str)
            Column 4: interaction_type (str) (Note: there are seven possible interaction types:chemical,reaction,ECrel,PPrel,GErel,PCrel,identical_id_connection)
            Column 5: interaction_id (str) (Note: only if present)
            Column 6: interaction_info (str) (e.g. reversible,irreversible,..   see KEGG KGML webpage for more info about futher details about the possible interaction types)
    """

    Total_network_table = pd.DataFrame(columns=['source_id','source_type','target_id','target_type','interaction_type','interaction_id','interaction_info']) #per pathway a dataframe is created. It is then added to the total network table.

    requests.adapters.DEFAULT_RETRIES = 5  # increase retries number
    s = requests.session()
    #s.keep_alive = False

    for pathway in KEGG_pathways_list:
        # connect to KEGG database
        time.sleep(0.34) #KEGG has a 3/s call limit. So wait 1/3 of a second before requesting information about the next pathway.
        api_url = f"https://rest.kegg.jp/get/{pathway}/kgml"
        response = s.get(api_url, timeout=5)
        pathway_xml = io.StringIO(str(response.text))

        # read out the kgml file
        kgml_info_list = kgml_reader(pathway_xml)

        # make a pandas dataframe that contains the network information of the specified pathway
        Pathway_network_table = network_builder(kgml_info_list,organism_code,KEGG_interaction_types,reverse_interaction_doubler)

        # merge new pathway dataframe with the dataframe which is going to contain all pathway info
        Total_network_table = pd.concat([Total_network_table,Pathway_network_table])

    #remove duplicates (duplicate rows)
    Total_network_table.drop_duplicates(inplace=True)

    #genes with variants of the same name (genename_reaction1, genename_reaction2,...) should be connected.
    #we will make a small dataframe containing these identical_id_connections
    network_source_ids = []
    network_source_types = []
    network_target_ids = []
    network_target_types = []
    network_interaction_types = []
    network_interaction_ids = []
    network_interaction_infos = []


    rows = Total_network_table.shape[0]
    reaction_genes_extended_name_set = set()
    relation_genes_short_name_set = set()
    for row in range(rows):
        # make list of all genes in reactions (format id: name_reactionid)
        source_name = str(Total_network_table.iloc[row,0])
        target_name = str(Total_network_table.iloc[row,2])
        interaction_type = Total_network_table.iloc[row,4]
        if interaction_type == "reaction":
            if "cpd" not in source_name:
                reaction_genes_extended_name_set.add(source_name)
            if "cpd" not in target_name:
                reaction_genes_extended_name_set.add(target_name)

        # make list of all genes in relations (format id: name)
        if (interaction_type != "reaction") and (interaction_type != "chemical"):
            if "cpd" not in source_name:
                relation_genes_short_name_set.add(source_name)
            if "cpd" not in target_name:
                relation_genes_short_name_set.add(target_name)


    reaction_genes_extended_name_list = list(reaction_genes_extended_name_set)
    reaction_genes_short_name_list = [extended_name.split("_")[0] for extended_name in reaction_genes_extended_name_list]

    # genes for which we have to make a connection between a node in a relation
    genes_part_of_relation_and_reaction_interactions_to_be_selfconnected_list = list(set(reaction_genes_short_name_list).intersection(relation_genes_short_name_set))

    # genes that have to be checked for name variants within reaction interactions
    genes_that_have_to_checked_for_indentical_variant_names = list(set(reaction_genes_short_name_list)-relation_genes_short_name_set)

    # We will first connect genes of the relation interactions with name variants part of reaction interactions
    # the genes in the reaction network are now still unconnected to the genes in the relation networks, due to the gene names in the reaction network being formatted as (genename_reactionid)
    if ("reaction" in KEGG_interaction_types) and ("ECrel" in KEGG_interaction_types or "PPrel" in KEGG_interaction_types or "GErel" in KEGG_interaction_types or "PCrel" in KEGG_interaction_types):
        #make the connections
        for gene_short in genes_part_of_relation_and_reaction_interactions_to_be_selfconnected_list:
            for gene_extended in reaction_genes_extended_name_list:
                if gene_short == gene_extended.split('_')[0]:
                    #define edge from short to extended id
                    network_source_ids.append(gene_short)
                    network_source_types.append("gene")
                    network_target_ids.append(gene_extended)
                    network_target_types.append("gene")
                    network_interaction_types.append("identical_id_connection")
                    network_interaction_ids.append("NA")
                    network_interaction_infos.append("NA")

    # Add identical_id_connection between genes in the reaction network that aren't linked through nodes from the relation networks
    for gene_short_reaction_name in genes_that_have_to_checked_for_indentical_variant_names:
        for gene_extended_reaction_name in reaction_genes_extended_name_list:
            if (reaction_genes_short_name_list.count(gene_short_reaction_name) > 1) and (gene_extended_reaction_name.split('_')[0] == gene_short_reaction_name):
                # define edge from short to extended id
                network_source_ids.append(gene_short_reaction_name)
                network_source_types.append("gene")
                network_target_ids.append(gene_extended_reaction_name)
                network_target_types.append("gene")
                network_interaction_types.append("identical_id_connection")
                network_interaction_ids.append("NA")
                network_interaction_infos.append("NA")


    #make dataframe
    identical_id_connection_network_dict = {"source_id": network_source_ids, "source_type": network_source_types,
                            "target_id": network_target_ids, "target_type": network_target_types,
                            "interaction_type": network_interaction_types,
                            "interaction_id": network_interaction_ids,
                            "interaction_info": network_interaction_infos}
    identical_id_connection_network_table = pd.DataFrame.from_dict(identical_id_connection_network_dict)

    #remove duplicate rows
    identical_id_connection_network_table.drop_duplicates(inplace=True)

    #merge with total dataframe
    Total_network_table = pd.concat([Total_network_table, identical_id_connection_network_table])

    #writing pandas dataframe to tsv
    Total_network_table.to_csv(path_outputfile,sep='\t',index=False)


def KEGG_network_constructor_from_string(KEGG_pathways_string: str,organism_code: str,selected_KEGG_interaction_types_string: str,path_outputfile: str, reverse_interaction_doubler = False):
    """
    This method will create tsv file containing network information of the KEGG pathways described in the KEGG pathways string.
        Note: this method assumes that the pathway list is totally correct (see method: KEGG_pathway_list_decoder).
        Note: this method assumes that the list of selected interaction types is totally correct (see method: KEGG_interacion_type_decoder).

    :param KEGG_pathways_string: a string containing KEGG pathway ids of a specific organism in the KEGG database separated by the delimiter '$'.
    :param organism_code: the code that KEGG uses for a specific organism (e.g. 'pae' for Pseudomonas aeruginosa PAO1)
    :param selected_KEGG_interaction_types_string: a string of all the interaction types that need to be included in the network table separated by the delimiter '$'.
        This network may contain the following types of interactions: chemical, reaction, ECrel, PPrel, GErel, PCrel
        note: chemical and reaction can not be selected together
        note: reaction and ECrel can not be selected together
    :param path_outputfile: path to the output file (includes file name) (e.g. C:/test/all_metabolic_reactions_pae.txt)
    :param reverse_interaction_doubler: if set to true, reversible reaction edges will be contained in the network table in two directions (A to B, B to A). False by default.
    :return: a tsv file containing info on KEGG pathways in a network format.
        Format network tsv file:
            Column 0: source_id (str)
            Column 1: source_type (str)
            Column 2: target_id (str)
            Column 3: target_type (str)
            Column 4: interaction_type (str) (Note: there are seven possible interaction types:chemical,reaction,ECrel,PPrel,GErel,PCrel,identical_id_connection)
            Column 5: interaction_id (str) (Note: only if present)
            Column 6: interaction_info (str) (e.g. reversible,irreversible,..   see KEGG KGML webpage for more info about futher details about the possible interaction types)
    """
    selected_KEGG_interaction_types_list = KEGG_interacion_type_decoder(selected_KEGG_interaction_types_string)
    KEGG_pathways_list = KEGG_pathway_list_decoder(KEGG_pathways_string,organism_code)
    KEGG_network_constructor_from_list(KEGG_pathways_list,organism_code,selected_KEGG_interaction_types_list,path_outputfile,reverse_interaction_doubler)


def network_merger(path_inputfile_network_1:str,path_inputfile_network_2:str,connect_identical_names:bool,output_directory_and_filename: str):
    """
    This method merges two networks, and then returns a tsv file containing the merged network table in the specified location.
        The first seven columns should be (in this order):
            0) source_id
            1) source_type
            2) target_id
            3) target_type
            4) interaction_type
            5) interaction_id
            6) interaction_info
    This method assumes that if you have multiple id that are near identical, they are connected via a short variant of that gene id using identical_id_connections (interacion type).
    """
    # INPUT VALIDATON --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #read tsv files containing the two networks
    network_1 = pd.read_csv(path_inputfile_network_1,sep="\t")
    network_2 = pd.read_csv(path_inputfile_network_2, sep="\t")

    #check format of first seven columns
    required_column_names = ["source_id","source_type","target_id","target_type","interaction_type","interaction_id","interaction_info"]

    network_1_column_names = list(network_1.columns)
    network_2_column_names = list(network_2.columns)
    for index in range(7):
        if network_1_column_names[index] != required_column_names[index] or network_2_column_names[index] != required_column_names[index]:
            raise Exception("Please check the format of the input network tsv files. The first seven columns should be named: source_id,source_type,target_id,target_type,interaction_type,interaction_id,interaction_info.")

    # MERGE TWO NETWORKS ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # first add two datafames together
    network_merged_pd = pd.concat(network_1,network_2)

    # check for edges with same source and target node











    # ADDITION OF CONNECTING EDGES FOR IDENTICAL ID VARIANTS (if requested) ---------------------------------------------------------------------------------------------------------------------------------------
    if connect_identical_names:

        # make list of all long gene ids in network 1
        genes_network1_sources_long_id_list = network_1.loc[(network_1['source_type'] == "gene") or (network_1['source_type'] == 'group_gene'), "source_id"]
        genes_network1_targets_long_id_list = network_1.loc[(network_1['target_type'] == "gene") or (network_1['target_type'] == 'group_gene'), "target_id"]
        genes_network1_long_id_set = set(genes_network1_sources_long_id_list + genes_network1_targets_long_id_list)
        genes_network1_long_id_list = list(genes_network1_long_id_set)

        # make list of all genes in network 2
        genes_network2_sources_long_id_list = network_2.loc[(network_2['source_type'] == "gene") or (network_2['source_type'] == 'group_gene'), "source_id"]
        genes_network2_targets_long_id_list = network_2.loc[(network_2['target_type'] == "gene") or (network_2['target_type'] == 'group_gene'), "target_id"]
        genes_network2_long_id_set = set(genes_network2_sources_long_id_list + genes_network2_targets_long_id_list)
        genes_network2_long_id_list = list(genes_network2_long_id_set)

        #make dict of short id to all long ids in network 1
        genes_network1_short_to_long_id_dict = {}

        for long_id in genes_network1_long_id_list:
            short_id = long_id.split('_')[0]
            if short_id not in genes_network1_short_to_long_id_dict:
                genes_network1_short_to_long_id_dict[short_id] = [long_id]
            else:
                list_long_ids = genes_network1_short_to_long_id_dict[short_id]
                list_long_ids.append(long_id)
                genes_network1_short_to_long_id_dict[short_id] = list_long_ids

        #make dict of short id to all long ids in network 2
        genes_network2_short_to_long_id_dict = {}

        for long_id in genes_network2_long_id_list:
            short_id = long_id.split('_')[0]
            if short_id not in genes_network2_short_to_long_id_dict:
                genes_network2_short_to_long_id_dict[short_id] = [long_id]
            else:
                list_long_ids = genes_network2_short_to_long_id_dict[short_id]
                list_long_ids.append(long_id)
                genes_network2_short_to_long_id_dict[short_id] = list_long_ids

        #make list for which we need to make connections
        network1_short_ids_set = set(genes_network1_short_to_long_id_dict.keys())
        network2_short_ids_set = set(genes_network2_short_to_long_id_dict.keys())
        overlapping_short_ids_set = network1_short_ids_set.intersection(network2_short_ids_set) #this still contains ids that match exactly between the networks

        #make list of ids that match exactly between the networks
        genes_network1_long_id_set = set(genes_network1_long_id_list)
        genes_network2_long_id_set = set(genes_network2_long_id_list)
        overlapping_long_ids_set = genes_network1_long_id_set.intersection(genes_network2_long_id_set)
        overlapping_long_ids_list = list(overlapping_long_ids_set)
        overlapping_exact_long_id_matches_short_ids_list =  [x.split('_')[0] for x in overlapping_long_ids_list]
        overlapping_exact_long_id_matches_short_ids_set = set(overlapping_exact_long_id_matches_short_ids_list)

        #remove genes ids that match exactly between the networks
        overlapping_short_ids_removed_exact_matches_set = overlapping_short_ids_set - overlapping_exact_long_id_matches_short_ids_set
        overlapping_short_ids_removed_exact_matches_list = list(overlapping_short_ids_removed_exact_matches_set)

        #make connections
        # we will make a small dataframe containing these identical_id_connections
        network_source_ids = []
        network_source_types = []
        network_target_ids = []
        network_target_types = []
        network_interaction_types = []
        network_interaction_ids = []
        network_interaction_infos = []

        for short_id in overlapping_short_ids_removed_exact_matches_list:
            connection_code = 0 #0: connect both networks to short new node|1: connect nodes in network 2 to short name that is already in network 1|2: connect nodes in network 1 to short name that is already in network 2

            if short_id in genes_network1_short_to_long_id_dict[short_id]: #if the short version of an id is already in network 1 or 2, than
                connection_code = 1
            else:
                if short_id in genes_network2_short_to_long_id_dict[short_id]:
                    connection_code = 2

            if ((connection_code == 1) or (connection_code == 0)):
                for network2_id in genes_network2_short_to_long_id_dict[short_id]:
                    network_source_ids.append(short_id)
                    network_source_types.append("gene")
                    network_target_ids.append(network2_id)
                    network_target_types.append("gene")
                    network_interaction_types.append("identical_id_connection")
                    network_interaction_ids.append("NA")
                    network_interaction_infos.append("NA")

            if ((connection_code == 2) or (connection_code == 0)):
                for network1_id in genes_network1_short_to_long_id_dict[short_id]:
                    network_source_ids.append(short_id)
                    network_source_types.append("gene")
                    network_target_ids.append(network1_id)
                    network_target_types.append("gene")
                    network_interaction_types.append("identical_id_connection")
                    network_interaction_ids.append("NA")
                    network_interaction_infos.append("NA")

    # make dataframe that contains all identical id connections
    identical_id_connection_network_dict = {"source_id": network_source_ids,
                                                    "source_type": network_source_types,
                                                    "target_id": network_target_ids,
                                                    "target_type": network_target_types,
                                                    "interaction_type": network_interaction_types,
                                                    "interaction_id": network_interaction_ids,
                                                    "interaction_info": network_interaction_infos}
    identical_id_connection_network_table = pd.DataFrame.from_dict(identical_id_connection_network_dict)

    #merge with total dataframe
    network_merged_pd = pd.concat([network_merged_pd, identical_id_connection_network_table])

    #remove duplicate rows
    network_merged_pd.drop_duplicates(inplace=True)

    #write to output tsv file
    network_merged_pd.to_csv(output_directory_and_filename,sep="\t",index=False)

#testing -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
metabolic_pathways_pae = KEGG_organism_all_metabolic_pathway_retriever("pae")

KEGG_network_constructor_from_string(metabolic_pathways_pae,"pae","reaction","C:/test/all_metabolic_reactions_pae_02072024.txt")


"""
    This method uses the list provided by KEGG_pathway_list_decoder to look up the relevant in the KEGG database using the KEGG API.
    During the course of this process,"""