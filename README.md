# MONIR: Multi-omics Network Impact-based Ranking

## Table of contents
1. [Introduction](#introduction)
2. [System requirements](#system-requirements)
3. [Installation](#installation)
4. [Basic installation verification tests](#basic-installation-verification-tests)
5. [Workflow](#workflow)
    1. [Network construction](#network-construction)
    2. [Impact analysis](#impact-analysis)
    4. [Subnetwork visualisation](#subnetwork-visualisation)
6. [Usage](#usage)
7. [Citing MONIR](#citing-monir)
8. [Contact](#contact)


## Introduction

MONIR is a versatile gene prioritisation tool that employs a network diffusion-like approach. 
It ranks genes based on the number of and proximity to observed omics changes in a network. 
To achieve this, it calculates an ‘impact score’ for each gene of interest using node weights after diffusion. 
MONIR supports a wide range of omics types, including transcriptomics, proteomics, phosphorylation, glycosylation, ubiquitination, and metabolomics. 
Notably, it simultaneously considers all available omics data, the network topology, and the specific types and directions of interactions between network nodes to calculate impact scores.


## System Requirements
- **Operating Systems**: The tool can run on Linux, macOS, and Windows.
- **Python Version**: MONIR requires Python 3.11 or higher.

## Installation
- MONIR is available on pip
- command: 'pip install monir'

## Basic installation verification tests
- command: 'monir-test-cli'
- command: 'monir-monir-test-methods'

## Workflow

A standard analysis involves three phases. 
In the first phase, a network should be constructed. 
Next, the network and available omics data are used to rank a set of genes-of-interest (GOI)/nodes-of-interest (NOI). 
Finally, subnetworks centred around high-scoring genes can be generated for downstream visualisation and inspection. 
In this section, we will give an overview of the methods for each phase and their input and output files.
More information, concerning method options, can be found in the help pages of each method via command line.


### network construction
- Available methods:
    - **KEGG_network_constructor**: creates a network from user-specified KEGG pathways, including only the selected interaction types.
    - **network_merger**: merges two user-specified networks
- Output (also input format for the network_merger method):
    - Network file
      - File format: TSV
      - File structure:
        - *Column 0:*
          - header: source_id
          - definition: the unique identifier of the node
          - format: compound identifiers should be represented as 'cpd:*compound*' (E.g. cpd:C00121).
            Gene and group_gene identifiers should be represented as '*gene*' or '*gene*_*reaction* (e.g. PA5546 or PA3375_rn:R10185).
        - *Column 1:* 
          - header: source_type
          - definition: the type of node
          - supported values: compound, gene, group_gene (as defined by KEGG)
        - *Column 2:* 
          - header: target_id
          - definition: the unique identifier of the node
          - format: compound identifiers should be represented as 'cpd:*compound*' (E.g. cpd:C00121).
            Gene and group_gene identifiers should be represented as '*gene*' or '*gene*_*reaction* (e.g. PA5546 or PA3375_rn:R10185).
        - *Column 3:* 
          - header: target_type
          - definition: the type of node
          - supported values: compound, gene, group_gene (as defined by KEGG)
        - *Column 4:* 
          - header: interaction_type
          - definition: the type of edge between the target and source node
          - supported values: PPrel, GErel, PCrel, ECrel, chemical, reaction, identical_id_connection
          - additional information: PPrel, GErel and PCrel represent protein-protein, gene expression and protein-compound interactions, respectively. ECrel, chemical and reaction are                  three different ways to represent metabolic reactions. ECrel interactions link enzymes that catalyze the neighbouring reactions. Chemical interactions link                                  reagents and products of a reaction with each other.  In contrast, the ‘reaction’ format links both compounds and enzymes of a reaction together. Notably, enzymes                           catalysing several reactions, are given a unique identifier per reaction (e.g. PA3375_rn:R10185). These are then connected to a node with the short variation of the identifier              (e.g. PA3375) via identical_id_connection type edges. This keeps the reactions separate in the network.
        - *Column 5:* 
          - header: interaction_id
          - definition: identifier of the interaction
          - additional information: this identifier is only unique for reaction type interactions
        - *Column 6:* 
          - header: interaction_info
          - definition: gives additional information on the interaction types.
          - supported values: For the interaction type specific impact calculation, this value should correspond to one or multiple of the possible relation name attributes in the KEGG                 Markup Language (e.g. expression, inhibition, phosphorylation,...) separated with a dollar sign or reversible/irreversibe in case of a metabolic reaction.

### impact analysis
- Available methods:
    - **node_impact_assessor:** calculates impact scores for all genes-of-interest based on provided omics data and network information using a network diffusion-like approach.
      This diffusion approach is explained in detail in the MONIR paper (see the section 'Citing MONIR').
    
- Input:
  - File containing nodes-of-interest
    - File format: TSV
    - File structure:
      - *Column 0:*
          - header: node_id
          - definition: the gene identifiers for which you want to calculate impact scores
          - additional information: if you list an identifier, but there is not a node with that identifier in the network, this method will also look for matches with longer                   identifiers in the network. If such a longer match is found, the impact score will be calculated for the node with the longer identifier. For example, if you provide PA3375 and no exact match is found, the impact score will be caculated for PA3375_rn:R10185. Importantly, the part before the underscore must match with the user-provided identifier.
            
  - File containing omics information:
    - File format: TSV
    - File structure:
      - *Column 0:* 
          - header: node_id
          - definition: the identifier of the node
          - format: compound identifiers should be represented as 'cpd:*compound*' (E.g. cpd:C00121).
            Gene and group_gene identifiers should be represented as '*gene*' or '*gene*_*reaction*' (e.g. PA5546).
          - additional information: if you list an identifier, but there is not a node with that identifier in the network, this method will also look for matches with longer                   identifiers in the network. If such a longer match is found, the provided omics information will be connected to the node with the longer identifier. Importantly, the part before the underscore in the identifier must match with the user-provided identifier.
      - *Column 1:*
          -  header: changed
          -  definition: indicates whether an omics change was registered for this gene/compoun
          -  supported values: True, False
      - *Column 2:* 
          - header: changed_omics_type
          - definition: the type of omics data that was measured
          - supported values: transcriptomics, proteomics, phosphorylation, glycosylation, ubiquitination, metabolomics
            
  - Network file:
      - Supported values: metabolic reactions should be formatted using the 'reaction' interaction type
      - Additional information: see network contruction > output
- Output:
    - File with impact analysis results
      - File format: TSV
      - File structure:
        - *Column 0:*
            - header: NOI_id
            - definition: node identifier as specified in the file containing the identifiers of the nodes-of-interest
        - *Column 1:*
            - header: NOI_network_id
            - definition: node indentifier as found in the network
        - *Column 2:*
            - header: total_impact_score
            - definition: see MONIR paper
        - *Column 3:*
            - header: hypothetical_maximum_total_impact_score
            - definition: the total impact score of a gene if all nodes in the network display changes in the user-measured omics layers 
        - *Column 4:*
            - topological_independent_total_impact_score
            - definition: the total impact score of a gene devided by the hypothetical maximum total impact score of the same gene
        - *Column 5:*
            - header: contributing_nodes
            - definition: nodes that contribute to the impact score of the gene-of-interest
            - additional information: the identifiers of the nodes are separated by a dollar sign
        - *Column 6:*
            - header: contributing_nodes_sub_scores
            - definition: the weights that each contributing node adds to the total impact score of the gene-of-interest 
            - additional information: the weights are separated with a dollar sign and are in the same order as the identifiers of the contributing nodes
        - *Column 7:*
            - header: contributing_nodes_max_level
            - definiton: the maximum number of propagation steps between the gene-of-interest and a contributing node
            - additional information: for reaction type interactions, one propagation step goes from one compound to the next compound (skips the gene). For the other interaction types,                  one propagation step corresponds to going to the direct neighbour (compound/gene).
        - *Column 8:*
            - header: contributing_nodes_level
            - definition: the number of propagation steps between the gene-of-interest and each contributing node
            - additional information: the levels are separated with a dollar sign and are in the same order as the identifiers of the contributing nodes
        - *Column 9:*
            - header: top_impacting_node_type
            - definition: the type of node which contributes the most to the total impact score
        - *Column 10:*
            - header: node_types
            - defintion: the types of nodes that have contributed to the total impact score
        - *Column 11:*
            - header: node_types_sub_scores
            - definition: the amount that each node type contributes to the total impact score
            - additional information: these amounts are separated with a dollar sign and are in the same order as the node types

- Additional information:
    - The impact analysis can be performed in an interaction-specific way. This involves checking whether an omics change might be 'informative' or not (see MONIR paper). We use the              following rules to assess this:
        - the layer that we have info about for a certain gene node should be able to have been affected by an upstream node. In other words:
            - have transcriptomics or proteomics info about a changing node that is a target of a GErel type interaction.
            - have methylation info about a changing node that is a target of a PPrel type interaction with methylation interaction info.
            - have ubiquitination info about a changing node that is a target of a PPrel type interaction with ubiquitination interaction info.
            - have glycosylation info about a changing node that is a target of a PPrel type interaction with glycosylation interaction info.
            - have phosphorylation info about a changing node that is a target of a PPrel type interaction with phosphorylation or dephosphorylation interaction info.
        - the layer that we have info about for a certain gene node should be able to affect a downstream node. In other words:
            - have transcriptomics or proteomics info about a changing node that is a source of a GErel type interaction.
            - have transcriptomics or proteomics info about a changing node that is a source of a PPrel type interaction with methylation interaction info.
            - have transcriptomics or proteomics info about a changing node that is a source of a PPrel type interaction with ubiquitination interaction info.
            - have transcriptomics or proteomics info about a changing node that is a source of a PPrel type interaction with glycosylation interaction info.
            - have transcriptomics or proteomics info about a changing node that is a source of a PPrel type interaction with phosphorylation or dephosphorylation interaction info.
        - we have metabolomics data for a compound node

### subnetwork visualisation
- Available methods:
    - **subnetwork_constructor:** makes a subnetwork (TSV file) for each gene-of-interest that passes a user-set impact threshold.
    - **node_table_id_and_type_extender:** makes a table (TSV file) containing a column that lists all the (network) node ids and a column that contains the short node ids (part before '_') and a column that contains all the node types (gene or compound). It can be used to extend the node table in cytoscape (this enables colouring based on type).
    RISK: separation of types is based on compound ids containing "cpd:"
    - **annotation_table_id_extender:** will add a column containing network node ids to your annotation table (TSV file). This table can then be used to extend the node table in cytoscape with the info in the annotation file.
- Networks and subnetworks can be imported to Cytoscape. 

## Usage
A usage case is described in the supplementary information of the MONIR paper (see section 'citing MONIR').

## Citing MONIR
- If you use MONIR in a publication, please cite:
**coming soon**

- Furthermore, please cite the databases used in the construction of your input network.
In case you used our in-built methods to construct KEGG networks, also reference:

  Kanehisa, M., Furumichi, M., Sato, Y., Matsuura, Y. and Ishiguro-Watanabe, M.; KEGG: biological systems database as a model of the real world. Nucleic Acids Res. 53, D672-D677 (2025)

  Kanehisa, M; Toward understanding the origin and evolution of cellular organisms. Protein Sci. 28, 1947-1951 (2019)

  Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000)

- If you used Cytoscape for visualisation, please cite:

  Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T. Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Research 2003 Nov; 13(11):2498-504

## Contact
In case of questions, please contact the Laboratory of Computational Systems Biology (KU Leuven).
