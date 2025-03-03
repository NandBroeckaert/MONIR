# MONIR: Multi-omics Network Impact-based Ranking

## Table of contents
1. [Introduction](#introduction)
2. [System requirements](#system-requirements)
3. [Installation](#installation)
4. [Workflow](#workflow)
    1. [Network construction](#network-construction)
    2. [Impact analysis](#impact-analysis)
    4. [Subnetwork visualisation](#subnetwork-visualisation)
5. [Usage](#usage)
6. [Citing MONIR](#citing-MONIR)
7. [Contact](#contact)


## Introduction

MONIR is a versatile gene prioritisation tool that employs a network diffusion-like approach. 
It ranks genes based on the number of and proximity to observed omics changes in a network. 
To achieve this, it calculates an ‘impact score’ for each gene of interest using node weights after diffusion. 
MONIR supports a wide range of omics types, including transcriptomics, proteomics, phosphorylation, glycosylation, ubiquitination, and metabolomics. 
Notably, it simultaneously considers all available omics data, the network topology, and the specific types and directions of interactions between network nodes to calculate impact scores. 


## System Requirements
- **Operating Systems**: The tool can run on Linux, macOS, and Windows.
- **Python Version**: MONIR requires Python 3.9 or higher.

## Installation
- MONIR is available on pip
- **insert steps**

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
- Output (also inputs for the network_merger method):
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
          - supported values: compound, gene, group_gene
        - *Column 2:* 
          - header: target_id
          - definition: the unique identifier of the node
          - format: 
        - *Column 3:* 
          - header: target_type
          - definition:
          - supported values: compound, gene, group_gene
        - *Column 4:* 
          - header: interaction_type
          - definition:
          - supported values:
        - *Column 5:* 
          - header: interaction_id
          - definition:
          - supported values:
        - *Column 6:* 
          - header: interaction_info
          - definition:
          - supported values:
      - Supported values:
        - *source_id and target_id:* Gene node identifiers are formatted as 
        - *source_type and target_type:* A node can be a gene, a gene_group or a compound (as defined by KEGG).
        - *interaction_type:* 'PPrel' (protein-protein), GErel (gene expression), PCrel (protein-compound), metabolic interactions and identical_id_connection 
        - *interaction_id:*
        - *interaction_info:* This gives additional information on the interaction types.
        - For the interaction type specific impact calculation, these values should correspond to one or multiple of the possible relation name attributes in the KEGG Markup Language (e.g. expression, inhibition, phosphorylation,...).
        - If multiple val
      - No

  - unique: indentical_id_interaction --> zeg hoe reaction interactions opgeslagen worden + hoe de id variants er uit zien


### impact analysis
- Available methods:
    - **node_impact_assessor:** calculates impact scores for all genes-of-interest based on provided omics data and network information using a network diffusion-like approach.
      This diffusion approach is explained in detail in the MONIR paper (see the section 'Citing MONIR').
- Input:
  - File containing nodes-of-interest
    - File format: TSV
    - File structure:
      - *Column 0:* node_id
  - File containing omics information:
    - File format: TSV
    - File structure:
      - *Column 0:* node_id
      - *Column 1:* changed
      - *Column 2:* changed_omics_type
    - Supported values:
  - Network file:
      - see network contruction>output
- Output:
    - File with impact analysis results
      - File format: TSV
      - File structure:
        - *Column 0:* NOI_id
        - *Column 1:* NOI_network_id
        - *Column 2:* total_impact_score
        - *Column 3:* hypothetical_maximum_total_impact_score
        - *Column 4:* topological_independent_total_impact_score
        - *Column 5:* contributing_nodes
        - *Column 6:* contributing_nodes_sub_scores
        - *Column 7:* contributing_nodes_max_level
        - *Column 8:* contributing_nodes_level
        - *Column 9:* top_impacting_node_type
        - *Column 10:* node_types
        - *Column 11:* node_types_sub_scores
      - Notes:

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
**insert paper**

- Furthermore, please cite the databases used in the construction of your input network.
In case you used our in-built methods to construct KEGG networks, also reference:

  Kanehisa, M., Furumichi, M., Sato, Y., Matsuura, Y. and Ishiguro-Watanabe, M.; KEGG: biological systems database as a model of the real world. Nucleic Acids Res. 53, D672-D677 (2025)

  Kanehisa, M; Toward understanding the origin and evolution of cellular organisms. Protein Sci. 28, 1947-1951 (2019)

  Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000)

- If you used Cytoscape for visualisation, please cite:

  Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T. Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Research 2003 Nov; 13(11):2498-504

## Contact
In case of questions, please contact the Laboratory of Computational Systems Biology (KU Leuven).
