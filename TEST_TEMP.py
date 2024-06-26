import pandas as pd
import requests
import time
import io
import xml.etree.ElementTree as ET


network_table = pd.read_csv("C:/test/pae00010_reaction_pae_14052024.txt", sep="\t")
print(network_table)

previous_nodes_id_reversible_KEGG_reactions = network_table.loc[(network_table["source_id"]=="PA1770_rn:R00199")&(network_table["interaction_info"]=="reversible"),]
print(previous_nodes_id_reversible_KEGG_reactions)

test_list = [['PA2275_rn:R00746', 'reaction', 'reversible'], ['PA2275_rn:R00746', 'reaction', 'reversible'], ['PA2275_rn:R00746', 'reaction', 'reversible'], ['PA2275_rn:R00746', 'reaction', 'reversible'], ['PA2275_rn:R00746', 'reaction', 'reversible']]
test_list2 =[[1,3,5],[1,5,3],[1,3,5],[8,7,6]]
unique_data = [list(x) for x in set(tuple(x) for x in test_list2)]

print(test_list)
print(unique_data)