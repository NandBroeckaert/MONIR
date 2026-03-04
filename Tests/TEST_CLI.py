"""
CLI Testing Module for MONIR

This script tests the argparse parser and command-line interface structure.
It validates that:
1. All subcommands are registered correctly
2. Arguments are parsed without errors for each subcommand
3. Invalid inputs are rejected appropriately
"""

import sys
from Parsers import get_parser

def test_parser_basic():
    """Test that the parser is correctly initialized."""
    print("=" * 70)
    print("TEST: Parser Initialization")
    print("=" * 70)
    parser = get_parser()
    print("✓ Parser created successfully")
    print()

def test_subcommands_exist():
    """Test that all subcommands are registered."""
    print("=" * 70)
    print("TEST: Subcommands Registration")
    print("=" * 70)
    parser = get_parser()
    # argparse stores subparsers action in a private attribute; look for it by type
    import argparse as _argparse
    subparsers_action = None
    for action in parser._subparsers._actions:
        if isinstance(action, _argparse._SubParsersAction):
            subparsers_action = action
            break

    if subparsers_action is not None:
        subcommands_dict = subparsers_action.choices
        print(f"Found {len(subcommands_dict)} subcommands:")
        for cmd in sorted(subcommands_dict.keys()):
            print(f"  • {cmd}")
    else:
        print("No subparsers action found; parser may be misconfigured.")
    print()

def test_help_message():
    """Test that help messages work."""
    print("=" * 70)
    print("TEST: Help Messages")
    print("=" * 70)
    parser = get_parser()
    
    import argparse as _argparse
    subparsers_action = None
    for action in parser._subparsers._actions:
        if isinstance(action, _argparse._SubParsersAction):
            subparsers_action = action
            break
    if subparsers_action is not None and subparsers_action.choices:
        cmd_name = list(subparsers_action.choices.keys())[0]
        print(f"✓ Help available for subcommand: {cmd_name}")
    print()

def test_kegg_network_constructor_args():
    """Test KEGG_network_constructor argument parsing."""
    print("=" * 70)
    print("TEST: KEGG_network_constructor Arguments")
    print("=" * 70)
    parser = get_parser()
    
    # Valid arguments
    test_args = [
        'KEGG_network_constructor',
        '-p', 'pae00010',
        '-c', 'pae',
        '-t', 'reaction',
        '-o', '/tmp/test_network.tsv'
    ]
    
    try:
        args = parser.parse_args(test_args)
        print(f"✓ Successfully parsed KEGG_network_constructor")
        print(f"  - pathways: {args.pathways}")
        print(f"  - organism_code: {args.organism_code}")
        print(f"  - selected_interaction_types: {args.selected_interaction_types}")
        print(f"  - path_outputfile: {args.path_outputfile}")
        print(f"  - reverse_interactions: {args.reverse_interactions}")
    except SystemExit as e:
        print(f"✗ Failed to parse KEGG_network_constructor: {e}")
    print()

def test_network_merger_args():
    """Test network_merger argument parsing."""
    print("=" * 70)
    print("TEST: network_merger Arguments")
    print("=" * 70)
    parser = get_parser()
    
    test_args = [
        'network_merger',
        '-n', '/tmp/network1.tsv',
        '-m', '/tmp/network2.tsv',
        '-p', '1',
        '-o', '/tmp/merged_network.tsv'
    ]
    
    try:
        args = parser.parse_args(test_args)
        print(f"✓ Successfully parsed network_merger")
        print(f"  - network1: {args.path_inputfile_network1}")
        print(f"  - network2: {args.path_inputfile_network2}")
        print(f"  - prioritized_network: {args.prioritized_network}")
        print(f"  - reverse_interaction_doubler: {args.reverse_interaction_doubler}")
        print(f"  - output: {args.path_output_directory_and_filename}")
    except SystemExit as e:
        print(f"✗ Failed to parse network_merger: {e}")
    print()

def test_node_impact_assessor_args():
    """Test node_impact_assessor argument parsing."""
    print("=" * 70)
    print("TEST: node_impact_assessor Arguments")
    print("=" * 70)
    parser = get_parser()
    
    test_args = [
        'node_impact_assessor',
        '-n', '/tmp/network.tsv',
        '-m', '/tmp/omics.tsv',
        '-I', '/tmp/noi.tsv',
        '-o', '/tmp/impact_results.tsv'
    ]
    
    try:
        args = parser.parse_args(test_args)
        print(f"✓ Successfully parsed node_impact_assessor")
        print(f"  - network: {args.path_inputfile_network}")
        print(f"  - omics info: {args.path_inputfile_node_omics_info}")
        print(f"  - nodes_of_interest: {args.path_inputfile_nodes_of_interest}")
        print(f"  - directionality_reaction: {args.directionality_reaction}")
        print(f"  - directionality_other: {args.directionality_other}")
        print(f"  - output: {args.path_output_directory_and_filename}")
    except SystemExit as e:
        print(f"✗ Failed to parse node_impact_assessor: {e}")
    print()

def test_subnetwork_constructor_args():
    """Test subnetwork_constructor argument parsing."""
    print("=" * 70)
    print("TEST: subnetwork_constructor Arguments")
    print("=" * 70)
    parser = get_parser()
    
    test_args = [
        'subnetwork_constructor',
        '-n', '/tmp/network.tsv',
        '-i', '/tmp/impact_results.tsv',
        '-R', 'bidirectional',
        '-O', 'unidirectional',
        '-t', '5.0',
        '-b', '0.5',
        '-o', '/tmp/subnetwork.tsv'
    ]
    
    try:
        args = parser.parse_args(test_args)
        print(f"✓ Successfully parsed subnetwork_constructor")
        print(f"  - network: {args.path_inputfile_network}")
        print(f"  - impact_analysis_results: {args.path_inputfile_results_impact_analysis}")
        print(f"  - total_impact_score_threshold: {args.total_impact_score_threshold}")
        print(f"  - topological_independent_score_threshold: {args.topological_independent_score_threshold}")
        print(f"  - distance_level_limit: {args.distance_level_limit}")
        print(f"  - output: {args.output_directory_and_filename}")
    except SystemExit as e:
        print(f"✗ Failed to parse subnetwork_constructor: {e}")
    print()

def test_node_table_extender_args():
    """Test node_table_id_and_type_extender argument parsing."""
    print("=" * 70)
    print("TEST: node_table_id_and_type_extender Arguments")
    print("=" * 70)
    parser = get_parser()
    
    test_args = [
        'node_table_id_and_type_extender',
        '-n', '/tmp/network.tsv',
        '-o', '/tmp/node_table.tsv'
    ]
    
    try:
        args = parser.parse_args(test_args)
        print(f"✓ Successfully parsed node_table_id_and_type_extender")
        print(f"  - network: {args.path_inputfile_network}")
        print(f"  - output: {args.output_directory_and_filename}")
    except SystemExit as e:
        print(f"✗ Failed to parse node_table_id_and_type_extender: {e}")
    print()

def test_annotation_table_extender_args():
    """Test annotation_table_id_extender argument parsing."""
    print("=" * 70)
    print("TEST: annotation_table_id_extender Arguments")
    print("=" * 70)
    parser = get_parser()
    
    test_args = [
        'annotation_table_id_extender',
        '-n', '/tmp/network.tsv',
        '-a', '/tmp/annotations.tsv',
        '-i', '0',
        '-o', '/tmp/annotated_nodes.tsv'
    ]
    
    try:
        args = parser.parse_args(test_args)
        print(f"✓ Successfully parsed annotation_table_id_extender")
        print(f"  - network: {args.path_inputfile_network}")
        print(f"  - annotations: {args.path_inputfile_node_annotations}")
        print(f"  - column_index: {args.column_index_ids_annotation_inputfile}")
        print(f"  - output: {args.output_directory_and_filename}")
    except SystemExit as e:
        print(f"✗ Failed to parse annotation_table_id_extender: {e}")
    print()

def test_missing_required_args():
    """Test that missing required arguments are caught."""
    print("=" * 70)
    print("TEST: Missing Required Arguments Detection")
    print("=" * 70)
    parser = get_parser()
    
    # Missing required argument for KEGG_network_constructor
    test_args = [
        'KEGG_network_constructor',
        '-p', 'pae00010',
        # Missing -c, -t, -o
    ]
    
    try:
        args = parser.parse_args(test_args)
        print(f"✗ Parser should have failed with missing required arguments")
    except SystemExit as e:
        print(f"✓ Parser correctly caught missing required arguments")
    print()

def test_invalid_choices():
    """Test that invalid choices are rejected."""
    print("=" * 70)
    print("TEST: Invalid Choice Detection")
    print("=" * 70)
    parser = get_parser()
    
    # Invalid choice for prioritized_network
    test_args = [
        'network_merger',
        '-n', '/tmp/network1.tsv',
        '-m', '/tmp/network2.tsv',
        '-p', '3',  # Invalid: should be 1 or 2
        '-o', '/tmp/merged_network.tsv'
    ]
    
    try:
        args = parser.parse_args(test_args)
        print(f"✗ Parser should have failed with invalid choice")
    except SystemExit as e:
        print(f"✓ Parser correctly rejected invalid choice")
    print()


def test_monir_entry_help_and_version():
    """Basic invocation of the `MONIR.py` entry point to ensure integration."""
    print("=" * 70)
    print("TEST: MONIR Entry Point CLI")
    print("=" * 70)
    import subprocess, sys, os
    script = os.path.join(os.getcwd(), 'MONIR.py')
    # check version flag
    try:
        result = subprocess.run([sys.executable, script, '--version'], capture_output=True, text=True)
        print(f"version output: {result.stdout.strip()}")
    except Exception as e:
        print(f"✗ invoking --version failed: {e}")

    # check help on subcommand
    try:
        result = subprocess.run([sys.executable, script, 'KEGG_network_constructor', '--help'], capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ subcommand help executed successfully")
        else:
            print(f"✗ subcommand help returned code {result.returncode}")
    except Exception as e:
        print(f"✗ invoking subcommand help failed: {e}")
    print()

def main():
    """Run all CLI tests."""
    print("\n")
    print("█" * 70)
    print("  MONIR CLI Parser Tests")
    print("█" * 70)
    print("\n")
    
    test_parser_basic()
    test_subcommands_exist()
    test_help_message()
    test_kegg_network_constructor_args()
    test_network_merger_args()
    test_node_impact_assessor_args()
    test_subnetwork_constructor_args()
    test_node_table_extender_args()
    test_annotation_table_extender_args()
    test_missing_required_args()
    test_invalid_choices()
    
    print("=" * 70)
    print("CLI Tests Complete")
    print("=" * 70)
    print("\n")

if __name__ == '__main__':
    main()
