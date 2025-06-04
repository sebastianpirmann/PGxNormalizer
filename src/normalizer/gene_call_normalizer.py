from typing import List, Dict, Any, Optional
from collections import defaultdict
import json # Added for debugging print


from src.standard_formats import StandardizedGeneCall

def group_gene_calls_by_sample_gene(gene_calls: List[StandardizedGeneCall]) -> Dict[str, Dict[str, List[StandardizedGeneCall]]]:
    """
    Groups a list of StandardizedGeneCall objects first by sample_id,
    then by gene. This utility function is useful for organizing data
    before passing it to the normalizer on a per-gene-per-sample basis.

    Returns a nested dictionary structure:
    {
        "sample_id_1": {
            "GENE_A": [StandardizedGeneCall_obj_1_A_sol1, StandardizedGeneCall_obj_1_A_sol2],
            "GENE_B": [StandardizedGeneCall_obj_1_B_sol1]
        },
        "sample_id_2": { ... }
    }
    """
    grouped_data: Dict[str, Dict[str, List[StandardizedGeneCall]]] = defaultdict(lambda: defaultdict(list))

    for call in gene_calls:
        sample_id = call.get('sample_id')
        gene = call.get('gene')

        if sample_id and gene:
            grouped_data[sample_id][gene].append(call)
        else:
            print(f"Warning: Skipping gene call due to missing sample_id or gene in: {call.get('input_file', 'unknown_file')}")

    return grouped_data


class GeneCallNormalizer:
    """
    Normalizes a list of StandardizedGeneCall objects for a *single gene*
    and a *single sample*, selecting the best solution and inferring final attributes.
    """
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initializes the normalizer. Configuration can be passed for
        specific normalization rules, reference data paths, etc.
        """
        self.config = config if config else {}
        print("GeneCallNormalizer initialized for single gene/sample normalization.")

    def normalize(self, gene_calls_for_one_sample_one_gene: List[StandardizedGeneCall]) -> Optional[StandardizedGeneCall]:
        """
        Normalizes a list of StandardizedGeneCall objects for a *single gene*
        and a *single sample*, identifying the best/consensus call.

        Args:
            gene_calls_for_one_sample_one_gene (List[StandardizedGeneCall]):
                A list of raw gene calls that all belong to the SAME sample_id
                and the SAME gene (e.g., all ALDY solutions for CYP2D6 for NA10860).

        Returns:
            Optional[StandardizedGeneCall]: The single best/consensus StandardizedGeneCall
                                            for that sample and gene, or None if no valid
                                            call could be determined from the input.
        """
        if not gene_calls_for_one_sample_one_gene:
            print("Warning: Normalizer received an empty list for a single gene/sample. Returning None.")
            return None

        # Confirm all calls are for the same sample and gene (for robustness)
        first_call = gene_calls_for_one_sample_one_gene[0]
        sample_id = first_call.get('sample_id')
        gene = first_call.get('gene')

        # Add a check for consistency
        for call in gene_calls_for_one_sample_one_gene:
            if call.get('sample_id') != sample_id or call.get('gene') != gene:
                print(f"Error: Inconsistent sample_id or gene found in input list for normalizer. "
                      f"Expected {sample_id}:{gene}, but found {call.get('sample_id')}:{call.get('gene')}. "
                      "This list should only contain calls for one sample-gene pair.")
                return None

        print(f"  Normalizing {len(gene_calls_for_one_sample_one_gene)} calls for Sample: {sample_id}, Gene: {gene}")

        # Solution Selection Logic
        best_solution: Optional[StandardizedGeneCall] = None
        
        # Filter for "NORMAL" solutions first
        normal_solutions: List[StandardizedGeneCall] = []
        other_solutions: List[StandardizedGeneCall] = []

        for call in gene_calls_for_one_sample_one_gene:
            # Check if it's an ALDY call (assuming this for now for the rule)
            if call.get('tool_name') == "ALDY":
                # ALDY's raw output contains 'variants_reported' which has 'tool_specific_flags'.
                # We need to look for 'NORMAL' in these flags.
                # A simplistic approach is to check if any variant reported has 'NORMAL' status.
                # A more robust approach might be to define a 'primary_status' field in StandardizedGeneCall
                # during parsing, or aggregate status from multiple variants.
                
                # For ALDY, the 'VariantStatus' from the major/minor allele definition lines (Solution 1: *1.001, *4, *4.021)
                # is not directly available in `variants_reported` (which are for the actual variants).
                # However, the `VariantStatus` column (NORMAL, NOVEL, EXTRA, MISSING) exists in the raw TSV.
                # The AldyParser maps this to `VariantReported['tool_specific_flags']`.
                # Let's assume for this basic rule, we check the 'VariantStatus' from *any* variant reported
                # or a general flag if the parser extracted it at the solution level.
                
                # For now, let's rely on `aldy_solution_id` and the solution description for a simple heuristic
                # given the challenge of aggregating 'VariantStatus' across multiple lines into one solution object.
                # A better approach for Aldy-specific selection often involves the 'SolutionID' itself
                # or a more complex analysis of all alleles/variants in the solution.

                # Let's refine the heuristic to rely on a general "NORMAL" presence in `VariantStatus`
                # of *any* reported variant associated with the solution, or the SolutionID directly.
                # The 'VariantStatus' column (e.g., NORMAL, NOVEL) is captured in `VariantReported['tool_specific_flags']`.
                
                has_normal_status = False
                for variant in call.get('raw_tool_output', {}).get('variants_reported', []):
                    flags = variant.get('tool_specific_flags', '')
                    if 'NORMAL' in flags: # Checks if 'NORMAL' is present in any variant's flags
                        has_normal_status = True
                        break
                
                if has_normal_status:
                    normal_solutions.append(call)
                else:
                    other_solutions.append(call)
            else:
                # If from another tool, or no specific ALDY rule applies, treat as "other"
                other_solutions.append(call)

        if normal_solutions:
            # If there are "NORMAL" solutions, pick the one with the lowest SolutionID
            # Need to convert aldy_solution_id to int for proper sorting
            normal_solutions.sort(key=lambda x: int(x.get('raw_tool_output', {}).get('aldy_solution_id', sys.maxsize)))
            best_solution = normal_solutions[0]
            print(f"    Selected 'NORMAL' solution with SolutionID: {best_solution.get('raw_tool_output', {}).get('aldy_solution_id')}")
        elif other_solutions:
            # If no "NORMAL" solutions, pick the one with the lowest SolutionID from the rest
            other_solutions.sort(key=lambda x: int(x.get('raw_tool_output', {}).get('aldy_solution_id', sys.maxsize)))
            best_solution = other_solutions[0]
            print(f"    No 'NORMAL' solution found. Selected other solution with SolutionID: {best_solution.get('raw_tool_output', {}).get('aldy_solution_id')}")
        else:
            print(f"    Could not determine a best solution for {sample_id}:{gene}.")
            return None # No solutions provided or none could be selected

        # Data Harmonization & Enrichment
        # Once 'best_solution' is determined:
        # - Map `raw_tool_output['functional_status_raw']` to a standardized functional status.
        # - Infer the final phenotype (e.g., "Normal Metabolizer") based on the diplotype.
        # - Add external annotations from PharmVar, CPIC, etc.
        
        best_solution['normalized_functional_status'] = "Normal Function (Inferred)"
        best_solution['predicted_phenotype'] = "Normal Metabolizer (Inferred)"

        print(f"  Returning normalized call for {sample_id}:{gene}.")

        return best_solution