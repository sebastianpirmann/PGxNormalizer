# pgx_normalizer/src/normalizer/gene_call_normalizer.py

from typing import List, Dict, Any, Optional
from collections import defaultdict # Useful for efficient grouping

# Import the StandardizedGeneCall TypedDict
from src.standard_formats import StandardizedGeneCall

class GeneCallNormalizer:
    """
    Normalizes a list of StandardizedGeneCall objects.
    This involves grouping calls by sample and gene, selecting a 'best' solution
    among multiple options (e.g., from different tools or multiple solutions
    from one tool like ALDY), and inferring final normalized attributes.
    """
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initializes the normalizer. Configuration can be passed for
        specific normalization rules, reference data paths, etc.
        """
        self.config = config if config else {}
        print("GeneCallNormalizer initialized.")

    def normalize(self, raw_gene_calls: List[StandardizedGeneCall]) -> List[StandardizedGeneCall]:
        """
        Processes a list of raw StandardizedGeneCall objects (potentially from multiple tools/solutions)
        and returns a refined, normalized list containing one definitive call per gene per sample.

        Args:
            raw_gene_calls (List[StandardizedGeneCall]): The aggregated list of raw gene calls
                                                          from all parsers.

        Returns:
            List[StandardizedGeneCall]: A list of normalized gene calls, where each entry
                                        represents the single best/consensus call for a gene
                                        for a given sample.
        """
        print(f"Normalizer received {len(raw_gene_calls)} raw gene calls for processing.")

        # Step 1: Group gene calls by sample_id and gene
        # This is the foundational step to consolidate all information for a specific gene-sample pair.
        grouped_calls = self._group_by_sample_gene(raw_gene_calls)

        normalized_results: List[StandardizedGeneCall] = []

        # --- Future Normalization & Consensus Logic Will Be Implemented Here ---
        # For each sample_id in the grouped data:
        #   For each gene within that sample:
        #     - 'gene_calls_for_gene' will be a list containing all StandardizedGeneCall objects
        #       for that specific gene and sample (e.g., Aldy Solution 1, Aldy Solution 2, etc.).
        #     - Here you'll implement logic to:
        #       a) **Select the best solution** from `gene_calls_for_gene` (e.g., by SolutionID, confidence, etc.).
        #          For Aldy, this means choosing which of its multiple solutions is the 'true' one.
        #       b) **Resolve conflicts** if multiple tools provided different calls for the same gene.
        #       c) **Infer additional attributes** (e.g., final copy number, functional status, phenotype)
        #          based on the selected/consensus call and potentially external knowledge bases.
        #       d) **Add clinical annotations.**
        # -----------------------------------------------------------------------

        for sample_id, genes_data in grouped_calls.items():
            print(f"\n  Normalizing for Sample: {sample_id}")
            for gene, gene_calls_for_gene in genes_data.items():
                print(f"    Processing Gene: {gene} (found {len(gene_calls_for_gene)} raw solutions/entries)")

                # --- Placeholder for actual Normalization Logic per (sample, gene) ---
                # For now, as a starting point, we'll just pick the first solution found.
                # In a real scenario, this is where complex rules would go.
                if gene_calls_for_gene:
                    best_solution = gene_calls_for_gene[0] # Very simplistic selection: just take the first one
                    # Example of adding a new normalized field (this would be based on real logic)
                    best_solution['normalized_status_placeholder'] = "Normalized and Selected"
                    normalized_results.append(best_solution)
                # -------------------------------------------------------------------

        print(f"\nNormalizer finished. Returning {len(normalized_results)} normalized gene calls (one per sample-gene pair).")
        return normalized_results

    def _group_by_sample_gene(self, gene_calls: List[StandardizedGeneCall]) -> Dict[str, Dict[str, List[StandardizedGeneCall]]]:
        """
        Helper method to group a list of StandardizedGeneCall objects first by sample_id,
        then by gene. This facilitates processing all information relevant to a
        specific gene in a given sample.

        Returns a nested dictionary structure:
        {
            "sample_id_1": {
                "GENE_A": [StandardizedGeneCall_obj_1_A_sol1, StandardizedGeneCall_obj_1_A_sol2],
                "GENE_B": [StandardizedGeneCall_obj_1_B_sol1]
            },
            "sample_id_2": {
                "GENE_A": [StandardizedGeneCall_obj_2_A_sol1],
                "GENE_C": [StandardizedGeneCall_obj_2_C_sol1, StandardizedGeneCall_obj_2_C_sol2]
            }
        }
        """
        grouped_data: Dict[str, Dict[str, List[StandardizedGeneCall]]] = defaultdict(lambda: defaultdict(list))

        for call in gene_calls:
            sample_id = call.get('sample_id')
            gene = call.get('gene')

            if sample_id and gene:
                grouped_data[sample_id][gene].append(call)
            else:
                # This should ideally be caught earlier or indicate malformed data
                print(f"Warning: Skipping gene call due to missing sample_id or gene in: {call.get('input_file', 'unknown_file')}")

        return grouped_data