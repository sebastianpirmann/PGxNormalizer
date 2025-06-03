import json
import os
from typing import Dict, List, Optional, Any

class PharmVarManager:
    """
    Manages access to processed PharmVar allele definitions (from raw_pharmvar_download)
    and provides manually curated functional annotations and phenotype mappings.
    """
    def __init__(self, processed_pharmvar_path: str = "src/data/pharmvar_processed.json"):
        """
        Initializes the PharmVarManager by loading the pre-processed PharmVar data (allele definitions).
        """
        self.processed_pharmvar_path = processed_pharmvar_path
        self.pharmvar_db: Dict[str, Any] = self._load_pharmvar_db()
        self.supported_genes: List[str] = list(self.pharmvar_db.get("genes", {}).keys()) # Adjusted for new JSON structure
        print(f"PharmVarManager initialized. Loaded allele definitions for {len(self.supported_genes)} genes.")

    # --- Manually Curated Mappings for Functionality and Phenotype ---
    # These are crucial and based on CPIC guidelines and PharmVar functional annotations.
    # This is a simplified example; a full implementation would cover all known alleles
    # and nuances, potentially loading from a separate configuration file for larger projects.
    _allele_functionality_map = {
        "CYP2D6": {
            "*1": "Normal Function", "*2": "Normal Function", "*10": "Decreased Function",
            "*4": "No Function", "*5": "No Function", "*6": "No Function",
            "*xN": "Increased Function", # For gene duplications like *1x2, *2xN
            "UNKNOWN": "Unknown" # Placeholder for alleles not found in this map
        },
        "CYP2C9": {
            "*1": "Normal Function", "*2": "Decreased Function", "*3": "No Function",
            "UNKNOWN": "Unknown"
        },
        "CYP2C19": {
            "*1": "Normal Function", "*2": "No Function", "*3": "No Function",
            "*17": "Increased Function",
            "UNKNOWN": "Unknown"
        },
        "CYP2B6": {
            "*1": "Normal Function", "*6": "Decreased Function", "*18": "No Function",
            "UNKNOWN": "Unknown"
        },
        "CYP3A5": {
            "*1": "Normal Function", "*3": "No Function",
            "UNKNOWN": "Unknown"
        },
        "TPMT": {
            "*1": "Normal Function", "*2": "Decreased Function", "*3A": "No Function",
            "*3B": "No Function", "*3C": "No Function",
            "UNKNOWN": "Unknown"
        },
        "DPYD": {
            "*1": "Normal Function", "*2A": "No Function", # Common alias for *2A
            "UNKNOWN": "Unknown"
        },
        "UGT1A1": {
            "*1": "Normal Function", "*6": "Decreased Function", "*28": "Decreased Function",
            "UNKNOWN": "Unknown"
        },
        "SLCO1B1": { # Note: SLCO1B1 function is more complex; simplified here
            "*1A": "Normal Function", "*5": "Decreased Function", "*15": "Decreased Function",
            "UNKNOWN": "Unknown"
        },
        "CYP3A4": { # Most common functional variants
            "*1": "Normal Function", "*22": "Decreased Function",
            "UNKNOWN": "Unknown"
        },
        "VKORC1": { # Not a metabolizer, but relevant for dosing. Assuming variant-level mapping if no star allele.
            "*1": "Normal Function", # Assuming *1 exists or represents WT
            "rs9923231": "Decreased Function", # Common variant for VKORC1
            "UNKNOWN": "Unknown"
        },
        "F5": { # Not a metabolizer gene, no star alleles. Need to handle as variant.
            "*1": "Normal Function", # Assuming *1 exists or represents WT
            "rs6025": "Increased Function", # Factor V Leiden
            "UNKNOWN": "Unknown"
        },
        "CYP1A2": { # Added for your example
            "*1": "Normal Function", "*1F": "Increased Function", "*1C": "Decreased Function",
            "UNKNOWN": "Unknown"
        }
        # Add other genes and their common alleles as needed
    }

    # Standardized Phenotype Abbreviations (CPIC-like)
    _phenotype_map = {
        "Normal Function": "NM",
        "Increased Function": "UM",
        "Decreased Function": "IM",
        "No Function": "PM",
        "Unknown": "Unknown",
        "Uncertain Function": "Indeterminate",
        # For non-metabolizer genes, you might use different terms for functionality
        "Normal": "Normal", # Used for VKORC1, F5 if their function is just "Normal"
        "Reduced": "Reduced",
        "Increased": "Increased",
        "Non-functional": "Non-functional",
    }


    def _load_pharmvar_db(self) -> Dict[str, Any]:
        """Loads the processed PharmVar allele definition data from the JSON file."""
        if not os.path.exists(self.processed_pharmvar_path):
            raise FileNotFoundError(
                f"Processed PharmVar DB not found at '{self.processed_pharmvar_path}'. "
                "Please run 'src/data/load_pharmvar.py' first."
            )
        try:
            with open(self.processed_pharmvar_path, 'r') as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Error decoding PharmVar JSON from '{self.processed_pharmvar_path}': {e}")
        except Exception as e:
            raise RuntimeError(f"Unexpected error loading PharmVar DB: {e}")

    def get_gene_info(self, gene_name: str) -> Optional[Dict[str, Any]]:
        """Returns all PharmVar allele definition information for a given gene."""
        # Adjusted for new JSON structure: access via "genes" key
        return self.pharmvar_db.get("genes", {}).get(gene_name)

    def get_normalized_allele(self, gene_name: str, raw_allele: str) -> Optional[str]:
        """
        Normalizes a raw allele string to its PharmVar standard nomenclature (e.g., '*1', '*4').
        Handles aliases defined in the processed PharmVar DB.
        Returns "UNKNOWN" if no mapping is found.
        """
        if not gene_name or not raw_allele:
            return "UNKNOWN" # Return "UNKNOWN" for invalid inputs

        gene_data = self.pharmvar_db.get("genes", {}).get(gene_name) # Adjusted for new JSON structure
        if not gene_data:
            return "UNKNOWN"

        definitions = gene_data.get("definitions", {})

        # 1. Try direct lookup for exact match (e.g., if tool output is already normalized or an alias)
        if raw_allele in definitions:
            if "maps_to" in definitions[raw_allele]:
                return definitions[raw_allele]["maps_to"] # Return target if it's an alias
            else:
                # It's a primary definition, but ensure it's in the * form
                normalized_form = raw_allele.replace(gene_name, '').replace('*', '').strip()
                if not normalized_form.startswith('*'):
                    normalized_form = '*' + normalized_form
                return normalized_form

        # 2. Try normalizing by removing gene prefix and checking for '*'
        # This handles cases like 'CYP2D6*4' becoming '*4'
        normalized_form = raw_allele.replace(gene_name, '').replace('*', '').strip()
        if not normalized_form.startswith('*'):
            normalized_form = '*' + normalized_form
        
        if normalized_form in definitions:
            if "maps_to" in definitions[normalized_form]:
                return definitions[normalized_form]["maps_to"]
            else:
                return normalized_form
        
        # 3. Handle specific non-star allele cases like rsIDs for VKORC1, F5 etc.
        # If the raw_allele matches a known variant in our *functionality map*,
        # we treat it as a normalized allele for functionality lookup.
        if gene_name in self._allele_functionality_map and \
           raw_allele in self._allele_functionality_map[gene_name]:
            return raw_allele # Return the raw allele as the 'normalized' form for these cases

        # If still not found, return "UNKNOWN"
        return "UNKNOWN"


    def get_allele_functionality(self, gene_name: str, normalized_allele: str) -> Optional[str]:
        """
        Returns the functional consequence (e.g., 'Normal Function', 'No Function')
        for a given normalized PharmVar allele using the internal mapping.
        Returns "Unknown" if no mapping is found.
        """
        if gene_name not in self._allele_functionality_map:
            print(f"Warning: No functionality mapping defined for gene '{gene_name}'. Defaulting to 'Unknown'.")
            return "Unknown"
        
        # Try to get direct allele functionality (e.g., *1, *4, or a specific rsID for non-star allele genes)
        functionality = self._allele_functionality_map[gene_name].get(normalized_allele)
        if functionality:
            return functionality
        
        # Handle allele groups like *xN for duplications (e.g., *1x2, *2xN)
        # This is a simplification; a robust solution might need a regex or more complex logic.
        if gene_name == "CYP2D6" and 'x' in normalized_allele:
            # Check for patterns like *1X2, *2XN, *1x3 etc.
            if 'x' in normalized_allele.lower():
                return self._allele_functionality_map[gene_name].get("*xN", "Unknown")

        # Fallback for alleles not explicitly defined in the map
        return self._allele_functionality_map[gene_name].get("UNKNOWN", "Unknown")


    def get_standard_phenotype(self, functionality: str) -> Optional[str]:
        """
        Maps a functional consequence (e.g., 'Normal Function') to a standardized
        phenotype abbreviation (e.g., 'NM', 'PM') using the internal mapping.
        Returns "Unknown" if no mapping is found.
        """
        # Use the generic phenotype map directly
        return self._phenotype_map.get(functionality, "Unknown")


# Example usage (for testing this module directly)
if __name__ == "__main__":
    # Ensure src/data/pharmvar_processed.json exists by running load_pharmvar.py first
    
    print("--- Testing PharmVarManager ---")
    try:
        manager = PharmVarManager()

        # Test CYP2D6
        print("\n--- Testing CYP2D6 ---")
        gene = "CYP2D6"
        print(f"Supported genes: {manager.supported_genes}")
        
        # Test normalization
        print(f"Normalized 'CYP2D6*1': {manager.get_normalized_allele(gene, 'CYP2D6*1')}") # Should be *1
        print(f"Normalized '*4.001': {manager.get_normalized_allele(gene, '*4.001')}")   # Should be *4.001
        print(f"Normalized '*5': {manager.get_normalized_allele(gene, '*5')}")       # Should be *5
        print(f"Normalized '2D6*10': {manager.get_normalized_allele(gene, '2D6*10')}") # Should be *10
        print(f"Normalized 'UNKNOWN_ALLELE': {manager.get_normalized_allele(gene, 'UNKNOWN_ALLELE')}") # Should be UNKNOWN
        print(f"Normalized 'CYP2D6*1X2': {manager.get_normalized_allele(gene, 'CYP2D6*1X2')}") # Should be *1X2

        # Test functionality and phenotype mapping
        norm_allele_1 = manager.get_normalized_allele(gene, 'CYP2D6*1')
        func_1 = manager.get_allele_functionality(gene, norm_allele_1)
        phen_1 = manager.get_standard_phenotype(func_1) # Removed gene param from get_standard_phenotype call
        print(f"Allele {norm_allele_1}: Functionality '{func_1}', Phenotype '{phen_1}'")

        norm_allele_4 = manager.get_normalized_allele(gene, '*4')
        func_4 = manager.get_allele_functionality(gene, norm_allele_4)
        phen_4 = manager.get_standard_phenotype(func_4) # Removed gene param
        print(f"Allele {norm_allele_4}: Functionality '{func_4}', Phenotype '{phen_4}'")

        norm_allele_dup = manager.get_normalized_allele(gene, 'CYP2D6*1X2') # Example for duplication
        func_dup = manager.get_allele_functionality(gene, norm_allele_dup)
        phen_dup = manager.get_standard_phenotype(func_dup) # Removed gene param
        print(f"Allele {norm_allele_dup}: Functionality '{func_dup}', Phenotype '{phen_dup}'")


        # Test TPMT
        print("\n--- Testing TPMT ---")
        gene = "TPMT"
        print(f"Normalized 'TPMT*2': {manager.get_normalized_allele(gene, 'TPMT*2')}")
        func_tpmt2 = manager.get_allele_functionality(gene, manager.get_normalized_allele(gene, 'TPMT*2'))
        print(f"Functionality of 'TPMT*2': {func_tpmt2}")
        print(f"Standard phenotype of 'No Function' for TPMT: {manager.get_standard_phenotype('No Function')}") # Removed gene param

        # Test VKORC1 (example for non-star allele gene)
        print("\n--- Testing VKORC1 ---")
        gene = "VKORC1"
        print(f"Normalized 'rs9923231': {manager.get_normalized_allele(gene, 'rs9923231')}")
        func_vkorc1 = manager.get_allele_functionality(gene, manager.get_normalized_allele(gene, 'rs9923231'))
        print(f"Functionality of 'rs9923231': {func_vkorc1}")
        print(f"Standard phenotype of 'Decreased Function' for VKORC1: {manager.get_standard_phenotype('Decreased Function')}") # Removed gene param

        # Test unknown gene
        print("\n--- Testing Unknown Gene ---")
        gene = "UNKNOWN_GENE"
        raw_allele = "*1"
        norm_allele = manager.get_normalized_allele(gene, raw_allele)
        func = manager.get_allele_functionality(gene, norm_allele)
        phen = manager.get_standard_phenotype(func)
        print(f"Gene '{gene}', Raw Allele '{raw_allele}': Normalized '{norm_allele}', Functionality '{func}', Phenotype '{phen}'")

    except FileNotFoundError as e:
        print(f"Setup Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")