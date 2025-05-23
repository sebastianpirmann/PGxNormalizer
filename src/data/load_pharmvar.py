import os
import json
import argparse
import pandas as pd
from typing import Dict, Any, Optional, List

def load_and_process_gene_pharmvar(gene_name: str, tsv_path: str) -> Optional[Dict[str, Any]]:
    """
    Loads and processes a single PharmVar haplotype definitions TSV file for a given gene.
    Extracts allele information (name and defining variants).
    """
    print(f"Processing PharmVar data for {gene_name} from {tsv_path}...")
    
    try:
        # We expect the TSV to be tab-separated and comments starting with '#'
        df = pd.read_csv(tsv_path, sep='\t', comment='#', dtype=str)
    except FileNotFoundError:
        print(f"Error: PharmVar TSV file not found for {gene_name} at {tsv_path}. Skipping this gene.")
        return None
    except Exception as e:
        print(f"Error reading TSV for {gene_name}: {e}. Skipping this gene.")
        return None

    gene_data = {"definitions": {}} # No phenotype_mapping here, will be external

    # Create a mapping for primary alleles and their aliases/sub-alleles
    for index, row in df.iterrows():
        allele_name_raw = str(row.get('Haplotype Name')).strip() # Corrected column name
        
        if not allele_name_raw:
            continue # Skip malformed rows

        # Normalize the allele name to the standard * format
        # Remove gene prefix if present (e.g., 'CYP2D6*1' -> '*1')
        normalized_allele = allele_name_raw.replace(gene_name, '').replace('*', '').strip()
        if not normalized_allele.startswith('*'):
            normalized_allele = '*' + normalized_allele
        
        # Store the definition, including variant details
        # Check if this is a primary definition (not just a single variant defining part of an allele)
        # We assume 'REFERENCE' or 'substitution'/'deletion' etc. in 'Type' column
        # and non-empty 'rsID'/'Variant Start' for actual variants.
        is_primary_definition = (
            str(row.get('rsID', '')).strip().lower() == 'reference' or
            (str(row.get('rsID', '')).strip() == '' and str(row.get('Variant Start', '')).strip() == '')
        )
        
        if normalized_allele not in gene_data["definitions"] or is_primary_definition:
            # Only overwrite if it's the first entry for this normalized allele,
            # or if it's explicitly a 'REFERENCE' definition (which we want as the primary)
            gene_data["definitions"][normalized_allele] = {
                "raw_pharmvar_name": allele_name_raw,
                "variant_details": [], # Will store a list of associated variants
                # Functionality will be added externally via PharmVarManager
            }
            if is_primary_definition: # Mark primary definition for later
                gene_data["definitions"][normalized_allele]["is_primary_definition"] = True


        # Add variant details to the allele definition if available
        # Haplotypes.tsv lists variants for an allele, so we collect them.
        rsid = str(row.get('rsID', '')).strip()
        variant_start = str(row.get('Variant Start', '')).strip()
        variant_type = str(row.get('Type', '')).strip()
        
        if rsid or variant_start: # If there's actual variant information
            variant_info = {
                "rsID": rsid,
                "reference_sequence": str(row.get('ReferenceSequence', '')).strip(),
                "variant_start": variant_start,
                "variant_stop": str(row.get('Variant Stop', '')).strip(),
                "reference_allele": str(row.get('Reference Allele', '')).strip(),
                "variant_allele": str(row.get('Variant Allele', '')).strip(),
                "type": variant_type
            }
            # Ensure 'variant_details' is a list
            if "variant_details" not in gene_data["definitions"][normalized_allele]:
                gene_data["definitions"][normalized_allele]["variant_details"] = []
            gene_data["definitions"][normalized_allele]["variant_details"].append(variant_info)

        # Add common aliases pointing to the normalized allele
        # This is where we account for different ways tools might name an allele.
        # This will need to be refined based on actual tool outputs.
        if allele_name_raw not in gene_data["definitions"]:
            gene_data["definitions"][allele_name_raw] = {"maps_to": normalized_allele}
        if f"{gene_name}{normalized_allele}" not in gene_data["definitions"]:
             gene_data["definitions"][f"{gene_name}{normalized_allele}"] = {"maps_to": normalized_allele}


    return gene_data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process raw PharmVar haplotype definitions into a consolidated JSON database."
    )
    parser.add_argument(
        "--ref-genome",
        type=str,
        default="GRCh37", # Default to GRCh37
        choices=["GRCh37", "GRCh38"],
        help="Specify the reference genome assembly to use (e.g., GRCh37, GRCh38)."
    )
    parser.add_argument(
        "--pharmvar-version",
        type=str,
        default="6.2.7",
        help="Specify the PharmVar database version (e.g., 6.2.7)."
    )
    parser.add_argument(
        "--raw-download-dir",
        type=str,
        default="raw_pharmvar_download", # Default subfolder name
        help="Name of the sub-directory within src/data where the raw PharmVar download resides."
    )

    args = parser.parse_args()

    # Extract arguments
    selected_ref_genome = args.ref_genome
    pharmvar_version = args.pharmvar_version
    raw_download_dir_name = args.raw_download_dir

    core_panel = [
        "CYP2D6", "CYP2B6", "CYP3A5", "TPMT", "CYP2C9", "CYP2C19",
        "DPYD", "UGT1A1", "SLCO1B1", "CYP3A4", "VKORC1", "F5", "CYP1A2"
    ]

    pharmvar_root_dir = os.path.join(os.path.dirname(__file__), "data")
    output_json_path = os.path.join(pharmvar_root_dir, "pharmvar_processed.json")

    all_pharmvar_data = {
        "metadata": {
            "pharmvar_version": pharmvar_version,
            "reference_genome": selected_ref_genome,
            "processed_on": pd.Timestamp.now().isoformat(),
            "core_panel_genes": core_panel
        },
        "genes": {}
    }

    for gene in core_panel:
        # Path now includes raw_download_dir_name, pharmvar_version, gene, and selected_ref_genome
        gene_specific_dir = os.path.join(
            pharmvar_root_dir,
            raw_download_dir_name,
            f"pharmvar-{pharmvar_version}", # Dynamic version in path
            gene,
            selected_ref_genome
        )

        if not os.path.isdir(gene_specific_dir):
            print(f"Warning: Gene specific directory not found for {gene} ({selected_ref_genome}) at {gene_specific_dir}. Skipping.")
            continue

        tsv_candidates = [f for f in os.listdir(gene_specific_dir) if f.startswith(f"{gene}.NC_") and f.endswith(".haplotypes.tsv")]

        if not tsv_candidates:
            print(f"Warning: No .haplotypes.tsv file found for {gene} in {gene_specific_dir}. Skipping.")
            continue

        tsv_filename = tsv_candidates[0]
        tsv_path = os.path.join(gene_specific_dir, tsv_filename)

        gene_processed_data = load_and_process_gene_pharmvar(gene, tsv_path)

        if gene_processed_data:
            all_pharmvar_data["genes"][gene] = gene_processed_data
        else:
            print(f"Warning: No PharmVar data successfully loaded for {gene}.")

    try:
        with open(output_json_path, "w") as f:
            json.dump(all_pharmvar_data, f, indent=4)
        print(f"\nSuccessfully processed PharmVar data for {len(all_pharmvar_data['genes'])} genes.")
        print(f"Aggregated data saved to: {output_json_path}")
        print(f"Reference Genome used: {selected_ref_genome}")
    except Exception as e:
        print(f"Error saving processed PharmVar data: {e}")