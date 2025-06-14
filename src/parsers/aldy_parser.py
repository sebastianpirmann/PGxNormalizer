import sys
import pandas as pd
import os
import json
import io   # To treat string as a file
import re   # For regular expressions, to parse solution description

from typing import List, Dict, Any, Optional, Union

# Import standardized formats
from src.standard_formats import StandardizedGeneCall, RawToolOutput, VariantReported, RawAlleleComponent
# Import the BaseParser for interface adherence
from src.parsers.base_parser import BaseParser


class AldyParser(BaseParser):
    """
    Parses output files generated by the ALDY genotyping tool
    and converts them into a standardized intermediate format.
    Handles multiple solutions per sample/gene.
    """
    def __init__(self):
        super().__init__()
        self.tool_name = "ALDY"
        # Define expected raw columns to help with validation if needed
        self.expected_raw_columns = [
            'Sample', 'Gene', 'SolutionID', 'Major', 'Minor',
            'AlleleCopyIdentifier', 'Allele', 'Location', 'VariantType',
            'Coverage', 'VariantFunctionalityRaw', 'dbSNP', 'KarolinskaCode', 'VariantStatus'
        ]
        print(f"Initialized {self.tool_name} Parser.")

    def _infer_reference_genome_from_location(self, location: str) -> str:
        """
        Heuristically infers the reference genome based on variant location.
        This is a simple guess; more robust solutions might require explicit
        reference genome input or more complex logic.
        """
        try:
            loc_int = int(location)
            # Example: A common range for human chromosomes (e.g., chr22 where CYP2D6 is)
            # on GRCh37/hg19. This is a very broad and simplistic heuristic.
            if 1_000_000 <= loc_int <= 250_000_000:
                 return "GRCh37" # Assuming GRCh37 based on typical PharmVar/CYP2D6 data
        except (ValueError, TypeError):
            pass # Not an integer location, or location is missing
        return "UNKNOWN"


    def _infer_copy_number_from_diplotype(self, diplotype_string: str) -> Optional[Union[int, float]]:
        """
        Heuristically infers a copy number from a diplotype string.
        This is a simplified approach and may not capture all complexities (e.g., Xn/Ym).
        """
        if not diplotype_string:
            return None
        # Clean the diplotype string to count distinct allele components
        # Example: CYP2D6*1/*4+*4.021 -> *1/*4+*4.021 -> 1, 4, 4.021 (3 copies)
        # Handle various common delimiters and remove gene prefixes
        cleaned_diplotype = diplotype_string.replace('CYP2D6', '').replace('CYP2C9', '').replace('CYP2C19', '') # TODO: Implement general solution
        # Split by '/', '+', and filter out empty strings
        alleles_list = [a.strip() for a in re.split(r'[/\+]', cleaned_diplotype) if a.strip()]
        return len(alleles_list)


    def parse(self, filepath: str) -> List[StandardizedGeneCall]:
        """
        Parses a single ALDY output file (TSV) into a list of StandardizedGeneCall objects.
        Each ALDY "Solution" block is converted into a separate StandardizedGeneCall.

        Args:
            aldy_output_filepath (str): Path to the ALDY output file.

        Returns:
            List[StandardizedGeneCall]: A list of dictionaries, each representing a gene call
                                   in the standardized format. Returns an empty list if parsing fails.
        """
        if not os.path.exists(filepath):
            print(f"Error: ALDY output file not found at {filepath}")
            return []

        print(f"Parsing ALDY output from: {filepath}")
        
        data_lines: List[str] = []
        column_names: List[str] = []
        # Maps SolutionID (str) to its descriptive text from the #Solution line
        solution_descriptions: Dict[str, str] = {}
        
        # Temporary variables to link a #Solution line to the first data row's SolutionID
        current_solution_id_group_from_header: Optional[str] = None
        current_solution_description_text: Optional[str] = None

        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line: # Skip empty lines
                        continue

                    if line.startswith('#Sample'):
                        # This is the actual header line with column names
                        column_names = [col.strip() for col in line[1:].split('\t')] # Remove '#' and split
                        # Basic validation: Check if parsed columns match our expected
                        if not all(col in column_names for col in self.expected_raw_columns):
                            print(f"Warning: ALDY file '{filepath}' header missing expected columns. "
                                  f"Expected: {self.expected_raw_columns}, Found: {column_names}")
                        continue # Skip to next line after finding header

                    if line.startswith('#Solution '):
                        # This line contains the overall solution description
                        match = re.match(r'#Solution (\d+): (.*)', line)
                        if match:
                            current_solution_id_group_from_header = match.group(1) # e.g., '1'
                            current_solution_description_text = match.group(2)    # e.g., '*1.001, *4, *4.021'
                        continue # Skip to next line

                    # If it's a data line (i.e., not a comment and we've seen the header)
                    if column_names: # Ensure we've parsed the header first
                        data_lines.append(line)
                        
                        # If we just read a solution description from a '#' line,
                        # associate it with the SolutionID found in this first data row of the new solution block.
                        if current_solution_id_group_from_header and current_solution_description_text:
                            parts = line.split('\t')
                            # Assuming SolutionID is the 3rd column (index 2) as per the example format
                            if len(parts) > 2:
                                solution_id_from_data_row = parts[2] 
                                if solution_id_from_data_row == current_solution_id_group_from_header:
                                    solution_descriptions[solution_id_from_data_row] = current_solution_description_text
                                    # Reset temporary vars for the next solution block
                                    current_solution_id_group_from_header = None
                                    current_solution_description_text = None
                    else:
                        print(f"Warning: Data line found before header in {filepath}: {line[:50]}...")

            if not column_names:
                print(f"Error: Could not find header line starting with '#Sample' in {filepath}")
                return []
            if not data_lines:
                print(f"Warning: No data rows found in {filepath} after header.")
                return []

            # Use StringIO to read the collected data lines into a DataFrame
            # `names=column_names` ensures our cleaned column names are used.
            df = pd.read_csv(io.StringIO("\n".join(data_lines)), sep='\t', names=column_names, dtype=str)

            # Add the SolutionDescription column by mapping SolutionID
            df['SolutionDescription'] = df['SolutionID'].map(solution_descriptions)
            df['SolutionDescription'].fillna("No Description Found", inplace=True) # Fallback for any unmapped

            # --- Now, process the DataFrame to populate StandardizedGeneCall objects ---
            parsed_gene_calls: List[StandardizedGeneCall] = []
            
            # Key step: Group by Sample, Gene, and SolutionID.
            # Each group corresponds to one StandardizedGeneCall in our output list.
            grouped_solutions = df.groupby(['Sample', 'Gene', 'SolutionID'])

            for (sample_id, gene, solution_id), group_df in grouped_solutions:
                # Get the first row of the group to extract common solution-level data
                first_row = group_df.iloc[0]

                # Map ALDY columns to StandardizedGeneCall top-level fields
                ref_genome = self._infer_reference_genome_from_location(
                    first_row.get('Location', '') if pd.notna(first_row.get('Location')) else ''
                )
                
                # --- Prepare data for raw_tool_output ---
                aldy_diplotype_string = first_row.get('Major', '') # ALDY's Major column maps to diplotype_string
                estimated_cn = self._infer_copy_number_from_diplotype(aldy_diplotype_string)
                aldy_alleles_in_solution_raw_string = first_row.get('Minor', '') # ALDY's Minor column

                variants_reported: List[VariantReported] = []
                aldy_alleles_parsed: List[RawAlleleComponent] = []

                # Iterate through each row within the current solution group to extract variants and allele components
                for _, row in group_df.iterrows():
                    # Populate RawAlleleComponent (from 'Allele' and 'AlleleCopyIdentifier')
                    current_allele = row.get('Allele')
                    current_allele_copy_id_str = row.get('AlleleCopyIdentifier')
                    current_allele_copy_id = int(current_allele_copy_id_str) if pd.notna(current_allele_copy_id_str) and current_allele_copy_id_str.isdigit() else None

                    if current_allele and current_allele_copy_id is not None:
                        # Ensure uniqueness for ald_alleles_parsed to avoid duplicate entries for the same allele copy
                        if {"raw_allele_name": current_allele, "allele_copy_id": current_allele_copy_id} not in aldy_alleles_parsed:
                             aldy_alleles_parsed.append({
                                "raw_allele_name": current_allele,
                                "allele_copy_id": current_allele_copy_id
                            })

                    # Populate VariantReported (if location data is present for a row)
                    if pd.notna(row.get('Location')) and row['Location'] != '':
                        variant_type_parts = row.get('VariantType', '').split('>')
                        ref_allele = variant_type_parts[0].strip() if len(variant_type_parts) > 0 else None
                        alt_allele = variant_type_parts[1].strip() if len(variant_type_parts) > 1 else None

                        # Combine tool-specific flags into a single string for tool_specific_flags
                        tool_flags = []
                        if pd.notna(row.get('VariantStatus')): tool_flags.append(row['VariantStatus'])
                        # NORMAL: variant is associated with the star-allele in the database and is found in the sample
                        # NOVEL: gene-disrupting (core) variant is NOT associated with the star-allele in the database, but is found in the sample (this indicates that Aldy found a novel major star-allele)
                        # EXTRA: neutral variant is NOT associated with the star-allele in the database, but is found in the sample (this indicates that Aldy found a novel minor star-allele)
                        # MISSING: neutral variant is associated with the star-allele in the database, but is NOT found in the sample (this also indicates that Aldy found a novel minor star-allele)
                        
                        if pd.notna(row.get('VariantFunctionalityRaw')): tool_flags.append(f"FUNC:{row['VariantFunctionalityRaw']}")
                        # DISRUPTING for gene-disrupting (core, functional) variants, and
                        # NEUTRAL for neutral (silent) variants
                        
                        if pd.notna(row.get('KarolinskaCode')): tool_flags.append(f"CODE:{row['KarolinskaCode']}")
                        tool_specific_flags_str = "|".join(tool_flags) if tool_flags else None
                        
                        variant_data: VariantReported = {
                            "rsid": row.get('dbSNP', None) if pd.notna(row.get('dbSNP')) else None,
                            "location": row.get('Location', None) if pd.notna(row.get('Location')) else None,
                            "ref_allele": ref_allele,
                            "alt_allele": alt_allele,
                            "quality_score": int(row['Coverage']) if pd.notna(row.get('Coverage')) and row['Coverage'].isdigit() else None,
                            "allele_assignment": row.get('Allele', None) if pd.notna(row.get('Allele')) else None,
                            "tool_specific_flags": tool_specific_flags_str,
                            # Other fields like genotype, zygosity, etc., are not directly available in ALDY raw output
                            # and would be None or inferred later if needed.
                        }
                        variants_reported.append(variant_data)
                
                # Sort aldy_alleles_parsed by allele_copy_id for consistent output
                aldy_alleles_parsed.sort(key=lambda x: x['allele_copy_id'] if isinstance(x['allele_copy_id'], int) else float('inf'))

                # Construct the RawToolOutput dictionary
                raw_output_data: RawToolOutput = {
                    "diplotype_string": aldy_diplotype_string,
                    "copy_number_raw": estimated_cn,
                    "comments_raw": first_row.get('SolutionDescription', None), # The extracted solution description
                    "variants_reported": variants_reported,
                    "aldy_solution_id": str(solution_id), # Ensure it's a string
                    "aldy_alleles_in_solution_raw_string": aldy_alleles_in_solution_raw_string,
                    "aldy_alleles_parsed": aldy_alleles_parsed,
                }
                
                # Construct the main StandardizedGeneCall dictionary
                gene_call: StandardizedGeneCall = {
                    "sample_id": sample_id,
                    "gene": gene,
                    "tool_name": self.tool_name,
                    "reference_genome": ref_genome,
                    "input_file": os.path.basename(filepath),
                    "raw_tool_output": raw_output_data
                }
                parsed_gene_calls.append(gene_call)

        except pd.errors.EmptyDataError:
            print(f"Warning: ALDY file {filepath} is empty or contains only comments.")
            return []
        except Exception as e:
            print(f"Error parsing ALDY file {filepath}: {e}")
            # print(f"Problematic DataFrame head:\n{df.head().to_string() if 'df' in locals() else 'DataFrame not created'}")
            return []

        return parsed_gene_calls

# Example Usage (for testing the parser directly) - unchanged
if __name__ == "__main__":
    dummy_aldy_output_content = """#Sample	Gene	SolutionID	Major	Minor	AlleleCopyIdentifier	Allele	Location	VariantType	Coverage	VariantFunctionalityRaw	dbSNP	KarolinskaCode	VariantStatus
#Solution 1: *1.001, *4, *4.021
NA10860	CYP2D6	1	CYP2D6*1/*4+*4.021	1.001;4;4.021	0	*1.001				
NA10860	CYP2D6	1	CYP2D6*1/*4+*4.021	1.001;4;4.021	1	*4	42522612	C>G	15	S486T	rs1135840	C1	NORMAL
NA10860	CYP2D6	1	CYP2D6*1/*4+*4.021	1.001;4;4.021	2	*4.021	42522612	A>T	18	DISRUPTING	rs1234567	C2	NOVEL
#Solution 2: *4, *4, *139.001
NA10860	CYP2D6	2	CYP2D6*4/*4+*139.001	4;139.001;4	0	*4	42522612	C>G	15	S486T	rs1135840	C3	NORMAL
NA10860	CYP2D6	2	CYP2D6*4/*4+*139.001	4;139.001;4	1	*4	42524946	C>T	32	splicing defect/169frameshift	rs3892097	C4	NORMAL
NA10860	CYP2D6	2	CYP2D6*4/*4+*139.001	4;139.001;4	2	*139.001	42525000	G>A	25	NEUTRAL	rs9876543	C5	EXTRA
"""
    dummy_file_path = "src/parsers/test_aldy_output.tsv"
    with open(dummy_file_path, "w") as f:
        f.write(dummy_aldy_output_content)
    print(f"Created dummy ALDY output file: {dummy_file_path}")

    parser = AldyParser()
    parsed_results = parser.parse(dummy_file_path)

    print("\n--- Parsed ALDY Results ---")
    if parsed_results:
        for i, entry in enumerate(parsed_results):
            print(f"\n--- Solution {i+1} ---")
            print(json.dumps(entry, indent=2))
    else:
        print("No data parsed or file is empty.")

    if os.path.exists(dummy_file_path):
        os.remove(dummy_file_path)
        print(f"\nCleaned up dummy ALDY output file: {dummy_file_path}")