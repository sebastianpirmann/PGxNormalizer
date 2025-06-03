# pgx_normalizer/src/standard_formats.py

from typing import List, Dict, Any, Optional, Union, TypedDict

# --- Helper TypedDicts for nested structures ---

class VariantReported(TypedDict, total=False):
    """
    Represents an individual variant (SNP/indel) reported by a genotyping tool.
    All fields are optional by default (total=False). Parsers should populate
    as much information as is available and relevant from the raw tool output.
    """
    rsid: Optional[str]                   # The dbSNP identifier (e.g., "rs1135840")
    location: Optional[str]               # Genomic coordinate (e.g., "chr22:42522612" or just "42522612")
    ref_allele: Optional[str]             # Reference allele (e.g., "C", "G")
    alt_allele: Optional[str]             # Alternative allele (e.g., "G", "A")
    genotype: Optional[str]               # The specific genotype at this variant locus (e.g., "A/G", "G/G")
    zygosity: Optional[str]               # Zygosity of the variant (e.g., "homozygous", "heterozygous", "hemizygous")
    quality_score: Optional[Union[int, float]] # A numerical score indicating quality or confidence (e.g., coverage depth, Phred score)
    allele_assignment: Optional[str]      # Indicates which specific raw allele/haplotype this variant defines/belongs to
                                          # (e.g., "haplotype1", "*4", "CYP2D6*4")
    tool_specific_flags: Optional[str]    # Any specific flags, notes, or raw functional annotations for this variant
                                          # directly from the tool's output (e.g., "NORMAL|FUNC:S486T", "DISRUPTING")


class StructuralVariantRaw(TypedDict, total=False):
    """
    Represents a raw structural variant (e.g., gene deletion, duplication, hybrid gene)
    as directly reported by a genotyping tool.
    """
    type: str           # The type of structural variant (e.g., "deletion", "duplication", "hybrid")
    description: str    # A descriptive string for the structural variant from the tool (e.g., "CYP2D6 gene deletion")
    location: str       # Genomic coordinates of the structural variant, if available (e.g., "chr22:start-end")
    tool_specific_id: str # Any internal ID or specific nomenclature the tool uses for this SV


class RawAlleleComponent(TypedDict):
    """
    A generic structure to capture individual allele components or raw allele names
    along with any associated copy identifiers, as might be parsed from tool-specific fields.
    This is particularly useful for tools like ALDY that delineate individual "allele copies".
    """
    raw_allele_name: str                 # The allele string as identified by the tool (e.g., "*1", "*4", "CYP2D6*1")
    allele_copy_id: Optional[Union[int, str]] # An identifier for the specific copy (e.g., 0, 1, 2 for distinct copies, or a string ID)


class RawToolOutput(TypedDict, total=False):
    """
    A flexible dictionary designed to store ALL raw, tool-specific output details
    that are relevant to the gene call but are not yet mapped to a universally
    standardized field (or might never be, if they are tool-unique).

    By setting `total=False`, all fields within this TypedDict are optional by default.
    However, the `diplotype_string` is considered REQUIRED by the parser to be present.
    """
    # --- Core Interpretation from the Tool (Highly recommended for parsers to provide) ---
    diplotype_string: str # REQUIRED: The primary diplotype call string as reported by the tool
                          # (e.g., "CYP2D6*1/*4", "*1x2/*4"). This is the most fundamental output.

    # --- Haplotype-level Raw Data (if tool distinguishes them) ---
    haplotype1_raw: Optional[str] # The raw string for the first inferred haplotype
    haplotype2_raw: Optional[str] # The raw string for the second inferred haplotype

    # --- Other Summary-level Raw Data from the Tool ---
    copy_number_raw: Optional[Union[int, float]] # Any reported gene/segmental copy number from the tool
    functional_status_raw: Optional[str]        # Tool's direct functional prediction (e.g., "Normal Function", "Decreased Function")
    phenotype_prediction_raw: Optional[str]     # Tool's direct phenotype prediction (e.g., "Normal Metabolizer", "UM")
    confidence_score_raw: Optional[Union[int, float]] # Any overall quality/confidence score for the call from the tool
    comments_raw: Optional[str]                 # General comments, notes, or supplementary text from the tool's output
                                                # (e.g., the "#Solution X:" description from ALDY)

    # --- Detailed Variant and Structural Information ---
    variants_reported: List[VariantReported]     # A list of individual variant calls (SNPs/indels) contributing to the diplotype
    structural_variants_raw: List[StructuralVariantRaw] # A list of structural variants reported by the tool

    # --- Tool-Specific Raw Data (Add as needed for each new parser) ---
    # These fields are crucial for capturing unique details from specific tools.
    # Example fields for ALDY parser (as discussed previously):
    aldy_solution_id: Optional[str]                 # ALDY's internal solution identifier
    aldy_alleles_in_solution_raw_string: Optional[str] # The raw string from ALDY's 'Minor' column (e.g., "1.001;4;4.021")
    aldy_alleles_parsed: List[RawAlleleComponent]   # A structured representation of individual alleles and copy IDs from ALDY

    # Example for other hypothetical tools (e.g., Stargazer):
    # stargazer_prediction_method: Optional[str]
    # stargazer_raw_gene_coverage_data: Optional[Dict[str, Any]]
    # ... any other fields unique to Stargazer's raw output ...


class StandardizedGeneCall(TypedDict):
    """
    The main standardized format for a single pharmacogene call from a single tool for a single sample.
    This TypedDict serves as the primary "contract" between all parsers and the PGx Normalizer.
    All fields in this top-level TypedDict are REQUIRED unless specified otherwise.
    """
    # --- Core Identifiers (REQUIRED for every StandardizedGeneCall) ---
    sample_id: str             # Unique identifier for the patient or sample (e.g., "NA10860")
    gene: str                  # The pharmacogene symbol (e.g., "CYP2D6", "CYP2C19", "TPMT")
    tool_name: str             # The name of the genotyping tool that generated this call (e.g., "ALDY", "Stargazer", "PharmCAT")
    reference_genome: str      # The genomic reference assembly used by the tool for this call (e.g., "GRCh37", "GRCh38")

    # --- Traceability and Raw Data (Important for Normalizer) ---
    input_file: Optional[str]  # The base name or path of the original input file processed by the parser (e.g., "aldy_output.tsv")

    # This crucial nested dictionary holds all the raw, tool-specific details.
    # It must always be present, although its contents are flexible.
    raw_tool_output: RawToolOutput