import abc # Import the Abstract Base Classes module
from typing import List # For type hinting lists

# Import standardized gene call format
from src.standard_formats import StandardizedGeneCall

class BaseParser(abc.ABC):
    """
    Abstract Base Class for all pharmacogenomics genotyping tool parsers.
    All concrete parser implementations must inherit from this class and
    implement the 'parse' method.
    """
    def __init__(self):
        # tool_name should be set by concrete parser implementations in their __init__
        self.tool_name: str = "GenericTool" 
        
    @abc.abstractmethod
    def parse(self, filepath: str) -> List[StandardizedGeneCall]:
        """
        Abstract method to parse a genotyping tool's output file.
        Concrete parser implementations must provide their specific logic here.

        Args:
            filepath (str): The path to the input file to be parsed.

        Returns:
            List[StandardizedGeneCall]: A list of StandardizedGeneCall dictionaries,
                                       each representing a gene call from the tool.
                                       Returns an empty list if parsing fails or no data.
        """
        pass # Abstract methods typically have 'pass' as their body