Dog DNA Project

# LOAD SEQUENCES


def load_sequences(database_file: str, mystery_file: str) -> Tuple[List[SeqRecord], SeqRecord]:
    """
    Load DNA sequences from FASTA files.

    Args:
        database_file (str): File containing known dog breed sequences
        mystery_file (str): File containing unknown sequence

    Returns:
        Tuple[List[SeqRecord], SeqRecord]:
            - List of known breed sequences
            - Mystery sequence
    """
    database = list(SeqIO.parse(database_file, "fasta"))
    mystery = next(SeqIO.parse(mystery_file, "fasta"))
    return database, mystery
 
