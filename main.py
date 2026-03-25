Dog DNA Project
from Bio import SeqIO, AlignIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from typing import List, Tuple
import subprocess
import shutil

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

# STEP 1: PAIRWISE RANKING


def rank_breeds(database: List[SeqRecord], mystery: SeqRecord) -> List[Tuple[str, float, float]]:
    """
    Rank all breeds by similarity to the mystery sequence using global alignment.

    Uses PairwiseAligner with simple scoring:
    - match = 1
    - mismatch = 0
    - no gap penalties

    Args:
        database: List of known breed sequences
        mystery: Unknown sequence

    Returns:
        List of tuples: (breed_id, score, percent_identity)
    """

    # Configure aligner (equivalent to old globalxx behaviour)
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0

    results: List[Tuple[str, float, float]] = []

    for breed in database:
        # Alignment score = number of matching positions
        score: float = aligner.score(mystery.seq, breed.seq)

        # Convert to approximate % identity
        max_len = max(len(mystery.seq), len(breed.seq))
        identity = (score / max_len) * 100

        results.append((breed.id, score, identity))

    # Sort highest score first (best match)
    results.sort(key=lambda x: x[1], reverse=True)

    # Print top results
    print("\n=== Pairwise Ranking ===")
    for i, (bid, score, pct) in enumerate(results[:10], 1):
        print(f"{i}. {bid:25} Score={score:.1f}  Identity={pct:.2f}%")

    best = results[0]
    print(f"\nBest match (pairwise): {best[0]} ({best[2]:.2f}%)")

    return results
    
def write_combined(database: List[SeqRecord], mystery: SeqRecord, output: str = "combined.fa") -> str:
    """
    Combine database and mystery sequences into a single FASTA file.

    Required for multiple sequence alignment.

    Args:
        database: Known sequences
        mystery: Unknown sequence
        output: Output FASTA file name

    Returns:
        str: Path to combined FASTA file
    """
    SeqIO.write(database + [mystery], output, "fasta")
    return output
    
def run_mafft(input_file: str, output_file: str = "aligned.fa"):
    """
    Run MAFFT multiple sequence alignment.
    """
    if not shutil.which("mafft"):
        raise RuntimeError("MAFFT is not installed or not in PATH")

    try:
        with open(output_file, "w") as output_handle:
            subprocess.run(
                ["mafft", "--auto", input_file],
                stdout=output_handle,
                stderr=subprocess.DEVNULL,
                check=True
            )

        print(f"\nAlignment complete: {output_file}")
        return output_file

    except subprocess.CalledProcessError as e:
        print("MAFFT failed:", e)
        return None
def build_tree(aligned_file: str):
    alignment = AlignIO.read(aligned_file, "fasta")

    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    print("\nPhylogenetic Tree:")
    Phylo.draw_ascii(tree)

    return tree
