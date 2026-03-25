#Dog DNA Project
from Bio import SeqIO, AlignIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from typing import List, Tuple
import subprocess
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
def rank_breeds(database_sequences: List[SeqRecord], mystery_sequence: SeqRecord) -> List[Tuple[str, float, float]]:
    """
    Compare mystery DNA to each database sequence and rank similarity.

    Returns:
        List of (sequence_id, score, percent_identity)
    """

    aligner = PairwiseAligner()
    aligner.mode = "global"

    # Scoring scheme (biologically realistic)
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    results: List[Tuple[str, float, float]] = []

    for breed in database_sequences:
        # Compute alignment score
        score: float = aligner.score(mystery_sequence.seq, breed.seq)

        # Convert to percent identity
        max_length: int = max(len(mystery_sequence.seq), len(breed.seq))
        identity: float = (score / max_length) * 100

        results.append((breed.id, score, identity))

    # Sort best match first
    results.sort(key=lambda x: x[1], reverse=True)

    # Print top 10 matches
    print("\n=== Pairwise Ranking ===")
    for i, (breed_id, score, identity) in enumerate(results[:10], 1):
        print(f"{i}. {breed_id} | Score={score:.2f} | Identity={identity:.2f}%")

    # Report best match
    best_match = results[0]
    print(f"\nBest match: {best_match[0]} ({best_match[2]:.2f}%)")

    return results

# STEP 2: MULTIPLE SEQUENCE ALIGNMENT(MAFFT)
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
    Run MAFFT to perform multiple sequence alignment.

    MAFFT must be installed and available in the system PATH.

    Args:
        input_file (str): FASTA file containing sequences to align.
        output_file (str): Output file for aligned sequences.

    Returns:
        str or None:
            Path to aligned file if successful, otherwise None.

    Raises:
        RuntimeError: If MAFFT is not installed or not found in PATH.
    """
    # Check if MAFFT is available on the system
    if not shutil.which("mafft"):
        raise RuntimeError("MAFFT is not installed or not in PATH")

    try:
        # Run MAFFT and write output to file
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
        
# STEP 3: PHYLOGENETIC TREE CONSTRUCTION
def build_tree(aligned_file: str, output_tree="tree.nwk"):
    """
    Construct a phylogenetic tree from aligned sequences.

    Uses:
    - Identity-based distance calculation
    - Neighbor-Joining (NJ) tree construction

    Args:
        aligned_file (str): Path to aligned FASTA file.
        output_tree (str): Output file for Newick tree.

    Returns:
        Tree: Constructed phylogenetic tree object.

    Side Effects:
        - Prints ASCII representation of the tree
        - Saves tree in Newick format
    """
    # Read aligned sequences
    alignment = AlignIO.read(aligned_file, "fasta")
    # Compute pairwise distances
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    # Build tree using Neighbor-Joining algorithm
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)
    # Display tree in console
    print("\nPhylogenetic Tree:")
    Phylo.draw_ascii(tree)

    # Save tree in Newick format
    Phylo.write(tree, output_tree, "newick")
    print(f"\nTree saved to: {output_tree}")

    return tree

def main():
    """
    Main pipeline execution.
    """
    database_file = "dog_breeds.fa"
    mystery_file = "mystery.fa"
    # Load input sequences
    database, mystery = load_sequences(database_file, mystery_file)
    # Rank similarity of mystery sequence
    rank_breeds(database, mystery)
    # Prepare combined FASTA for alignment
    combined_file = write_combined(database, mystery)
    # Run multiple sequence alignment
    aligned_file = run_mafft(combined_file)
    #Build phylogenetic tree
    if aligned_file is not None:
        build_tree(aligned_file, "tree.nwk")
# Entry point of the script
if __name__ == "__main__":
    main()



