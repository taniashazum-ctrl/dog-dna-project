# Dog DNA Project

## Dog Breed Identification & Phylogenetic Analysis

## Overview

This project identifies the closest dog breed to an unknown DNA sequence and visualizes evolutionary relationships using phylogenetic analysis. It combines pairwise alignment for similarity ranking with multiple sequence alignment and phylogenetic tree construction to provide both quantitative and evolutionary insights.

## Key Features

* Compares an unknown (“mystery”) DNA sequence against a database of known dog breeds
* Ranks breeds by similarity using global pairwise alignment
* Performs multiple sequence alignment using MAFFT
* Builds a phylogenetic tree using the Neighbor-Joining method
* Outputs similarity scores and a tree in Newick format

## Project Structure

* `dog_breeds.fa` – database of known dog breed DNA sequences
* `mystery.fa` – unknown sequence to classify
* `main.py` – main pipeline script
* `tree.nwk` – generated phylogenetic tree

## Methodology

### 1. Pairwise Alignment

Each database sequence is aligned with the mystery sequence using global alignment. A biologically informed scoring system (match, mismatch, and gap penalties) is used to compute similarity scores and percentage identity.

### 2. Multiple Sequence Alignment

All sequences (database + mystery) are combined and aligned using **MAFFT**.

**Note:**
MAFFT was used instead of MUSCLE because MUSCLE had significantly longer runtime, making it inefficient for this dataset.

### 3. Phylogenetic Tree Construction

* A distance matrix is computed from the aligned sequences
* A Neighbor-Joining (NJ) algorithm is used to construct the tree
* The tree is displayed in ASCII format and saved in Newick format

## Requirements

* Python 3.x
* Biopython
* MAFFT (installed and available in PATH)

## Output

* Ranked list of closest breeds with scores and identity percentages
* Best matching breed
* Aligned sequences file (`aligned.fa`)
* Phylogenetic tree (`tree.nwk`)

## Summary

This project implements a bioinformatics pipeline for DNA-based breed identification. It combines sequence similarity analysis with phylogenetic methods to identify the closest match and reveal evolutionary relationships.
