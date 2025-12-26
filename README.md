# SIMD Sequence Alignment

An optimized C++ implementation of Global, Local, and Semi-Global sequence alignment algorithms using SIMD (Single Instruction, Multiple Data) intrinsics. This project explores performance acceleration for bioinformatics workloads, comparing raw C++ SIMD implementations against the SEQ programming language.

## 🧬 Supported Algorithms

This repository implements the following alignment strategies with SIMD acceleration:

- **Global Alignment (Needleman-Wunsch):** Optimal alignment across the entire length of two sequences.
- **Local Alignment (Smith-Waterman):** Finds the most similar regions between two sequences.
- **Semi-Global Alignment:** Useful for finding a sequence within another (overlap alignment).

## 🚀 Key Improvements

- **SIMD Parallelization:** Utilizes instruction-level parallelism to process multiple data points in a single clock cycle.
- **Performance Comparison:** Includes data comparing execution times with the **SEQ** programming language.
- **Synthetic Data Generation:** A Python-based generator to create large-scale FASTA/TXT sequences for stress testing.

## 📂 Project Structure

```text
├── main.cpp           # Core SIMD implementation of alignment algorithms
├── generator.py       # Python script to generate random DNA/Protein sequences
├── Report.pdf         # Detailed analysis and performance evaluation
├── sequences.fa       # Sample input sequences in FASTA format
├── c++_results.txt    # Benchmark outputs for the C++ implementation
└── seq_results.txt    # Benchmark outputs for the SEQ language comparison
```
