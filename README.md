# SIMD Sequence Alignment

An optimized C++ implementation of Global, Local, and Semi-Global sequence alignment algorithms using SIMD (Single Instruction, Multiple Data) intrinsics. This project explores performance acceleration for bioinformatics workloads, comparing raw C++ SIMD implementations against the SEQ programming language.

## Supported Algorithms

This repository implements the following alignment strategies with SIMD acceleration:

- **Global Alignment (Needleman-Wunsch):** Optimal alignment across the entire length of two sequences.
- **Local Alignment (Smith-Waterman):** Finds the most similar regions between two sequences.
- **Semi-Global Alignment:** Useful for finding a sequence within another (overlap alignment).

## Key Improvements

- **SIMD Parallelization:** Utilizes instruction-level parallelism to process multiple data points in a single clock cycle.
- **Performance Comparison:** Includes data comparing execution times with the **SEQ** programming language.
- **Synthetic Data Generation:** A Python-based generator to create large-scale FASTA/TXT sequences for stress testing.

## Project Structure

```
├── README.md
├── Data/
│   ├── sequences.fa       # Sample input sequences in FASTA format
│   └── sequences.txt      # Sample input sequences in plain text
├── docs/
│   ├── c++_results.txt    # Benchmark outputs for the C++ implementation
│   ├── cpu.txt            # CPU specifications used for benchmarking
│   └── seq_results.txt    # Benchmark outputs for the SEQ language comparison
└── src/
    ├── generator.py       # Python script to generate random DNA sequences
    ├── main.cpp           # Core SIMD implementation of alignment algorithms
    └── main.seq           # SEQ implementation for comparison
```

## Prerequisites

- C++ compiler with SIMD support (e.g., g++ with AVX and SSE flags)
- Python 3 for data generation
- SEQ programming language for comparison benchmarks

## Installation and Setup

1. **Clone the repository:**

   ```bash
   git clone https://github.com/yourusername/SIMD-coding.git
   cd SIMD-coding
   ```

2. **Generate test data (optional, pre-generated data is included):**
   ```bash
   cd src
   python3 generator.py
   mv sequences.txt ../Data/
   mv sequences.fa ../Data/
   ```

## Building and Running

1. **Compile the C++ code:**

   ```bash
   cd src
   g++ -O3 -mavx -msse4.2 -std=c++11 main.cpp -o alignment
   ```

2. **Run the C++ benchmarks:**

   ```bash
   ./alignment
   ```

   This will generate `c++_results.txt` in the `docs/` directory with timing results.

3. **Run the SEQ benchmarks (requires SEQ installed):**
   ```bash
   seq main.seq
   ```
   Results will be printed to stdout and can be redirected to `docs/seq_results.txt`.

## Results

Benchmark results are stored in the `docs/` folder:

- `c++_results.txt`: Timing for scalar, SSE, and AVX implementations with and without CIGAR string generation.
- `seq_results.txt`: Timing for SEQ implementations.
- `cpu.txt`: CPU details for the benchmarking machine.

Example C++ results (times in milliseconds for 10 runs averaged):

- Global AVX: ~5729 ms (no CIGAR), ~6032 ms (with CIGAR)
- Local AVX: ~7685 ms (no CIGAR), ~7714 ms (with CIGAR)
- Semi-Global AVX: ~6039 ms (no CIGAR), ~6079 ms (with CIGAR)

SEQ results (times in seconds):

- Scalar without CIGAR: 0.304817 s
- Scalar with CIGAR: 0.36378 s

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## License

This project is licensed under the MIT License.
