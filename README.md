# TFHE Leveled Homomorphic Evaluation

This repository contains the implementation and experimental code for "Refined TFHE Leveled Homomorphic Evaluation and Its Application" submitted to ACM CCS.

## Abstract

This artifact provides two complementary implementations demonstrating our refined TFHE leveled homomorphic evaluation techniques. The repository showcases novel approaches to improving the efficiency of fully homomorphic encryption through optimized computational flows and advanced mathematical transforms.

## Paper Contributions

This artifact demonstrates the following key contributions:

- **Enhanced NTT-based LHE scheme**: Significantly improves the efficiency of Number Theoretic Transform-based Leveled Homomorphic Evaluation by optimizing the computational flow and reducing overhead.

- **Novel FFT-based LHE scheme**: Introduces a new Leveled Homomorphic Evaluation construction based on Fast Fourier Transform, achieving superior performance compared to existing approaches.

- **Transciphering application**: Provides comprehensive evaluation of transciphering tasks using our LHE schemes, demonstrating real-world performance improvements.

- **Integer-input LHE extension**: Extends the FFT-based LHE scheme from binary input operations to support integer inputs, broadening its practical applicability.

## Repository Structure

The implementation is organized into two main directories, each focusing on specific optimizations and mathematical approaches:

### 📁 `ntt_based_lhe/`

Contains the Number Theoretic Transform (NTT) based implementation of our leveled homomorphic evaluation scheme. This implementation corresponds to **Contribution 1** and focuses on optimizing traditional NTT-based approaches through refined computational flows.

[🔗 See detailed build and run instructions](./ntt_based_lhe/README.md)

### 📁 `fft_based_lhe/`

Contains the Fast Fourier Transform (FFT) based implementation of our leveled homomorphic evaluation scheme. This implementation encompasses **Contributions 2, 3, and 4**, representing our most significant algorithmic innovations.

[🔗 See detailed build and run instructions](./fft_based_lhe/README.md)

## Getting Started

Each implementation directory contains comprehensive documentation, build instructions, and usage examples. We recommend starting with the specific README files in each subdirectory to understand the requirements and setup procedures for your intended use case.
