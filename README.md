# Merge peak iteratively

## Introduction

This script is aimed to merge multiple peak sets iteratively. This method is first proposed in this paper: [The chromatin accessibility landscape of primary human cancers](https://www.science.org/doi/10.1126/science.aav1898).

Most codes are borrowed from [SnapATAC2](https://scverse.org/SnapATAC2).

## Usage

Compile
```
cargo build --release
```
Usage
```
merge_peaks --chrom-sizes <YOUR_CHROM_SIZE> --output <OUT_BED> <INPUT1> <INPUT2> <INPUT3>  --half-width 250
```

