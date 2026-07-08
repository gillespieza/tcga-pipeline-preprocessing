#!/usr/bin/env python3
"""
tcga_fetch.py — Entry point for the TCGA Data Preprocessing Pipeline.

Usage:
    python tcga_fetch.py                          # interactive mode
    python tcga_fetch.py --study kirc_tcga_pan_can_atlas_2018 --data clinical rnaseq
    python tcga_fetch.py --list-studies
    python tcga_fetch.py --help
"""

from tcga.cli import main

if __name__ == "__main__":
    main()
