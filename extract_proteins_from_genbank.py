#!/usr/bin/env python3
"""
extract_proteins_from_genbank.py

Extracts all protein sequences from GenBank files (.gb or .gbff)
and appends the assembly name and strain in the FASTA headers.

Usage:
    python extract_proteins_from_genbank.py input.gbff output.faa
    python extract_proteins_from_genbank.py *.gbff proteins.faa
"""

import sys
import os
from Bio import SeqIO

def get_assembly_name(record, file_path):
    """Try to determine an assembly name."""
    if record.annotations.get("assembly"):
        return record.annotations["assembly"]
    if record.name:
        return record.name
    if record.id:
        return record.id
    return os.path.splitext(os.path.basename(file_path))[0]

def get_strain_name(record):
    """Extract strain name from the source feature or annotations."""
    strain = None
    for feature in record.features:
        if feature.type == "source":
            quals = feature.qualifiers
            if "strain" in quals:
                strain = quals["strain"][0]
                break
            elif "isolate" in quals:
                strain = quals["isolate"][0]
                break
    if not strain:
        # Try organism annotation
        strain = record.annotations.get("organism", None)
    if not strain:
        strain = "unknown_strain"
    return strain.replace(" ", "_")  # safer for FASTA headers

def extract_proteins_from_genbank(gb_file):
    """Extract CDS protein sequences from a GenBank file."""
    proteins = []
    for record in SeqIO.parse(gb_file, "genbank"):
        assembly_name = get_assembly_name(record, gb_file)
        strain = get_strain_name(record)
        for feature in record.features:
            if feature.type == "CDS":
                qualifiers = feature.qualifiers
                if "translation" not in qualifiers:
                    continue
                protein_seq = qualifiers["translation"][0]
                locus_tag = qualifiers.get("locus_tag", ["unknown"])[0]
                product = qualifiers.get("product", ["hypothetical protein"])[0]
                header = f">{assembly_name}|{strain}|{locus_tag} {product}"
                proteins.append((header, protein_seq))
    return proteins

def main(genbank_files, output_file):
    with open(output_file, "w") as out:
        for gb_file in genbank_files:
            for header, seq in extract_proteins_from_genbank(gb_file):
                out.write(f"{header}\n{seq}\n")

    print(f"✅ Extracted protein sequences saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("Usage: python extract_proteins_from_genbank.py input.gbff [more.gbff ...] output.faa")

    *input_files, output_file = sys.argv[1:]
    main(input_files, output_file)
