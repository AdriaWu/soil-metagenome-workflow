#!/usr/bin/env python3
"""
Annotate Contigs with Taxonomy

Description: Assign taxonomy to contigs based on MAG mapping results
Input: 
    - Mapped contigs CSV (contig_id, MAG_ID)
    - SMAG taxonomy CSV (MAG_ID, taxonomic ranks)
    - Optional: Plant genome MAG lists for host assignment
Output: Annotated contig table with full taxonomy
"""

import pandas as pd
import sys
from pathlib import Path

# ---- I/O Paths ----
# Input files
MAPPED_CONTIGS_CSV = "<INPUT_MAPPED_CONTIGS_CSV>"
SMAG_TAXONOMY_CSV = "<INPUT_SMAG_TAXONOMY_CSV>"
OUTPUT_DIR = "<OUTPUT_ANNOTATION_DIR>"

# Plant genome MAG lists (optional)
PLANT_MAG_FILES = {
    "barley": "<BARLEY_CONTIGS_TXT>",
    "wheat": "<WHEAT_CONTIGS_TXT>",
    "canola": "<CANOLA_CONTIGS_TXT>",
    "carinata": "<CARINATA_CONTIGS_TXT>",
    "chickpea": "<CHICKPEA_CONTIGS_TXT>",
    "medicago": "<MEDICAGO_CONTIGS_TXT>",
    "eucalyptus": "<EUCALYPTUS_CONTIGS_TXT>"
}

# Plant taxonomies
PLANT_TAXONOMIES = {
    "barley": ["d__Eukaryote", "p__Tracheophyta", "c__Liliopsida", "o__Poales", 
               "f__Poaceae", "g__Hordeum", "s__H.vulgare"],
    "wheat": ["d__Eukaryote", "p__Tracheophyta", "c__Liliopsida", "o__Poales",
              "f__Poaceae", "g__Triticum", "s__T.aestivum"],
    "canola": ["d__Eukaryote", "p__Tracheophyta", "c__Magnoliopsida", "o__Brassicales",
               "f__Brassicaceae", "g__Brassica", "s__B.napus"],
    "carinata": ["d__Eukaryote", "p__Tracheophyta", "c__Magnoliopsida", "o__Brassicales",
                 "f__Brassicaceae", "g__Brassica", "s__B.carinata"],
    "chickpea": ["d__Eukaryote", "p__Tracheophyta", "c__Magnoliopsida", "o__Fabales",
                 "f__Fabaceae", "g__Cicer", "s__C.arietinum"],
    "medicago": ["d__Eukaryote", "p__Tracheophyta", "c__Magnoliopsida", "o__Fabales",
                 "f__Fabaceae", "g__Medicago", "s__M.truncatula"],
    "eucalyptus": ["d__Eukaryote", "p__Tracheophyta", "c__Magnoliopsida", "o__Myrtales",
                   "f__Myrtaceae", "g__Eucalyptus", "s__E.regnans"]
}

def clean_taxonomy(tax_string):
    """Remove GTDB prefixes from taxonomy strings"""
    if pd.isna(tax_string) or tax_string == "" or tax_string in ["g__", "s__"]:
        return "Unclassified"
    return tax_string.replace("d__", "").replace("p__", "").replace("c__", "") \
                     .replace("o__", "").replace("f__", "").replace("g__", "") \
                     .replace("s__", "")

def read_mag_ids(file_path):
    """Read MAG IDs from text file"""
    if not Path(file_path).exists():
        print(f"Warning: File not found: {file_path}")
        return []
    
    with open(file_path, 'r') as f:
        mag_ids = [line.strip() for line in f if line.strip()]
    
    return mag_ids

def assign_plant_taxonomy(contig_df, plant_mag_files, plant_taxonomies):
    """Assign taxonomy to contigs based on plant genome matches"""
    result_df = contig_df.copy()
    
    # Initialize taxonomy columns
    tax_ranks = ["classification", "phylum", "class", "order", "family", "genus", "species"]
    for rank in tax_ranks:
        if rank not in result_df.columns:
            result_df[rank] = None
    
    # Track statistics
    match_stats = {}
    
    # Process each plant species
    for plant_name, mag_file in plant_mag_files.items():
        print(f"Processing {plant_name}...")
        
        mag_ids = read_mag_ids(mag_file)
        
        if not mag_ids:
            print(f"  No MAG IDs found for {plant_name}")
            continue
        
        print(f"  Loaded {len(mag_ids)} MAG IDs")
        
        # Find matching contigs
        matches = result_df['MAG_ID'].isin(mag_ids)
        n_matches = matches.sum()
        
        if n_matches > 0:
            # Assign taxonomy
            tax = plant_taxonomies[plant_name]
            result_df.loc[matches, tax_ranks] = tax
            
            print(f"  Assigned taxonomy to {n_matches} contigs")
            match_stats[plant_name] = n_matches
        else:
            print(f"  No matching contigs found")
            match_stats[plant_name] = 0
    
    # Print summary
    print("\n=== Summary ===")
    print(f"Total contigs: {len(result_df)}")
    print(f"Contigs with assigned taxonomy: {result_df['classification'].notna().sum()}")
    print(f"Contigs without taxonomy: {result_df['classification'].isna().sum()}")
    
    print("\nMatches by plant species:")
    for plant, count in match_stats.items():
        print(f"  {plant}: {count}")
    
    return result_df

def main():
    print("=" * 60)
    print("Contig Annotation Pipeline")
    print("=" * 60)
    
    # Create output directory
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    
    # Load mapped contigs
    print("\nLoading mapped contigs...")
    mapped_contigs = pd.read_csv(MAPPED_CONTIGS_CSV)
    print(f"Loaded {len(mapped_contigs)} mapped contigs")
    
    # Load SMAG taxonomy
    print("\nLoading SMAG taxonomy...")
    smag_taxonomy = pd.read_csv(SMAG_TAXONOMY_CSV)
    
    # Clean taxonomy columns
    tax_cols = ["phylum", "class", "order", "family", "genus", "species"]
    for col in tax_cols:
        if col in smag_taxonomy.columns:
            smag_taxonomy[col] = smag_taxonomy[col].apply(clean_taxonomy)
    
    # Handle unclassified
    for col in ["order", "family"]:
        if col in smag_taxonomy.columns:
            smag_taxonomy[col] = smag_taxonomy[col].replace("", "Unclassified")
            smag_taxonomy[col] = smag_taxonomy[col].replace("unclassified", "Unclassified")
    
    print(f"Loaded taxonomy for {len(smag_taxonomy)} SMAGs")
    
    # Merge contig data with SMAG taxonomy
    print("\nMerging contig data with SMAG taxonomy...")
    annotated_contigs = mapped_contigs.merge(
        smag_taxonomy,
        left_on='MAG_ID',
        right_on='user_genome',
        how='left'
    )
    
    print(f"Merged {len(annotated_contigs)} contigs with taxonomy")
    
    # Assign plant taxonomy (optional)
    print("\nAssigning plant taxonomy...")
    annotated_contigs = assign_plant_taxonomy(
        annotated_contigs,
        PLANT_MAG_FILES,
        PLANT_TAXONOMIES
    )
    
    # Save annotated contigs
    output_file = Path(OUTPUT_DIR) / "annotated_contigs.csv"
    annotated_contigs.to_csv(output_file, index=False)
    print(f"\nSaved annotated contigs to: {output_file}")
    
    # Generate taxonomy summary
    print("\nGenerating taxonomy summary...")
    
    # Count by phylum
    if 'phylum' in annotated_contigs.columns:
        phylum_counts = annotated_contigs['phylum'].value_counts()
        phylum_file = Path(OUTPUT_DIR) / "taxonomy_summary_phylum.csv"
        phylum_counts.to_csv(phylum_file, header=['count'])
        print(f"Saved phylum summary to: {phylum_file}")
    
    # Count by family
    if 'family' in annotated_contigs.columns:
        family_counts = annotated_contigs['family'].value_counts()
        family_file = Path(OUTPUT_DIR) / "taxonomy_summary_family.csv"
        family_counts.to_csv(family_file, header=['count'])
        print(f"Saved family summary to: {family_file}")
    
    print("\nAnnotation complete!")

if __name__ == "__main__":
    main()