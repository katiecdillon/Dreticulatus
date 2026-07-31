#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from collections import defaultdict
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import multiprocessing as mp
import gzip

matplotlib.use("Agg")

te_classes = ["LINE", "SINE", "LTR", "DNA", "RC", "Satellite"]

# --------------------------------------------------
# Mapping: assembly → species
# --------------------------------------------------
def load_species_mapping(mapping_file):
    mapping = {}
    with open(mapping_file) as f:
        next(f)  # skip header
        for line in f:
            assembly, genus, species = line.strip().split("\t")
            mapping[assembly] = f"{genus} {species}"
    return mapping


# --------------------------------------------------
# Extract TE occupancy (bp)
# --------------------------------------------------
def extract_occupancy(file_path, assembly_to_species):
    class_bp = defaultdict(int)

    assembly = os.path.basename(file_path).replace(".out.gz", "")
    species = assembly_to_species.get(assembly, "Unknown")

    with gzip.open(file_path, "rt") as f:
        for _ in range(3):  # skip header lines
            next(f)
        for line in f:
            cols = line.split()
            if len(cols) > 10:
                clsfam = cols[10]
                cls = clsfam.split("/")[0]

                start = int(cols[5])
                end = int(cols[6])
                length = end - start + 1

                if cls in te_classes:
                    class_bp[cls] += length

    return assembly, species, class_bp


# --------------------------------------------------
# Build matrix
# --------------------------------------------------
def build_matrix(files, assembly_to_species, num_proc):
    with mp.Pool(num_proc) as pool:
        results = pool.starmap(
            extract_occupancy,
            [(fp, assembly_to_species) for fp in files]
        )

    matrix = pd.DataFrame()
    species_labels = {}

    for assembly, species, class_bp in results:
        species_labels[assembly] = species
        matrix = pd.concat(
            [matrix, pd.DataFrame(class_bp, index=[assembly])]
        ).fillna(0)

    return matrix, species_labels


# --------------------------------------------------
# PCA + scatter + biplot
# --------------------------------------------------
def run_pca(matrix, species_labels, outprefix, ylabel):
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(matrix)

    df = pd.DataFrame(pcs, index=matrix.index, columns=["PC1", "PC2"])
    df["Species"] = df.index.map(species_labels)

    species = sorted(df["Species"].unique())
    palette = dict(zip(species, sns.color_palette("tab10", len(species))))
    markers = dict(zip(species, ["o", "s", "D", "^", "v"]))

    df.to_csv(f"{outprefix}_scores.tsv", sep="\t")

    # ---------- PCA scatter ----------
    plt.figure(figsize=(9, 7))
    sns.scatterplot(
        data=df, x="PC1", y="PC2",
        hue="Species", style="Species",
        palette=palette, markers=markers,
        s=140, edgecolor="black"
    )
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title(ylabel)
    plt.legend(bbox_to_anchor=(1.02, 0.5), loc="center left")
    plt.tight_layout()
    plt.savefig(f"{outprefix}_PCA.png")
    plt.close()

    # ---------- BIPLOT ----------
    plt.figure(figsize=(9, 7))
    sns.scatterplot(
        data=df, x="PC1", y="PC2",
        hue="Species", style="Species",
        palette=palette, markers=markers,
        s=140, edgecolor="black"
    )

    # Draw arrows for TE classes
    for i, cls in enumerate(matrix.columns):
        plt.arrow(
            0, 0,
            pca.components_[0, i] * df["PC1"].max(),
            pca.components_[1, i] * df["PC2"].max(),
            color="black", alpha=0.6,
            head_width=0.02, length_includes_head=True
        )
        plt.text(
            pca.components_[0, i] * df["PC1"].max() * 1.1,
            pca.components_[1, i] * df["PC2"].max() * 1.1,
            cls, ha="center", va="center"
        )

    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title(ylabel + " (Biplot)")
    plt.legend(bbox_to_anchor=(1.02, 0.5), loc="center left")
    plt.tight_layout()
    plt.savefig(f"{outprefix}_BIPLOT.png")
    plt.close()


# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--mapping_file", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--num_proc", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    assembly_to_species = load_species_mapping(args.mapping_file)

    files = [
        os.path.join(args.out_dir, f)
        for f in os.listdir(args.out_dir)
        if f.endswith(".out.gz")
    ]

    matrix, species_labels = build_matrix(
        files, assembly_to_species, args.num_proc
    )

    matrix.to_csv(
        f"{args.output_dir}/class_occupancy_matrix.tsv", sep="\t"
    )

    run_pca(
        matrix,
        species_labels,
        f"{args.output_dir}/class_occupancy",
        "TE Occupancy (bp)"
    )


if __name__ == "__main__":
    main()
