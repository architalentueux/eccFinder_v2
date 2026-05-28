import pandas as pd
from pybedtools import BedTool
import re
import os

def parse_gtf_attributes(attr_string, key):
    """Extraire une valeur (gene_name, gene_id, etc.) depuis un champ GTF."""
    match = re.search(fr'{key} "([^"]+)"', attr_string)
    return match.group(1) if match else None


def intersect_ecc(
    ecc_file,
    gene_annot,
    te_annot,
    alu_annot,
    output_file,
    overlap_fraction=0.2
):
    # Charger ecc (supposé TSV ou CSV)
    ecc = pd.read_csv(ecc_file, sep=None, engine="python")

    # Convertir en BedTool
    ecc_bed = BedTool.from_dataframe(ecc.iloc[:, :7])

    # --- GENE ---
    gene_bed = BedTool(gene_annot)
    gene_inter = ecc_bed.intersect(gene_bed, wa=True, wb=True, f=overlap_fraction)

    gene_dict = {}
    geneid_dict = {}

    for entry in gene_inter:
        fields = entry.fields
        key = tuple(fields[:7])
        info = " ".join(fields[7:])

        gene_name = parse_gtf_attributes(info, "gene_name")
        gene_id = parse_gtf_attributes(info, "gene_id")

        if gene_name and "Unclassified" not in gene_name:
            gene_dict[key] = gene_name
            geneid_dict[key] = gene_id

    # --- TE ---
    te_bed = BedTool(te_annot)
    te_inter = ecc_bed.intersect(te_bed, wa=True, wb=True, f=overlap_fraction)

    te_dict = {}
    for entry in te_inter:
        fields = entry.fields
        key = tuple(fields[:7])
        info = " ".join(fields[7:])
        te_name = parse_gtf_attributes(info, "gene_id")
        te_dict[key] = te_name

    # --- ALU ---
    alu_bed = BedTool(alu_annot)
    alu_inter = ecc_bed.intersect(alu_bed, wa=True, wb=True, f=overlap_fraction)

    alu_dict = {}
    for entry in alu_inter:
        fields = entry.fields
        key = tuple(fields[:7])
        info = " ".join(fields[7:])
        alu_name = parse_gtf_attributes(info, "gene_id")
        alu_dict[key] = alu_name

    # --- CONSTRUCTION SORTIE ---
    output_rows = []

    for _, row in ecc.iterrows():
        key = tuple(map(str, row.iloc[:7]))

        chrom = row.iloc[0]
        mito = "YES" if chrom in ["chrM", "MT", "chrMT"] else "NO"

        output_rows.append([
            *row.iloc[:7],
            row.iloc[7] if len(row) > 7 else "",  # ID si présent
            gene_dict.get(key, ""),
            te_dict.get(key, ""),
            alu_dict.get(key, ""),
            mito
        ])

    columns = [
        "chrom","chromStart","chromEnd","Nb_reads",
        "Nb_hits","Peak_size","Pvalue","ID",
        "Gene","TE","ALU","Mito"
    ]

    df_out = pd.DataFrame(output_rows, columns=columns)
    df_out.to_csv(output_file, sep="\t", index=False)

    print(f"Fichier généré : {output_file}")

######################### parse gff ####################################
def parse_attributes(attributes):
    """
    Parse automatiquement GTF ou GFF (RepeatMasker)
    """

    # --- cas GTF classique ---
    gene_name = re.search(r'gene_name "([^"]+)"', attributes)
    if gene_name:
        return gene_name.group(1)

    gene_id = re.search(r'gene_id "([^"]+)"', attributes)
    if gene_id:
        return gene_id.group(1)

    # --- cas GFF RepeatMasker ---
    motif = re.search(r'Motif:([^"\s]+)', attributes)
    if motif:
        return motif.group(1)

    # fallback : ID=
    gff_id = re.search(r'ID=([^;]+)', attributes)
    if gff_id:
        return gff_id.group(1)

    return ""


def intersect_ecc_flexible(
    ecc_file,
    annot_file,
    output_file,
    overlap_fraction=0.2
):
    # --- ECC ---
    ecc = pd.read_csv(ecc_file, sep=None, engine="python")
    ecc_bed = BedTool.from_dataframe(ecc.iloc[:, :7])

    # --- ANNOT ---
    annot_bed = BedTool(annot_file)
    inter = ecc_bed.intersect(annot_bed, wa=True, wb=True, f=overlap_fraction)

    annot_dict = {}

    for entry in inter:
        fields = entry.fields
        key = tuple(fields[:7])

        attributes = " ".join(fields[8:])  # colonne attributs

        value = parse_attributes(attributes)

        if value and "Unclassified" not in value:
            annot_dict.setdefault(key, set()).add(value)

    # --- SORTIE ---
    output_rows = []

    for _, row in ecc.iterrows():
        key = tuple(map(str, row.iloc[:7]))

        chrom = row.iloc[0]
        mito = "YES" if chrom in ["chrM", "MT", "chrMT"] else "NO"

        values = ";".join(annot_dict.get(key, []))

        output_rows.append([
            *row.iloc[:7],
            #row.iloc[7] if len(row) > 7 else "",
            values,
            #mito
        ])

    columns = [
        "chrom","chromStart","chromEnd","Nb_reads",
        "Nb_hits","Peak_size","Pvalue",
        "Annotation"
    ]

    df_out = pd.DataFrame(output_rows, columns=columns)
    df_out.to_csv(output_file, sep="\t", index=False)

    print(f"Fichier généré : {output_file}")

#intersect_ecc(
#    ecc_file="./OUTPUTACOBIOM_45_SA/ecc.ont_fusion.csv",
#    gene_annot=os.path.expanduser("~/Documents/annotations/gencode.v38.annotation.bed"),
#    te_annot=os.path.expanduser("~/Documents/annotations/hg38_rmsk_TE.sorted.gtf"),
#    alu_annot=os.path.expanduser("~/Documents/annotations/hg38_rmsk_TE.sorted_AluONLY.gtf"),
#    output_file="final_output.tsv"
#)
#bash intersect_all.sh final_circles_OUTPUTACOBIOM_17D_SA.csv ~/Documents/annotations/gencode.v38.annotation.bed ~/Documents/annotations/hg38_rmsk_TE.sorted.gtf ~/Documents/annotations/hg38_rmsk_TE.sorted_AluONLY.gtf  OUT_GEN_TE_ALU_MIT
