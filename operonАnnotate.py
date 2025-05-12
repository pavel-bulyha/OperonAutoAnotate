from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import operondbSearch


def get_gene_dict(record):
    """
    Creates a dictionary that maps the gene name (converted to lowercase and with
    underscores removed) to its feature from the GenBank record.
    Features of type "gene" or "CDS" are considered. First, it looks for the
    "gene" qualifier, then "locus_tag", and finally "note".
    """
    gene_dict = {}
    for feature in record.features:
        # Note: feature.type can be exactly "gene" – using a check with 'in' or equality
        if feature.type == "gene" or feature.type == "CDS":
            gene_name = None
            if "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
            elif "locus_tag" in feature.qualifiers:
                gene_name = feature.qualifiers["locus_tag"][0]
            elif "note" in feature.qualifiers:
                gene_name = feature.qualifiers["note"][0]
            if gene_name:
                key = gene_name.lower().replace("_", "")
                if key not in gene_dict:
                    gene_dict[key] = feature
    return gene_dict

def filter_operon_features(features):
    """
    Accepts a list of candidate annotations (SeqFeature objects with type "operon"),
    each containing coordinates and qualifiers.

    For each candidate, the following are computed:
      - start, end, and the operon length (end - start)
      - a flag (has_info): True if at least one of the "definition" or "name" qualifiers
        contains non-empty text
      - order — the index of the candidate in the original list

    Then:
      1. Operons with has_info == True remain unfiltered.
      2. For operons with has_info == False, filtering is applied:
         If operon A is completely contained within operon B
         (B.start <= A.start and A.end <= B.end) and the length of B is strictly greater
         than the length of A, then A is excluded.
         If the lengths are equal, the one with the smaller index (i.e., encountered earlier) is kept.

    Returns the filtered list of SeqFeature objects.
    """
    features_with_meta = []
    for i, feat in enumerate(features):
        start = int(feat.location.start)
        end = int(feat.location.end)
        length = end - start
        qualifiers = feat.qualifiers
        has_info = False
        # Set has_info to True if at least one of the "definition" or "name" contains non-empty text
        if "definition" in qualifiers and qualifiers["definition"].strip():
            has_info = True
        if "name" in qualifiers and qualifiers["name"].strip():
            has_info = True
        features_with_meta.append({
            "index": i,
            "start": start,
            "end": end,
            "length": length,
            "has_info": has_info,
            "feature": feat
        })

    # Split into two groups:
    groupA = [f for f in features_with_meta if f["has_info"]]
    groupB = [f for f in features_with_meta if not f["has_info"]]

    filteredB = []
    for candidate in groupB:
        keep = True
        for other in groupA + groupB:
            if other["index"] == candidate["index"]:
                continue
            # Check if candidate is completely contained within other
            if other["start"] <= candidate["start"] and candidate["end"] <= other["end"]:
                # If the length of other is strictly greater or, if equal, the other appeared earlier,
                # then exclude candidate
                if (other["length"] > candidate["length"]) or (
                        other["length"] == candidate["length"] and other["index"] < candidate["index"]):
                    keep = False
                    break
        if keep:
            filteredB.append(candidate)

    final = groupA + filteredB
    # Sort by the original order (index) to preserve the priority of the first occurrence
    final_sorted = sorted(final, key=lambda x: x["index"])
    return [item["feature"] for item in final_sorted]

def annotate_operons(genbank_file, output_file, include_all=True):
    """
    Reads a GenBank file, extracts the organism name, finds operons using the
    function fetch_operon_data, determines the operon boundaries using the first
    and last gene, and adds a new operon feature.

    The operon boundaries are defined as:
      start = the minimum coordinate of the first and last gene;
      end = the maximum coordinate of the first and last gene.
    The strand of the operon is taken into account.
    The new feature receives the following qualifiers: operon_id, name, definition,
    first_gene, last_gene, and note.
    """
    # Read the GenBank file (assumes a single record)
    record = SeqIO.read(genbank_file, "genbank")

    organism = record.annotations.get("organism", "").strip()
    if not organism:
        print("Organism name not found in the GenBank file.")
        sys.exit(1)
    print(f"Detected organism: {organism}")

    # Use the operon search function (imported from operondbSearch)
    df_operons = operondbSearch.fetch_operon_data(organism, include_all)
    print(f"Found operons: {len(df_operons)}")

    if df_operons.empty:
        print("No operons found, exiting.")
        sys.exit(0)

    gene_dict = get_gene_dict(record)

    new_operon_features = []
    for idx, row in df_operons.iterrows():
        genes_list = row["Genes"]
        if not genes_list:
            print(f"Skipping row {idx}: empty list of genes.")
            continue
        # Take the first and last gene from the list
        first_gene = genes_list[0]
        last_gene = genes_list[-1]
        # Convert names to lowercase and remove underscores for matching
        first_gene_key = first_gene.lower().replace("_", "")
        last_gene_key = last_gene.lower().replace("_", "")

        if first_gene_key not in gene_dict:
            print(f"Gene {first_gene} not found in GenBank for the operon, row {idx}.")
            continue
        if last_gene_key not in gene_dict:
            print(f"Gene {last_gene} not found in GenBank for the operon, row {idx}.")
            continue

        feature_first = gene_dict[first_gene_key]
        feature_last = gene_dict[last_gene_key]

        start = min(int(feature_first.location.start), int(feature_last.location.start))
        end = max(int(feature_first.location.end), int(feature_last.location.end))
        strand = feature_first.location.strand  # Assume the operon is unidirectional

        qualifiers = {}
        qualifiers["operon_id"] = row["Operon ID"]
        # If there is data in the Name/Definition columns, add them
        if "Name" in row and row["Name"]:
            qualifiers["name"] = row["Name"]
        if "Definition" in row and row["Definition"]:
            qualifiers["definition"] = row["Definition"]
        qualifiers["first_gene"] = first_gene
        qualifiers["last_gene"] = last_gene
        qualifiers["note"] = f"Operon boundaries: {start}-{end}; strand: {strand}"

        new_feature = SeqFeature(FeatureLocation(start, end, strand=strand), type="operon", qualifiers=qualifiers)
        new_operon_features.append(new_feature)

    if not new_operon_features:
        print("Failed to create any annotations for operons.")
        sys.exit(0)

    # Filter operons according to the specified hierarchy:
    filtered_features = filter_operon_features(new_operon_features)
    print(f"After filtering, {len(filtered_features)} operons remain.")

    # Add the filtered operon features to the existing GenBank record
    record.features.extend(filtered_features)

    SeqIO.write(record, output_file, "genbank")
    print(f"Annotated GenBank record saved to: {output_file}")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description="Annotate operons in a GenBank file using data fetched from https://operondb.jp/search."
    )
    parser.add_argument(
        "input_file",
        help="Path to the input GenBank file (e.g., Haemophilus_influenzae.gb)."
    )
    parser.add_argument(
        "-o", "--output",
        default="annotated.gb",
        help="Path for the output annotated GenBank file (default: annotated.gb)."
    )
    parser.add_argument(
        "--include_all",
        dest="include_all",
        action="store_true",
        help="Include all rows from OperonDB (default behavior)."
    )
    parser.add_argument(
        "--no_include_all",
        dest="include_all",
        action="store_false",
        help="Include only rows with a non-empty Definition."
    )

    parser.set_defaults(include_all=True)
    args = parser.parse_args()
    annotate_operons(args.input_file, args.output, include_all=args.include_all)