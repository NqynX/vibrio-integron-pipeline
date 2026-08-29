#!/usr/bin/env python3"""Multi-genome integron comparative visualization using pyGenomeViz.Designed for Google Colab — clone the repo, install pygenomeviz, then run this script.Genomes included:  - V. cholerae N16961 (superintegron on chrII)  - V. cholerae O395 (attC-only region on chrII, no integrase)  - V. parahaemolyticus RIMD2210633 (superintegron on chrI)  - A. fischeri ES114 (small integron)  - V. vulnificus CMCP6: SKIPPED (no integrons)"""import osimport reimport pandas as pdfrom pathlib import Pathfrom collections import defaultdictimport pygenomeviz as pgv# ============================================================# CONFIGURATION# ============================================================REPO_DIR = Path("vibrio-integron-pipeline")IF_DIR = REPO_DIR / "results" / "integron_detection" / "integronfinder"CASSETTE_DIR = REPO_DIR / "results" / "cassette_reconstruction"ANNOT_DIR = REPO_DIR / "results" / "annotation"OUTPUT_DIR = Path("integron_figures")OUTPUT_DIR.mkdir(exist_ok=True)# Genome metadata: directory name -> display name + expected integron chromosome#GENOMES = {
    "V_cholerae_N16961": {
        "display": "V. cholerae N16961",
        "chrom": "NC_002506.1",  # chrII — superintegron
    },
    "V_cholerae_O395": {
        "display": "V. cholerae O395",
        "chrom": "NC_009456.1",  # chrII — attC-only region
    },
    "V_parahaemolyticus": {
        "display": "V. parahaemolyticus RIMD",
        "chrom": "NC_004603.1",  # chrI — superintegron (unusual)
    },
    "A_fischeri": {
        "display": "A. fischeri ES114",
        "chrom": None,  # auto-detect from .integrons file
    },
}

# BLASTx TSV filenames (N16961 has no genome prefix)
BLASTX_FILES = {
    "V_cholerae_N16961": "cassettes_blastx_swissprot.tsv",
    "V_cholerae_O395": "V_cholerae_O395_blastx_swissprot.tsv",
    "V_parahaemolyticus": "V_parahaemolyticus_blastx_swissprot.tsv",
    "A_fischeri": "A_fischeri_blastx_swissprot.tsv",
}

# Colour palette (colourblind-safe, Paul Tol)
CATEGORY_COLORS = {
    "no_hit":           "#BBBBBB",  # light gray — no BLASTx hit
    "hypothetical":      "#DDDDDD",  # very light gray
    "toxin_antitoxin":   "#0077BB",  # blue
    "amr":               "#CC3311",  # red
    "metabolic":         "#009988",  # teal
    "transport":         "#33BBEE",  # cyan
    "regulatory":        "#EE7733",  # orange
    "mobile_element":    "#EE3377",  # pink
    "other_known":       "#AA4499",  # purple
}

INTEGRASE_COLOR = "#222222"   # dark — integrase gene
ATTC_COLOR = "#FFD700"        # gold — attC site markers

# ============================================================# HELPER FUNCTIONS# ============================================================


def find_integrons_file(genome_dir):
    """Find the .integrons TSV file inside the IntegronFinder results subdir."""
    results_subdir = list(genome_dir.glob("Results_Integron_Finder_*"))
    if not results_subdir:
        return None
    integrons_files = list(results_subdir[0].glob("*.integrons"))
    return integrons_files[0] if integrons_files else None


def parse_integrons_file(filepath):
    """Parse IntegronFinder .integrons TSV.

    Returns dict with:
      - replicons: {replicon_id: {attCs: [(start,end,strand),...], integrase: (start,end,strand) or None}}
      - has_complete_integron: bool
    """
    replicons = defaultdict(lambda: {"attcs": [], "integrase": None})
    has_complete = False
    if filepath is None:
        return replicons, has_complete

    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue
            replicon = cols[0]
            feat_type = cols[2]
            try:
                start = int(cols[3])
                end = int(cols[4])
                strand = int(cols[5]) if len(cols) > 5 else 1
            except (ValueError, IndexError):
                continue

            if feat_type == "complete_integron":
                has_complete = True
            elif feat_type == "intI":
                replicons[replicon]["integrase"] = (start, end, strand)
            elif feat_type.startswith("attc"):
                replicons[replicon]["attcs"].append((start, end, strand))

    # Sort attC sites by position
    for rep_data in replicons.values():
        rep_data["attcs"].sort(key=lambda x: x[0])

    return dict(replicons), has_complete


def parse_cassette_tsv(filepath):
    """Parse cassette reconstruction TSV.

    Expects columns with start/end info. Returns list of dicts:
      [{cassette_id, start, end, length}, ...]
    """
    if filepath is None or not filepath.exists():
        return []
    df = pd.read_csv(filepath, sep="\t")
    if df.empty:
        return []
    # Auto-detect column names (case-insensitive partial match)
    col_map = {c.lower(): c for c in df.columns}
    id_col = col_map.get("cassette_id", col_map.get("cassette", col_map.get("id", None)))
    start_col = col_map.get("start", col_map.get("cassette_start", None))
    end_col = col_map.get("end", col_map.get("cassette_end", None))
    len_col = col_map.get("length", col_map.get("cassette_length", col_map.get("size", None)))

    if start_col is None or end_col is None:
        print(f"  WARNING: Could not detect start/end columns in {filepath.name}")
        print(f"  Columns found: {list(df.columns)}")
        return []

    cassettes = []
    for _, row in df.iterrows():
        c = {
            "start": int(row[start_col]),
            "end": int(row[end_col]),
        }
        if id_col:
            c["cassette_id"] = str(row[id_col])
        else:
            c["cassette_id"] = f"cassette_{len(cassettes)+1}"
        if len_col:
            c["length"] = int(row[len_col])
        else:
            c["length"] = c["end"] - c["start"]
        cassettes.append(c)
    return cassettes


def parse_blastx_tsv(filepath):
    """Parse BLASTx outfmt 6 TSV.

    Returns dict: {query_id: {pident, bitscore, stitle, category}}
    Only keeps the top hit (best bitscore) per query.
    """
    if filepath is None or not filepath.exists():
        return {}
    df = pd.read_csv(filepath, sep="\t", header=None,
                      names=["qseqid", "sseqid", "pident", "length", "mismatch",
                              "gapopen", "qstart", "qend", "sstart", "send",
                              "evalue", "bitscore", "stitle"])
    if df.empty:
        return {}
    # Keep top hit per query (highest bitscore)
    best = df.loc[df.groupby("qseqid")["bitscore"].idxmax()]
    result = {}
    for _, row in best.iterrows():
        qid = str(row["qseqid"])
        stitle = str(row["stitle"]) if pd.notna(row["stitle"]) else ""
        category = categorize_annotation(stitle)
        result[qid] = {
            "pident": float(row["pident"]),
            "bitscore": float(row["bitscore"]),
            "stitle": stitle,
            "category": category,
        }
    return result


def categorize_annotation(stitle):
    """Categorize a Swiss-Prot protein description into a functional group."""
    s = stitle.lower()
    if not s or s == "nan":
        return "no_hit"
    # AMR
    if any(kw in s for kw in ["resistance", "antibiotic", "chloramphenicol", "beta-lactamase",
                                "catb", "kanamycin", "tetracycline", "sulfonamide"]):
        return "amr"
    # Toxin-antitoxin
    if any(kw in s for kw in ["toxin", "antitoxin", "rele", "relb", "ccda", "ccdb",
                                "higa", "higb", "pare", "pard", "mazf", "maze",
                                "chpk", "vagc", "vagd", "pin", "doc", "phd", "yoe"]):
        return "toxin_antitoxin"
    # Transport
    if any(kw in s for kw in ["transport", "permease", "abc transporter", "porin",
                                "efflux", "importer", "exporter", "symporter",
                                "antiporter", "channel"]):
        return "transport"
    # Regulatory
    if any(kw in s for kw in ["regulator", "transcriptional", "repressor", "activator",
                                "sigma", "response regulator", "two-component",
                                "transcription factor", "lytic repressor"]):
        return "regulatory"
    # Metabolic
    if any(kw in s for kw in ["kinase", "phosphatase", "synthase", "synthetase",
                                "dehydrogenase", "oxidoreductase", "hydrolase",
                                "transferase", "isomerase", "lyase", "ligase",
                                "protease", "peptidase", "nuclease", "polymerase",
                                "methyltransferase", "decarboxylase", "reductase",
                                "oxygenase", "cyclase", "phosphorylase"]):
        return "metabolic"
    # Mobile elements
    if any(kw in s for kw in ["transposase", "integrase", "recombinase", "insertion",
                                "mobilization", "conjugation", "plasmid"]):
        return "mobile_element"
    # Hypothetical / uncharacterised
    if any(kw in s for kw in ["hypothetical", "uncharacterized", "uncharacterised",
                                "duf", "unknown protein", "putative",
                                "predicted protein", "ygg", "yjba", "yci"]):
        return "hypothetical"
    return "other_known"


def get_replicon_id(genome_name, replicons):
    """Determine which replicon contains the integron for a given genome."""
    meta = GENOMES.get(genome_name, {})
    expected_chrom = meta.get("chrom")
    if expected_chrom and expected_chrom in replicons:
        return expected_chrom
    # Auto-detect: pick replicon with the most attC sites
    best_rep = None
    best_count = 0
    for rep_id, data in replicons.items():
        n_atcc = len(data["attcs"])
        if n_atcc > best_count:
            best_count = n_atcc
            best_rep = rep_id
    return best_rep


# ============================================================# MAIN: BUILD THE FIGURE# ============================================================

def main():
    print("=" * 60)
    print("pyGenomeViz Multi-Genome Integron Comparison")
    print("=" * 60)

    # ---- Step 1: Parse all data per genome ----
    genome_data = {}  # {genome_name: {replicon, attcs, integrase, cassettes, blastx, has_complete}}

    for genome_name, meta in GENOMES.items():
        print(f"\n--- {meta['display']} ---")
        genome_dir = IF_DIR / genome_name
        integrons_file = find_integrons_file(genome_dir)
        if integrons_file is None:
            print(f"  No .integrons file found, skipping.")
            continue
        replicons, has_complete = parse_integrons_file(integrons_file)
        if not replicons:
            print(f"  No integron data found, skipping.")
            continue
        rep_id = get_replicon_id(genome_name, replicons)
        rep_data = replicons[rep_id]
        n_atcc = len(rep_data["attcs"])
        print(f"  Replicon: {rep_id}")
        print(f"  attC sites: {n_atcc}")
        print(f"  Integrase: {'Yes' if rep_data['integrase'] else 'No'}")
        print(f"  Complete integron: {has_complete}")

        # Parse cassette TSV
        cassette_file = CASSETTE_DIR / f"{genome_name}_cassettes.tsv"
        cassettes = parse_cassette_tsv(cassette_file)
        print(f"  Cassettes: {len(cassettes)}")

        # Parse BLASTx results
        blastx_file = ANNOT_DIR / BLASTX_FILES.get(genome_name, f"{genome_name}_blastx_swissprot.tsv")
        blastx = parse_blastx_tsv(blastx_file)
        n_annotated = sum(1 for v in blastx.values() if v["category"] not in ("no_hit", "hypothetical"))
        print(f"  BLASTx hits (non-hypothetical): {n_annotated}/{len(cassettes)}")

        genome_data[genome_name] = {
            "replicon": rep_id,
            "attcs": rep_data["attcs"],
            "integrase": rep_data["integrase"],
            "cassettes": cassettes,
            "blastx": blastx,
            "has_complete": has_complete,
            "display": meta["display"],
        }

    if not genome_data:
        print("\nNo genome data found. Check file paths.")
        return

    # ---- Step 2: Create Genomeviz figure ----
    # Calculate track sizes (integron region span + padding)
    PADDING = 5000  # 5 kb padding on each side
    track_info = {}
    for gname, gdata in genome_data.items():
        attcs = gdata["attcs"]
        if not attcs:
            continue
        region_start = attcs[0][0]
        region_end = attcs[-1][1]
        # Include integrase if it's upstream
        if gdata["integrase"]:
            int_start, int_end, _ = gdata["integrase"]
            region_start = min(region_start, int_start)
            region_end = max(region_end, int_end)
        track_start = max(0, region_start - PADDING)
        track_end = region_end + PADDING
        track_info[gname] = {
            "size": track_end - track_start,
            "offset": track_start,  # offset to convert genome coords to track coords
            "region_start": region_start,
            "region_end": region_end,
        }
        print(f"\n{gdata['display']}: track region {track_start:,}-{track_end:,} bp ({(track_end-track_start)/1000:.1f} kb)")

    # Create figure
    n_tracks = len(track_info)
    gv = pgv.Genomeviz(
        fig_width=16,
        fig_track_height=1.0,
        track_align_type="center",
        link_curve=False,
        feature_track_ratio=0.30,
        scale_style="full",  # show scale bar for each track
    )

    # ---- Step 3: Add tracks and features ----
    for gname, tinfo in track_info.items():
        gdata = genome_data[gname]
        offset = tinfo["offset"]
        display_name = gdata["display"]
        n_cass = len(gdata["cassettes"])
        n_atcc = len(gdata["attcs"])
        label = f"{display_name} ({n_cass} cassettes, {n_atcc} attCs)"
        if not gdata["has_complete"]:
            label += " [attC-only]"
        print(f"\nAdding track: {label}")

        track = gv.add_track(
            track_name=gname,
            size=tinfo["size"],
            labelsize=11,
        )

        # Add integrase gene
        if gdata["integrase"]:
            int_s, int_e, int_strand = gdata["integrase"]
            track.add_feature(
                start=int_s - offset,
                end=int_e - offset,
                strand=int_strand,
                label="intI",
                facecolor=INTEGRASE_COLOR,
                edgecolor="black",
                linewidth=0.8,
                labelsize=9,
                plotstyle="bigarrow",
            )
            print(f"  + intI at {int_s:,}-{int_e:,} bp")

        # Add cassettes coloured by BLASTx annotation
        blastx = gdata["blastx"]
        for i, cass in enumerate(gdata["cassettes"]):
            # Try to match cassette to BLASTx result
            cass_id = cass.get("cassette_id", f"cassette_{i+1}")
            # Try various ID matching strategies
            blx_hit = None
            for query_id in blastx:
                if cass_id in query_id or query_id in cass_id:
                    blx_hit = blastx[query_id]
                    break
            if blx_hit is None:
                # Try numeric match
                for query_id in blastx:
                    if cass_id.split("_")[-1] in query_id.split("_")[-1]:
                        blx_hit = blastx[query_id]
                        break
            category = blx_hit["category"] if blx_hit else "no_hit"
            color = CATEGORY_COLORS[category]
            track.add_feature(
                start=cass["start"] - offset,
                end=cass["end"] - offset,
                strand=1,
                facecolor=color,
                edgecolor="#333333",
                linewidth=0.3,
                labelsize=6,
                plotstyle="box",
            )

        # Add attC site markers (small features)
        for attc_s, attc_e, attc_strand in gdata["attcs"]:
            mid = (attc_s + attc_e) // 2
            marker_half = 50  # small marker
            track.add_feature(
                start=mid - offset - marker_half,
                end=mid - offset + marker_half,
                strand=attc_strand,
                facecolor=ATTC_COLOR,
                edgecolor="#B8860B",
                linewidth=0.2,
                plotstyle="box",
                labelsize=0,
                alpha=0.7,
            )

    # ---- Step 4: Save figure ----
    print(f"\nSaving figures to {OUTPUT_DIR}/...")
    gv.savefig(str(OUTPUT_DIR / "integron_comparison.png"), dpi=300)
    gv.savefig(str(OUTPUT_DIR / "integron_comparison.svg"))
    gv.savefig(str(OUTPUT_DIR / "integron_comparison.pdf"))
    print("Done! Files saved:")
    print(f"  {OUTPUT_DIR / 'integron_comparison.png'}")
    print(f"  {OUTPUT_DIR / 'integron_comparison.svg'}")
    print(f"  {OUTPUT_DIR / 'integron_comparison.pdf'}")

    # ---- Step 5: Print legend ----
    print(f"\n{'='*60}")
    print("COLOUR LEGEND")
    print(f"{'='*60}")
    for cat, color in CATEGORY_COLORS.items():
        label = cat.replace("_", " ").title()
        print(f"  {color}  {label}")
    print(f"  {INTEGRASE_COLOR}  Integrase (intI)")
    print(f"  {ATTC_COLOR}  attC recombination site")


if __name__ == "__main__":
    main()
