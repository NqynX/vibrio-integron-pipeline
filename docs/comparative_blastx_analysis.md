# Comparative BLASTx Analysis of Integron Cassette Genes

## Summary Statistics

| Species | Total cassettes | With hits | % Annotated | High (>=80%) | Medium (50-80%) | Low (<50%) |
|---|---|---|---|---|---|---|
| V. cholerae N16961 (chrII) | 59 | 21 | 35.6% | 11 | 33 | 58 |
| V. cholerae O395 (chrII) | 50 | 20 | 40.0% | 9 | 4 | 7 |
| V. parahaemolyticus RIMD (chrI) | 44 | 16 | 36.4% | 0 | 3 | 13 |
| A. fischeri ES114 (chrII) | 3 | 1 | 33.3% | 0 | 0 | 1 |
| V. vulnificus CMCP6 | 0 | 0 | -- | -- | -- | -- |

Note: N16961 BLASTx used default max_target_seqs (multiple hits per cassette), so identity distribution counts all hits. Other genomes used max_target_seqs=1 (best hit per cassette).

## Overall: 156 cassettes across 4 genomes with integrons; 58 cassettes (37.2%) had Swiss-Prot BLASTx hits

## Key Functional Categories

### V. cholerae N16961 (21 unique cassettes with hits)
- AMR: chloramphenicol acetyltransferase CatB2 (cassette_008, ~67%)
- Toxin-antitoxin: ParD/ParE (cassette_025, 100%/46%), HigA-1/HigB-1 (cassette_031, 100%/100%), RelE (cassette_047, 56%)
- V. cholerae-specific: VC_A0308 (100%), VC_A0337 carboxypeptidase (99-100%), VC_A0354 (100%), VC_A0395 isomerase (100%)
- Outer membrane lipoprotein Blc (cassettes 021, 037, 046; ~99%)
- Metabolic: sulfate-binding protein Sbp (cassette_050, ~70%), glyoxalase domain (cassettes 018, 056; ~61-66%)
- N-acetyltransferases (cassette_043, ~27-32%)
- Phosphinothricin N-acetyltransferase (cassette_030, ~41%)

### V. cholerae O395 (20 cassettes with hits)
- Many N16961 homologs: Blc (3 copies, 99-100%), VC_A0337 (2 copies, 97-100%), ParD (100%), HigA-1 (100%), VC_A0354 (100%), VC_A0395 (100%)
- Toxin-antitoxin: RelE-like toxin (56%), TacT toxin (31%)
- Metabolic: sulfate-binding protein (66%), glyoxalase (66%)
- Transposase IS200 (41%)
- Phosphinothricin N-acetyltransferase (41%)

### V. parahaemolyticus RIMD (16 cassettes with hits)
- YoeB toxin (cassette_024, 65%)
- Most hits low identity (13/16 < 50%); many to non-bacterial proteins (likely spurious)
- Repeated: cold shock protein (2x), HIT-like protein (2x), IRC4-like (2x)

### A. fischeri ES114 (1 cassette with hit)
- Chloramphenicol 3-O phosphotransferase (cassette_002, 31.5%) - potential AMR relevance but low identity

## RGI AMR Detection
- N16961: 1 AMR gene (catB9)
- All other genomes: 0 AMR genes detected
- Total across 156 cassettes: 1 AMR gene (0.6%)
