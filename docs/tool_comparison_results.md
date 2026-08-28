# Tool Comparison: IntegronFinder 2.0.2 vs HattCI 1.0b

## Benchmark Genome
- *Vibrio cholerae* N16961 chromosome II (NC_002506.1, 1,072,315 bp)
- Contains a well-characterised class 4 superintegron

## IntegronFinder 2.0.2 Results
- Complete integrons: 1 (12 attC sites + integrase intI)
- CALINs: 4 (52 attC sites)
- **Total attC sites: 64**
- Output: structured TSV + GFF + summary

## HattCI 1.0b Results
- Raw hits: 184 (183 with Vscore > 0)
- Vscore distribution: median = 2.82, range = 0.03–10.09
  - Vscore < 5: 152 (83.1%)
  - Vscore 5–10: 29 (15.8%)
  - Vscore >= 10: 2 (1.1%)
- No integrase detection, no integron/CALIN classification

## Positional Agreement
- Sites detected by both tools at exactly matching coordinates (0 bp): **55**
- IntegronFinder unique sites: 9 (across 4 integrons/CALINs, e-values 2.8e-4 to 3.7e-1)
- HattCI unique sites: 128 (majority likely false positives)

## Vscore Threshold Analysis (HattCI vs IntegronFinder as reference)

| Vscore > | HattCI hits | Recall | Precision | F1    |
|----------|-------------|--------|-----------|-------|
| 0        | 183         | 0.859  | 0.301     | 0.445 |
| 1        | 150         | 0.781  | 0.333     | 0.467 |
| 2        | 109         | 0.641  | 0.376     | 0.474 |
| 3        | 85          | 0.500  | 0.376     | 0.430 |
| 5        | 31          | 0.234  | 0.484     | 0.316 |
| 7.5      | 9           | 0.094  | 0.667     | 0.164 |

- Best F1: 0.474 at Vscore > 2 (still substantially below IF completeness)
- Recommended >7.5 threshold: recall of only 0.094, unsuitable for Vibrio
- No HattCI threshold achieved F1 > 0.50

## Summary Paragraph (for thesis)

IntegronFinder 2.0.2 and HattCI 1.0b were benchmarked on the *Vibrio cholerae* N16961 chromosome II (NC_002506.1, 1,072,315 bp), which harbours a well-characterised class 4 superintegron. IntegronFinder identified a total of 64 attC sites, distributed across one complete integron containing 12 attC cassettes and an associated integrase (intI), plus four CALIN (Cassettes Associated with Lacking Integron integrase) clusters comprising 52 additional attC sites. HattCI, using its default detection mode without Viterbi score filtering, returned 184 raw attC hits (183 with Vscore > 0). Direct positional comparison revealed that 55 attC sites (85.9% of IntegronFinder detections) were identified by both tools at exactly matching coordinates (0 bp start-position difference), providing a high-confidence consensus set. IntegronFinder uniquely detected 9 additional attC sites across four integron/CALIN regions (e-values ranging from 2.8e-4 to 3.7e-1) that HattCI failed to recover at any threshold, while HattCI reported 128 sites not detected by IntegronFinder, the majority of which are likely false positives given their low Viterbi scores. The Vscore distribution of HattCI detections was heavily skewed toward low confidence values (median = 2.82, range = 0.03–10.09), with 83.1% of hits (152/183) scoring below 5.0 and only 1.1% (2/183) exceeding 10.0. A systematic threshold analysis demonstrated that HattCI's recommended quality cutoff of Vscore > 7.5 retained only 9 hits, yielding a recall of just 0.094 against IntegronFinder's detections and an F1 score of 0.164. The optimal F1 score of 0.474 was achieved at Vscore > 2.0 (109 hits, recall = 0.641, precision = 0.376), though this still fell substantially below IntegronFinder's apparent completeness. No HattCI Vscore threshold achieved an F1 score exceeding 0.50, indicating that the tool's probabilistic model is poorly calibrated for the Vibrio superintegron context. Beyond attC detection, IntegronFinder provides additional capabilities essential for downstream cassette array analysis, including integrase (intI) identification, integron completeness classification (complete, CALIN, or incomplete), and structured output in standard TSV and GFF formats. Based on these results, IntegronFinder 2.0.2 was selected as the primary detection tool for all subsequent pipeline stages.
