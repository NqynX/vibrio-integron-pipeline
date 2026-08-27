# I-VIP Installation and Compatibility Assessment

## Repository
- URL: https://github.com/caozhichongchong/I-VIP
- Cloned: 27 August 2026

## Issues Encountered

### 1. Wrong repository URL in literature
- The Zhang et al. (2018) paper does not provide a GitHub URL.
- The repo is at `caozhichongchong/I-VIP`, not `AnnieZhang2018/I-VIP`.

### 2. Python 2 code despite "Python 3" release label
- The 2021 release (I-VIP.2021release.zip) claims Python 3 compatibility.
- Both the git repo and the release zip contain Python 2 syntax (print statements).
- Required `2to3 -w` conversion to fix.

### 3. Biopython API incompatibility
- `feature.location.start.position` → `feature.location.start` (Biopython >= 1.78).
- Fixed with sed replacement in GbffParser.py.

### 4. Protein translation extraction failure
- GbffParser.py fails to extract CDS translations from GenBank files.
- Produces "N no amni acids found!" warnings (typo in original code).
- Results in malformed output file, causing IndexError in I-VIP.py line 145.
- Not fixed — would require further debugging of Biopython qualification access patterns.

### 5. Scope limitation
- I-VIP is designed specifically for **class 1 integrons** only.
- Uses IntI1-specific BLAST database and class 1-specific attC covariance model.
- Not applicable to chromosomal (class 3/4) integrons in Vibrio.

## Conclusion
I-VIP excluded from quantitative benchmarking due to:
1. Non-reproducible installation requiring multiple patches.
2. Incomplete Python 3 / modern Biopython compatibility.
3. Biological scope limited to class 1 integrons, not relevant for Vibrio chromosomal integron detection.
