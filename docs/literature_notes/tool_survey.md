# Integron Bioinformatics Tool Survey

| Tool | Purpose | Language | Last Updated | Installable? | Benchmarked? | Notes |
|------|---------|----------|-------------|-------------|-------------|-------|
| IntegronFinder 2 | intI + attC + integron classification | Python/R | 2022 | ✅ Yes | ✅ Yes | v2.0.2, primary tool |
| HattCI | attC detection (generalised HMM) | C | 2016 | ✅ Yes | ✅ Yes | v1.0b, compiled from source |
| I-VIP | class 1 integron ID + visualisation | Python | 2021 | ⚠️ Partial | ❌ No | Python 2/3 issues, class 1 only. See ivip_installation_notes.md |
| XXR | attC detection in Vibrio superintegrons | ? | 2003 | ❌ No | ❌ No | Not found in any public repository |
| IntegronFinder v1 | Original integron finder | Python | 2016 | ❌ No | ❌ No | Superseded by v2 |
| ACID | attC detection | ? | ? | ❌ No | ❌ No | Not reproducibly installable |
| ATTACCA/RAC | attC detection | ? | ? | ❌ No | ❌ No | Not reproducibly installable |

## Summary
Two tools (IntegronFinder 2 and HattCI) are included in quantitative benchmarking.
Five tools were excluded due to: unavailable source code (XXR, ACID, ATTACCA/RAC),
supersession (IntegronFinder v1), or reproducibility/scope limitations (I-VIP).
