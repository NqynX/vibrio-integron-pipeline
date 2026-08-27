# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.3.2 (Nov 2020); http://hmmer.org/
# Copyright (C) 2020 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /scratch/user/uqcngu19/envs/integronfinder/lib/python3.9/site-packages/integron_finder/data/Models/integron_integrase.hmm
# target sequence database:        /scratch/user/uqcngu19/vibrio-integron-pipeline/results/integronfinder/reference_test/Results_Integron_Finder_GCF_000006745.1_ASM674v1_genomic/tmp_NC_002506.1/NC_002506.1.prt
# output directed to file:         /scratch/user/uqcngu19/vibrio-integron-pipeline/results/integronfinder/reference_test/Results_Integron_Finder_GCF_000006745.1_ASM674v1_genomic/tmp_NC_002506.1/NC_002506.1_intI.res
# per-seq hits tabular output:     /scratch/user/uqcngu19/vibrio-integron-pipeline/results/integronfinder/reference_test/Results_Integron_Finder_GCF_000006745.1_ASM674v1_genomic/tmp_NC_002506.1/NC_002506.1_intI_table.res
# model-specific thresholding:     GA cutoffs
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       intI_Cterm  [M=59]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence        Description
    ------- ------ -----    ------- ------ -----   ---- --  --------        -----------
    9.5e-27   88.4   1.8    1.7e-26   87.5   1.8    1.5  1  NC_002506.1_269  # 309750 # 310712 # -1 # ID=1_269;partial=00


Domain annotation for each sequence (and alignments):
>> NC_002506.1_269  # 309750 # 310712 # -1 # ID=1_269;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   87.5   1.8   1.7e-29   1.7e-26       2      59 .]     187     245 ..     186     245 .. 0.97

  Alignments for each domain:
  == domain 1  score: 87.5 bits;  conditional E-value: 1.7e-29
       intI_Cterm   2 lhekDl.aegyggVyLPnaLarKYPnaakelaWqylFPsaklsvdprsgelrRHHldes 59 
                      ++++Dl +++yggV+LP+aL++KYPna+ e++W+ylFPs +ls dp+s+++rRHH++e+
  NC_002506.1_269 187 YYDRDLhQKNYGGVWLPTALKEKYPNAPYEFRWHYLFPSFQLSLDPESDVMRRHHMNET 245
                      8*******************************************************996 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                            1  (59 nodes)
Target sequences:                       1008  (301046 residues searched)
Passed MSV filter:                        35  (0.0347222); expected 20.2 (0.02)
Passed bias filter:                       29  (0.0287698); expected 20.2 (0.02)
Passed Vit filter:                         3  (0.00297619); expected 1.0 (0.001)
Passed Fwd filter:                         1  (0.000992063); expected 0.0 (1e-05)
Initial search space (Z):               1008  [actual number of targets]
Domain search space  (domZ):               1  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: 5441.98
//
[ok]
