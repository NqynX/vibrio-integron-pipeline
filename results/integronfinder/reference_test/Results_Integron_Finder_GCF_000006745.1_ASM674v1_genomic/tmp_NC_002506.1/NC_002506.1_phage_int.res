# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.3.2 (Nov 2020); http://hmmer.org/
# Copyright (C) 2020 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /scratch/user/uqcngu19/envs/integronfinder/lib/python3.9/site-packages/integron_finder/data/Models/phage-int.hmm
# target sequence database:        /scratch/user/uqcngu19/vibrio-integron-pipeline/results/integronfinder/reference_test/Results_Integron_Finder_GCF_000006745.1_ASM674v1_genomic/tmp_NC_002506.1/NC_002506.1.prt
# output directed to file:         /scratch/user/uqcngu19/vibrio-integron-pipeline/results/integronfinder/reference_test/Results_Integron_Finder_GCF_000006745.1_ASM674v1_genomic/tmp_NC_002506.1/NC_002506.1_phage_int.res
# per-seq hits tabular output:     /scratch/user/uqcngu19/vibrio-integron-pipeline/results/integronfinder/reference_test/Results_Integron_Finder_GCF_000006745.1_ASM674v1_genomic/tmp_NC_002506.1/NC_002506.1_phage_int_table.res
# model-specific thresholding:     GA cutoffs
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       Phage_integrase  [M=173]
Accession:   PF00589.16
Description: Phage integrase family
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence        Description
    ------- ------ -----    ------- ------ -----   ---- --  --------        -----------
    4.9e-47  156.0   0.2    7.3e-47  155.4   0.2    1.3  1  NC_002506.1_269  # 309750 # 310712 # -1 # ID=1_269;partial=00
    6.8e-15   51.3   0.0    1.3e-14   50.4   0.0    1.3  1  NC_002506.1_259  # 300473 # 301624 # 1 # ID=1_259;partial=00;


Domain annotation for each sequence (and alignments):
>> NC_002506.1_269  # 309750 # 310712 # -1 # ID=1_269;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  155.4   0.2   1.4e-49   7.3e-47       2     171 ..     105     307 ..     104     309 .. 0.96

  Alignments for each domain:
  == domain 1  score: 155.4 bits;  conditional E-value: 1.4e-49
                      HHHHHHHHHHHHCCCTHHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHHHHH.HHTTSTTS CS
  Phage_integrase   2 vLtedeverllaaleeslsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkeilsdrkkeaeere 95 
                      vLt de++rll+ ++     +++l ++ll++ glR  E+++lrv+d+d+d g ir+   ++K++k+rtv+l +el++ lke+++ + k + +r+
  NC_002506.1_269 105 VLTRDEIRRLLEIVDP----KHQLPIKLLYGSGLRLMECMRLRVQDIDFDYGAIRI--WQGKGGKNRTVTLAKELYPHLKEQIALA-KRYYDRD 191
                      89**************....************************************..**************************99.4455555 PP

                      .............................BSSBEC...........TSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HH CS
  Phage_integrase  96 .............................llfvsk...........rgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlk 149
                                                   +lf+s            r+++++++ +++a++r ++eagiek +t HtLRhsfat+lle G+d++
  NC_002506.1_269 192 lhqknyggvwlptalkekypnapyefrwhYLFPSFqlsldpesdvmRRHHMNETVLQKAVRRSAQEAGIEKTVTCHTLRHSFATHLLEVGADIR 285
                      66678999************************************************************************************** PP

                      HHHHH----SHHHHHHHHCCSH CS
  Phage_integrase 150 vvqkllGHssisttkiYthvak 171
                      +vq++lGH++++tt+iYthv +
  NC_002506.1_269 286 TVQEQLGHTDVKTTQIYTHVLD 307
                      *******************987 PP

>> NC_002506.1_259  # 300473 # 301624 # 1 # ID=1_259;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;g
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   50.4   0.0   2.5e-17   1.3e-14       5     166 ..     207     375 ..     204     379 .. 0.91

  Alignments for each domain:
  == domain 1  score: 50.4 bits;  conditional E-value: 2.5e-17
                      HHHHHHHHHCCCT.HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHHHHH.HHTTSTTSBS CS
  Phage_integrase   5 edeverllaalee.slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkeilsdrkkeaeerell 97 
                      e+ v++l +a+++ s + +++ +++l+l  + R +El+  +  d+dl+++++ vp+++ K +k+   ++ +++ ++ ++i+ ++   ++++++ 
  NC_002506.1_259 207 EQGVKALWKAIDDiSIHESNKNFLRLMLIFANRSNELRLAKKADFDLEKRVWTVPEENNKVRKKQGGAIRRAIPPLAEKIIMEQFAIWPNHTMM 300
                      566899999999999999**************************************************************************** PP

                      SBEC...TSSB..HHHHHHHHHHHHHHTT--..CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHH.HHH CS
  Phage_integrase  98 fvsk...rgkplsdstvnrafkravkeagie..keltpHtLRhsfatallesGvdlkvvqkllGHssistt.kiY 166
                      f+     +++p+s++      ++++ + +    ++ t H++R+++      +Gv+ +v + +lGH+    + ++Y
  NC_002506.1_259 301 FPPVnteQDRPMSANVPVDFGSKLADRIEELgfPRTTNHDMRRTARNIWESMGVPYHVAETMLGHKVHTGVqSHY 375
                      *98888999*****99888888888883333667*******************************9776665565 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                            1  (173 nodes)
Target sequences:                       1008  (301046 residues searched)
Passed MSV filter:                        34  (0.0337302); expected 20.2 (0.02)
Passed bias filter:                       26  (0.0257937); expected 20.2 (0.02)
Passed Vit filter:                         6  (0.00595238); expected 1.0 (0.001)
Passed Fwd filter:                         3  (0.00297619); expected 0.0 (1e-05)
Initial search space (Z):               1008  [actual number of targets]
Domain search space  (domZ):               2  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: 7777.02
//
[ok]
