# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.3.2 (Nov 2020); http://hmmer.org/
# Copyright (C) 2020 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  /scratch/user/uqcngu19/envs/integronfinder/lib/python3.9/site-packages/integron_finder/data/Models/phage-int.hmm
# target sequence database:        /scratch/user/uqcngu19/vibrio-integron-pipeline/results/integronfinder/reference_test/Results_Integron_Finder_GCF_000006745.1_ASM674v1_genomic/tmp_NC_002505.1/NC_002505.1.prt
# output directed to file:         /scratch/user/uqcngu19/vibrio-integron-pipeline/results/integronfinder/reference_test/Results_Integron_Finder_GCF_000006745.1_ASM674v1_genomic/tmp_NC_002505.1/NC_002505.1_phage_int.res
# per-seq hits tabular output:     /scratch/user/uqcngu19/vibrio-integron-pipeline/results/integronfinder/reference_test/Results_Integron_Finder_GCF_000006745.1_ASM674v1_genomic/tmp_NC_002505.1/NC_002505.1_phage_int_table.res
# model-specific thresholding:     GA cutoffs
# number of worker threads:        1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       Phage_integrase  [M=173]
Accession:   PF00589.16
Description: Phage integrase family
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence         Description
    ------- ------ -----    ------- ------ -----   ---- --  --------         -----------
    4.1e-52  173.8   0.1    6.8e-52  173.1   0.1    1.3  1  NC_002505.1_2259  # 2593375 # 2594283 # -1 # ID=1_2259;partia
      1e-45  153.0   1.8    1.2e-45  152.7   0.5    1.7  1  NC_002505.1_112   # 122005 # 122940 # 1 # ID=1_112;partial=00
    1.4e-28   97.2   0.2    2.4e-28   96.4   0.2    1.4  1  NC_002505.1_474   # 548780 # 550021 # -1 # ID=1_474;partial=0
    1.5e-28   97.1   0.0    2.8e-28   96.2   0.0    1.5  1  NC_002505.1_784   # 912856 # 914124 # -1 # ID=1_784;partial=0
    2.4e-27   93.2   0.3    4.5e-27   92.3   0.3    1.5  1  NC_002505.1_1641  # 1896092 # 1897327 # 1 # ID=1_1641;partial
    5.2e-14   49.8   0.0    6.3e-14   49.5   0.0    1.1  1  NC_002505.1_1642  # 1897423 # 1897878 # 1 # ID=1_1642;partial
    6.6e-11   39.7   0.9    1.3e-10   38.7   0.9    1.5  1  NC_002505.1_164   # 188166 # 189368 # -1 # ID=1_164;partial=0


Domain annotation for each sequence (and alignments):
>> NC_002505.1_2259  # 2593375 # 2594283 # -1 # ID=1_2259;partial=00;start_type=GTG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10b
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  173.1   0.1   1.8e-54   6.8e-52       3     172 ..     117     289 ..     115     290 .. 0.98

  Alignments for each domain:
  == domain 1  score: 173.1 bits;  conditional E-value: 1.8e-54
                       HHHHHHHHHHHCCCT..HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHHHHH....HHT CS
   Phage_integrase   3 Ltedeverllaalee..slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkeilsdr...kke 90 
                       L e +ve+ll a +   +l++rd++++ell++TglR+ El+sl++++++l++g++rv   ++K++ker vp++e+++e++++ l++       
  NC_002505.1_2259 117 LSEAQVEALLSAPDPqsPLELRDKAMLELLYATGLRVTELVSLTMENMSLRQGVVRV---MGKGGKERLVPMGENAIEWIETFLQQGrslLLG 206
                       8899**********999****************************************...***************************999999 PP

                       TSTTSBSSBECTSSB..HHHHHHHHHHHHHHTT--.CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHHHHHHCCSHH CS
   Phage_integrase  91 aeerellfvskrgkplsdstvnrafkravkeagie.keltpHtLRhsfatallesGvdlkvvqkllGHssisttkiYthvake 172
                       ++ ++++f+s rg++++++t+++++k+++  agi+ ++l+pH LRh fat+ll+ G+dl+vvq llGHs++stt+iYthva e
  NC_002505.1_2259 207 EQTSDIVFPSSRGQQMTRQTFWHRIKHYAVIAGIDvEKLSPHVLRHAFATHLLNYGADLRVVQMLLGHSDLSTTQIYTHVATE 289
                       ********************************************************************************976 PP

>> NC_002505.1_112  # 122005 # 122940 # 1 # ID=1_112;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  152.7   0.5   3.3e-48   1.2e-45       2     172 ..     122     291 ..     121     292 .. 0.97

  Alignments for each domain:
  == domain 1  score: 152.7 bits;  conditional E-value: 3.3e-48
                      HHHHHHHHHHHHCCCT.HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-HHHHHHHHHHHHH..HHTTST CS
  Phage_integrase   2 vLtedeverllaalee.slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvplseellevlkeilsdr.kkeaee 93 
                       L+ de+ +ll+  ++ +lsirdr+++el+++ glR +El+s+++kd++l +g+irv    +K++ker+v +  ++ e++ ++l+ r + +++ 
  NC_002505.1_112 122 NLDVDEMAQLLEVTDDdPLSIRDRAIMELMYGAGLRLAELVSIDIKDVNLSEGEIRV---IGKGNKERKVWFAGQAQEWVGKWLKLRsQLADSA 212
                      6999*****************************************************...**************************98888999 PP

                      TSBSSBECTSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHHHHHHCCSHH CS
  Phage_integrase  94 rellfvskrgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssisttkiYthvake 172
                      +++lfvsk g ++s++ v++++ ++++++++  +++pH+LRhsfat++les  +l++vq+llGH++i+tt+iYth++ +
  NC_002505.1_112 213 ETALFVSKLGTRISHRSVQKRMAEWGQKQAVASHISPHKLRHSFATHMLESSNNLRAVQELLGHENIATTQIYTHLDFQ 291
                      99*************************************************************************9987 PP

>> NC_002505.1_474  # 548780 # 550021 # -1 # ID=1_474;partial=00;start_type=ATG;rbs_motif=GGxGG;rbs_spacer=5-10bp;gc_con
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   96.4   0.2   6.5e-31   2.4e-28       3     170 ..     214     377 ..     212     380 .. 0.92

  Alignments for each domain:
  == domain 1  score: 96.4 bits;  conditional E-value: 6.5e-31
                      HHHHHHHHHHHCCCT.HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-...HHHHHHHHHHHHH.HHTTS CS
  Phage_integrase   3 Ltedeverllaalee.slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvpls...eellevlkeilsdrkkeae 92 
                      ++++ +++l++a+ + +  + +r+l+e++l T  R  E ++ r++d+d++n+++ +pa+++K+k+ +t+pl+    +lle +k i s+r     
  NC_002505.1_474 214 IEPERLPELMQAIANaNITLATRCLLEWQLHTMTRPIESATARWQDIDFKNKVWVIPAERMKMKRPHTIPLTeqtLALLEIMKPISSHR----- 302
                      78999**********99999****************************************************66666666666666666..... PP

                      TTSBSSBEC..TSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHH.HHHHCCS CS
  Phage_integrase  93 erellfvsk..rgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssistt.kiYthva 170
                        e++f+s   ++++++++ +n a+kr++++     +l +H LR  ++t l e+G++++ ++  l H + + + k+Y+++ 
  NC_002505.1_474 303 --EYIFPSNknPKSHVNSQSANMALKRMGYK----GQLVSHGLRALASTTLNEQGFNPDIIEAALAHVDKNEVrKAYNRAE 377
                      ..*******9999******************....5688**************************************9876 PP

>> NC_002505.1_784  # 912856 # 914124 # -1 # ID=1_784;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   96.2   0.0   7.7e-31   2.8e-28       2     171 ..     213     378 ..     212     380 .. 0.92

  Alignments for each domain:
  == domain 1  score: 96.2 bits;  conditional E-value: 7.7e-31
                      HHHHHHHHHHHHCCCT.HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-...HHHHHHHHHHHHH.HHTT CS
  Phage_integrase   2 vLtedeverllaalee.slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvpls...eellevlkeilsdrkkea 91 
                      +Lt+ e+++l+ a+ + s +  +r+l+e++l T  R +E+ + r+++++++++++ +pa+++K+++e+++pl+    +llev+k i ++     
  NC_002505.1_784 213 ALTPAELPELMSAIANaSIKRTTRCLLEWQLHTMTRPAEASGARWDEINWEEKVWTIPAERMKKRREHRIPLTeqmLALLEVMKPISGH----- 301
                      79**************99999****************************************************5544555555555555..... PP

                      STTSBSSBEC..TSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHH.HHHHCCSH CS
  Phage_integrase  92 eerellfvsk..rgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssistt.kiYthvak 171
                        r+++f+s   ++kp +++t+n a+kr++ +     +l +H LR  ++t l e+G+d++ v++ l H + + + ++Y++ ++
  NC_002505.1_784 302 --RDFIFPSDraPKKPCNSQTANMALKRMGFA----GRLVSHGLRSLASTTLNEQGFDPDLVESALAHVDDNQVrSAYNRTDY 378
                      ..5*******999****************999....6789**************************************99876 PP

>> NC_002505.1_1641  # 1896092 # 1897327 # 1 # ID=1_1641;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   92.3   0.3   1.2e-29   4.5e-27       1     171 [.     209     375 ..     209     377 .. 0.90

  Alignments for each domain:
  == domain 1  score: 92.3 bits;  conditional E-value: 1.2e-29
                       -HHHHHHHHHHHHCCCT.HHHHHHHHHHHHHHHHT--HHHHHC-BGGGEECTTEEEEEE..CCSSSCCEEEEE-...HHHHHHHHHHHHH.HH CS
   Phage_integrase   1 kvLtedeverllaalee.slsirdrllvellleTglRisEllslrvkdldldngtirvparetKtkkertvpls...eellevlkeilsdrkk 89 
                       k+L+++e+  l+ ++ + + +  +r+l+e++l T +R +E+ + r++++d+ n+ +++p++++K+++e+ vpl+    ++le++k i ++r  
  NC_002505.1_1641 209 KALQPHEMHDLIRTVATaNIQHVTRFLIEWQLHTMVRPNEASGARWEEIDMVNKLWIIPKERMKMNREHVVPLTaqtLAILEAIKPISGHR-- 299
                       57999***********97777778999***********************************************77777778888877777.. PP

                       TTSTTSBSSBEC..TSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HHHHHHH----SHHHH.HHHHCCSH CS
   Phage_integrase  90 eaeerellfvsk..rgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlkvvqkllGHssistt.kiYthvak 171
                            e++f+s   ++ p +++t+n+a+ r++ +    ++ t H LR  ++t l e+G+ ++v++  l H++ + + k+Y++ ++
  NC_002505.1_1641 300 -----EFIFPSSrnPKVPTDSETANKALGRMGFK----DRTTAHGLRALASTTLNEQGFEPDVIEAALAHTDKNQIrKAYNRTDY 375
                       .....*******6667777888888888877777....6899**************************************99876 PP

>> NC_002505.1_1642  # 1897423 # 1897878 # 1 # ID=1_1642;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_con
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   49.5   0.0   1.7e-16   6.3e-14      62     171 ..       1     105 [.       1     107 [. 0.87

  Alignments for each domain:
  == domain 1  score: 49.5 bits;  conditional E-value: 1.7e-16
                       CSSSCCEEEEE-...HHHHHHHHHHHHH.HHTTSTTSBSSBEC..TSSB..HHHHHHHHHHHHHHTT--CC-HHHHHHHHHHHHHHH----HH CS
   Phage_integrase  62 tKtkkertvpls...eellevlkeilsdrkkeaeerellfvsk..rgkplsdstvnrafkravkeagiekeltpHtLRhsfatallesGvdlk 149
                       +K+++e+++pl     +llev+k i ++       r+++f+s   ++kp +++t+n a+kr++ +     +l +H LR  ++t l e+G++++
  NC_002505.1_1642   1 MKKRREHRIPLKeqmLALLEVMKPISGH-------RDFIFPSDrdPKKPCNSQTANMALKRMGFA----GRLVSHGLRSLASTTLNEQGFNPD 82 
                       69**********5544555555555555.......5*******9999***************999....6789******************** PP

                       HHHHH----SHHHH.HHHHCCSH CS
   Phage_integrase 150 vvqkllGHssistt.kiYthvak 171
                        v+  l H + + + ++Y++ ++
  NC_002505.1_1642  83 LVEAALAHVDDNQVrSAYNRTDY 105
                       ******************99876 PP

>> NC_002505.1_164  # 188166 # 189368 # -1 # ID=1_164;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   38.7   0.9   3.6e-13   1.3e-10      21     169 ..     202     377 ..     184     381 .. 0.78

  Alignments for each domain:
  == domain 1  score: 38.7 bits;  conditional E-value: 3.6e-13
                      HHHHHHHH..HHHHHT--HHHHHC-BGGGEECTTEEEEEE..........CCSSSCCEEEEE-HHHHHHHHHHHHH.......HHTTSTTSBSS CS
  Phage_integrase  21 irdrllve..llleTglRisEllslrvkdldldngtirvp........aretKtkkertvplseellevlkeilsdr......kkeaeerellf 98 
                      ++d+++++  l+ + glR++E+++   + +++ n++i+          + +tK +k r+v ++++l+  + ++   +      kk++  +++l+
  NC_002505.1_164 202 CSDEFIIHqlLQIQSGLRVEEACTFPFSIVEMPNPHIHRYeveigihnGVHTKFNKTRKVEIPNQLMRKMYDYSVSErrlkreKKTDGVNKTLL 295
                      66777754116778***************99999988877899**999999************7777666666544447777788888899*** PP

                      BECTSSB..HHHHHHHHHHHHHHTT--....CC-HHHHHHHHHHHHHHH----.......HHHHHHH----SHHHHHHHHCC CS
  Phage_integrase  99 vskrgkplsdstvnrafkravkeagie....keltpHtLRhsfatallesGvd.......lkvvqkllGHssisttkiYthv 169
                      ++  g+pl+++++++ f+r+  + + +       ++H+LR ++ t  l+s +d       l  ++  +GH++ +tt +Y + 
  NC_002505.1_164 296 LNNLGNPLCSNNIQQHFRRLRHHIQNKhnivFSHRTHDLRATYGTYRLDSLLDhlpvgdaLALIMGWMGHKDDKTTWKYLRY 377
                      ********************99944444566999***********8777755433333337889999*********999765 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                            1  (173 nodes)
Target sequences:                       2597  (863360 residues searched)
Passed MSV filter:                        80  (0.0308048); expected 51.9 (0.02)
Passed bias filter:                       56  (0.0215633); expected 51.9 (0.02)
Passed Vit filter:                        16  (0.00616095); expected 2.6 (0.001)
Passed Fwd filter:                         9  (0.00346554); expected 0.0 (1e-05)
Initial search space (Z):               2597  [actual number of targets]
Domain search space  (domZ):               7  [number of targets reported over threshold]
# CPU time: 0.01u 0.01s 00:00:00.02 Elapsed: 00:00:00.06
# Mc/sec: 2138.24
//
[ok]
