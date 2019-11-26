# RiboMethSeqPipeline
You will find here the different RiboMethSeq pipelines used by the NGS Core Facility "EpiRNA-Seq" from Nancy, France (https://umsibslor.univ-lorraine.fr/en/facility/epitranscriptomics-sequencing-epirna-seq).

!!! Be aware !!!
All the scripts has been used in a local computer, with relatives paths for the inputs and outputs. 
If you want to use these scripts on your own, please be careful of the inputs and directories names along the script.

# Pipeline workflow #
Trimming (Unix script) --> Alignment (Unix Script) --> Fusion_5prime3prime_Counts.R --> Standard_RiboMethSeq_script_5prime3prime_Score2.R/RiboMethSeq_Performance_assay_5prime3prime_Score2.R                                                                                                    
If you want to use the old pipeline using the 'UCount5prime' files :
Standard_RiboMethSeq_script_5primeOnly_score6.R/RiboMethSeq_Performance_assay_5primeOnly_Score6.R
