Raw and archived input files for the PTB V1-V7 workflow.

Current R workflow:

1. Run:
   git/codes/00 - cleaning-pre-processing and descriptive/make_V1_V7_k_averaged.R

2. Main raw inputs used by that script:
   - UIC_QUS_8.8.22_C425updated.xlsx
       Raw V1/V2 QUS acquisitions. The script reads the "per-image (SLD)" sheet.
   - postpartum_data.csv
       Raw postpartum V3-V7 QUS acquisitions.

3. Clinical scaffold:
   - ../V1_V7_wideformat_k1.csv
       Current participant-level scaffold used to preserve the downstream
       clinical/outcome columns and participant set.

4. Main outputs written to git/data:
   - V1_V7_wideformat_k{k}.csv
   - V1_V7_longformat_k{k}.csv

Notes:

- k is the number of first V1/V2 QUS acquisition rows averaged within each
  Participant_ID and Visit_ID.
- V3-V7 postpartum QUS values are averaged over all available acquisitions by
  default, matching the previous R merge behavior.
- Old generated intermediates QUS_means_k.csv and fulldata_k.csv were removed
  from this folder because the new R script generates the final V1-V7 files
  directly.
- Set PTB_WRITE_QUS_MEANS=true before running the script only if you need the
  QUS_means_k.csv diagnostic intermediates regenerated.
