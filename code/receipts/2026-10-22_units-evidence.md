# Units of the DBS, EEG, EMG and ECG channels stay n/a: no source states a scale

**Claim**: `units` stays n/a for the 603 DBS, EEG, EMG and ECG rows of the 44 `channels.tsv` (312 DBS, 273 EEG, 14 EMG,
4 ECG), and the README explains why no unit is asserted.

**Source**
- Computed for PR head `8ba48227260579488819b8d4061d7abdba377e7e` plus the revert to n/a on this branch. Files read from NEMAR
  nm000446 `origin/main` `57f0dd99c42dbb5c516e565401197aaf005c28b7`: all 44 `*_ieeg.vhdr`, `*_ieeg.eeg` and `*_channels.tsv`.
  The sha256 of each of these 132 files equals its SHA256E git-annex key.
- Release: MRC BNDU, doi:10.5287/bodleian:mzJ7YwXvo, read as the 61 `sourcedata/bndu/*/*.mat` (the NEMAR copies, which differ
  only in the header date and `FileName`/`FilePath`, see `code/b3w12_convert_onoff.py`). They were fetched from the NEMAR S3
  remote by annex key; sha256 = key for 61 of 61.
- Paper: Wiest et al. 2023, eLife 12:e82467 (PMC10005762), Methods, and its Source code 1 (`elife-82467-code1-v2.zip`,
  sha256 `61001c11a86a51b217e597dc183c3d7332ca81db1431f64b158299e2bfdada4d`; also doi:10.5287/bodleian:rJ7jyjX97).
- Converter: `code/b3w12_convert_onoff.py` l. 215-235 writes the release samples unchanged, with unit `n/a` in every `.vhdr`
  channel entry and in `channels.tsv`.

**Source check per channel type**

| type | rows | release .mat | paper | authors' code |
|---|---|---|---|---|
| DBS (STN LFP) | 312 | no unit, scale or gain field (fields below). `MATRIX_DBS.metadata` is a table with columns signal_index, hemisphere_id, medication | names the amplifiers (TMSi Porti at 2048 Hz, SAGA32 at 4096 Hz); no unit, gain or calibration | every spectrum is divided by its mean (`spec = spec./mean(spec)`, Med_OFFs.m l. 71, Med_ONs.m l. 66, Fig3.m l. 71 and 129; Methods: "divided by the mean power from 1 to 90 Hz"). The µV²/Hz label (Fig2.m l. 15) is on the exponent axis. Fig. 2B/3A therefore set no scale |
| EEG | 273 | same | patients' EEG not described | not used |
| EMG | 14 | same | not described | not used |
| ECG | 4 | same | not described (ECG is mentioned only for the rats) | not used |

**Method**: read the unit field of each `[Channel Infos]` entry of the `.vhdr`. Reader: MNE 1.13.2 `read_raw_brainvision`;
`get_data()` is divided by the stored float32 numbers over the first second (non-zero samples) to check the scale MNE
applies. For every release `.mat`, list the fields of its variable (h5py for v7.3, scipy for v5) and flag any field name
that matches unit, scale, gain, calib, offset or resol.

**Result**

| recording | DBS / EEG / EMG / ECG channels | .vhdr unit field | MNE get_data / stored (first second) |
|---|---|---|---|
| sub-G10_ses-medoff_task-stim_acq-left | 2 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G10_ses-medoff_task-stim_acq-right | 3 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G10_ses-medon_task-stim_acq-left | 2 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G10_ses-medon_task-stim_acq-right | 2 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G23_ses-medoff_task-rest_acq-left | 4 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G23_ses-medoff_task-stim_acq-right | 4 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G23_ses-medon_task-rest_acq-left | 3 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G23_ses-medon_task-stim_acq-right | 4 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G24_ses-medoff_task-stim_acq-left | 4 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G24_ses-medoff_task-stim_acq-right | 4 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G24_ses-medon_task-stim_acq-left | 4 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G24_ses-medon_task-stim_acq-right | 4 / 8 / 0 / 0 | n/a | 1.0 |
| sub-G27_ses-medoff_task-unspecified_acq-left | 1 / 0 / 0 / 0 | n/a | 1.0 |
| sub-G27_ses-medon_task-stim_acq-left | 24 / 6 / 0 / 0 | n/a | 1.0 |
| sub-G27_ses-medon_task-stim_acq-right | 24 / 6 / 0 / 0 | n/a | 1.0 |
| sub-G28_ses-medoff_task-rest_acq-bilateral | 4 / 6 / 2 / 0 | n/a | 1.0 |
| sub-G28_ses-medon_task-stim_acq-left | 8 / 6 / 0 / 0 | n/a | 1.0 |
| sub-G28_ses-medon_task-stim_acq-right | 8 / 6 / 0 / 0 | n/a | 1.0 |
| sub-G30_ses-medoff_task-rest_acq-bilateral | 4 / 6 / 2 / 0 | n/a | 1.0 |
| sub-G30_ses-medon_task-stim_acq-left | 12 / 6 / 0 / 0 | n/a | 1.0 |
| sub-G30_ses-medon_task-stim_acq-right | 12 / 6 / 0 / 0 | n/a | 1.0 |
| sub-G31_ses-medoff_task-rest_acq-bilateral | 4 / 6 / 2 / 0 | n/a | 1.0 |
| sub-G31_ses-medon_task-stim_acq-left | 9 / 6 / 0 / 0 | n/a | 1.0 |
| sub-G31_ses-medon_task-stim_acq-right | 9 / 6 / 0 / 0 | n/a | 1.0 |
| sub-G32_ses-medoff_task-rest_acq-bilateral | 4 / 4 / 2 / 0 | n/a | 1.0 |
| sub-G32_ses-medon_task-stim_acq-left | 17 / 6 / 0 / 0 | n/a | 1.0 |
| sub-G32_ses-medon_task-stim_acq-right | 17 / 6 / 0 / 0 | n/a | 1.0 |
| sub-G33_ses-medoff_task-rest_acq-bilateral | 8 / 7 / 0 / 1 | n/a | 1.0 |
| sub-G33_ses-medon_task-rest_acq-bilateral | 8 / 7 / 0 / 1 | n/a | 1.0 |
| sub-G34_ses-medoff_task-rest_acq-bilateral | 8 / 7 / 0 / 1 | n/a | 1.0 |
| sub-G34_ses-medon_task-rest_acq-bilateral | 8 / 7 / 0 / 1 | n/a | 1.0 |
| sub-K11_ses-medoff_task-rest_acq-right | 2 / 2 / 2 / 0 | n/a | 1.0 |
| sub-K11_ses-medon_task-stim_acq-right | 4 / 4 / 0 / 0 | n/a | 1.0 |
| sub-K6_ses-medoff_task-rest_acq-bilateral | 8 / 4 / 0 / 0 | n/a | 1.0 |
| sub-K6_ses-medon_task-rest_acq-bilateral | 8 / 4 / 0 / 0 | n/a | 1.0 |
| sub-K7_ses-medoff_task-rest_acq-left | 6 / 6 / 1 / 0 | n/a | 1.0 |
| sub-K7_ses-medon_task-stim_acq-left | 2 / 6 / 1 / 0 | n/a | 1.0 |
| sub-K8_ses-medoff_task-rest_acq-bilateral | 4 / 6 / 2 / 0 | n/a | 1.0 |
| sub-K8_ses-medon_task-stim_acq-left | 8 / 6 / 0 / 0 | n/a | 1.0 |
| sub-K8_ses-medon_task-stim_acq-right | 8 / 6 / 0 / 0 | n/a | 1.0 |
| sub-XG37_ses-medoff_task-rest_acq-left | 4 / 6 / 0 / 0 | n/a | 1.0 |
| sub-XG37_ses-medon_task-stim_acq-left | 13 / 6 / 0 / 0 | n/a | 1.0 |
| sub-XG39_ses-medoff_task-rest_acq-right | 4 / 5 / 0 / 0 | n/a | 1.0 |
| sub-XG39_ses-medon_task-stim_acq-right | 12 / 6 / 0 / 0 | n/a | 1.0 |

| type | recordings | channels.tsv rows | .vhdr unit field |
|---|---|---|---|
| DBS | 44 | 312 | n/a |
| EEG | 43 | 273 | n/a |
| EMG | 8 | 14 | n/a |
| ECG | 4 | 4 | n/a |

Release .mat files read: 61
- 1 x MATRIX_DBS: fs, metadata, signal_base, signal_dbs
- 34 x SmrData: EvData, EvTits, FileName, FilePath, Fs, WvData, WvTits
- 20 x SmrData: FileName, FilePath, Fs, WvData, WvTits
- 6 x SmrData: Fs, WvData, WvTits
Fields matching unit/scale/gain/calib/offset/resol: none

The header (`n/a`), the reader (scale 1) and `channels.tsv` (`n/a`) agree, and none of the 61 release files carries a
unit or calibration. µV is not asserted for any type.

**Coverage**: 44 of 44 recordings (all 603 rows the PR touched) and 61 of 61 release `.mat` files.

Reproduce: `python code/receipts/2026-10-22_units-evidence.py <nm000446 root with sourcedata/bndu>` (Python 3.12.15,
MNE 1.13.2, numpy 2.5.4, scipy 1.18.1, pandas 3.1.0, h5py 3.16.0).
