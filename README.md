[![DOI](https://img.shields.io/badge/DOI-10.82901%2Fnemar.nm000446-blue)](https://doi.org/10.82901/nemar.nm000446)

# STN LFPs from awake Parkinson's patients ON and OFF levodopa, with and without DBS (Wiest et al. 2023)

Local field potentials from the subthalamic nucleus (STN), recorded together with scalp EEG and accelerometers, from patients
with Parkinson's disease 3-6 days after DBS surgery while the leads were externalised. Recordings were made OFF and ON
dopaminergic medication, at rest and during stimulation sessions. The data were analysed in Wiest et al. (2023), who showed
that the aperiodic exponent of STN LFPs rises with levodopa and with 130 Hz STN DBS.

## Source
- Wiest C, Torrecillos F, Pogosyan A, Baig F, Pereira E, Morgante F, Ashkan K, Tan H (2022). STN local field potential
  recordings from awake patients with Parkinson's, ON and OFF meds, and during 130 Hz DBS. Oxford University Research
  Archive / MRC BNDU data sharing platform. doi:10.5287/bodleian:mzJ7YwXvo (https://data.mrc.ox.ac.uk/stn-lfp-on-off-and-dbs)
- Wiest C, Torrecillos F, Pogosyan A, Bange M, Muthuraman M, Groppa S, Hulse N, Hasegawa H, Ashkan K, Baig F, Morgante F,
  Pereira EA, Mallet N, Magill PJ, Brown P, Sharott A, Tan H (2023). The aperiodic exponent of subthalamic field potentials
  reflects excitation/inhibition balance in Parkinsonism. eLife 12:e82467. doi:10.7554/eLife.82467
- Analysis code: doi:10.5287/bodleian:rJ7jyjX97
- Funding (paper): Medical Research Council (MC_UU_00003/2, MC_UU_00003/5, MC_UU_00003/6, MR/V00655X/1, MR/P012272/1); National Institute for Health Research (NIHR) Oxford Biomedical Research Centre; Rosetrees Trust; Parkinson's UK (G-0806).

## Licence
CC BY-SA 4.0, as stated on the record page and on the download form. The release zips contain only .mat files (no licence
or README file), so nothing more restrictive applies.

## Ethics (Wiest et al. 2023)
Approved by the Health Research Authority UK and the National Research Ethics Service local Research Ethics Committee
(IRAS: 46576) and by the local ethics committee at the University of Mainz (837.208.17 (11042)). Written informed consent,
and consent to publish, was obtained before surgery.

## Recordings (methods from Wiest et al. 2023)
- Leads: Medtronic 3389 or SenSight, Boston Scientific Vercise DB-2202 or DB-2201, or Abbott Infinity 6172ANS (the paper
  does not say which lead each patient had). Directional contacts were fused to ring contacts. Ground on the non-dominant forearm.
- Amplifiers: TMSi Porti at 2048 Hz or TMSi SAGA32 at 4096 Hz. Exported from Spike2 to MATLAB by the authors. Each
  file's sampling rate is in the sidecar; the model name is inferred from that rate as the paper describes.
- Stimulation: custom in-house neurostimulator, symmetric biphasic constant-current pulses (60 us, negative phase first),
  130 Hz at the middle contacts; a self-adhesive electrode on the back was the stimulation reference.
- The paper analysed 60 s rest segments per medication state (30 hemispheres, 17 patients) and 60 s segments without and with
  DBS (26 hemispheres, 17 patients).

## Contents
- 16 patients (`sub-<release code>`, e.g. sub-G10, sub-K6, sub-XG37), sessions `ses-medoff` and `ses-medon`.
  44 continuous recordings, 18.21 h in total, at 2048 or 4096 Hz. Durations range from 120 s to 96 min.
  - `task-rest` (17 recordings): rest recordings.
  - `task-stim` (26 recordings): sessions in which stimulation was delivered for parts of the recording. The original
    Spike2 file names refer to ERNA (evoked resonant neural activity) and to continuous DBS. Periods without stimulation are
    included. The stimulation amplitude trace (Amp*/Am*) is included where the release has it.
  - `task-unspecified` (1 recording(s): sub-G27 ses-medoff acq-left): single-channel exports without the original file name. The release does
    not say what the patient was doing.
  - `acq-left` / `acq-right` follow the hemisphere in the release file name (`<code>_leSTN_*` / `<code>_riSTN_*`).
    `acq-bilateral`: the release has the same file (byte-identical) under both hemisphere names. These are the rest files
    that hold both leads, and they are converted once. See `source_files` in scans.tsv.
- Every channel of each file is kept with its release label (channels.tsv):
  - DBS: STN lead channels (L../R../Le../Lv../Rv..). The digits are contact numbers. The release does not document
    monopolar vs bipolar montage.
  - EEG: scalp electrodes. ECG and EMG where the release labels them so.
  - MISC: the online-derived copies of channels (prefix F_/F-, S_/S-), the accelerometers (Ac*/Acc*), the stimulation
    amplitude (Amp*/Am*), and channels the release does not document (Ch1/Ch2, Rs*, Ax*, *Fex/*Ffl, Code, Ph/Pho, Snd,
    Virtual, ...). Labels that occur more than once in a file (e.g. 'Virtual') get a suffix _2, _3, ...
  - Channels that are constant in the release are marked `bad`.
- Units are not stated in the release, so channels.tsv and the .vhdr headers list them as n/a, and MNE reads the stored
  numbers with scale 1. No source gives a unit or a calibration for the DBS, EEG, EMG or ECG channels. The release .mat files
  have no unit, scale or gain field. Wiest et al. (2023) name the amplifiers (TMSi Porti, SAGA32) for the LFPs only and do not
  describe the patients' EEG, EMG or ECG. Their spectra (Fig. 2B, 3A) are divided by the mean 1-90 Hz power, so they do not set
  a scale either. The evidence is in code/receipts/2026-10-22_units-evidence.md. Treat the values as numbers on an
  undocumented scale: compare waveforms, spectral shapes and relative amplitudes within a channel, not absolute amplitudes.
  The values are the release's float64 numbers.
  Every sample is exactly representable in float32, so the BrainVision IEEE_FLOAT_32 files hold the values unchanged. This
  was checked per channel by the converter (code/conversion_report.json).
- events.tsv: Spike2 marker times (EvData/EvTits, e.g. 'Mark') where the release has them. They are undocumented keyboard or
  marker events, kept as `spike2_marker`.
- electrodes.tsv: the DBS channels, with x, y, z = n/a. The paper localised the leads with Lead-DBS, but no coordinates are
  released.
- Recordings whose release file name gives one hemisphere while the internal Spike2 file name refers to the other: sub-G10 ses-medon task-stim acq-right. We keep the release naming; see `spike2_name_hemisphere_tag` in code/conversion_report.json.
- participants.tsv: release code, diagnosis and recruiting hospital (from the original site folder name). Age and sex are not
  published per patient, so they are n/a.

## Not converted: 60 s DBS analysis epochs
`sourcedata/bndu/Human-DBS/MATRIX_DBS.mat` (from the 'Human DBS' zip) holds 26 single-channel excerpts per
condition: `signal_base` (no stimulation) and `signal_dbs` (130 Hz DBS), 58-60 s each, at 2048/4096 Hz,
plus a MATLAB table. They are the 60 s segments the authors cut for the DBS analysis, not complete recordings, so they are
kept unchanged in sourcedata and not converted to BIDS recordings.

## Not converted: single-channel files with author-processed values
These release files hold one channel each, with no file name or path. Their values are not plain amplifier samples: they are
float64 numbers on a one-third grid, consistent with averages of fused directional contacts, and they are not representable in
float32. They are kept unchanged in sourcedata (only the header date day set to 01) and are not presented as raw recordings:
- G25_leSTN_OFF.mat: channel(s) L24, 448 s at 2048 Hz (author-processed values (not float32-representable; e.g. averages of fused directional contacts))
- G25_leSTN_ON.mat: channel(s) L24, 402 s at 2048 Hz (author-processed values (not float32-representable; e.g. averages of fused directional contacts))
- G25_riSTN_OFF.mat: channel(s) R24, 184 s at 2048 Hz (author-processed values (not float32-representable; e.g. averages of fused directional contacts))
- G25_riSTN_ON.mat: channel(s) R24, 201 s at 2048 Hz (author-processed values (not float32-representable; e.g. averages of fused directional contacts))
- G27_riSTN_OFF.mat: channel(s) R24, 896 s at 4096 Hz (author-processed values (not float32-representable; e.g. averages of fused directional contacts))

## sourcedata/bndu
All 61 .mat files of the release, in folders named after their zips. Changes, made for privacy only:
- The MAT header 'Created on' date has its day set to 01.
- `SmrData.FileName` (the original Spike2 file name, which held dates and hospital study codes) is replaced by the release
  file name.
- `SmrData.FilePath` (the original folder, which held dates, site folder names and researchers' user names) is blanked.
In MATLAB v7.3 files the edits are made in place, and the rest of the file is byte-identical. MATLAB v5 files were re-saved
(compressed v5) and every variable was checked to be equal after reloading. Each file was checked for the removed strings.
SHA-256 of the original and the edited files: code/conversion_report.json.

## Privacy
The release has no names in the signals. The identifying metadata in the exports (file names and paths with recording dates,
hospital study numbers, site folders and staff user names) were removed as described above. MAT-file creation dates keep
only the year and month.

## Data-quality notes (fix-20261030, Oct 2026)

- sub-K6 `ses-medoff` (rest, 157 s): the DBS and EEG channels carry large slow excursions of the baseline. Per 1-s window,
  the RMS around the file mean is 30-150 stored units (median per DBS channel) against about 8 in `ses-medon`, with 10-s
  averages up to 318 (R2) in the first 60 s and again rising over the last 40 s; the within-second peak-to-peak range
  (50-110) is much less affected, so the excursions are slow drifts or transients rather than broadband noise. The
  release does not annotate them. High-pass filter or exclude these stretches before analysis; they are kept unchanged.
- sub-G24: channels that are constant for the whole release file are already marked `bad` with "constant value 0 in all
  samples of the release": FC3/FC4 (both `acq-left` files), R13 (`acq-right`, both sessions), L13 (`ses-medon` `acq-left`),
  the accelerometer axis AcZ in `ses-medoff` `acq-right`, and the unused stimulator inputs Amp2/Ch2.
