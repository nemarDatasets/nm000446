"""d3-mrc-stn-onoff (Oxford MRC BNDU / ORA doi:10.5287/bodleian:mzJ7YwXvo; Wiest, Torrecillos, Pogosyan, Baig, Pereira,
Morgante, Ashkan, Tan 2022; Wiest et al. 2023 eLife 12:e82467) -> BIDS-iEEG (raw).

Release (Cambium node 43, links 43/0-43/3; probe Jobs ieeg-b3w12-onoff-probe*-*, ieeg-b3w12-onoff-probe3-1008122146):
- 'Human- Meds- OFF', 'ON- part 1', 'ON- part 2' zips: 60 .mat files named <patient>_<le|ri>STN_<OFF|ON>.mat, each one Spike2
  MATLAB export 'SmrData' (WvData channels x samples, float64; WvTits titles; Fs 2048 Hz = TMSi Porti or 4096 Hz = TMSi SAGA32
  per the paper; EvData/EvTits Spike2 marker times in s; FileName/FilePath = original .smr name and the researcher's folder).
  50 files are MATLAB v7.3 (HDF5), 10 are MATLAB v5 (scipy). 11 le/ri pairs are byte-identical (bilateral rest files released
  under both hemisphere names): converted once, acq-bilateral. Every sample is exactly representable in float32 (probe3), so
  BrainVision IEEE_FLOAT_32 keeps the values unchanged; the converter re-checks this per channel (rt_all).
- 'Human DBS' zip: DBS/MATRIX_DBS.mat = struct with 26 single-channel 60 s excerpts (signal_base / signal_dbs, 2048/4096 Hz) cut
  by the authors for the DBS analysis, plus a MATLAB table. These are 60 s analysis epochs, not complete recordings: kept in
  sourcedata only (day of the MAT header date set to 01), not converted.
sourcedata/bndu: all 61 release .mat files. Edits: MAT header 'Created on' day -> 01; SmrData.FileName -> '<release name>.smr',
SmrData.FilePath -> blank (they held dates, hospital study codes, site folder names and user names); v7.3 edits are in place
(same length, rest of the file byte-identical), v5 files are re-written with scipy (compressed v5, every variable verified
equal on reload). The binary of every sourcedata file is searched for the removed strings (UTF-16LE and ASCII): must be absent.
Usage: python b3w12_convert_onoff.py <work_dir> <bids_root>
"""
import glob, hashlib, json, os, re, shutil, sys
import h5py
import numpy as np
import scipy.io as sio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b3w3_common import wjson, write_vhdr, wtsv, scrub_mat_header, sha256_file, add_electrodes_na  # noqa: E402

W, OUT = sys.argv[1], sys.argv[2]
EX = os.path.join(W, "extracted")
SRC = sorted(glob.glob(os.path.join(EX, "*Meds*", "*", "*.mat")))
DBSM = glob.glob(os.path.join(EX, "*Human DBS", "DBS", "MATRIX_DBS.mat"))
assert len(SRC) == 60 and len(DBSM) == 1, (len(SRC), DBSM)
DOI = "10.5287/bodleian:mzJ7YwXvo"
PAGE = "https://data.mrc.ox.ac.uk/stn-lfp-on-off-and-dbs"
EEG = {"Fz", "FCz", "Cz", "Oz", "Pz", "F3", "F4", "FC3", "FC4", "C3", "C4", "CP3", "CP4", "CPz", "P3", "P4", "O1", "O2"}
NAMERX = re.compile(r"^(X?[GK]\d+)_(?:ERNA_)?(le|ri)STN_(OFF|ON)\.mat$")
# static tokens (site folders / cloud-drive names); user names in the paths are checked through the runtime path components
BAD_TOKENS = [b"StGeorges", b"St. George", b"Kings", b"Nexus365", b"OneDrive", b"Recordings_20", b"Users\\"]


def ctype(n):
    if n in EEG:
        return "EEG", "scalp EEG electrode %s (recorded together with the STN LFPs; montage not documented in the release)" % n
    if re.fullmatch(r"(L|R|Le|Ri)v?\d+", n):
        side = "left" if n.startswith("L") else "right"
        return "DBS", ("%s STN DBS-lead channel '%s' as labelled in the release (digits = contact numbers; the release does not "
                       "document monopolar vs bipolar montage; Wiest et al. 2023 analysed bipolar channels from the contacts adjacent to "
                       "the stimulation contact, directional contacts fused to a ring)" % (side, n))
    m = re.fullmatch(r"([FS])[_-](.+)", n)
    if m:
        return "MISC", "online-derived copy of channel %s (prefix '%s' as released; computation not documented in the release)" % (m.group(2), m.group(1))
    if re.fullmatch(r"Ac(c|\d|l|r)?[XYZxyz][LRFH]?", n):
        return "MISC", "accelerometer axis (label '%s' as released; placement and unit not documented)" % n
    if n == "ECG":
        return "ECG", "ECG (label as released)"
    if re.fullmatch(r"EMG\w*", n):
        return "EMG", "EMG (label '%s' as released; muscle not documented)" % n
    if re.fullmatch(r"Amp?\d", n):
        return "MISC", "stimulation amplitude trace (label as released; scale/unit not documented; values 0 when no stimulation)"
    return "MISC", "release channel '%s', not documented in the release" % n


def h5str(ds):
    if "MATLAB_empty" in ds.attrs:
        return ""
    return "".join(chr(int(c)) for c in np.asarray(ds[()]).ravel(order="F"))


def load(fn):
    """-> dict(titles, X (nch x ns float64), fs, ev=[(title, [times])], filename, filepath, fmt)"""
    hdr = open(fn, "rb").read(128)
    if b"MATLAB 7.3" in hdr:
        with h5py.File(fn, "r") as f:
            g = f["SmrData"]
            tits = [h5str(f[r]) for r in np.asarray(g["WvTits"][()]).ravel(order="F")]
            X = np.asarray(g["WvData"][()], dtype=np.float64)
            fs = float(np.asarray(g["Fs"][()]).ravel()[0])
            ev = []
            if "EvData" in g:
                et = [h5str(f[r]) for r in np.asarray(g["EvTits"][()]).ravel(order="F")]
                for t, r in zip(et, np.asarray(g["EvData"][()]).ravel(order="F")):
                    d = f[r]
                    ev.append((t, [] if "MATLAB_empty" in d.attrs else np.asarray(d[()], dtype=float).ravel().tolist()))
            fnm, fp = h5str(g["FileName"]), h5str(g["FilePath"])
        fmt = "v7.3"
    else:
        m = sio.loadmat(fn, squeeze_me=True, struct_as_record=False)
        s = m["SmrData"]
        tits = [str(t) for t in np.atleast_1d(s.WvTits)]
        X = np.asarray(s.WvData, dtype=np.float64)
        fs = float(np.atleast_1d(s.Fs)[0])
        ev = []
        if hasattr(s, "EvData"):
            et = [str(t) for t in np.atleast_1d(s.EvTits)]
            ed = s.EvData
            ed = [ed] if (isinstance(ed, np.ndarray) and ed.dtype != object) or np.isscalar(ed) else list(np.atleast_1d(ed))
            for t, d in zip(et, ed):
                ev.append((t, np.atleast_1d(np.asarray(d, dtype=float)).ravel().tolist()))
        fnm, fp = str(getattr(s, "FileName", "")), str(getattr(s, "FilePath", ""))
        fmt = "v5"
    if X.ndim == 1:
        X = X[None, :]
    if X.shape[0] != len(tits) and X.shape[1] == len(tits):
        X = X.T
    assert X.shape[0] == len(tits), (fn, X.shape, len(tits))
    return dict(titles=tits, X=X, fs=fs, ev=ev, filename=fnm, filepath=fp, fmt=fmt)


def forbidden(b, strings):
    hits = []
    for s in strings:
        if not s:
            continue
        for enc in ("utf-16-le", "ascii"):
            try:
                pat = s.encode(enc)
            except UnicodeEncodeError:
                continue
            if len(pat) >= 6 and pat in b:
                hits.append((s, enc))
    for t in BAD_TOKENS:
        for pat in (t, t.decode().encode("utf-16-le")):
            if pat in b:
                hits.append((t.decode(), "token"))
    return hits


def scrub_v73(src, dst, newname):
    b = open(src, "rb").read()
    nb, ch = scrub_mat_header(b)
    open(dst, "wb").write(nb)
    with h5py.File(dst, "r+") as f:
        g = f["SmrData"]
        for key, val in (("FileName", newname), ("FilePath", "")):
            ds = g[key]
            if "MATLAB_empty" in ds.attrs:
                continue
            n = ds.size
            v = (val[:n]).ljust(n, " ")
            ds[...] = np.array([ord(c) for c in v], dtype=ds.dtype).reshape(ds.shape, order="F")
    return ch


def scrub_v5(src, dst, newname):
    raw = open(src, "rb").read(116)
    m = sio.loadmat(src, squeeze_me=False, struct_as_record=True, mat_dtype=False)
    s = m["SmrData"]
    if "FileName" in s.dtype.names:
        s["FileName"][0, 0] = np.array([newname])
    if "FilePath" in s.dtype.names:
        s["FilePath"][0, 0] = np.array([""], dtype="<U1")
    sio.savemat(dst, {"SmrData": s}, format="5", do_compression=True, long_field_names=True)
    nh, ch = scrub_mat_header(raw + b"\0" * 12)
    hb = nh[:116]
    with open(dst, "r+b") as f:
        f.seek(0)
        f.write(hb)
    # verify every field equal on reload
    a = sio.loadmat(src, squeeze_me=False, struct_as_record=True)["SmrData"]
    c = sio.loadmat(dst, squeeze_me=False, struct_as_record=True)["SmrData"]
    ok = True
    for k in a.dtype.names:
        if k in ("FileName", "FilePath"):
            continue
        x, y = a[k][0, 0], c[k][0, 0]
        if x.dtype == object:
            ok &= x.shape == y.shape and all(np.array_equal(np.asarray(p), np.asarray(q)) for p, q in zip(x.ravel(), y.ravel()))
        else:
            ok &= np.array_equal(x, y) and x.dtype == y.dtype
    return ch, bool(ok)


os.makedirs(OUT, exist_ok=True)
SD = os.path.join(OUT, "sourcedata", "bndu")
rep = {"recordings": [], "sourcedata": [], "excluded": []}
by_sha = {}
for fn in SRC:
    by_sha.setdefault(sha256_file(fn), []).append(fn)
groups = sorted(by_sha.values(), key=lambda g: os.path.basename(g[0]))
parts = {}
scans = {}
for grp in groups:
    names = sorted(os.path.basename(g) for g in grp)
    ms = [NAMERX.match(n) for n in names]
    assert all(ms), names
    pats = {m.group(1) for m in ms}
    meds = {m.group(3) for m in ms}
    assert len(pats) == 1 and len(meds) == 1, names
    pat, med = pats.pop(), meds.pop()
    hemis = sorted({m.group(2) for m in ms})
    acq = "bilateral" if len(hemis) == 2 else {"le": "left", "ri": "right"}[hemis[0]]
    d = load(grp[0])
    X, lab0, fs = d["X"], d["titles"], d["fs"]
    ns = X.shape[1]
    lab, seen = [], {}
    for n in lab0:
        n = n.strip() or "unnamed"
        seen[n] = seen.get(n, 0) + 1
        lab.append(n if seen[n] == 1 else f"{n}_{seen[n]}")
    task = ("rest" if "rest" in d["filename"].lower() else "stim") if d["filename"].strip() else "unspecified"
    if ns / fs < 60:
        rep["excluded"].append(dict(files=names, reason="shorter than 60 s", duration_s=ns / fs))
        continue
    nonexact = [lab0[i] for i in range(len(lab0)) if not np.array_equal(X[i].astype(np.float32).astype(np.float64), X[i])]
    if nonexact:
        rep["excluded"].append(dict(files=names, reason="author-processed values (not float32-representable; e.g. averages of fused directional contacts) - kept in sourcedata only, not relabelled as raw",
                                    channels=lab0, nonexact_channels=nonexact, duration_s=ns / fs, sfreq=fs, fmt=d["fmt"]))
        print("EXCLUDED (processed values)", names, lab0, round(ns / fs, 1), "s", flush=True)
        continue
    sub = "sub-" + pat
    ses = "ses-med" + med.lower()
    stem = f"{sub}_{ses}_task-{task}_acq-{acq}"
    D = os.path.join(OUT, sub, ses, "ieeg")
    os.makedirs(D, exist_ok=True)
    assert not os.path.exists(os.path.join(D, stem + "_ieeg.eeg")), stem
    X32 = np.ascontiguousarray(X.T.astype("<f4"))
    X32.tofile(os.path.join(D, stem + "_ieeg.eeg"))
    del X32
    write_vhdr(os.path.join(D, stem + "_ieeg"), len(lab), fs, lab, ["n/a"] * len(lab),
               comment="b3w12_convert_onoff.py: SmrData.WvData of %s, values unchanged (all exactly representable in float32; unit not stated in the release)" % " = ".join(names))
    back = np.fromfile(os.path.join(D, stem + "_ieeg.eeg"), dtype="<f4").reshape(ns, len(lab))
    exact = [bool(np.array_equal(back[:, i].astype(np.float64), X[i])) for i in range(len(lab))]
    del back
    types = [ctype(n) for n in lab]
    rows = []
    for i, (n, (t, desc)) in enumerate(zip(lab, types)):
        if n != lab0[i]:
            t, desc = ctype(lab0[i])
            types[i] = (t, desc)
            desc = "release label '%s' occurs more than once in this file; renamed %s. " % (lab0[i], n) + desc
        flat = bool(np.all(X[i] == X[i, 0]))
        grpn = ("left STN" if n.startswith("L") else "right STN") if t == "DBS" else "n/a"
        rows.append([n, t, "n/a", "n/a", "n/a", fs, grpn, "bad" if flat else "good",
                     ("constant value %g in all samples of the release" % X[i, 0]) if flat else "n/a", desc])
    wtsv(os.path.join(D, stem + "_channels.tsv"), ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "group", "status", "status_description", "description"], rows)
    ev = []
    for t, times in d["ev"]:
        for x in times:
            if np.isfinite(x) and 0 <= x <= ns / fs:
                ev.append([round(float(x), 6), 0.0, "spike2_marker", t, int(round(x * fs))])
    ev.sort(key=lambda r: r[0])
    wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "trial_type", "marker_channel", "sample"], ev)
    ncnt = {k: sum(1 for tt, _ in types if tt == k) for k in ("DBS", "EEG", "EMG", "ECG", "TRIG", "MISC")}
    has_amp = any(re.fullmatch(r"Amp?\d", n) and np.any(X[i] != 0) for i, n in enumerate(lab))
    stim = task == "stim"
    js = {
        "TaskName": task,
        "TaskDescription": ("Rest recording: patient awake and at rest (Wiest et al. 2023 analysed 60 s of rest from each medication state)." if task == "rest" else
                            "Single-channel export without the original file name; the release does not say what the patient was doing." if task == "unspecified" else
                            "Recording during a stimulation session with the leads externalised; the original Spike2 file name refers to ERNA "
                            "(evoked resonant neural activity) / continuous DBS testing. Stimulation was delivered during parts of the recording; "
                            "periods without stimulation are included. The stimulation amplitude trace is included where the release has it."),
        "SamplingFrequency": fs, "PowerLineFrequency": 50,
        "SoftwareFilters": "n/a",
        "HardwareFilters": "n/a",
        "iEEGReference": "not documented in the release; channel labels as released (Wiest et al. 2023: bipolar STN recordings from the contacts adjacent to the stimulation contact)",
        "iEEGGround": "non-dominant forearm (Wiest et al. 2023)",
        "Manufacturer": "TMS International",
        "ManufacturersModelName": "Porti" if fs == 2048 else ("SAGA32" if fs == 4096 else "n/a"),
        "SoftwareVersions": "Spike2 (CED) export to MATLAB (SmrData struct)",
        "iEEGPlacementScheme": "STN DBS leads, externalised 3-6 days after surgery (Wiest et al. 2023); hemisphere label from the release file name: %s" % acq,
        "RecordingType": "continuous", "RecordingDuration": round(ns / fs, 4),
        "DBSChannelCount": ncnt["DBS"], "EEGChannelCount": ncnt["EEG"], "EMGChannelCount": ncnt["EMG"], "ECGChannelCount": ncnt["ECG"], "TriggerChannelCount": 0,
        "MiscChannelCount": ncnt["MISC"], "ECOGChannelCount": 0, "SEEGChannelCount": 0,
    }
    if task != "unspecified":
        js["ElectricalStimulation"] = bool(stim)
    if stim:
        js["ElectricalStimulationParameters"] = ("Not documented per recording in the release. Wiest et al. 2023: custom-built in-house neurostimulator, symmetric "
                                                 "constant-current biphasic pulses (60 us, negative phase first), 130 Hz DBS at the middle contact(s), current "
                                                 "increased in 0.5 mA steps; a self-adhesive electrode on the back served as stimulation reference."
                                                 + (" A non-zero stimulation amplitude trace (Amp*/Am*) is present in this file." if has_amp else ""))
    wjson(os.path.join(D, stem + "_ieeg.json"), js)
    scans.setdefault((sub, ses), []).append([f"ieeg/{stem}_ieeg.vhdr", "n/a", " ".join(names), d["fmt"]])
    p = parts.setdefault(sub, {"code": pat, "site": "n/a", "files": []})
    fpl = d["filepath"].lower()
    if "george" in fpl:
        p["site"] = "St George's University Hospital, London"
    elif "kings" in fpl:
        p["site"] = "King's College Hospital, London"
    p["files"] += names
    hemi_in_name = re.findall(r"STN(le|ri)|_(L|R)\d", d["filename"])
    rep["recordings"].append(dict(src=names, fmt=d["fmt"], sub=sub, ses=ses, task=task, acq=acq, stem=stem, n_ch=len(lab), labels=lab, release_labels=lab0, sfreq=fs,
                                  n_samples=int(ns), duration_s=ns / fs, n_events=len(ev), types=ncnt, has_stim_amp=bool(has_amp),
                                  exact_float32=dict(zip(lab, exact)), roundtrip_ok=all(exact),
                                  spike2_name_hemisphere_tag=[h[0] or h[1] for h in hemi_in_name]))
    print(stem, len(lab), ns, fs, round(ns / fs, 1), "s ev", len(ev), "rt", all(exact), flush=True)
    del X, d

# scans + participants
for (sub, ses), rows in scans.items():
    wtsv(os.path.join(OUT, sub, ses, f"{sub}_{ses}_scans.tsv"), ["filename", "acq_time", "source_files", "source_format"], rows)
    wjson(os.path.join(OUT, sub, ses, f"{sub}_{ses}_scans.json"), {
        "source_files": {"Description": "Release .mat file(s) this recording was converted from (two names = byte-identical files released under both hemisphere names)"},
        "source_format": {"Description": "MAT-file version of the release file (v7.3 = HDF5, v5)"}})
MISM = []
for r in rep["recordings"]:
    sides = {"left" if t in ("le", "L") else "right" for t in r["spike2_name_hemisphere_tag"]}
    r["spike2_name_hemisphere_mismatch"] = bool(r["acq"] in ("left", "right") and sides and r["acq"] not in sides)
    if r["spike2_name_hemisphere_mismatch"]:
        MISM.append(f"{r['sub']} {r['ses']} task-{r['task']} acq-{r['acq']}")
MISM_TXT = ("- Recordings whose release file name gives one hemisphere while the internal Spike2 file name refers to the other: "
            + "; ".join(MISM) + ". We keep the release naming; see `spike2_name_hemisphere_tag` in code/conversion_report.json.") if MISM else ""
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "release_code", "diagnosis", "recruitment_site", "age", "sex"],
     [[s, parts[s]["code"], "Parkinson's disease", parts[s]["site"], "n/a", "n/a"] for s in sorted(parts)])
wjson(os.path.join(OUT, "participants.json"), {
    "release_code": {"Description": "Patient code used in the release file names (G = St George's, K = King's College Hospital series; X prefix as released, meaning not documented)"},
    "diagnosis": {"Description": "Idiopathic Parkinson's disease (UK PD Society Brain Bank criteria; Wiest et al. 2023)"},
    "recruitment_site": {"Description": "Recruiting hospital, from the site folder name in the original Spike2 export path (path removed from the shared files)"},
    "age": {"Description": "Not published per patient", "Units": "years"},
    "sex": {"Description": "Not published per patient"}})

# sourcedata (all 61 release files, scrubbed)
os.makedirs(SD, exist_ok=True)
for fn in SRC + DBSM:
    zipdir = os.path.basename(os.path.dirname(os.path.dirname(fn))).replace("aperiodic exponent of subthalamic field potentials- ", "")
    zd = re.sub(r"[^A-Za-z0-9]+", "-", zipdir).strip("-")
    dst = os.path.join(SD, zd, os.path.basename(fn))
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    base = os.path.basename(fn)
    rec = dict(file=f"{zd}/{base}", src_sha256=sha256_file(fn))
    if base == "MATRIX_DBS.mat":
        b = open(fn, "rb").read()
        nb, ch = scrub_mat_header(b)
        open(dst, "wb").write(nb)
        with h5py.File(fn, "r") as f:
            M = f["MATRIX_DBS"]
            fsv = np.asarray(M["fs"][()]).ravel().tolist()
            lens = {k: [int(f[r].size) for r in np.asarray(M[k][()]).ravel(order="F")] for k in ("signal_base", "signal_dbs")}
        rec.update(header_changed=ch, edits="MAT header day -> 01", content_identical=True,
                   matrix_dbs=dict(n=len(fsv), fs=fsv, base_s=[l / s for l, s in zip(lens["signal_base"], fsv)], dbs_s=[l / s for l, s in zip(lens["signal_dbs"], fsv)]))
        old = []
    else:
        d = load(fn)
        old = [d["filename"], d["filepath"]] + [x for x in re.split(r"[\\/]", d["filepath"]) if len(x) >= 6]
        newname = base[:-4] + ".smr"
        if d["fmt"] == "v7.3":
            ch = scrub_v73(fn, dst, newname)
            e = load(dst)
            same = bool(np.array_equal(e["X"], d["X"]) and e["titles"] == d["titles"] and e["fs"] == d["fs"] and e["ev"] == d["ev"])
            ok = same
        else:
            ch, ok = scrub_v5(fn, dst, newname)
            e = load(dst)
            same = bool(np.array_equal(e["X"], d["X"]) and e["titles"] == d["titles"] and e["fs"] == d["fs"])
            ok = ok and same
        rec.update(fmt=d["fmt"], header_changed=ch, edits="MAT header day -> 01; SmrData.FileName -> release name; SmrData.FilePath blanked",
                   content_identical=bool(ok), filename_new=e["filename"].strip(), filepath_new=e["filepath"].strip())
        del d, e
    bb = open(dst, "rb").read()
    hits = forbidden(bb, [s for s in old if s and s.strip() and s not in (base[:-4] + ".smr",)])
    if base == "MATRIX_DBS.mat":
        hits = [h for h in hits if h[1] == "token"]
    rec.update(sourcedata_sha256=hashlib.sha256(bb).hexdigest(), residual_identifier_hits=[list(h) for h in hits])
    del bb
    rep["sourcedata"].append(rec)
    print("sourcedata", rec["file"], rec.get("fmt", ""), "ok", rec["content_identical"], "hits", hits, flush=True)
DBSINFO = next(r["matrix_dbs"] for r in rep["sourcedata"] if "matrix_dbs" in r)

# electrodes (no coordinates released)
SPREF = "none: no electrode coordinates or localisation are included in the source release (x, y, z = n/a)"
add_electrodes_na(OUT, "No contact coordinates are released (Wiest et al. 2023 localised leads with Lead-DBS, but the coordinates are not shared). "
                  "Electrode rows are the DBS channels as labelled in the release, x/y/z n/a.",
                  group_fn=lambda n: "left STN" if n.startswith("L") else "right STN", types=("DBS",), spatial_ref=SPREF)

os.makedirs(os.path.join(OUT, "code"), exist_ok=True)
for fn in (os.path.abspath(__file__), os.path.join(os.path.dirname(os.path.abspath(__file__)), "b3w3_common.py")):
    shutil.copy(fn, os.path.join(OUT, "code", os.path.basename(fn)))
json.dump(rep, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1)

R_ = rep["recordings"]
hrs = sum(r["duration_s"] for r in R_) / 3600
nsub = len(parts)
nrest = sum(1 for r in R_ if r["task"] == "rest"); nstim = sum(1 for r in R_ if r["task"] == "stim"); nuns = len(R_) - nrest - nstim
EXC = rep["excluded"]
EXC_TXT = "\n".join("- %s: channel(s) %s, %.0f s at %g Hz (%s)" % (" = ".join(e["files"]), ", ".join(e.get("channels", [])), e["duration_s"], e.get("sfreq", 0), e["reason"].split(" - ")[0]) for e in EXC)
UNS = ", ".join(f"{r['sub']} {r['ses']} acq-{r['acq']}" for r in R_ if r["task"] == "unspecified")
durs = sorted(r["duration_s"] for r in R_)
FUND = ["Medical Research Council (MC_UU_00003/2, MC_UU_00003/5, MC_UU_00003/6, MR/V00655X/1, MR/P012272/1)",
        "National Institute for Health Research (NIHR) Oxford Biomedical Research Centre", "Rosetrees Trust", "Parkinson's UK (G-0806)"]
wjson(os.path.join(OUT, "dataset_description.json"), {
    "Name": "STN local field potentials from awake Parkinson's disease patients ON and OFF medication, with and without DBS (Wiest et al. 2023)",
    "BIDSVersion": "1.10.0", "DatasetType": "raw", "License": "CC-BY-SA-4.0",
    "Authors": ["Christoph Wiest", "Flavie Torrecillos", "Alek Pogosyan", "Fahd Baig", "Erlick A. Pereira", "Francesca Morgante", "Keyoumars Ashkan", "Huiling Tan"],
    "Acknowledgements": "Data collected and shared by the MRC Brain Network Dynamics Unit, University of Oxford (Tan group), with St George's University Hospital and King's College Hospital, London.",
    "HowToAcknowledge": "Cite the dataset (Wiest C, Torrecillos F, Pogosyan A, Baig F, Pereira E, Morgante F, Ashkan K, Tan H. STN local field potential recordings from awake patients with Parkinson's, ON and OFF meds, and during 130 Hz DBS. Oxford University Research Archive, 2022, doi:10.5287/bodleian:mzJ7YwXvo) and Wiest C, et al. (2023) The aperiodic exponent of subthalamic field potentials reflects excitation/inhibition balance in Parkinsonism. eLife 12:e82467, doi:10.7554/eLife.82467.",
    "Funding": FUND,
    "EthicsApprovals": ["Health Research Authority UK / National Research Ethics Service local Research Ethics Committee (IRAS: 46576)", "Local ethics committee, University of Mainz (837.208.17 (11042))"],
    "ReferencesAndLinks": ["https://doi.org/" + DOI, "https://doi.org/10.7554/eLife.82467", PAGE, "https://doi.org/10.5287/bodleian:rJ7jyjX97"],
    "Keywords": ["subthalamic nucleus", "local field potential", "aperiodic exponent", "excitation/inhibition balance", "levodopa", "deep brain stimulation", "Parkinson's disease", "EEG"],
    "SourceDatasets": [{"DOI": "doi:" + DOI, "URL": PAGE, "Version": "MRC BNDU Cambium release (4 zips: Human- Meds- OFF, Human- Meds- ON- part 1/2, Human DBS)"}],
    "GeneratedBy": [{"Name": "b3w12_convert_onoff.py", "Description": "Spike2 MATLAB exports (SmrData) repackaged to BrainVision float32 with values unchanged (bit-exact); Spike2 marker times as events; code in code/."}],
})
open(os.path.join(OUT, "README"), "w").write(f"""# STN LFPs from awake Parkinson's patients ON and OFF levodopa, with and without DBS (Wiest et al. 2023)

Local field potentials from the subthalamic nucleus (STN), recorded together with scalp EEG and accelerometers, from patients
with Parkinson's disease 3-6 days after DBS surgery while the leads were externalised. Recordings were made OFF and ON
dopaminergic medication, at rest and during stimulation sessions. The data were analysed in Wiest et al. (2023), who showed
that the aperiodic exponent of STN LFPs rises with levodopa and with 130 Hz STN DBS.

## Source
- Wiest C, Torrecillos F, Pogosyan A, Baig F, Pereira E, Morgante F, Ashkan K, Tan H (2022). STN local field potential
  recordings from awake patients with Parkinson's, ON and OFF meds, and during 130 Hz DBS. Oxford University Research
  Archive / MRC BNDU data sharing platform. doi:{DOI} ({PAGE})
- Wiest C, Torrecillos F, Pogosyan A, Bange M, Muthuraman M, Groppa S, Hulse N, Hasegawa H, Ashkan K, Baig F, Morgante F,
  Pereira EA, Mallet N, Magill PJ, Brown P, Sharott A, Tan H (2023). The aperiodic exponent of subthalamic field potentials
  reflects excitation/inhibition balance in Parkinsonism. eLife 12:e82467. doi:10.7554/eLife.82467
- Analysis code: doi:10.5287/bodleian:rJ7jyjX97
- Funding (paper): {'; '.join(FUND)}.

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
- {nsub} patients (`sub-<release code>`, e.g. sub-G10, sub-K6, sub-XG37), sessions `ses-medoff` and `ses-medon`.
  {len(R_)} continuous recordings, {hrs:.2f} h in total, at 2048 or 4096 Hz. Durations range from {durs[0]:.0f} s to {durs[-1] / 60:.0f} min.
  - `task-rest` ({nrest} recordings): rest recordings.
  - `task-stim` ({nstim} recordings): sessions in which stimulation was delivered for parts of the recording. The original
    Spike2 file names refer to ERNA (evoked resonant neural activity) and to continuous DBS. Periods without stimulation are
    included. The stimulation amplitude trace (Amp*/Am*) is included where the release has it.
  - `task-unspecified` ({nuns} recording(s): {UNS}): single-channel exports without the original file name. The release does
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
- Units are not stated in the release, so channels.tsv lists them as n/a. The values are the release's float64 numbers.
  Every sample is exactly representable in float32, so the BrainVision IEEE_FLOAT_32 files hold the values unchanged. This
  was checked per channel by the converter (code/conversion_report.json).
- events.tsv: Spike2 marker times (EvData/EvTits, e.g. 'Mark') where the release has them. They are undocumented keyboard or
  marker events, kept as `spike2_marker`.
- electrodes.tsv: the DBS channels, with x, y, z = n/a. The paper localised the leads with Lead-DBS, but no coordinates are
  released.
{MISM_TXT}
- participants.tsv: release code, diagnosis and recruiting hospital (from the original site folder name). Age and sex are not
  published per patient, so they are n/a.

## Not converted: 60 s DBS analysis epochs
`sourcedata/bndu/Human-DBS/MATRIX_DBS.mat` (from the 'Human DBS' zip) holds {DBSINFO['n']} single-channel excerpts per
condition: `signal_base` (no stimulation) and `signal_dbs` (130 Hz DBS), {min(DBSINFO['base_s']):.0f}-{max(DBSINFO['base_s']):.0f} s each, at 2048/4096 Hz,
plus a MATLAB table. They are the 60 s segments the authors cut for the DBS analysis, not complete recordings, so they are
kept unchanged in sourcedata and not converted to BIDS recordings.

## Not converted: single-channel files with author-processed values
These release files hold one channel each, with no file name or path. Their values are not plain amplifier samples: they are
float64 numbers on a one-third grid, consistent with averages of fused directional contacts, and they are not representable in
float32. They are kept unchanged in sourcedata (only the header date day set to 01) and are not presented as raw recordings:
{EXC_TXT}

## sourcedata/bndu
All {len(rep['sourcedata'])} .mat files of the release, in folders named after their zips. Changes, made for privacy only:
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
""")
open(os.path.join(OUT, "CHANGES"), "w").write("1.0.0 2026-10-08\n  - Conversion of MRC BNDU / ORA doi:10.5287/bodleian:mzJ7YwXvo to BIDS-iEEG for NEMAR.\n")
open(os.path.join(OUT, ".bidsignore"), "w").write("code/\n")
wjson(os.path.join(OUT, "task-unspecified_events.json"), {"trial_type": {"Description": "Event type", "Levels": {"spike2_marker": "Spike2 marker (EvData) as released; meaning not documented"}},
                                                          "marker_channel": {"Description": "Spike2 marker channel title (EvTits)"}, "sample": {"Description": "0-based sample index"}})
wjson(os.path.join(OUT, "task-stim_events.json"), {"trial_type": {"Description": "Event type", "Levels": {"spike2_marker": "Spike2 marker (EvData) as released; meaning not documented"}},
                                                   "marker_channel": {"Description": "Spike2 marker channel title (EvTits)"}, "sample": {"Description": "0-based sample index"}})
wjson(os.path.join(OUT, "task-rest_events.json"), {"trial_type": {"Description": "Event type", "Levels": {"spike2_marker": "Spike2 marker (EvData) as released; meaning not documented"}},
                                                   "marker_channel": {"Description": "Spike2 marker channel title (EvTits)"}, "sample": {"Description": "0-based sample index"}})
sd_ok = all(s["content_identical"] for s in rep["sourcedata"]) and not any(s["residual_identifier_hits"] for s in rep["sourcedata"])
print(json.dumps({"recordings": len(R_), "subjects": nsub, "hours": hrs, "fs": sorted({r["sfreq"] for r in R_}), "rest": nrest, "stim": nstim,
                  "min_s": durs[0], "max_s": durs[-1], "excluded": rep["excluded"], "sourcedata_ok": sd_ok,
                  "rt_all": all(r["roundtrip_ok"] for r in R_) and sd_ok}))
print("CONVERT_DONE")
