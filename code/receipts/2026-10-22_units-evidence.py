import argparse
import collections
import re
import sys
from pathlib import Path

import h5py
import mne
import numpy as np
import pandas as pd
import scipy
import scipy.io

TYPES = ("DBS", "EEG", "EMG", "ECG")
SCALE_FIELD = re.compile(r"unit|scale|gain|calib|offset|resol", re.IGNORECASE)


def header_units(vhdr: Path) -> list[str]:
    units, section = [], None
    for line in vhdr.read_text(encoding="utf-8").splitlines():
        if line.startswith("["):
            section = line.strip()
        elif section == "[Channel Infos]" and line.startswith("Ch") and "=" in line:
            fields = line.split("=", 1)[1].split(",")
            units.append(fields[3] if len(fields) > 3 else "")
    return units


def mat_fields(path: Path) -> tuple[str, list[str]]:
    try:
        with h5py.File(path, "r") as f:
            var = next(k for k in f if not k.startswith("#"))
            return var, sorted(f[var].keys())
    except OSError:
        m = scipy.io.loadmat(path, squeeze_me=True, struct_as_record=False)
        var = next(k for k in m if not k.startswith("__"))
        return var, sorted(m[var]._fieldnames)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path, help="nm000446 dataset root (with sourcedata/bndu)")
    root = parser.parse_args().root

    print("| recording | DBS / EEG / EMG / ECG channels | .vhdr unit field | MNE get_data / stored (first second) |")
    print("|---|---|---|---|")
    header = collections.defaultdict(set)
    n_rows, n_recordings = collections.Counter(), collections.Counter()
    for tsv in sorted(root.glob("sub-*/ses-*/ieeg/*_channels.tsv")):
        stem = str(tsv).removesuffix("_channels.tsv")
        ch = pd.read_csv(tsv, sep="\t")
        picked = ch.index[ch.type.isin(TYPES)]
        units = header_units(Path(f"{stem}_ieeg.vhdr"))
        raw = mne.io.read_raw_brainvision(f"{stem}_ieeg.vhdr", preload=False, verbose="error")
        n = int(raw.info["sfreq"])
        stored = np.fromfile(f"{stem}_ieeg.eeg", dtype="<f4", count=n * len(ch)).reshape(n, len(ch)).T[picked]
        read = raw.get_data(picks=list(picked), stop=n)
        nonzero = stored != 0
        ratio = sorted({round(float(r), 6) for r in read[nonzero] / stored[nonzero]})
        for t in TYPES:
            idx = ch.index[ch.type == t]
            n_rows[t] += len(idx)
            n_recordings[t] += len(idx) > 0
            header[t] |= {units[i] for i in idx}
        print(f"| {Path(stem).name} | {' / '.join(str((ch.type == t).sum()) for t in TYPES)} "
              f"| {', '.join(sorted({units[i] for i in picked}))} | {', '.join(map(str, ratio))} |")

    print("\n| type | recordings | channels.tsv rows | .vhdr unit field |")
    print("|---|---|---|---|")
    for t in TYPES:
        print(f"| {t} | {n_recordings[t]} | {n_rows[t]} | {', '.join(sorted(header[t]))} |")

    mats = sorted(root.glob("sourcedata/bndu/*/*.mat"))
    fields = collections.Counter()
    scale_fields = set()
    for mat in mats:
        var, names = mat_fields(mat)
        fields[f"{var}: {', '.join(names)}"] += 1
        scale_fields |= {f for f in names if SCALE_FIELD.search(f)}
    print(f"\nRelease .mat files read: {len(mats)}")
    for k, v in fields.items():
        print(f"- {v} x {k}")
    print(f"Fields matching unit/scale/gain/calib/offset/resol: {sorted(scale_fields) or 'none'}")
    print(f"\nPython {sys.version.split()[0]}, MNE {mne.__version__}, numpy {np.__version__}, "
          f"scipy {scipy.__version__}, pandas {pd.__version__}, h5py {h5py.__version__}")


if __name__ == "__main__":
    main()
