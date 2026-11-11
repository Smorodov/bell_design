#!/usr/bin/env python3
"""Compare WAV level and modal-energy metadata without external dependencies."""
import argparse
import csv
import math
import struct
from pathlib import Path


def read_wav(path):
    data = path.read_bytes()
    if data[:4] != b"RIFF" or data[8:12] != b"WAVE":
        raise ValueError(f"{path}: not a RIFF/WAVE file")
    pos, fmt, payload = 12, None, None
    while pos + 8 <= len(data):
        tag = data[pos:pos+4]
        size = struct.unpack_from("<I", data, pos+4)[0]
        chunk = data[pos+8:pos+8+size]
        if tag == b"fmt ":
            fmt = struct.unpack_from("<HHIIHH", chunk, 0)
        elif tag == b"data":
            payload = chunk
        pos += 8 + size + (size & 1)
    if fmt is None or payload is None:
        raise ValueError(f"{path}: missing fmt or data chunk")
    code, channels, rate, _, align, bits = fmt
    if code == 3 and bits == 32:
        count = len(payload) // 4
        samples = struct.unpack("<" + "f"*count, payload)
    elif code == 1 and bits == 16:
        count = len(payload) // 2
        samples = tuple(x / (32768.0 if x < 0 else 32767.0)
                        for x in struct.unpack("<" + "h"*count, payload))
    else:
        raise ValueError(f"{path}: unsupported WAV format tag={code}, bits={bits}")
    if not samples or any(not math.isfinite(x) for x in samples):
        raise ValueError(f"{path}: empty or non-finite samples")
    peak = max(abs(x) for x in samples)
    rms = math.sqrt(sum(x*x for x in samples) / len(samples))
    return {"sample_rate": rate, "channels": channels, "format": f"{code}/{bits}",
            "frames": len(samples)//channels, "duration_s": len(samples)//channels/rate,
            "peak": peak, "rms": rms, "sum_squares": sum(x*x for x in samples)}


def metadata(path):
    result = {}
    meta = path.with_name(path.stem + "_metadata.txt")
    if meta.exists():
        for line in meta.read_text(encoding="utf-8", errors="replace").splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                result[key.strip()] = value.strip()
    return result


def mode_energy(path):
    modes = path.with_name(path.stem + "_modes.csv")
    if not modes.exists():
        return None, 0
    total, count = 0.0, 0
    with modes.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            try:
                total += float(row.get("modal_energy_at_separation_j", "0"))
                count += 1
            except (TypeError, ValueError):
                pass
    return total, count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wav", nargs="+", type=Path)
    parser.add_argument("--csv", type=Path, help="also write a CSV summary")
    args = parser.parse_args()
    rows = []
    for path in args.wav:
        audio = read_wav(path)
        meta = metadata(path)
        energy, modes = mode_energy(path)
        row = {"file": str(path), **audio,
               "normalization": meta.get("Output normalization", "unknown"),
               "initial_energy_j": meta.get("Initial clapper energy J", meta.get("Initial striker energy J", "")),
               "bell_energy_j": meta.get("Bell modal energy at separation J", meta.get("Bell modal energy at release J", "")),
               "transfer_fraction": meta.get("Bell energy transfer fraction", ""),
               "max_mode_hz": meta.get("Maximum available modal frequency Hz", meta.get("Maximum returned modal frequency Hz", "")),
               "band_covered": meta.get("Synthesis band covered", meta.get("Synthesis frequency band covered", "")),
               "modal_csv_energy_sum_j": energy if energy is not None else "",
               "modal_csv_rows": modes}
        rows.append(row)
        print(f"{path}: {audio['format']} {audio['sample_rate']} Hz, {audio['channels']} ch, "
              f"peak={audio['peak']:.7g}, RMS={audio['rms']:.7g}, duration={audio['duration_s']:.3f}s, "
              f"transfer={row['transfer_fraction'] or 'n/a'}, max-mode={row['max_mode_hz'] or 'n/a'} Hz")
    if args.csv:
        keys = list(rows[0])
        with args.csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=keys)
            writer.writeheader(); writer.writerows(rows)
        print(f"Wrote {args.csv}")

if __name__ == "__main__":
    main()
