#!/usr/bin/env python3
"""Create geometrically similar round-bell JSON fixtures for scale audio audits."""
import argparse
import copy
import json
from pathlib import Path

SCALES = (0.2, 0.5, 1.0, 2.0)
DIMENSIONS = ("H", "R_mouth", "t_bottom", "t_top")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="bells/best_geometry.json")
    parser.add_argument("--output-dir", default="research/scale_audio_audit/cases")
    args = parser.parse_args()
    source_path = Path(args.input)
    source = json.loads(source_path.read_text(encoding="utf-8"))
    if "geometry" not in source:
        raise SystemExit("input JSON has no 'geometry' object")
    base = source["geometry"]
    missing = [key for key in DIMENSIONS if key not in base]
    if missing:
        raise SystemExit("geometry is missing dimensions: " + ", ".join(missing))
    if abs(float(base.get("ellipticity", 1.0)) - 1.0) > 1e-8:
        raise SystemExit("scale audit fixtures require a circular bell (ellipticity=1)")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for scale in SCALES:
        data = copy.deepcopy(source)
        for key in DIMENSIONS:
            data["geometry"][key] = float(base[key]) * scale
        path = output_dir / f"geometry_{scale:g}.json"
        path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        manifest.append({"scale": scale, "geometry": str(path)})
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Created {len(manifest)} circular-bell scale cases in {output_dir}")
    print("Scale laws for a gravity-similar physical clapper:")
    print("  mass = 5*s^3 kg; COM = 0.46*s m; tip length = 0.52*s m; tip radius = 0.035*s m")
    print("  initial angular velocity = 3.2/sqrt(s) rad/s; contact step and horizon scale as sqrt(s)")

if __name__ == "__main__":
    main()
