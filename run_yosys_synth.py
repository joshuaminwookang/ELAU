#!/usr/bin/env python3
"""Helper script to run Yosys synthesis with generic technology mapping."""

import argparse
import shlex
import shutil
import subprocess
from pathlib import Path


def run_synth(sources, top=None, out_json="synth.json", flatten=False):
    """Invoke Yosys to synthesize ``top`` from ``sources`` and write ``out_json``.

    If ``top`` is ``None``, use the basename of the first source file.
    """

    if shutil.which("yosys") is None:
        raise FileNotFoundError("yosys binary not found in PATH")

    if top is None:
        top = Path(sources[0]).stem

    cmd_script = [
        "read_verilog " + " ".join(shlex.quote(str(Path(src))) for src in sources),
        f"synth -top {shlex.quote(top)}" + (" -flatten" if flatten else ""),
        f"write_json {shlex.quote(out_json)}",
    ]

    subprocess.run(["yosys", "-q", "-p", "; ".join(cmd_script)], check=True)


def main():
    p = argparse.ArgumentParser(
        description="Run Yosys synthesis with generic technology mapping",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("sources", nargs="+", help="Verilog source files")
    p.add_argument(
        "--top",
        help="Top module name (default: basename of first source)",
    )
    p.add_argument(
        "--out",
        default="synth.json",
        help="Output JSON file",
    )
    p.add_argument(
        "--flatten",
        action="store_true",
        help="Flatten design before mapping",
    )
    args = p.parse_args()

    run_synth(args.sources, args.top, args.out, args.flatten)


if __name__ == '__main__':
    main()
