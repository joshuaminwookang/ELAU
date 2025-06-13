#!/usr/bin/env python3
"""Helper script to run Yosys synthesis with generic technology mapping.

Optionally runs ABC with a custom sequence of logic optimizations. The command
string is passed directly to ``abc`` using ``-c``. Refer to the ABC
documentation for available optimizations (e.g. ``-b`` to balance the AIG).
"""

import argparse
import random
import shlex
import shutil
import subprocess
from pathlib import Path

# Common ABC logic optimization commands. Keys are descriptive names while the
# values are the commands passed to ABC.
ABC_OPTIMIZATIONS = {
    "balance": "balance",
    "rewrite": "rewrite",
    "refactor": "refactor",
    "resub": "resub",
    "dc2": "dc2",
}


def random_abc_sequence(length=3):
    """Return a random sequence of ABC optimization commands.

    Parameters
    ----------
    length: int
        Number of random commands to choose.
    """

    cmds = random.choices(list(ABC_OPTIMIZATIONS.values()), k=length)
    return "; ".join(cmds)


def run_synth(sources, top=None, out_json="synth.json", flatten=False, abc_cmds=None):
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
    ]

    if abc_cmds:
        cmd_script.append(f"abc -c {shlex.quote(abc_cmds)}")

    cmd_script.append(f"write_json {shlex.quote(out_json)}")

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
    p.add_argument(
        "--abc-cmds",
        help="ABC optimization command sequence to apply",
    )
    p.add_argument(
        "--random-abc",
        type=int,
        metavar="N",
        help="Apply N random ABC optimization steps from a built-in pool",
    )
    args = p.parse_args()

    abc_cmds = args.abc_cmds
    if abc_cmds is None and args.random_abc:
        abc_cmds = random_abc_sequence(args.random_abc)

    run_synth(args.sources, args.top, args.out, args.flatten, abc_cmds)


if __name__ == '__main__':
    main()
