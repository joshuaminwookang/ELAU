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
import os
import re
from pathlib import Path

# Common ABC logic optimization commands. Keys are descriptive names while the
# values are the commands passed to ABC.
# ... existing code ...

ABC_OPTIMIZATIONS = [
    '&st', 
    '&synch2', 
    '&dc2', 
    '&syn2', 
    '&sweep', 
    '&scorr', 
    '&b',
    '&b -d', 
    '&dch', 
    '&dch -f', 
    '&syn3',
    '&syn4', 
]

# ... rest of code ...


def random_abc_sequence(length=3):
    """Return a random sequence of ABC optimization commands.

    Parameters
    ----------
    length: int
        Number of random commands to choose.
    """

    cmds = random.choices(ABC_OPTIMIZATIONS, k=length)
    return cmds


def run_synth(source, out_dir="./", flatten=False, abc_cmds: list[str]=[]):
    """Invoke Yosys to synthesize ``top`` from ``sources`` and write ``out_json``.

    If ``top`` is ``None``, use the basename of the first source file.
    """

    if shutil.which("yosys") is None:
        raise FileNotFoundError("yosys binary not found in PATH")
    
    if Path(source).is_dir():
        sources = list(Path(source).glob("*.sv"))
    else:
        sources = [source]
    if len(sources) == 0:
        raise FileNotFoundError("No source files found")

    tmp_yosys_tcl = Path("./tmp_yosys_script.tcl")
    for source in sources:
        if not Path(source).is_file():
            raise FileNotFoundError(f"Source file {source} not found")
        top = Path(sources[0]).stem

        cmd_script = [
            "read_verilog -sv " + f"{str(source)}",
            # f"synth -top {shlex.quote(top)}" + (" -flatten" if flatten else ""),
            f"synth -top {shlex.quote(top)} -run coarse",
            "opt -fast -full", "memory_map", "opt -full", "techmap", "opt -fast",
        ]
        if abc_cmds:
            abc_cmds_str = "strash; &get -n; &fraig -x; &put; scorr; dc2; dretime; strash; &get -n; "
            abc_cmds_str += "; ".join(abc_cmds) + "; &nf -v; &put; print_stats"
            with open(tmp_yosys_tcl, "w") as f:
                f.write(abc_cmds_str)
            cmd_script.append(f"abc -script {str(tmp_yosys_tcl)}")
        else:
            cmd_script.append(f"abc -fast")
        out_verilog = Path(out_dir) / f"{top}.synth.v"
        cmd_script.append(f"write_verilog {str(out_verilog)}")

        result = subprocess.run(
            ["yosys", "-p", "; ".join(cmd_script)],
            check=True,
            stdout=subprocess.PIPE,
            text=True
        )
        
        # Now you can access the stdout from the result
        stdout_output = result.stdout
        # Extract area using regex by finding "area =" followed by numbers
        area_match = re.search(r'area\s*=\s*(\d+\.\d+)', stdout_output)
        parsed_area = float(area_match.group(1)) if area_match else None
        
        # Extract delay using regex by finding "delay =" followed by numbers
        delay_match = re.search(r'delay\s*=\s*(\d+\.\d+)', stdout_output)
        parsed_delay = float(delay_match.group(1)) if delay_match else None
        import pdb;  pdb.set_trace()

def main():
    p = argparse.ArgumentParser(
        description="Run Yosys synthesis with generic technology mapping",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--source",
        "-i",
        help="input source dir or path",
        default="src-flattened"
    )
    p.add_argument(
        "--out",
        "-o",
        default="synth-out",
        help="Output file dir",
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
        default=3,
        help="Apply N random ABC optimization steps from a built-in pool",
    )
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    abc_cmds = args.abc_cmds
    if abc_cmds is None and args.random_abc:
        abc_cmds = random_abc_sequence(args.random_abc)

    run_synth(args.source, args.out, args.flatten, abc_cmds)


if __name__ == '__main__':
    main()
