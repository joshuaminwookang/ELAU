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
import re
import numpy as np
import subprocess
import os
from pathlib import Path

# Common ABC logic optimization commands. Keys are descriptive names while the
# values are the commands passed to ABC.
# ... existing code ...

ABC_OPTIMIZATIONS = [
    '&st', 
    '&synch2', 
    '&dc2', 
    '&syn2', 
    # '&sweep', 
    # '&scorr', 
    '&b',
    '&b -d', 
    '&dch', 
    '&dch -f', 
    '&syn3',
]

# ... rest of code ...


def random_abc_sequence(avg_length=3):
    """Return a random sequence of ABC optimization commands.

    Parameters
    ----------
    length: int
        Number of random commands to choose.
    """    
    k = np.random.geometric(1/avg_length)
    return random.choices(ABC_OPTIMIZATIONS, k=k)

def run_synth(src_file_path: str, num_runs, mean_abc_seq_len=5):
    """Invoke Yosys to synthesize ``top`` from ``sources`` and write ``out_json``.

    If ``top`` is ``None``, use the basename of the first source file.
    """

    if shutil.which("yosys") is None:
        raise FileNotFoundError("yosys binary not found in PATH")
    
    assert not Path(src_file_path).is_dir()
    if not Path(src_file_path).is_file():
        raise FileNotFoundError(f"Source file {src_file_path} not found")

    results: list[dict] = []

    top = Path(src_file_path).stem.split("_")[0]
    tmp_yosys_tcl = Path("./tmp_yosys_script.tcl")
    for _ in range(num_runs):
        abc_cmd = random_abc_sequence(avg_length=mean_abc_seq_len)
        cmd_script = [
            "read_verilog -sv " + f"{str(src_file_path)}",
            # f"synth -top {shlex.quote(top)}" + (" -flatten" if flatten else ""),
            f"synth -top {shlex.quote(top)} -run coarse",
            "opt -fast -full", "memory_map", "opt -full", "techmap", "opt -fast",
        ]
        abc_cmds_str = "strash; &get -n; &fraig -x; &put; scorr; dretime; strash; &get -n; "
        # abc_cmds_str = "strash; &get -n; &fraig -x; "
        abc_cmds_str += "; ".join(abc_cmd) + "; &put; map; print_stats"
        with open(tmp_yosys_tcl, "w") as f:
            f.write(abc_cmds_str)
        cmd_script.append(f"abc -script {str(tmp_yosys_tcl)}")

        stdout_output=""
        # try:
        print(f"abc_cmds_str:\n{abc_cmds_str}")
        result = subprocess.run(
            ["yosys", "-p", "; ".join(cmd_script)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        import pdb; pdb.set_trace()
        # print(f"Yosys command completed successfully: {result.stdout}")
    
        # Now you can access the stdout from the result
        stdout_output = result.stdout
        # Extract area using regex by finding "area =" followed by numbers
        area_match = re.search(r'area\s*=\s*(\d+\.\d+)', stdout_output)
        parsed_area = float(area_match.group(1)) if area_match else None
        
        # Extract delay using regex by finding "delay =" followed by numbers
        delay_match = re.search(r'delay\s*=\s*(\d+\.\d+)', stdout_output)
        parsed_delay = float(delay_match.group(1)) if delay_match else None

        print(f"Parsed area: {parsed_area}")
        print(f"Parsed delay: {parsed_delay}")
            
        # except:
        #     print("Error during yosys run")
        #     parsed_area = -1
        #     parsed_delay = -1
        results.append({
            "abc_cmd_full": abc_cmds_str,
            "abc_cmd": "; ".join(abc_cmd),
            "area": parsed_area,
            "delay": parsed_delay,
            "cmd_script": cmd_script,
            "stdout_output": stdout_output
        })
    return results

def run_all(source, out_dir="./", num_runs_per_design=10, mean_abc_seq_len=5):
    """Invoke Yosys to synthesize ``top`` from ``sources`` and write ``out_json``.

    If ``top`` is ``None``, use the basename of the first source file.
    """

    if shutil.which("yosys") is None:
        raise FileNotFoundError("yosys binary not found in PATH")
    
    assert Path(source).is_dir()
    sources = list(Path(source).glob("*.sv"))
    if len(sources) == 0:
        raise FileNotFoundError("No source files found")
    if not Path(out_dir).is_dir():
        Path(out_dir).mkdir(exist_ok=True)
    
    for source in sources:
        top = Path(source).stem
        print(f"Running synthesis for {top}...")
        results = run_synth(
            src_file_path=source,
            num_runs=num_runs_per_design,
            mean_abc_seq_len=mean_abc_seq_len,
        )

        # save as .txt file for each entry in results
        for i, result in enumerate(results):
            with open(f"{out_dir}/{top}_run{i}.txt", "w") as f:
                f.write(str(result))

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
        "--num-runs-per-design",
        "-n",
        help="Number of runs per design",
        default=10,
        type=int,
    )
    p.add_argument(
        "--seed",
        default=0,
        type=int
    )
    args = p.parse_args()


    random.seed(args.seed)
    np.random.seed(args.seed)
    
    os.makedirs(args.out, exist_ok=True)
    run_all(args.source, args.out, args.num_runs_per_design)


if __name__ == '__main__':
    main()
