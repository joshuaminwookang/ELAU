import re
from pathlib import Path

SRC = Path('src')
DST = Path('src-flattened')
DST.mkdir(exist_ok=True)

# Manual replacement for lau_pkg contents to avoid package use
LOG2FLOOR_FUNC = """
function automatic integer log2floor;
    input integer n;
    integer m;
    integer p;
    begin
        m = -1;
        p = 1;
        while (p <= n) begin
            m = m + 1;
            p = p * 2;
        end
        log2floor = m;
    end
endfunction
"""

# Gather module/package definitions
module_defs = {}
module_deps = {}
for path in SRC.glob('*.sv'):
    text = path.read_text()
    # modules
    for m in re.finditer(r'(?s)\bmodule\s+([A-Za-z_][A-Za-z0-9_]*)\b(.*?)endmodule', text):
        name = m.group(1)
        body = 'module ' + name + m.group(2) + 'endmodule'
        module_defs[name] = body
    # packages
    for m in re.finditer(r'(?s)\bpackage\s+([A-Za-z_][A-Za-z0-9_]*)\b(.*?)endpackage', text):
        name = m.group(1)
        body = 'package ' + name + m.group(2) + 'endpackage'
        module_defs[name] = body

module_names = list(module_defs.keys())

# Simple comment stripper
_comment_re = re.compile(r'//.*?$|/\*.*?\*/', re.S | re.M)

def strip_comments(s: str) -> str:
    return re.sub(_comment_re, '', s)

# Determine dependencies for each module
def find_deps(body: str):
    body_nc = strip_comments(body)
    deps = set()
    for name in module_names:
        pattern = re.compile(r'\b' + re.escape(name) + r'\b\s*(?:#\s*\(|\()', re.S)
        for m in pattern.finditer(body_nc):
            before = body_nc[max(0, m.start()-20):m.start()]
            if re.search(r'module\s+' + re.escape(name), before):
                continue
            deps.add(name)
            break
    return deps

for name, body in module_defs.items():
    module_deps[name] = find_deps(body)

# Recursive gather
def gather(name, seen):
    content = []
    for dep in module_deps.get(name, []):
        if dep not in seen:
            seen.add(dep)
            content.extend(gather(dep, seen))
            content.append(module_defs[dep])
    return content

for path in SRC.glob('*.sv'):
    text = path.read_text()
    mods = [m.group(1) for m in re.finditer(r'\bmodule\s+([A-Za-z_][A-Za-z0-9_]*)', text)]
    seen = set(mods)
    dep_contents = []
    for mname in mods:
        dep_contents.extend(gather(mname, seen))
    out = [text]
    out.extend(dep_contents)
    base_text = '\n\n'.join(out)

    def apply_common(txt: str) -> str:
        txt = re.sub(r'(?s)package\s+lau_pkg\s*;.*?endpackage', '', txt)
        txt = txt.replace('lau_pkg::log2floor', 'log2floor')
        txt = txt.replace('lau_pkg::FAST', '2')
        txt = txt.replace('lau_pkg::MEDIUM', '1')
        txt = txt.replace('lau_pkg::SLOW', '0')
        txt = txt.replace('lau_pkg::speed_e', 'int')
        return txt

    for val, suff in [(2, 'fast'), (1, 'medium'), (0, 'slow')]:
        txt = apply_common(base_text)
        txt = re.sub(r'(parameter\s+int\s+speed\s*=)\s*\d+', f"\\1 {val}", txt)
        (DST / f"{path.stem}_{suff}{path.suffix}").write_text(
            LOG2FLOOR_FUNC + '\n\n' + txt
        )
