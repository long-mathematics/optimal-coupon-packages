#!/usr/bin/env python3
"""Reject unresolved references, LaTeX/package warnings, and overfull boxes."""
from pathlib import Path
import re
import sys

log = Path(sys.argv[1]).read_text(errors='replace')
problems = [line for line in log.splitlines() if re.search(
    r'(^!|LaTeX Warning:|Package .* Warning:|Overfull \\[hv]box)', line)]
if problems:
    raise SystemExit('\n'.join(problems))
if not re.search(r'Output written on .*\.pdf', log):
    raise SystemExit('No successful PDF output recorded in the LaTeX log.')
print('PASS: PDF compiled; no LaTeX/package warnings or overfull boxes')
