"""
Regenerate the tutorial reference outputs used by ``tests/test_tutorials.py``.

    python tests/update_tutorial_references.py            # all tutorials
    python tests/update_tutorial_references.py choke      # only ids containing "choke"

Review ``git diff tests/tutorial_reference`` afterwards: every changed number
is a changed tutorial result.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _tutorial_golden import (  # noqa: E402
    GUI_FILES,
    SCRIPTS,
    gui_ref_path,
    run_all,
    run_gui,
    run_script,
    sanity_problems,
    script_ref_path,
    tutorial_id,
)


def main(pattern: str = "") -> int:
    scripts = [p for p in SCRIPTS if pattern in tutorial_id(p)]
    guis = [p for p in GUI_FILES if pattern in tutorial_id(p)]
    failed = 0
    with tempfile.TemporaryDirectory() as tmp:
        for path, proc in run_all(scripts, run_script, Path(tmp)).items():
            problems = [f"exit code {proc.returncode}"] if proc.returncode else sanity_problems(proc.stdout)
            if problems:
                print(f"FAILED  {tutorial_id(path)}: {'; '.join(problems)}\n{proc.stderr[-2000:]}")
                failed += 1
                continue
            ref = script_ref_path(path)
            ref.parent.mkdir(parents=True, exist_ok=True)
            changed = not ref.exists() or ref.read_text(encoding="utf-8") != proc.stdout
            ref.write_text(proc.stdout, encoding="utf-8")
            print(f"{'updated' if changed else 'same   '} scripts/{ref.name}")
        for path, proc in run_all(guis, run_gui, Path(tmp)).items():
            if proc.returncode:
                print(f"FAILED  {tutorial_id(path)} (GUI)\n{proc.stderr[-2000:]}")
                failed += 1
                continue
            data = json.loads(proc.stdout)
            text = json.dumps(data, indent=1, sort_keys=True) + "\n"
            ref = gui_ref_path(path)
            ref.parent.mkdir(parents=True, exist_ok=True)
            changed = not ref.exists() or ref.read_text(encoding="utf-8") != text
            ref.write_text(text, encoding="utf-8")
            print(f"{'updated' if changed else 'same   '} gui/{ref.name}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else ""))
