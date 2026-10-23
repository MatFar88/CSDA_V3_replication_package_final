"""Reporting preflight. Standard library only; it never installs anything."""
from __future__ import annotations
import argparse
import importlib
import importlib.metadata
import json
from pathlib import Path
import platform
import re
import struct
import sys
import sysconfig

MODULES = {"numpy": "numpy", "pandas": "pandas", "scipy": "scipy",
           "matplotlib": "matplotlib", "h5py": "h5py"}

def assess(requirements: Path, base_only: bool = False,
           require_venv: bool = False) -> dict:
    result = {"python": platform.python_version(),
              "executable": sys.executable,
              "implementation": platform.python_implementation(),
              "bits": struct.calcsize("P") * 8,
              "in_virtual_environment": sys.prefix != sys.base_prefix,
              "packages": {}, "errors": []}
    if sys.version_info[:2] != (3, 13):
        result["errors"].append("The tested reporting environment requires Python 3.13.x.")
    if result["implementation"] != "CPython" or result["bits"] != 64:
        result["errors"].append("The reporting setup targets 64-bit CPython.")
    if sysconfig.get_config_var("Py_GIL_DISABLED"):
        result["errors"].append("The reference environment is not a free-threaded CPython build.")
    if require_venv and not result["in_virtual_environment"]:
        result["errors"].append("Refusing to install packages outside a virtual environment.")
    if not base_only and not result["errors"]:
        for line in requirements.read_text(encoding="utf-8-sig").splitlines():
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            m = re.fullmatch(r"([\w.-]+)==([\w.+-]+)", line)
            if m is None:
                result["errors"].append("Unsupported requirements entry: " + line)
                continue
            name, expected = m.groups()
            try:
                actual = importlib.metadata.version(name)
                importlib.import_module(MODULES.get(name, name.replace("-", "_")))
                result["packages"][name] = actual
                if actual != expected:
                    result["errors"].append(f"{name}: expected {expected}, found {actual}.")
            except Exception as exc:
                result["errors"].append(f"{name}: {type(exc).__name__}: {exc}")
    result["ready"] = not result["errors"]
    return result

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base-only", action="store_true")
    p.add_argument("--require-venv", action="store_true")
    p.add_argument("--requirements", type=Path,
                   default=Path(__file__).resolve().parents[1] / "requirements.txt")
    a = p.parse_args(argv)
    try:
        result = assess(a.requirements, a.base_only, a.require_venv)
    except Exception as exc:
        result = {"ready": False, "errors": [f"{type(exc).__name__}: {exc}"]}
    print(json.dumps(result, indent=2))
    return 0 if result["ready"] else 1

if __name__ == "__main__":
    sys.exit(main())
