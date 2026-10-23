# Windows x64 bootstrap, V1.1

## Contract

- Explicit SETUP only: downloads require a console confirmation or `Confirm=true`.
- Reuse a healthy package .venv; do not reinstall a working environment.
- Otherwise download the pinned uv 0.12.20 Windows x64 wheel from PyPI's JSON API.
- Validate HTTPS host and SHA-256 from PyPI metadata; extract only uv.exe to a fixed path.
- Request managed CPython 3.13 with `--no-bin --no-registry`.
- Store tool, Python and cache in .runtime; do not modify global PATH or registry.
- Install requirements.txt only into the .venv interpreter and require binary wheels.
- Run a stdlib preflight including import/version checks, then the original CHECK.
- A failure stops setup before REPORT. No setup command calls SELFTEST or FULL.
- Do not automatically delete an existing incompatible or incomplete .venv.

uv is an external setup utility by Astral; its managed CPython comes from
python-build-standalone, not from the python.org Windows installer. Neither
uv nor a Python runtime binary is bundled in this archive. The first setup
needs a permitted internet connection and permission to execute user-local software.
This is not a fully offline or zero-dependency executable distribution.

Pinned scientific Python dependencies are unchanged. The exact 3.13 patch downloaded
is chosen from the pinned uv release's known downloads and recorded by the preflight
and setup log. Previously configured 64-bit CPython 3.13 environments may be reused.
No compatibility with 3.12, 3.14, 32-bit Python or free-threaded builds is asserted.

PowerShell is started with a process-local ExecutionPolicy argument, without
Set-ExecutionPolicy or any persistent security-policy change. Group Policy and
institutional application restrictions remain in force; no TLS verification is disabled.

Documentation used for this implementation:
- https://docs.astral.sh/uv/guides/install-python/
- https://docs.astral.sh/uv/concepts/python-versions/
- https://docs.astral.sh/uv/reference/cli/
- https://docs.astral.sh/uv/reference/environment/
- https://docs.pypi.org/api/json/

Validation scope: Python preflight/unit tests and byte-preservation checks are local.
Windows PowerShell, the uv Windows executable download and the native no-Python
installation branch have not been run in the Linux build environment. A clean
Windows x64 test remains necessary before calling this installer production-tested.

Do not move .venv between machines or locations. Distribute the source ZIP, then
run setup at its destination. Existing code, data, outputs and MATLAB workspaces
are not removed by this setup.
