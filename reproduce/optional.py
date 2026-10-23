"""Resolver for the separately archived FINAL32 fixed-parameter benchmark."""
from pathlib import Path
import json
from .common import ReplicationError

def resolve_inputs(root, fixed_results=None):
    config=json.loads((root/'config/result_inputs.json').read_text())
    if fixed_results is None and config.get('fixed32'):
        fixed_results=root/config['fixed32']
    if fixed_results is not None:
        fixed_results=Path(fixed_results)
        if not fixed_results.is_file():
            raise ReplicationError(f'Missing configured fixed32 input: {fixed_results}')
    return fixed_results
