"""Local I/O used by the original Sachs numerical implementation."""
from pathlib import Path
import json
ROOT = Path(__file__).resolve().parent
def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False), encoding='utf-8')
