"""Shared helper for writing a single timestamped log file per update-cycle run."""

from datetime import datetime
from pathlib import Path


def start_log(log_dir, name="update_run"):
    """Create a new run log file (one per run, timestamped) and return its path."""
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    log_file.touch()
    return log_file


def log(log_file, message):
    """Append a timestamped message to the run log and print it."""
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {message}"
    print(line)
    with open(log_file, "a") as f:
        f.write(line + "\n")
