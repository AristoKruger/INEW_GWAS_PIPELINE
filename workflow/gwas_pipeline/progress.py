from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from .utils import ensure_parent


class ProgressLog:
    """Record per-run pipeline progress so resuming is straightforward."""

    def __init__(
        self,
        path: Path,
        *,
        session: str,
        dry_run: bool,
        run_order: List[str],
        started_at: Optional[datetime] = None,
    ) -> None:
        self.path = path
        self.latest_path = self.path.parent / "latest.json"
        self.session = session
        self.dry_run = dry_run
        self.started_at = started_at or datetime.now()
        self._active_starts: Dict[str, datetime] = {}
        self.data: Dict[str, object] = {
            "session": session,
            "dry_run": dry_run,
            "started_at": self._fmt(self.started_at),
            "status": "running",
            "run_order": list(run_order),
            "steps": [],
        }
        self._write()

    @staticmethod
    def _fmt(dt: datetime) -> str:
        return dt.isoformat(timespec="seconds")

    def _write(self) -> None:
        ensure_parent(self.path)
        tmp_path = self.path.parent / f"{self.path.name}.tmp"
        tmp_path.write_text(json.dumps(self.data, indent=2), encoding="utf-8")
        tmp_path.replace(self.path)

        latest_payload: Dict[str, object] = {
            "session": self.session,
            "progress_file": self.path.name,
            "status": self.data["status"],
            "started_at": self.data["started_at"],
        }
        if "finished_at" in self.data:
            latest_payload["finished_at"] = self.data["finished_at"]
        if "last_completed_step" in self.data:
            latest_payload["last_completed_step"] = self.data["last_completed_step"]

        ensure_parent(self.latest_path)
        tmp_latest = self.latest_path.parent / f"{self.latest_path.name}.tmp"
        tmp_latest.write_text(json.dumps(latest_payload, indent=2), encoding="utf-8")
        tmp_latest.replace(self.latest_path)

    def plan_step(self, step: str) -> None:
        stamp = datetime.now()
        self.data["steps"].append(
            {
                "step": step,
                "status": "planned",
                "noted_at": self._fmt(stamp),
            }
        )
        self._write()

    def start_step(self, step: str) -> None:
        started = datetime.now()
        self._active_starts[step] = started
        self.data["steps"].append(
            {
                "step": step,
                "status": "running",
                "started_at": self._fmt(started),
            }
        )
        self._write()

    def _find_step_entry(self, step: str) -> Dict[str, object]:
        for entry in reversed(self.data["steps"]):
            if entry.get("step") == step and entry.get("status") in {"running", "completed", "failed"}:
                return entry
        raise KeyError(f"No step entry recorded for {step!r}")

    def complete_step(
        self,
        step: str,
        *,
        status: str = "completed",
        message: Optional[str] = None,
    ) -> None:
        finished = datetime.now()
        entry = self._find_step_entry(step)
        entry["status"] = status
        entry["finished_at"] = self._fmt(finished)
        start_time = self._active_starts.pop(step, None)
        if start_time is not None:
            entry["duration_seconds"] = round((finished - start_time).total_seconds(), 2)
        if message:
            entry["message"] = message
        if status == "completed":
            self.data["last_completed_step"] = step
        self._write()

    def fail_step(self, step: str, message: Optional[str] = None) -> None:
        self.complete_step(step, status="failed", message=message)

    def finish(self, *, status: str, message: Optional[str] = None) -> None:
        finished = datetime.now()
        self.data["status"] = status
        self.data["finished_at"] = self._fmt(finished)
        if message:
            self.data["message"] = message
        # Fail any steps that were still marked as running.
        while self._active_starts:
            pending_step, started = self._active_starts.popitem()
            pending_entry = self._find_step_entry(pending_step)
            pending_entry["status"] = "failed"
            pending_entry["finished_at"] = self._fmt(finished)
            pending_entry["duration_seconds"] = round((finished - started).total_seconds(), 2)
            if message:
                pending_entry.setdefault("message", message)
        self._write()
