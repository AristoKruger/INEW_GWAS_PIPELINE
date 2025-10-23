from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

import yaml


class ConfigError(RuntimeError):
    """Raised when the pipeline configuration is invalid."""


@dataclass(frozen=True)
class PipelineConfig:
    data: Dict[str, Any]
    root: Path
    config_path: Path

    @classmethod
    def load(cls, path: Path | str) -> "PipelineConfig":
        config_path = Path(path).expanduser().resolve()
        if not config_path.exists():
            raise ConfigError(f"Config file not found: {config_path}")

        try:
            data = yaml.safe_load(config_path.read_text()) or {}
        except yaml.YAMLError as exc:  # pragma: no cover - pure defensive
            raise ConfigError(f"Failed to parse YAML config {config_path}: {exc}") from exc

        if not isinstance(data, dict):
            raise ConfigError("Top-level YAML structure must be a mapping")

        root_ref = data.get("project_root", ".")
        root = (config_path.parent / root_ref).resolve()
        return cls(data=data, root=root, config_path=config_path)

    def get(self, *keys: str, default: Any = None) -> Any:
        node: Any = self.data
        for key in keys:
            if isinstance(node, dict) and key in node:
                node = node[key]
            else:
                return default
        return node

    def require(self, *keys: str) -> Any:
        value = self.get(*keys)
        if value is None:
            dotted = ".".join(keys)
            raise ConfigError(f"Missing required config key: {dotted}")
        return value

    def has(self, *keys: str) -> bool:
        sentinel = object()
        return self.get(*keys, default=sentinel) is not sentinel

    def resolve_path(self, value: Optional[str], create_parent: bool = False) -> Optional[Path]:
        if value in (None, ""):
            return None
        path = Path(value)
        if not path.is_absolute():
            path = (self.root / path).resolve()
        if create_parent:
            path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def path(self, *keys: str, create_parent: bool = False, default: Optional[str] = None) -> Optional[Path]:
        raw = self.get(*keys, default=default)
        return self.resolve_path(raw, create_parent=create_parent)

    def list_steps(self) -> Iterable[str]:
        steps = self.get("steps", default={})
        if not isinstance(steps, dict):
            raise ConfigError("steps must be a mapping of step name to configuration block")
        return steps.keys()


load_config = PipelineConfig.load
