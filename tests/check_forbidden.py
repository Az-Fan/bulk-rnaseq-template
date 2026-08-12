from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = {
    r"renv::": "renv namespace",
    r"BiocManager::install": "BiocManager installation",
    r"install\.packages\s*\(": "install.packages",
    r"remotes::install_": "remotes installation",
    r"/usr/(local/)?bin/Rscript": "system Rscript path",
    r"system2\(\s*[\"']Rscript[\"']": "runtime Rscript bypass",
}
SKIP = {"pixi.lock", "check_forbidden.py", "AGENTS.md"}

def main():
    violations = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.name in SKIP or ".pixi" in path.parts:
            continue
        if "internal/migration_artifacts" in path.relative_to(ROOT).as_posix():
            continue
        # Fixed third-party source is audited at the adapter boundary. It is not
        # part of the downstream runtime and is never called for DE/enrichment.
        if "internal/upstream/omics-pipelines" in path.relative_to(ROOT).as_posix():
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for pattern, label in FORBIDDEN.items():
            if re.search(pattern, text):
                violations.append(f"{path.relative_to(ROOT)}: {label}")
    public_orchestrator = (ROOT / "workflow/run.R").read_text(encoding="utf-8")
    if re.search(r"internal[\\/]engines", public_orchestrator):
        violations.append("workflow/run.R: public workflow bypasses the seven modules")
    if violations:
        raise SystemExit("Forbidden environment paths found:\n" + "\n".join(violations))

if __name__ == "__main__":
    main()
