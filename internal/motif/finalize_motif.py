#!/usr/bin/env python3
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

path = Path(sys.argv[1])
manifest = json.loads(path.read_text(encoding="utf-8"))
if manifest.get("status") not in {"extended_complete", "comprehensive_v2_complete"}:
    raise SystemExit("Motif finalization requires an extended_complete or comprehensive_v2_complete run manifest")
manifest["status"] = "complete"
manifest["motif"] = {
    "status": "complete",
    "mode": "counts-only promoter motif",
    "promoter_definition": "-1000 to +100 bp from GENCODE v49 TSS",
    "assembly": "hg38",
    "known_motif_database": "JASPAR 2026 CORE vertebrates",
    "de_novo_method": "STREME",
    "known_motif_method": "AME",
    "interpretation": "Motif enrichment does not establish direct binding or causality.",
}
manifest["completed_at"] = datetime.now(timezone.utc).isoformat()
path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
