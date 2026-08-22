from pathlib import Path

import yaml

from internal.lib.core import import_resource_pack, validate_resources


def test_resource_import_is_checksummed(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "custom_gene_sets.tsv").write_text("set\tgene\nA\tG1\n", encoding="utf-8")
    registry = tmp_path / "resources/registry.yml"
    rows = import_resource_pack(source, tmp_path / "resources/data/test", registry)
    assert len(rows) == 1
    assert len(rows[0]["sha256"]) == 64
    assert len(validate_resources(registry)) == 1


def test_all_network_sync_declarations_pin_sha256():
    source_dir = Path(__file__).resolve().parents[1] / "resources" / "sources"
    for manifest_path in sorted(source_dir.glob("*.yml")):
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
        for resource in manifest.get("resources", []):
            digest = resource.get("sha256", "")
            assert len(digest) == 64, (
                f"Network resource must pin SHA256: {manifest_path}:{resource.get('resource_id')}"
            )
            int(digest, 16)
