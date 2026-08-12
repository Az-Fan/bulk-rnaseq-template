from pathlib import Path

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
