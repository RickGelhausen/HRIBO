"""Cross-condition reports keep detection, differential calls, and NA distinct."""

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import generate_diffex_summary as summary


def test_detection_is_replicate_aware_and_zero_depth_is_uncertain():
    labels = ["RNA-WT-1", "RNA-WT-2", "RNA-WT-3"]
    totals = dict.fromkeys(labels, 1_000_000)
    values = dict(zip(labels, [12, 14, 0], strict=True))

    call = summary.detection_call(values, labels, totals, 1, 10, 2)
    assert call["state"] == "detected"
    assert call["passing_replicates"] == 2

    values["RNA-WT-2"] = 0
    assert summary.detection_call(values, labels, totals, 1, 10, 2)["state"] == "uncertain"

    values["RNA-WT-1"] = 0
    assert summary.detection_call(values, labels, totals, 1, 10, 2)["state"] == "not_detected"

    totals["RNA-WT-2"] = 0
    totals["RNA-WT-3"] = 0
    assert summary.detection_call(values, labels, totals, 1, 10, 2)["state"] == "uncertain"


def test_effect_call_requires_both_thresholds_and_distinguishes_missing():
    result = {"RNA_log2FC": "2", "RNA_pvalue_adjusted": "0.01"}
    assert summary.effect_call(result, "RNA", 0.05, 1)["state"] == "up"
    result["RNA_log2FC"] = "-2"
    assert summary.effect_call(result, "RNA", 0.05, 1)["state"] == "down"
    result["RNA_log2FC"] = "0.5"
    assert summary.effect_call(result, "RNA", 0.05, 1)["state"] == "not_significant"
    result["RNA_pvalue_adjusted"] = "NA"
    assert summary.effect_call(result, "RNA", 0.05, 1)["state"] == "not_tested"
    assert summary.effect_call(None, "TE", 0.05, 1)["state"] == "not_tested"


def test_count_reader_rejects_duplicate_features(tmp_path):
    path = tmp_path / "counts.csv"
    path.write_text("Identifier,RNA-WT-1\nchr:1-3:+,10\nchr:1-3:+,11\n")
    with pytest.raises(ValueError, match="duplicate feature ID"):
        summary.read_counts(path)


def test_count_reader_ignores_other_mapped_library_types(tmp_path):
    path = tmp_path / "counts.csv"
    path.write_text(
        "Identifier,RNA-WT-1,RIBO-WT-1,TIS-WT-1,RNATIS-WT-1\n"
        "chr:1-3:+,10,20,30,40\n"
    )

    counts, labels, totals = summary.read_counts(path)

    assert labels == ["RNA-WT-1", "RIBO-WT-1"]
    assert counts["chr:1-3:+"] == {"RNA-WT-1": 10, "RIBO-WT-1": 20}
    assert totals == {"RNA-WT-1": 10, "RIBO-WT-1": 20}


def test_explicit_locus_tag_beats_gene_id_fallback_in_either_row_order(tmp_path):
    gene = "chr\tHRIBO\tgene\t1\t90\t.\t+\t.\tID=gene;locus_tag=parent_tag\n"
    cds = "chr\tHRIBO\tCDS\t1\t90\t.\t+\t.\tID=cds;gene_id=fallback_id\n"

    for index, records in enumerate(((gene, cds), (cds, gene))):
        annotation = tmp_path / f"annotation-{index}.gff"
        annotation.write_text("##gff-version 3\n" + "".join(records))
        coordinates = summary.read_annotation(annotation)
        assert summary._metadata("chr:1-90:+", coordinates)["locus_tag"] == (
            "parent_tag"
        )


def _write_fixture(tmp_path):
    count_path = tmp_path / "counts.csv"
    samples = [
        f"{assay}-{condition}-{replicate}"
        for assay in ("RNA", "RIBO")
        for condition in ("WT", "Mut")
        for replicate in (1, 2)
    ]
    ids = ["chr:1-90:+", "chr:101-190:-", "chr:201-290:+"]
    data = [
        [20, 20, 0, 0, 20, 20, 0, 0],
        [0, 0, 20, 20, 0, 0, 20, 20],
        [20, 0, 20, 0, 20, 0, 20, 0],
    ]
    with count_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Identifier", *samples])
        for gene, values in zip(ids, data, strict=True):
            writer.writerow([gene, *values])

    annotation = tmp_path / "annotation.gff"
    annotation.write_text(
        "#hribo-gff-read-counts-v1\tRNA-WT-1\n"
        "chr\tHRIBO\tCDS\t1\t90\t.\t+\t.\tID=A;Name=geneA;locus_tag=b0001\t10\n"
        "chr\tHRIBO\tCDS\t101\t190\t.\t-\t.\tID=B;Name=geneB\t10\n"
        "chr\tHRIBO\tCDS\t201\t290\t.\t+\t.\tID=C;Name=geneC;gene_id=gidC\t10\n"
    )
    deltate = tmp_path / "deltate.csv"
    with deltate.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "gene_id", "contrast", "RNA_log2FC", "RNA_pvalue_adjusted",
            "RIBO_log2FC", "RIBO_pvalue_adjusted", "TE_log2FC",
            "TE_pvalue_adjusted",
        ])
        writer.writeheader()
        writer.writerow({
            "gene_id": ids[0], "contrast": "deltate_Mut-WT",
            "RNA_log2FC": -2, "RNA_pvalue_adjusted": 0.01,
            "RIBO_log2FC": -3, "RIBO_pvalue_adjusted": 0.01,
            "TE_log2FC": -1, "TE_pvalue_adjusted": 0.2,
        })
        writer.writerow({
            "gene_id": ids[1], "contrast": "deltate_Mut-WT",
            "RNA_log2FC": 2, "RNA_pvalue_adjusted": 0.01,
            "RIBO_log2FC": 3, "RIBO_pvalue_adjusted": 0.01,
            "TE_log2FC": 1.5, "TE_pvalue_adjusted": 0.01,
        })
    xtail = tmp_path / "xtail.csv"
    xtail.write_text(
        "gene_id,contrast,log2FC_TE_final,pvalue_adjusted\n"
        f"{ids[1]},xtail_Mut-WT,1.7,0.02\n"
    )
    riborex = tmp_path / "riborex.csv"
    riborex.write_text(
        "gene_id,contrast,log2FC,pvalue_adjusted\n"
        f"{ids[1]},riborex_Mut-WT,1.6,0.03\n"
    )
    return count_path, annotation, xtail, riborex, deltate, ids


def _read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_end_to_end_report_and_browser_tracks_are_consistent(tmp_path):
    counts, annotation, xtail, riborex, deltate, ids = _write_fixture(tmp_path)
    out = tmp_path / "report"
    summary.main([
        "--counts", str(counts), "--annotation", str(annotation),
        "--xtail", str(xtail), "--riborex", str(riborex),
        "--deltate", str(deltate), "--contrasts", "Mut-WT",
        "--output_dir", str(out),
    ])

    condition_rows = _read_tsv(out / "condition_matrix.tsv")
    assert len(condition_rows) == 12
    assert list(condition_rows[0])[:8] == [
        "feature_id", "locus_tag", "genome", "start", "end", "strand",
        "feature_type", "name",
    ]
    keyed = {(r["feature_id"], r["condition"], r["assay"]): r for r in condition_rows}
    assert keyed[(ids[0], "WT", "RNA")]["state"] == "detected"
    assert keyed[(ids[0], "WT", "RNA")]["locus_tag"] == "b0001"
    assert keyed[(ids[0], "WT", "RNA")]["name"] == "geneA"
    assert keyed[(ids[1], "WT", "RNA")]["locus_tag"] == ""
    assert keyed[(ids[2], "WT", "RNA")]["locus_tag"] == "gidC"
    assert keyed[(ids[0], "Mut", "RNA")]["state"] == "not_detected"
    assert keyed[(ids[2], "WT", "RNA")]["state"] == "uncertain"
    assert "RNA-WT-1=20" in keyed[(ids[0], "WT", "RNA")]["sample_counts"]

    contrast_rows = _read_tsv(out / "contrast_matrix.tsv")
    assert len(contrast_rows) == 9
    assert list(contrast_rows[0])[:8] == [
        "feature_id", "locus_tag", "genome", "start", "end", "strand",
        "feature_type", "name",
    ]
    effect = {(r["feature_id"], r["assay"]): r for r in contrast_rows}
    assert effect[(ids[0], "RNA")]["state"] == "down"
    assert effect[(ids[0], "RNA")]["locus_tag"] == "b0001"
    assert effect[(ids[0], "TE")]["state"] == "not_significant"
    assert effect[(ids[1], "TE")]["state"] == "up"
    assert effect[(ids[1], "TE")]["xtail_te_log2fc"] == "1.7"
    assert effect[(ids[2], "RNA")]["state"] == "not_tested"

    tracks = _read_tsv(out / "browser_tracks.tsv")
    assert len(tracks) == 10
    for track in tracks:
        assert (out / track["path"]).is_file()
    manifest = json.loads((out / "browser" / "tracks_manifest.json").read_text())
    assert manifest["unmapped_feature_ids"] == []
    assert any(t["kind"] == "condition" and t["feature_count"] == 1 for t in manifest["tracks"])
    page = (out / "condition_overview.html").read_text()
    assert "RNA detection" in page and "TE change" in page
    assert "No directional call" in page
    assert "condition_matrix.tsv" in page
    assert "condition_overview.xlsx" in page
    assert (out / "condition_overview.xlsx").is_file()
    assert '"aliases":["geneA","b0001","A"]' in page
    assert '"id":"chr:1-90:+","locus_tag":"b0001","name":"geneA"' in page
    assert '"id":"chr:101-190:-","locus_tag":"","name":"geneB"' in page
    assert '"riborex":true' in page
    assert "RiboRex:" in page


def test_report_runs_without_optional_riborex(tmp_path):
    counts, annotation, xtail, _, deltate, _ = _write_fixture(tmp_path)
    out = tmp_path / "report without riborex"
    summary.main([
        "--counts", str(counts), "--annotation", str(annotation),
        "--xtail", str(xtail), "--deltate", str(deltate),
        "--contrasts", "Mut-WT", "--output_dir", str(out),
    ])

    contrast_rows = _read_tsv(out / "contrast_matrix.tsv")
    assert contrast_rows
    assert all(row["riborex_te_log2fc"] == "" for row in contrast_rows)
    assert all(row["riborex_te_padj"] == "" for row in contrast_rows)
    page = (out / "condition_overview.html").read_text()
    assert '"riborex":false' in page
    assert "RiboRex:" not in page


def test_output_is_deterministic(tmp_path):
    counts, annotation, xtail, riborex, deltate, _ = _write_fixture(tmp_path)
    arguments = [
        "--counts", str(counts), "--annotation", str(annotation),
        "--xtail", str(xtail), "--riborex", str(riborex),
        "--deltate", str(deltate), "--contrasts", "Mut-WT",
    ]
    first, second = tmp_path / "first", tmp_path / "second"
    summary.main([*arguments, "--output_dir", str(first)])
    summary.main([*arguments, "--output_dir", str(second)])
    for filename in ("condition_matrix.tsv", "contrast_matrix.tsv", "browser_tracks.tsv", "condition_overview.html"):
        assert (first / filename).read_bytes() == (second / filename).read_bytes()


def test_html_sorting_is_state_aware_stable_and_precedes_pagination(tmp_path):
    page_path = tmp_path / "condition_overview.html"
    settings = SimpleNamespace(
        min_cpm=1, min_count=10, min_replicates=2,
        padj_cutoff=0.05, log2fc_cutoff=1,
    )

    summary.write_html(
        page_path, [], {}, ["WT", "Mut"], ["Mut-WT"], settings,
    )
    page = page_path.read_text()

    assert "sortHeader('Locus tag','locus_tag','identity locus-tag')" in page
    assert "sortHeader('Identifier','identifier','identity identifier')" in page
    assert "sortHeader('Feature'" not in page
    assert "columns.append(sortHeader(col,cellSortKey(assay,kind,col)))" in page
    assert "th.setAttribute('aria-sort'" in page
    assert "indicator.className='sort-indicator'" in page
    assert "?'▲':'▼'):'↕'" in page
    assert "sortDirection=sortDirection==='asc'?'desc':'asc'" in page
    assert "const stateRank={condition:" in page
    assert "leftCall.mean_cpm" in page
    assert "Math.abs(leftCall.log2fc)" in page
    assert "sortKey==='locus_tag'||sortKey==='identifier'" in page
    assert "left.row.locus_tag" in page
    assert "left.row.id" in page
    assert "`${row.locus_tag} ${row.id} ${row.name}" in page
    assert "@media (max-width:600px)" in page
    assert "Display label:" in page
    assert "return order||left.index-right.index" in page
    assert page.index("sortRows(chooseRows(") < page.index("selected.slice(")


def test_annotation_names_cannot_break_out_of_json_script(tmp_path):
    counts, annotation, xtail, riborex, deltate, _ = _write_fixture(tmp_path)
    annotation.write_text(annotation.read_text().replace(
        "Name=geneA",
        "Name=__SUPPLEMENTARY_METHODS____RIBOREX_DETAIL__"
        "%3C%2FScRiPt%3E%3Cscript%3Ealert(1)%3C%2Fscript%3E",
    ))
    out = tmp_path / "report"
    summary.main([
        "--counts", str(counts), "--annotation", str(annotation),
        "--xtail", str(xtail), "--riborex", str(riborex),
        "--deltate", str(deltate), "--contrasts", "Mut-WT",
        "--output_dir", str(out),
    ])

    page = (out / "condition_overview.html").read_text()
    assert "</ScRiPt>" not in page
    assert "\\u003c/ScRiPt>" in page
    assert "__SUPPLEMENTARY_METHODS____RIBOREX_DETAIL__\\u003c/ScRiPt>" in page


def test_five_treatments_against_wildtype_share_one_condition_axis():
    treatments = ["A", "B", "C", "D", "E"]
    contrasts = [f"{condition}-WT" for condition in treatments]
    labels = [
        f"{assay}-{condition}-{replicate}"
        for condition in ["WT", *treatments]
        for assay in ("RNA", "RIBO")
        for replicate in (1, 2)
    ]
    feature_id = "chr:1-90:+"
    counts = {feature_id: dict.fromkeys(labels, 20)}
    totals = dict.fromkeys(labels, 100)
    coordinates = {feature_id: {
        "seqid": "chr", "start": 1, "end": 90, "strand": "+",
        "type": "CDS", "name": "geneA",
    }}
    settings = SimpleNamespace(
        min_cpm=1, min_count=10, min_replicates=2,
        padj_cutoff=0.05, log2fc_cutoff=1,
    )

    _, conditions, condition_table, contrast_table = summary.build_summary(
        counts, labels, totals, coordinates, contrasts,
        {"deltate": {}, "xtail": {}, "riborex": {}}, settings,
    )

    assert conditions == ["WT", *treatments]
    assert len(condition_table) == 6 * 2
    assert len(contrast_table) == 5 * 3
    assert all(row["state"] == "not_tested" for row in contrast_table)
