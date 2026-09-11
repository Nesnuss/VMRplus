"""Tests for the HMMs_repository organisation step (vmrplus.repository).

These touch the disk (via ``tmp_path``) but need no network and no external
binaries: the HMM files are stand-ins, since the step only ever moves them.
"""

import csv

import pytest

from vmrplus import repository

LIBRARY_COLUMNS = ["vHMM_ID", "Taxon", "Genus", "Protein", "Model type"]


def _make_output_dir(tmp_path, rows, create_files=True):
    """Build a miniature VMR+ output dir: hmm_library.csv + hmms/*.hmm."""
    hmms = tmp_path / "hmms"
    hmms.mkdir()

    with open(tmp_path / "hmm_library.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=LIBRARY_COLUMNS, delimiter=";")
        writer.writeheader()
        for hmm_id, model_type in rows:
            writer.writerow({
                "vHMM_ID": hmm_id,
                "Taxon": "Anelloviridae",
                "Genus": "Betatorquevirus",
                "Protein": "capsid",
                "Model type": model_type,
            })
            if create_files:
                (hmms / f"{hmm_id}.hmm").write_text(f"NAME  {hmm_id}\n//\n")
    return tmp_path


def _names(directory):
    return sorted(p.name for p in directory.iterdir())


# ── layout mapping ─────────────────────────────────────────────────────────
class TestModelLayout:
    def test_conservative(self):
        assert repository._model_layout("conservative") == ("conservative", "_C")

    def test_discriminative(self):
        assert repository._model_layout("discriminative") == ("discriminative", "_D")

    @pytest.mark.parametrize("raw", ["  Conservative ", "DISCRIMINATIVE"])
    def test_case_and_whitespace_tolerant(self, raw):
        assert repository._model_layout(raw) is not None

    @pytest.mark.parametrize("raw", ["", "structural", None])
    def test_unrecognised_returns_none(self, raw):
        assert repository._model_layout(raw) is None


# ── the relocation itself ──────────────────────────────────────────────────
class TestBuildHmmsRepository:
    def test_files_are_split_and_suffixed(self, tmp_path):
        _make_output_dir(tmp_path, [
            ("vHMM_1", "conservative"),
            ("vHMM_2", "discriminative"),
            ("vHMM_3", "discriminative"),
        ])
        counts = repository.build_hmms_repository(tmp_path)

        repo = tmp_path / "HMMs_repository"
        assert _names(repo / "conservative") == ["vHMM_1_C.hmm"]
        assert _names(repo / "discriminative") == ["vHMM_2_D.hmm", "vHMM_3_D.hmm"]
        assert (counts["conservative"], counts["discriminative"], counts["total"]) == (1, 2, 3)

    def test_replaces_the_flat_hmms_directory(self, tmp_path):
        _make_output_dir(tmp_path, [("vHMM_1", "conservative")])
        repository.build_hmms_repository(tmp_path)
        # The repository is the only HMM directory left.
        assert not (tmp_path / "hmms").exists()

    def test_file_contents_survive_the_move(self, tmp_path):
        _make_output_dir(tmp_path, [("vHMM_7", "discriminative")])
        repository.build_hmms_repository(tmp_path)
        moved = tmp_path / "HMMs_repository" / "discriminative" / "vHMM_7_D.hmm"
        assert moved.read_text() == "NAME  vHMM_7\n//\n"

    def test_both_subdirectories_exist_even_when_one_is_empty(self, tmp_path):
        _make_output_dir(tmp_path, [("vHMM_1", "conservative")])
        repository.build_hmms_repository(tmp_path)
        repo = tmp_path / "HMMs_repository"
        assert (repo / "conservative").is_dir()
        assert (repo / "discriminative").is_dir()

    def test_copy_mode_leaves_hmms_in_place(self, tmp_path):
        _make_output_dir(tmp_path, [("vHMM_1", "conservative")])
        repository.build_hmms_repository(tmp_path, move=False)
        assert (tmp_path / "hmms" / "vHMM_1.hmm").is_file()
        assert (tmp_path / "HMMs_repository" / "conservative" / "vHMM_1_C.hmm").is_file()

    def test_rerunning_is_idempotent(self, tmp_path):
        _make_output_dir(tmp_path, [("vHMM_1", "conservative"), ("vHMM_2", "discriminative")])
        first = repository.build_hmms_repository(tmp_path)
        second = repository.build_hmms_repository(tmp_path)

        assert second["total"] == first["total"] == 2
        assert second["already_present"] == 2
        assert second["missing"] == 0


# ── degraded inputs ────────────────────────────────────────────────────────
class TestDegradedInputs:
    def test_missing_library_csv_is_a_no_op(self, tmp_path):
        (tmp_path / "hmms").mkdir()
        counts = repository.build_hmms_repository(tmp_path)
        assert counts["total"] == 0
        assert not (tmp_path / "HMMs_repository").exists()

    def test_library_without_the_required_columns_aborts(self, tmp_path):
        (tmp_path / "hmms").mkdir()
        (tmp_path / "hmm_library.csv").write_text("Taxon;Genus\nAnelloviridae;Beta\n")
        counts = repository.build_hmms_repository(tmp_path)
        assert counts["total"] == 0

    def test_unknown_model_type_is_counted_not_placed(self, tmp_path):
        _make_output_dir(tmp_path, [("vHMM_1", "conservative"), ("vHMM_2", "structural")])
        counts = repository.build_hmms_repository(tmp_path)

        assert counts["unknown_type"] == 1
        assert counts["total"] == 1
        # Its file is untouched, so nothing is silently lost.
        assert (tmp_path / "hmms" / "vHMM_2.hmm").is_file()

    def test_listed_but_absent_file_is_reported(self, tmp_path):
        _make_output_dir(tmp_path, [("vHMM_1", "conservative")], create_files=False)
        counts = repository.build_hmms_repository(tmp_path)
        assert counts["missing"] == 1
        assert counts["total"] == 0

    def test_unlisted_leftovers_are_kept_not_deleted(self, tmp_path):
        _make_output_dir(tmp_path, [("vHMM_1", "conservative")])
        (tmp_path / "hmms" / "stray.hmm").write_text("NAME  stray\n//\n")

        repository.build_hmms_repository(tmp_path)

        # hmms/ survives precisely because deleting an unlisted model would
        # lose data the library cannot account for.
        assert (tmp_path / "hmms" / "stray.hmm").is_file()


# ── the tally written to the log ───────────────────────────────────────────
class TestSummaryReport:
    def test_summary_reports_each_subdirectory_and_the_total(self, tmp_path, caplog):
        _make_output_dir(tmp_path, [
            ("vHMM_1", "conservative"),
            ("vHMM_2", "discriminative"),
            ("vHMM_3", "discriminative"),
        ])
        with caplog.at_level("INFO"):
            repository.build_hmms_repository(tmp_path)

        log = caplog.text
        assert "conservative   :       1 model(s)" in log
        assert "discriminative :       2 model(s)" in log
        assert "TOTAL          :       3 model(s)" in log

    def test_summary_warns_about_missing_files(self, tmp_path, caplog):
        _make_output_dir(tmp_path, [("vHMM_1", "conservative")], create_files=False)
        with caplog.at_level("INFO"):
            repository.build_hmms_repository(tmp_path)
        assert "1 file(s) listed in the library were not found" in caplog.text
