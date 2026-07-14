"""Unit tests for the pure (no network / no disk) helper functions of the
VMR+ pipeline.

Why these tests exist
----------------------
The pipeline as a whole depends on NCBI (network) and external binaries
(BLAST+, MAFFT, tabajara.pl), so it cannot be exercised end-to-end in CI.
The functions covered here are *pure*: given the same input they always
return the same output, with no side effects. That makes them cheap, fast
and deterministic to test -- the ideal first safety net for a legacy script.

Loading the script as a module
------------------------------
The pipeline lives in ``VMR+_1.7.14.py``. Two obstacles:

1. The filename contains ``+`` and ``.``, so it is NOT a valid Python
   identifier and cannot be imported with a normal ``import`` statement.
   We load it by path via ``importlib``.
2. The script calls ``parser.parse_args()`` at module top-level (line ~2406).
   On import that would parse pytest's OWN argv and crash. We temporarily
   swap ``sys.argv`` for a clean list so ``parse_args()`` falls back to its
   defaults. The ``if __name__ == '__main__'`` block does not run because we
   give the module a custom name (``vmr_pipeline``), not ``"__main__"``.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

# Path to the pipeline script (one level up from tests/).
SCRIPT = Path(__file__).resolve().parent.parent / "VMR+_1.7.14.py"


def _load_pipeline_module():
    """Load VMR+_1.7.14.py as a module named 'vmr_pipeline', side-effect free."""
    saved_argv = sys.argv
    sys.argv = ["VMR+_1.7.14.py"]  # clean argv -> parse_args() uses defaults
    try:
        spec = importlib.util.spec_from_file_location("vmr_pipeline", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.argv = saved_argv  # always restore, even if loading fails
    return module


# Loaded once at collection time; reused by every test below.
vmr = _load_pipeline_module()


# ── _safe_path_name ────────────────────────────────────────────────────────
class TestSafePathName:
    def test_replaces_slash(self):
        assert vmr._safe_path_name("Ser/Thr") == "Ser_Thr"

    def test_collapses_internal_whitespace(self):
        # 'RNA pol  II' has two spaces -> collapsed to a single underscore.
        assert vmr._safe_path_name("RNA pol  II") == "RNA_pol_II"

    def test_strips_leading_and_trailing_underscores(self):
        # Unsafe chars at the edges must not leave dangling underscores.
        assert vmr._safe_path_name("/marker/") == "marker"

    def test_replaces_every_unsafe_character(self):
        # Each of \\ : * ? " < > | must become an underscore.
        assert vmr._safe_path_name('a\\b:c*d?e"f<g>h|i') == "a_b_c_d_e_f_g_h_i"


# ── _shorten_accession_filename ────────────────────────────────────────────
class TestShortenAccessionFilename:
    def test_short_name_is_unchanged(self):
        assert vmr._shorten_accession_filename("AY225133") == "AY225133"

    def test_long_name_is_truncated_and_suffixed(self):
        raw = "X" * 300
        result = vmr._shorten_accession_filename(raw, limit=100)
        assert result == "X" * 100 + "_others_genbankID"

    def test_is_deterministic(self):
        raw = "Z" * 250
        assert vmr._shorten_accession_filename(raw) == vmr._shorten_accession_filename(raw)


# ── _num_duplicates ────────────────────────────────────────────────────────
class TestNumDuplicates:
    # DUPLICATES_BY_COUNT = {1: 7, 2: 6, 3: 3}: pad small groups so the total
    # exceeds 5 (1+7=8, 2+6=8, 3+3=6) for MAFFT/tabajara.
    @pytest.mark.parametrize("n,expected", [(1, 7), (2, 6), (3, 3)])
    def test_known_counts(self, n, expected):
        assert vmr._num_duplicates(n) == expected

    @pytest.mark.parametrize("n", [0, 4, 5, 6, 100])
    def test_counts_needing_no_padding_return_zero(self, n):
        assert vmr._num_duplicates(n) == 0


# ── build_tabajara_args ────────────────────────────────────────────────────
class TestBuildTabajaraArgs:
    def test_flattens_dict_to_cli_list(self):
        # Insertion order is preserved (dicts are ordered in Python 3.7+).
        args = vmr.build_tabajara_args({"t": "0.5", "p": "50", "m": "c"})
        assert args == ["-t", "0.5", "-p", "50", "-m", "c"]

    def test_exclude_keys_are_skipped(self):
        args = vmr.build_tabajara_args({"i": "in", "o": "out", "m": "c"},
                                       exclude_keys={"i", "o"})
        assert args == ["-m", "c"]

    def test_empty_dict_yields_empty_list(self):
        assert vmr.build_tabajara_args({}) == []


# ── detect_cli_config_conflict ─────────────────────────────────────────────
class TestDetectCliConfigConflict:
    def test_no_config_flag_is_a_noop(self):
        # When -c/-config was not used, the function must return without exiting.
        assert vmr.detect_cli_config_conflict(["-i", "x", "-o", "y"],
                                              config_flag_used=False) is None

    def test_conflict_aborts_with_systemexit(self):
        # Using -c together with a config-covered flag (-i) must exit(1).
        with pytest.raises(SystemExit) as exc:
            vmr.detect_cli_config_conflict(["-c", "cfg.ini", "-i", "in.xlsx"],
                                           config_flag_used=True)
        assert exc.value.code == 1

    def test_config_flag_without_conflict_is_allowed(self):
        # -c plus only non-covered tokens (e.g. -thread) must NOT exit.
        assert vmr.detect_cli_config_conflict(["-c", "cfg.ini", "-thread", "4"],
                                              config_flag_used=True) is None
