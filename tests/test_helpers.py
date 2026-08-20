"""Unit tests for the pure (no network / no disk) helper functions of the
VMR+ pipeline.

Why these tests exist
----------------------
The pipeline as a whole depends on NCBI (network) and external binaries
(BLAST+, MAFFT, tabajara.pl), so it cannot be exercised end-to-end in CI.
The functions covered here are *pure*: given the same input they always
return the same output, with no side effects. That makes them cheap, fast
and deterministic to test -- the ideal first safety net for a legacy script.

Where the helpers live
----------------------
After the modularisation the helpers moved into the importable ``vmrplus``
package; the ``VMRplus.py`` entry point is now just a thin CLI wrapper around
it. They are gathered below under a single ``vmr`` namespace so the
assertions -- unchanged from before the refactor -- keep reading ``vmr.<fn>``.
"""

import socket
import types
import urllib.error

import pandas as pd
import pytest

from vmrplus import config, markers, ncbi, paths, pipeline

vmr = types.SimpleNamespace(
    _safe_path_name=paths._safe_path_name,
    _shorten_accession_filename=paths._shorten_accession_filename,
    build_tabajara_args=config.build_tabajara_args,
    detect_cli_config_conflict=config.detect_cli_config_conflict,
    count_fasta_sequences=markers.count_fasta_sequences,
)


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


# ── count_fasta_sequences ──────────────────────────────────────────────────
class TestCountFastaSequences:
    def test_counts_headers(self, tmp_path):
        f = tmp_path / "m.fasta"
        f.write_text(">a\nAAAA\n>b\nCCCC\n>c\nGGGG\n")
        assert vmr.count_fasta_sequences(f) == 3

    def test_single_sequence(self, tmp_path):
        # The single-sequence case that triggers the special tabajara params.
        f = tmp_path / "m.fasta"
        f.write_text(">only\nMKLV\n")
        assert vmr.count_fasta_sequences(f) == 1

    def test_missing_file_returns_zero(self, tmp_path):
        assert vmr.count_fasta_sequences(tmp_path / "nope.fasta") == 0


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


# ── _as_integer_columns ────────────────────────────────────────────────────
class TestAsIntegerColumns:
    """min_length/max_length are whole numbers, but a single empty cell in the
    terms table floats the column, so they rendered as '500.0' in the CSV and
    the .xlsx."""

    def _frame(self):
        return pd.DataFrame({
            'min_length': [500.0, 800.0, None],
            'max_length': [800.0, 1200.0, 700.0],
            'Protein_ID': ['BAA93606.1', '', 'WCS65251.1'],
        })

    def test_whole_floats_become_integers(self):
        df = pipeline._as_integer_columns(self._frame())
        assert df['min_length'].tolist()[:2] == [500, 800]
        assert str(df['min_length'].dtype) == 'Int64'

    def test_missing_values_stay_missing(self):
        # Nullable Int64, not plain int: an empty cell must not become 0.
        df = pipeline._as_integer_columns(self._frame())
        assert pd.isna(df['min_length'].iloc[2])

    def test_other_columns_untouched(self):
        df = pipeline._as_integer_columns(self._frame())
        assert df['Protein_ID'].tolist() == ['BAA93606.1', '', 'WCS65251.1']

    def test_absent_column_is_ignored(self):
        df = pipeline._as_integer_columns(pd.DataFrame({'Family': ['Anelloviridae']}))
        assert df['Family'].tolist() == ['Anelloviridae']


# ── _normalize_accession ───────────────────────────────────────────────────
class TestNormalizeAccession:
    """The VMR cell for VMR0013777 is 'PP173676.' -- a trailing period NCBI
    rejects outright, so that genome was never downloaded under -gb."""

    def test_strips_trailing_period(self):
        assert ncbi._normalize_accession("PP173676.") == "PP173676"

    def test_preserves_version_suffix(self):
        # The dot of a real version is followed by digits and must survive.
        assert ncbi._normalize_accession("MK861116.1") == "MK861116.1"
        assert ncbi._normalize_accession("FJ373894.1") == "FJ373894.1"

    def test_strips_surrounding_whitespace_and_quotes(self):
        assert ncbi._normalize_accession('  "AY702913"  ') == "AY702913"

    def test_strips_trailing_separators(self):
        assert ncbi._normalize_accession("MT138234,") == "MT138234"
        assert ncbi._normalize_accession("MT138234;") == "MT138234"

    def test_plain_accession_is_unchanged(self):
        assert ncbi._normalize_accession("MT138234") == "MT138234"


# ── _retry_delay ───────────────────────────────────────────────────────────
class TestRetryDelay:
    def test_grows_exponentially(self):
        # Compare the jitter-free targets: 5, 10, 20, 40 ...
        for attempt in range(4):
            delay = ncbi._retry_delay(attempt, base=5.0, cap=60.0)
            target = 5.0 * (2 ** attempt)
            assert target / 2 <= delay <= target

    def test_respects_the_cap(self):
        assert ncbi._retry_delay(20, base=5.0, cap=60.0) <= 60.0

    def test_is_jittered(self):
        # Identical attempts must not return the same delay every time, or
        # throttled threads retry in lockstep and re-trigger the 429 burst.
        delays = {ncbi._retry_delay(2) for _ in range(50)}
        assert len(delays) > 1


# ── _is_retryable ──────────────────────────────────────────────────────────
class TestIsRetryable:
    """429 and 5xx are transient; any other 4xx means the request itself is
    wrong, so repeating it verbatim only burns the attempt budget."""

    def _http(self, code):
        return urllib.error.HTTPError("u", code, "err", {}, None)

    @pytest.mark.parametrize("code", [429, 500, 502, 503, 504])
    def test_transient_statuses_are_retried(self, code):
        assert ncbi._is_retryable(self._http(code)) is True

    @pytest.mark.parametrize("code", [400, 401, 403, 404])
    def test_client_errors_are_not_retried(self, code):
        assert ncbi._is_retryable(self._http(code)) is False

    def test_non_http_errors_are_retried(self):
        # Timeouts and truncated reads carry no status and are transient.
        assert ncbi._is_retryable(socket.timeout()) is True
        assert ncbi._is_retryable(urllib.error.URLError("boom")) is True
