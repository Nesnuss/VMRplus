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

import configparser
import socket
import types
import urllib.error

import pandas as pd
import pytest

from vmrplus import config, hyperlinks, markers, ncbi, paths, pipeline

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


# ── resolve_genome_taxon ───────────────────────────────────────────────────
class TestResolveGenomeTaxon:
    """The taxon a '-gb yes' genome is filed under. VMR rows are not always
    complete: some carry only a Subfamily, and 1612 rows of MSL41 carry no
    Family at all. Those are grouped by the nearest rank they do carry
    instead of being dumped together into 'unclassified/'."""

    NAN = float("nan")

    def test_family_wins(self):
        row = {"Family": "Ackermannviridae", "Subfamily": "Weiservirinae",
               "Class": "Caudoviricetes", "Genus": "Agtrevirus"}
        assert ncbi.resolve_genome_taxon(row) == "Ackermannviridae"

    def test_subfamily_replaces_missing_family(self):
        row = {"Family": self.NAN, "Subfamily": "Weiservirinae",
               "Class": "Caudoviricetes", "Genus": "Unicornvirus"}
        assert ncbi.resolve_genome_taxon(row) == "Weiservirinae"

    def test_walks_up_when_family_and_subfamily_are_missing(self):
        # No family/subfamily, but the record is still a Caudoviricetes.
        row = {"Family": self.NAN, "Subfamily": self.NAN,
               "Class": "Caudoviricetes", "Realm": "Duplodnaviria",
               "Genus": "Actaeavirus"}
        assert ncbi.resolve_genome_taxon(row) == "Caudoviricetes"

    def test_nearest_rank_above_wins_over_a_broader_one(self):
        row = {"Family": self.NAN, "Suborder": "Tevenvirineae",
               "Order": "Caudovirales", "Class": "Caudoviricetes"}
        assert ncbi.resolve_genome_taxon(row) == "Tevenvirineae"

    def test_genus_is_the_last_resort(self):
        # Nothing above the family level either -> fall back below it.
        row = {"Family": self.NAN, "Subfamily": self.NAN, "Class": self.NAN,
               "Genus": "Actaeavirus"}
        assert ncbi.resolve_genome_taxon(row) == "Actaeavirus"

    def test_species_when_genus_is_missing_too(self):
        row = {"Family": self.NAN, "Genus": self.NAN,
               "Species": "Chlorella virus XW01"}
        assert ncbi.resolve_genome_taxon(row) == "Chlorella virus XW01"

    def test_empty_row_is_unclassified(self):
        assert ncbi.resolve_genome_taxon({}) == "unclassified"

    def test_literal_nan_string_is_not_a_name(self):
        # The .xlsx -> CSV -> DataFrame round-trip can surface 'nan' as text.
        row = {"Family": "nan", "Subfamily": "nan", "Class": "Caudoviricetes"}
        assert ncbi.resolve_genome_taxon(row) == "Caudoviricetes"

    def test_surrounding_whitespace_is_stripped(self):
        assert ncbi.resolve_genome_taxon({"Family": "  Ackermannviridae "}) \
            == "Ackermannviridae"

    def test_accepts_a_pandas_row(self):
        row = pd.DataFrame([{"Family": None, "Subfamily": None,
                             "Class": "Caudoviricetes"}]).iloc[0]
        assert ncbi.resolve_genome_taxon(row) == "Caudoviricetes"


# ── _genome_header_tag ─────────────────────────────────────────────────────
class TestGenomeHeaderTag:
    """The tag prepended to every '-gb yes' genome header: the resolved taxon
    paired with the genus, so the saved sequence id carries its taxonomy."""

    NAN = float("nan")

    def test_family_and_genus(self):
        assert ncbi._genome_header_tag(
            "Ackermannviridae", "Agtrevirus") == "Ackermannviridae_Agtrevirus"

    def test_resolved_upper_rank_replaces_unclassified(self):
        # A row with no family is tagged with the rank it does carry.
        assert ncbi._genome_header_tag(
            "Caudoviricetes", "Actaeavirus") == "Caudoviricetes_Actaeavirus"

    def test_genus_only_row_is_not_repeated(self):
        # resolve_genome_taxon() fell back to the genus itself.
        assert ncbi._genome_header_tag(
            "Actaeavirus", "Actaeavirus") == "Actaeavirus"

    def test_missing_genus_falls_back(self):
        assert ncbi._genome_header_tag(
            "Ackermannviridae", self.NAN) == "Ackermannviridae_unclassified"

    def test_unresolved_taxon_with_a_genus(self):
        assert ncbi._genome_header_tag(
            "unclassified", "Actaeavirus") == "unclassified_Actaeavirus"

    def test_nothing_at_all_yields_no_prefix(self):
        assert ncbi._genome_header_tag(self.NAN, self.NAN) == ""
        assert ncbi._genome_header_tag("unclassified", self.NAN) == ""

    def test_literal_nan_string_is_not_a_name(self):
        # The .xlsx -> CSV -> DataFrame round-trip can surface 'nan' as text.
        assert ncbi._genome_header_tag(
            "nan", "Actaeavirus") == "unclassified_Actaeavirus"

    def test_whitespace_is_collapsed(self):
        # A space would truncate the FASTA id at the tag.
        assert ncbi._genome_header_tag(
            "Some Family", "Some Genus") == "Some_Family_Some_Genus"


# ── excel_hyperlink_formula ────────────────────────────────────────────────
class TestExcelHyperlinkFormula:
    """The final table links through '=HYPERLINK()' formulas rather than real
    cell hyperlinks: Excel caps a worksheet at 65,530 of the latter and, once
    past it, opens the workbook in repair mode and discards every link. Four
    linked columns per row reach that cap at ~16k rows."""

    NAN = float("nan")

    def test_text_display_is_quoted(self):
        assert hyperlinks.excel_hyperlink_formula(
            "https://www.ncbi.nlm.nih.gov/nuccore/FJ373894", "FJ373894"
        ) == '=HYPERLINK("https://www.ncbi.nlm.nih.gov/nuccore/FJ373894","FJ373894")'

    def test_numeric_display_stays_a_number(self):
        # Quoting it would left-align the id and break numeric sorting.
        assert hyperlinks.excel_hyperlink_formula(
            "https://ictv.global/id/20240001", 20240001
        ) == '=HYPERLINK("https://ictv.global/id/20240001",20240001)'

    def test_integral_float_loses_its_trailing_zero(self):
        assert hyperlinks.excel_hyperlink_formula("https://x/1", 42.0) \
            == '=HYPERLINK("https://x/1",42)'

    def test_quotes_in_the_display_text_are_doubled(self):
        # A bare quote would terminate the string literal and corrupt the cell.
        assert hyperlinks.excel_hyperlink_formula("https://x/1", 'Seg "L"') \
            == '=HYPERLINK("https://x/1","Seg ""L""")'

    def test_blank_display_falls_back_to_the_url(self):
        assert hyperlinks.excel_hyperlink_formula("https://x/1", self.NAN) \
            == '=HYPERLINK("https://x/1","https://x/1")'

    @pytest.mark.parametrize("url", [None, "", "   ", "nan", float("nan")])
    def test_no_url_means_no_formula(self, url):
        assert hyperlinks.excel_hyperlink_formula(url, "FJ373894") is None

    def test_surrounding_whitespace_is_stripped(self):
        assert hyperlinks.excel_hyperlink_formula("  https://x/1 ", "  A  ") \
            == '=HYPERLINK("https://x/1","A")'

    def test_overlong_formula_is_refused(self):
        # Excel silently corrupts a cell whose formula exceeds 8192 chars.
        url = "https://x/" + "a" * hyperlinks.EXCEL_FORMULA_MAX_LEN
        assert hyperlinks.excel_hyperlink_formula(url, "A") is None

    def test_formula_at_the_limit_is_accepted(self):
        url = "https://x/" + "a" * 8000
        formula = hyperlinks.excel_hyperlink_formula(url, "A")
        assert formula is not None
        assert len(formula) <= hyperlinks.EXCEL_FORMULA_MAX_LEN


# ── TEMPLATE_CONFIG / format_tabajara_section ──────────────────────────────
class TestTemplateConfigTabajaraSections:
    """The template's [tabajara_*] sections are rendered from the defaults
    dicts, so editing TABAJARA_CON_DEFAULTS / TABAJARA_DIS_DEFAULTS cannot
    leave the generated config file advertising stale values."""

    def _parsed(self):
        parser = configparser.ConfigParser(inline_comment_prefixes=(';', '#'))
        parser.read_string(config.render_config_template())
        return parser

    def test_con_section_matches_defaults(self):
        assert dict(self._parsed().items('tabajara_con')) == config.TABAJARA_CON_DEFAULTS

    def test_dis_section_matches_defaults(self):
        assert dict(self._parsed().items('tabajara_dis')) == config.TABAJARA_DIS_DEFAULTS

    def test_extra_default_shows_up_in_the_section(self):
        section = config.format_tabajara_section({'m': 'c', 'zz': '7'})
        assert section.splitlines()[-1] == 'zz = 7'

    def test_known_keys_are_commented_and_ordered(self):
        section = config.format_tabajara_section({'t': '0.5', 'm': 'c'})
        assert section == (
            "; Mode ('c' conservative / 'd' discriminative)\n"
            "m = c\n"
            "; Threshold\n"
            "t = 0.5"
        )

    def test_single_seq_params_are_written_commented_out(self):
        rendered = config.render_config_template()
        for key, value in config.TABAJARA_SINGLE_SEQ_PARAMS.items():
            assert f"\n;{key} = {value}\n" in rendered

    def test_single_seq_params_stay_inactive_when_parsed(self):
        # Listed for reference only: they must not leak into either section.
        parsed = self._parsed()
        for section in ('tabajara_con', 'tabajara_dis'):
            keys = dict(parsed.items(section)).keys()
            assert not keys & config.TABAJARA_SINGLE_SEQ_PARAMS.keys()
