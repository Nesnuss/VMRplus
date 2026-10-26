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

from vmrplus import annotation, config, hyperlinks, markers, ncbi, paths, pipeline, segments

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


# ── segments.py (seg_library.csv, -gb headers) ───────────────────────────
WALLAL_CELL = ("Seg1: KF234259; Seg2: KF234260; Seg3: KF234261; Seg4: KF234262; "
               "Seg5: KF234263; Seg6: KF234264; Seg7: KF234265; Seg8: KF234266; "
               "Seg9: KF234267; Seg10: KF234268")


class TestParseAccessionField:
    def test_labelled_segments(self):
        pairs = segments.parse_accession_field(WALLAL_CELL)
        assert pairs[0] == ("Seg1", "KF234259")
        assert pairs[-1] == ("Seg10", "KF234268")
        assert len(pairs) == 10

    def test_single_code_is_seg1(self):
        assert segments.parse_accession_field("PP467602") == [("Seg1", "PP467602")]

    def test_unlabelled_list_uses_position(self):
        assert segments.parse_accession_field("KT732816; KT732815; KT732817") == [
            ("Seg1", "KT732816"), ("Seg2", "KT732815"), ("Seg3", "KT732817")]

    def test_mixed_labels_use_position_for_the_unlabelled(self):
        assert segments.parse_accession_field("A: EU623082; EU623083; EU623084") == [
            ("A", "EU623082"), ("Seg2", "EU623083"), ("Seg3", "EU623084")]

    def test_label_formats(self):
        assert segments.parse_accession_field(
            "L: OR365538 ; M: OR365537; RNA3:BK064252; DNA-S.1: MF926436") == [
            ("L", "OR365538"), ("M", "OR365537"), ("RNA3", "BK064252"),
            ("DNA-S.1", "MF926436")]

    def test_labels_are_not_taken_for_accessions(self):
        # 'dsRNA1' looks like an accession (letters + digits, 6 chars).
        assert segments.parse_accession_field("dsRNA1: LC350277") == [
            ("dsRNA1", "LC350277")]

    def test_coordinates_and_trailing_dot_are_dropped(self):
        assert segments.parse_accession_field("AE006468 (2844298.2877981)") == [
            ("Seg1", "AE006468")]
        assert segments.parse_accession_field("PP173676.") == [("Seg1", "PP173676")]

    def test_empty_cell(self):
        assert segments.parse_accession_field(float("nan")) == []
        assert segments.parse_accession_field("") == []


def _wallal_table():
    return pd.DataFrame([
        {"Family": None, "Subfamily": None, "Class": "Caudoviricetes",
         "Genus": "Actaeavirus", "Species": "Actaeavirus x",
         "Virus GENBANK accession": "OR420743"},
        {"Family": "Sedoreoviridae", "Subfamily": None, "Genus": "Orbivirus",
         "Species": "Orbivirus betamitchellense",
         "Virus GENBANK accession": WALLAL_CELL},
        {"Family": None, "Subfamily": "Weiservirinae", "Genus": "Unicornvirus",
         "Species": "Unicornvirus y", "Virus GENBANK accession": None},
    ])


class TestBuildSegLibrary:
    def test_one_row_per_code_sharing_the_sample_id(self):
        entries = segments.build_seg_library(_wallal_table())
        wallal = [e for e in entries if e["vSeg_ID"] == "vSeg2"]
        assert len(wallal) == 10
        assert wallal[0] == {
            "n": 2, "vSeg_ID": "vSeg2", "Segment": "Seg1", "Family": "Sedoreoviridae",
            "Genus": "Orbivirus", "Species": "Orbivirus betamitchellense",
            "GenBank": "KF234259"}

    def test_family_falls_back_to_nearest_rank(self):
        entries = segments.build_seg_library(_wallal_table())
        assert entries[0]["Family"] == "Caudoviricetes"
        assert entries[0]["Segment"] == "Seg1"

    def test_rows_without_codes_add_no_lines(self):
        entries = segments.build_seg_library(_wallal_table())
        assert {e["vSeg_ID"] for e in entries} == {"vSeg1", "vSeg2"}


def test_write_seg_library_writes_csv_and_xlsx(tmp_path):
    entries = segments.build_seg_library(_wallal_table())
    csv_path = segments.write_seg_library(entries, tmp_path)
    xlsx = pd.read_excel(tmp_path / "seg_library.xlsx", dtype=str)
    assert list(xlsx.columns) == list(segments.SEG_COLUMNS)
    assert len(xlsx) == len(entries) == len(csv_path.read_text().splitlines()) - 1
    assert xlsx.iloc[1].tolist() == ["vSeg2", "Seg1", "Sedoreoviridae", "Orbivirus",
                                     "Orbivirus betamitchellense", "KF234259"]


class TestSegHeaders:
    ENTRY = {"vSeg_ID": "vSeg1", "Segment": "Seg1", "Family": "Sedoreoviridae",
             "Genus": "Orbivirus", "GenBank": "KF234259"}
    TITLE = "KF234259.1 Wallal virus isolate 927 segment 1, complete sequence"
    EXPECTED = ("Sedoreoviridae_Orbivirus_vSeg1_Seg1_KF234259.1 "
                "Wallal virus isolate 927 segment 1, complete sequence")

    def test_header_keeps_version_and_description(self):
        _, description = segments.seg_header(self.ENTRY, "KF234259.1", self.TITLE)
        assert description == self.EXPECTED

    def test_names_with_spaces_are_underscored(self):
        entry = dict(self.ENTRY, Family="Some Family", Genus="Some Genus")
        new_id, _ = segments.seg_header(entry, "KF234259.1", self.TITLE)
        assert new_id == "Some_Family_Some_Genus_vSeg1_Seg1_KF234259.1"

    def test_rewrite_is_idempotent_and_upgrades_old_headers(self, tmp_path):
        lookup = {"KF234259": self.ENTRY}
        fasta = tmp_path / "VMR0000001_KF234259.fasta"
        fasta.write_text(f">Sedoreoviridae_Orbivirus_{self.TITLE}\nACGT\nAC\n")
        matched, unmatched = segments.rewrite_fasta_headers(fasta, lookup)
        assert matched == {"KF234259"} and unmatched == []
        once = fasta.read_text()
        assert once == f">{self.EXPECTED}\nACGT\nAC\n"
        segments.rewrite_fasta_headers(fasta, lookup)
        assert fasta.read_text() == once

    def test_unknown_record_is_left_alone(self, tmp_path):
        fasta = tmp_path / "x.fasta"
        fasta.write_text(">ZZ999999.1 other\nACGT\n")
        matched, unmatched = segments.rewrite_fasta_headers(
            fasta, {"KF234259": self.ENTRY})
        assert matched == set() and unmatched == ["ZZ999999.1"]
        assert fasta.read_text() == ">ZZ999999.1 other\nACGT\n"


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

    def test_single_seq_params_are_part_of_the_defaults(self):
        for defaults in (config.TABAJARA_CON_DEFAULTS, config.TABAJARA_DIS_DEFAULTS):
            assert config.TABAJARA_SINGLE_SEQ_PARAMS.items() <= defaults.items()

    def test_single_seq_params_are_active_in_both_sections(self):
        parsed = self._parsed()
        for section in ('tabajara_con', 'tabajara_dis'):
            items = dict(parsed.items(section)).items()
            assert config.TABAJARA_SINGLE_SEQ_PARAMS.items() <= items


# ── vmrplus/annotation.py ──────────────────────────────────────────────────

_HMM_HEADER = (
    "HMMER3/f [3.4 | Aug 2023]\n"
    "NAME  vHMM_61\n"
    "LENG  33\n"
    "STATS LOCAL FORWARD   -4.8521  0.71975\n"
    "CUTOFF SCORE\t38.08\n"
    "HMM          A        C\n"
    "  COMPO   2.65995  3.53488\n"
    "//\n"
)

_CON_FIELDS = [("TAXON", "Adamaviridae"), ("TXID", "3152106"),
               ("PROTEIN", "Replication-associated")]

_DIS_FIELDS = [("TAXON", "Poxviridae"), ("TXID", "10240"), ("PROTEIN", "Rna polymerase"),
               None,
               ("TAXON", "Avipoxvirus"), ("TXID", "10260"), ("PROTEIN", "Rna polymerase")]


class TestAnnotateHmmText:
    def test_inserts_after_cutoff_score(self):
        out = annotation.annotate_hmm_text(_HMM_HEADER, _CON_FIELDS)
        assert ("CUTOFF SCORE\t38.08\n"
                "TAXON   Adamaviridae\n"
                "TXID    3152106\n"
                "PROTEIN Replication-associated\n"
                "\n"
                "HMM ") in out

    def test_falls_back_to_stats_local_forward(self):
        text = _HMM_HEADER.replace("CUTOFF SCORE\t38.08\n", "")
        out = annotation.annotate_hmm_text(text, _CON_FIELDS)
        assert "0.71975\nTAXON   Adamaviridae\n" in out

    def test_no_anchor_returns_none(self):
        text = "NAME  x\nHMM  A\n//\n"
        assert annotation.annotate_hmm_text(text, _CON_FIELDS) is None

    def test_rest_of_file_untouched(self):
        out = annotation.annotate_hmm_text(_HMM_HEADER, _CON_FIELDS)
        for tag, _ in _CON_FIELDS:
            out = "".join(l for l in out.splitlines(True) if not l.startswith(tag))
        assert out == _HMM_HEADER.replace("HMM ", "\nHMM ", 1)

    def test_idempotent_and_updates(self):
        once = annotation.annotate_hmm_text(_HMM_HEADER, _CON_FIELDS)
        assert annotation.annotate_hmm_text(once, _CON_FIELDS) == once
        changed = [("TAXON", "X"), ("GENUS", "G"), ("TXID", "1"), ("PROTEIN", "P")]
        out = annotation.annotate_hmm_text(once, changed)
        assert out.count("TAXON") == 1
        assert "TAXON   X\nGENUS   G\nTXID    1\nPROTEIN P\n\nHMM " in out

    def test_desc_mode(self):
        out = annotation.annotate_hmm_text(_HMM_HEADER, _CON_FIELDS, mode="desc")
        assert ("NAME  vHMM_61\n"
                "DESC  TAXON=Adamaviridae;TXID=3152106;PROTEIN=Replication-associated\n"
                ) in out
        assert annotation.annotate_hmm_text(out, _CON_FIELDS, mode="desc") == out

    def test_discriminative_two_blocks(self):
        out = annotation.annotate_hmm_text(_HMM_HEADER, _DIS_FIELDS)
        assert ("CUTOFF SCORE\t38.08\n"
                "TAXON   Poxviridae\n"
                "TXID    10240\n"
                "PROTEIN Rna polymerase\n"
                "\n"
                "TAXON   Avipoxvirus\n"
                "TXID    10260\n"
                "PROTEIN Rna polymerase\n"
                "\n"
                "HMM ") in out
        assert annotation.annotate_hmm_text(out, _DIS_FIELDS) == out

    def test_old_genus_layout_is_rewritten(self):
        old = [("TAXON", "Poxviridae"), ("GENUS", "Avipoxvirus"),
               ("TXID", "10240"), ("PROTEIN", "Rna polymerase")]
        once = annotation.annotate_hmm_text(_HMM_HEADER, old)
        out = annotation.annotate_hmm_text(once, _DIS_FIELDS)
        assert "GENUS" not in out
        assert out == annotation.annotate_hmm_text(_HMM_HEADER, _DIS_FIELDS)

    def test_discriminative_desc_mode(self):
        out = annotation.annotate_hmm_text(_HMM_HEADER, _DIS_FIELDS, mode="desc")
        assert ("DESC  TAXON=Poxviridae;TXID=10240;PROTEIN=Rna polymerase|"
                "TAXON=Avipoxvirus;TXID=10260;PROTEIN=Rna polymerase\n") in out
        assert annotation.annotate_hmm_text(out, _DIS_FIELDS, mode="desc") == out


class TestBuildMetadata:
    TAXIDS = {"Adamaviridae": "3152106"}

    def test_conservative_has_no_genus(self):
        row = {"Taxon": " Adamaviridae ", "Genus": "Family-wide models",
               "Protein": "Rep", "Model type": "conservative"}
        assert annotation.build_metadata(row, self.TAXIDS) == [
            ("TAXON", "Adamaviridae"), ("TXID", "3152106"), ("PROTEIN", "Rep"),
            ("RANGE", "Conservative"), ("TYPE", "Short")]

    def test_discriminative_has_genus_block(self):
        row = {"Taxon": "Adamaviridae", "Genus": "Andorivirus", "Genus TxID": "3152140.0",
               "Protein": "Capsid", "Model type": "discriminative"}
        assert annotation.build_metadata(row, self.TAXIDS) == [
            ("TAXON", "Adamaviridae"), ("TXID", "3152106"), ("PROTEIN", "Capsid"),
            ("RANGE", "Discriminative"), ("TYPE", "Short"),
            None,
            ("TAXON", "Andorivirus"), ("TXID", "3152140"), ("PROTEIN", "Capsid")]

    def test_type_follows_full_length_per_model_type(self):
        row = {"Taxon": "Adamaviridae", "Genus": "Andorivirus", "Genus TxID": "3152140",
               "Protein": "Capsid", "Model type": "Discriminative"}
        fields = annotation.build_metadata(
            row, self.TAXIDS, full_length={"conservative": "no", "discriminative": "yes"})
        assert fields[3:5] == [("RANGE", "Discriminative"), ("TYPE", "Full-length")]
        # RANGE/TYPE only in the first block
        assert [t for t, _ in fields[6:]] == ["TAXON", "TXID", "PROTEIN"]

    def test_genus_taxid_falls_back_to_resolver(self):
        row = {"Taxon": "Adamaviridae", "Genus": "Andorivirus",
               "Protein": "Capsid", "Model type": "discriminative"}
        fields = annotation.build_metadata(row, self.TAXIDS, resolver=lambda g: "42")
        assert fields[-2] == ("TXID", "42")

    def test_genus_without_taxid_keeps_block(self):
        row = {"Taxon": "Adamaviridae", "Genus": "Andorivirus", "Genus TxID": "",
               "Protein": "Capsid", "Model type": "discriminative"}
        fields = annotation.build_metadata(row, self.TAXIDS, resolver=lambda g: "")
        assert fields[5:] == [None, ("TAXON", "Andorivirus"), ("PROTEIN", "Capsid")]

    def test_missing_taxid_keeps_block_without_txid(self):
        row = {"Taxon": "Unknownviridae", "Model type": "conservative"}
        assert annotation.build_metadata(row, self.TAXIDS) == [
            ("TAXON", "Unknownviridae"), ("PROTEIN", ""),
            ("RANGE", "Conservative"), ("TYPE", "Short")]


def test_model_fields_defaults_to_short():
    assert annotation.model_fields("conservative") == [
        ("RANGE", "Conservative"), ("TYPE", "Short")]
    assert annotation.model_fields("discriminative", {"discriminative": "yes"}) == [
        ("RANGE", "Discriminative"), ("TYPE", "Full-length")]


class TestResolveFullLength:
    def test_default_is_no(self):
        assert config.resolve_full_length(None) == "no"
        assert config.resolve_full_length(dict(config.TABAJARA_CON_DEFAULTS)) == "no"

    def test_fl_case_insensitive(self):
        assert config.resolve_full_length({"fl": " YES "}) == "yes"

    def test_fl_wins_over_full_length(self):
        assert config.resolve_full_length({"fl": "no", "full_length": "yes"}) == "no"
        assert config.resolve_full_length({"full_length": "yes"}) == "yes"

    def test_invalid_value_falls_back_to_no(self):
        assert config.resolve_full_length({"fl": "maybe"}) == "no"


_RANGE_FIELDS = _CON_FIELDS + [("RANGE", "Conservative"), ("TYPE", "Short")]


def test_annotate_hmm_text_range_type_idempotent_and_upgrades_old():
    old = annotation.annotate_hmm_text(_HMM_HEADER, _CON_FIELDS)
    out = annotation.annotate_hmm_text(old, _RANGE_FIELDS)
    assert ("PROTEIN Replication-associated\n"
            "RANGE   Conservative\n"
            "TYPE    Short\n"
            "\n"
            "HMM ") in out
    assert out == annotation.annotate_hmm_text(_HMM_HEADER, _RANGE_FIELDS)
    assert annotation.annotate_hmm_text(out, _RANGE_FIELDS) == out
    desc = annotation.annotate_hmm_text(out, _RANGE_FIELDS, mode="desc")
    assert ";PROTEIN=Replication-associated;RANGE=Conservative;TYPE=Short\n" in desc


def test_annotate_repository_writes_range_type_for_every_model(tmp_path, monkeypatch):
    monkeypatch.setattr(annotation, "_ncbi_genus_resolver", lambda: (lambda name: ""))
    monkeypatch.setattr(annotation, "_hmmer_accepts", lambda paths: None)
    repo = tmp_path / "HMMs_repository"
    for sub, name in (("conservative", "vHMM_1_C"), ("discriminative", "vHMM_2_D"),
                      ("discriminative", "vHMM_9_D"), ("conservative", "vHMM_3_C")):
        (repo / sub).mkdir(parents=True, exist_ok=True)
        (repo / sub / f"{name}.hmm").write_text(_HMM_HEADER)
    (tmp_path / "hmm_library.csv").write_text(
        "vHMM_ID;Taxon;Genus;Protein;Model type;Genus TxID\n"
        "vHMM_1;Adamaviridae;Family-wide models;Rep;conservative;\n"
        "vHMM_2;Adamaviridae;Andorivirus;Rep;discriminative;3152140\n"
        "vHMM_3;Unknownviridae;Family-wide models;Rep;conservative;\n")
    terms = tmp_path / "terms.csv"
    terms.write_text("tax_id;Name\n3152106;Adamaviridae\n")

    counts = annotation.annotate_hmms_repository(
        tmp_path, terms, full_length={"discriminative": "yes"})
    assert counts["annotated"] == 4
    assert counts["no_library_row"] == 1 and counts["no_txid"] == 1
    texts = {p.stem: p.read_text() for p in repo.glob("*/*.hmm")}
    for stem, text in texts.items():
        assert text.count("RANGE ") == 1 and text.count("TYPE ") == 1, stem
    assert "RANGE   Conservative\nTYPE    Short\n" in texts["vHMM_1_C"]
    assert "RANGE   Discriminative\nTYPE    Full-length\n" in texts["vHMM_2_D"]
    assert "CUTOFF SCORE\t38.08\nRANGE   Discriminative\nTYPE    Full-length\n\n" \
        in texts["vHMM_9_D"]
    assert "TXID" not in texts["vHMM_3_C"]


def test_library_keys_strip_model_suffix():
    assert annotation._library_keys("vHMM_61_C") == ["vHMM_61_C", "vHMM_61"]
    assert annotation._library_keys("vHMM_2_D") == ["vHMM_2_D", "vHMM_2"]


def test_load_taxids_trims_and_drops_float_suffix(tmp_path):
    terms = tmp_path / "terms.csv"
    terms.write_text("tax_id;Name\n3152106.0; Adamaviridae \n;Emptyviridae\n")
    assert annotation.load_taxids(terms) == {"Adamaviridae": "3152106"}
