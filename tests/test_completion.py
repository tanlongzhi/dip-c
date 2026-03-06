"""Tests for dip_c.completion."""

from dip_c.completion import SUBCOMMANDS, generate_completion_script


class TestGenerateCompletionScript:
    def test_contains_complete_directive(self):
        script = generate_completion_script()
        assert "complete -o default -F _dip_c_completions dip-c" in script

    def test_contains_subcommands(self):
        script = generate_completion_script()
        for cmd in ["color", "vis", "bincon", "con3", "force"]:
            assert cmd in script

    def test_contains_data_filenames(self):
        script = generate_completion_script()
        assert "hg19.chr.len" in script
        assert "hg19.cpg.20k.txt" in script

    def test_contains_data_file_flags(self):
        script = generate_completion_script()
        assert "-c" in script
        assert "-n" in script

    def test_subcommands_list_matches_cli(self):
        # Sanity check that key commands are in the list
        assert "seg" in SUBCOMMANDS
        assert "data-path" in SUBCOMMANDS
        assert "force" in SUBCOMMANDS
