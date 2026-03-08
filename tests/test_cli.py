"""Tests for dip_c.cli."""

import os
import sys

import pytest

from dip_c.cli import main


class TestCliUsage:
    def test_no_args_prints_usage(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["dip-c"])
        assert main() == 1
        captured = capsys.readouterr()
        assert "Usage: dip-c <command>" in captured.err
        assert "Commands:" in captured.err

    def test_usage_lists_subcommands(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["dip-c"])
        main()
        captured = capsys.readouterr()
        for cmd in ["seg", "con", "color", "vis", "bincon", "data-path"]:
            assert cmd in captured.err


class TestCliCompletion:
    def test_completion_prints_script(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["dip-c", "--completion"])
        assert main() == 0
        captured = capsys.readouterr()
        assert "complete -o default -F _dip_c_completions dip-c" in captured.out

    def test_completion_no_stderr(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["dip-c", "--completion"])
        main()
        captured = capsys.readouterr()
        assert captured.err == ""


class TestCliDataPath:
    def test_data_path_prints_path(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["dip-c", "data-path"])
        assert main() == 0
        captured = capsys.readouterr()
        path = captured.out.strip()
        assert path.endswith("data")
        assert os.path.isdir(path)


class TestCliUnknownCommand:
    def test_unknown_command_returns_error(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["dip-c", "not-a-command"])
        assert main() == 1
        captured = capsys.readouterr()
        assert "unknown command" in captured.err


ALL_COMMANDS = [
    "seg", "con", "dedup", "reg", "clean", "bincon", "impute", "impute3",
    "clean3", "reg3", "mkcon", "color", "color2", "mgcolor", "vis", "info",
    "ard", "cv", "con3", "dist", "align", "pd", "rg", "pos", "tad", "exp",
    "force",
]


class TestCliDispatch:
    @pytest.mark.parametrize("cmd", ALL_COMMANDS)
    def test_dispatch_returns_usage_when_no_args(self, monkeypatch, capsys, cmd):
        """Every command should print usage and return 1 when called with no arguments."""
        monkeypatch.setattr(sys, "argv", ["dip-c", cmd])
        ret = main()
        assert ret == 1
        captured = capsys.readouterr()
        assert "Usage:" in captured.err

    def test_successful_command_prints_timing(self, monkeypatch, capsys, tmp_path):
        """A successful command should print timing info to stderr."""
        con = tmp_path / "one.con"
        con.write_text("1,100,0\t1,200,0\n")
        monkeypatch.setattr(sys, "argv", ["dip-c", "info", str(con)])
        ret = main()
        assert ret == 0
        captured = capsys.readouterr()
        assert "finished in" in captured.err
