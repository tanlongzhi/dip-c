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
