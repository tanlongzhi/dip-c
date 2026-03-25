"""Tests for dip_c.cli."""

import importlib
import os
import runpy
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
    "force", "pairs2con",
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


# Step 2a: _version.py fallback coverage
class TestVersionFallbacks:
    def test_package_not_found_fallback(self, monkeypatch):
        """PackageNotFoundError → __version__ == 'unknown'."""
        from importlib.metadata import PackageNotFoundError

        def fake_version(name):
            raise PackageNotFoundError(name)

        monkeypatch.setattr("importlib.metadata.version", fake_version)
        import dip_c._version
        importlib.reload(dip_c._version)
        assert dip_c._version.__version__ == "unknown"
        # Restore
        importlib.reload(dip_c._version)

    def test_import_error_fallback(self, monkeypatch):
        """ImportError on importlib.metadata → __version__ == 'unknown'."""
        import builtins
        real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "importlib.metadata":
                raise ImportError("no metadata")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        import dip_c._version
        importlib.reload(dip_c._version)
        assert dip_c._version.__version__ == "unknown"
        # Restore
        monkeypatch.undo()
        importlib.reload(dip_c._version)


# Step 2b: cli.py __main__ block
class TestCliMain:
    def test_main_block_triggers_sys_exit(self, monkeypatch):
        """if __name__ == '__main__': sys.exit(main()) raises SystemExit."""
        monkeypatch.setattr(sys, "argv", ["dip-c"])
        with pytest.raises(SystemExit) as exc_info:
            runpy.run_module("dip_c.cli", run_name="__main__")
        assert exc_info.value.code == 1


# Step 3: Getopt error paths — every command returns 1 on invalid option
GETOPT_COMMANDS = [
    ("seg", "dip_c.commands.seg", "seg"),
    ("con", "dip_c.commands.con", "con"),
    ("dedup", "dip_c.commands.dedup", "dedup"),
    ("reg", "dip_c.commands.reg", "reg"),
    ("clean", "dip_c.commands.clean", "clean"),
    ("bincon", "dip_c.commands.bincon", "bincon"),
    ("impute", "dip_c.commands.impute", "impute"),
    ("impute3", "dip_c.commands.impute3", "impute3"),
    ("clean3", "dip_c.commands.clean3", "clean3"),
    ("reg3", "dip_c.commands.reg3", "reg3"),
    ("mkcon", "dip_c.commands.mkcon", "mkcon"),
    ("color", "dip_c.commands.color", "color"),
    ("color2", "dip_c.commands.color2", "color2"),
    ("mgcolor", "dip_c.commands.mgcolor", "mgcolor"),
    ("vis", "dip_c.commands.vis", "vis"),
    ("info", "dip_c.commands.info", "info"),
    ("ard", "dip_c.commands.ard", "ard"),
    ("cv", "dip_c.commands.cv", "cv"),
    ("con3", "dip_c.commands.con3", "con3"),
    ("dist", "dip_c.commands.dist", "dist"),
    ("align", "dip_c.commands.align", "align"),
    ("pd", "dip_c.commands.pd", "pd"),
    ("rg", "dip_c.commands.rg", "rg"),
    ("pos", "dip_c.commands.pos", "pos"),
    ("tad", "dip_c.commands.tad", "tad"),
    ("exp", "dip_c.commands.exp", "exp"),
    ("force", "dip_c.commands.force", "force"),
]


class TestGetoptError:
    @pytest.mark.parametrize("cmd,module,func", GETOPT_COMMANDS)
    def test_invalid_option_returns_error(self, cmd, module, func, tmp_path):
        """Passing --zzz should trigger getopt.GetoptError and return 1."""
        mod = importlib.import_module(module)
        fn = getattr(mod, func)
        # Provide a dummy file argument so we don't fail on "no args"
        dummy = tmp_path / "dummy"
        dummy.write_text("")
        ret = fn([cmd, "--zzz", str(dummy)])
        assert ret == 1
