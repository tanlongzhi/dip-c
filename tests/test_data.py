"""Tests for dip_c.data helpers."""

import os
import tempfile

import pytest

from dip_c.data import list_bundled_data_files, resolve_data_file


class TestListBundledDataFiles:
    def test_returns_sorted_list(self):
        files = list_bundled_data_files()
        assert files == sorted(files)

    def test_no_gz_suffix(self):
        for name in list_bundled_data_files():
            assert not name.endswith(".gz")

    def test_known_files_present(self):
        files = list_bundled_data_files()
        assert "hg19.chr.len" in files
        assert "hg19.cpg.20k.txt" in files
        assert "mm10.chr.len" in files


class TestResolveDataFile:
    def test_opens_bundled_file_by_basename(self):
        f = resolve_data_file("hg19.chr.len")
        try:
            content = f.read()
            assert len(content) > 0
        finally:
            f.close()

    def test_basename_fallback_strips_directory(self):
        f = resolve_data_file("color/hg19.chr.len")
        try:
            content = f.read()
            assert len(content) > 0
        finally:
            f.close()

    def test_local_file_takes_priority(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tmp:
            tmp.write("local content\n")
            tmp_path = tmp.name
        try:
            f = resolve_data_file(tmp_path)
            assert f.read() == "local content\n"
            f.close()
        finally:
            os.unlink(tmp_path)

    def test_nonexistent_raises_file_not_found(self):
        with pytest.raises(FileNotFoundError, match="nonexistent_file_xyz"):
            resolve_data_file("nonexistent_file_xyz")

    def test_error_message_lists_available_files(self):
        with pytest.raises(FileNotFoundError, match="Available bundled files"):
            resolve_data_file("no_such_file")

    def test_resolve_plain_bundled_file(self, tmp_path, monkeypatch):
        """Line 45: resolve_data_file returns plain file from COLOR_DIR."""
        import dip_c.data as data_mod
        # Create a temp dir with a plain (non-.gz) file
        plain = tmp_path / "test_plain.txt"
        plain.write_text("plain content\n")
        monkeypatch.setattr(data_mod, "COLOR_DIR", str(tmp_path))
        f = resolve_data_file("test_plain.txt")
        assert f.read() == "plain content\n"
        f.close()


class TestListBundledNonGz:
    def test_non_gz_file_listed(self, tmp_path, monkeypatch):
        """Line 17: non-.gz files in COLOR_DIR are listed as-is."""
        import dip_c.data as data_mod
        # Create a temp dir with a non-.gz file
        (tmp_path / "plain_file.txt").write_text("")
        (tmp_path / "compressed.txt.gz").write_bytes(b"")
        monkeypatch.setattr(data_mod, "COLOR_DIR", str(tmp_path))
        files = list_bundled_data_files()
        assert "plain_file.txt" in files
        assert "compressed.txt" in files
