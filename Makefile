.PHONY: all deps deps-rmsd deps-pdbx clean-deps install install-dev

PYTHON ?= python3
VENDOR_DIR := $(CURDIR)/vendor
PDBX_URL := https://mmcif.wwpdb.org/docs/sw-examples/python/src/pdbx.tar.gz
NUMPY_VERSION := numpy<2
SCIPY_VERSION := scipy<1.12

all: deps

# --- Modern installation (recommended) ---

install:
	"$(PYTHON)" -m pip install .

install-dev:
	"$(PYTHON)" -m pip install -e ".[dev,seg]"

# --- Legacy vendored installation (for backward compat) ---

deps: deps-rmsd deps-pdbx

deps-rmsd:
	@mkdir -p "$(VENDOR_DIR)"
	"$(PYTHON)" -m pip install --target "$(VENDOR_DIR)" "$(NUMPY_VERSION)" "$(SCIPY_VERSION)" rmsd

deps-pdbx:
	@mkdir -p "$(VENDOR_DIR)"
	@curl -L "$(PDBX_URL)" | tar -xz -C "$(VENDOR_DIR)"
	@2to3 -w -n "$(VENDOR_DIR)/pdbx" >/dev/null 2>&1 || true
	@echo "Applying PdbxReader escape sequence patch..."
	@if patch -p0 -d "$(VENDOR_DIR)" < "$(CURDIR)/patch/pdbx_reader_escape_sequences.patch" 2>/dev/null; then \
		echo "✓ Patch applied successfully"; \
	elif patch -p0 -d "$(VENDOR_DIR)" -R --dry-run < "$(CURDIR)/patch/pdbx_reader_escape_sequences.patch" >/dev/null 2>&1; then \
		echo "✓ Patch already applied"; \
	else \
		echo "✗ Patch failed to apply"; \
		exit 1; \
	fi

clean-deps:
	@rm -rf "$(VENDOR_DIR)"
