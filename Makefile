.PHONY: all deps deps-rmsd deps-pdbx clean-deps

PYTHON ?= python3
VENDOR_DIR := $(CURDIR)/vendor
PDBX_URL := https://mmcif.wwpdb.org/docs/sw-examples/python/src/pdbx.tar.gz
NUMPY_VERSION := numpy<2
SCIPY_VERSION := scipy<1.12

all: deps

deps: deps-rmsd deps-pdbx

deps-rmsd:
	@mkdir -p "$(VENDOR_DIR)"
	"$(PYTHON)" -m pip install --target "$(VENDOR_DIR)" "$(NUMPY_VERSION)" "$(SCIPY_VERSION)" rmsd

deps-pdbx:
	@mkdir -p "$(VENDOR_DIR)"
	@curl -L "$(PDBX_URL)" | tar -xz -C "$(VENDOR_DIR)"
	@2to3 -w -n "$(VENDOR_DIR)/pdbx" >/dev/null 2>&1 || true

clean-deps:
	@rm -rf "$(VENDOR_DIR)"
