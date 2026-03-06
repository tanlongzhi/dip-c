.PHONY: install install-dev

PYTHON ?= python3

install:
	"$(PYTHON)" -m pip install .

install-dev:
	"$(PYTHON)" -m pip install -e ".[dev,seg]"
