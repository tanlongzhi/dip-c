# Patch Directory

This directory contains modifications for external dependencies used by Dip-C.

## File Types

### 1. Replacement Source Files

Complete modified versions of source files from external projects. These are **manually copied** to replace the original files in their respective dependencies.

- **`nuc_dynamics.py`** - Modified version of [nuc_dynamics](https://github.com/TheLaueLab/nuc_dynamics) that:
  - Changes the backbone energy function
  - Skips removal of isolated contacts
  - Outputs in 3D Genome (3DG) format instead of PDB format (which has a 99,999-atom limit)
  - **Usage**: Manually replace the original `nuc_dynamics.py` after downloading nuc_dynamics (see [README_old.md](../README_old.md))

- **`dyn_util.pyx`** - Cython utilities for nuc_dynamics
  - **Usage**: Part of the nuc_dynamics patching process

- **`trim.c`** - Modified version of LIANTI's trim utility
  - **Usage**: Manually replace LIANTI's `trim.c` and recompile (see [README_old.md](../README_old.md))

### 2. Patch Files (`.patch`)

Standard unified diff patches applied automatically during the build process using the `patch` command.

- **`pdbx_reader_escape_sequences.patch`** - Fixes Python 3.12+ SyntaxWarnings in the vendored PdbxReader library
  - Adds raw string prefixes (`r`) to regex patterns with escape sequences (`\S`, `\s`)
  - **Applied**: Automatically during `make deps-pdbx` after `2to3` conversion
  - **Format**: Unified diff (`-p0` level, applied in `vendor/` directory)

## How Patching Works

### Automatic Patches (`.patch` files)

Applied via the Makefile during dependency installation:

```makefile
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
```

### Manual Replacement Files

Must be manually copied to their target locations:

1. Download the external dependency
2. Replace the specified files with versions from this `patch/` directory
3. Recompile if necessary

See [README_old.md](../README_old.md) for detailed instructions on patching nuc_dynamics and LIANTI.

## Creating New Patches

### For Automatic Patches

1. Make changes to the file in `vendor/` directory
2. Generate a unified diff:
   ```bash
   diff -u original.py modified.py > patch/description.patch
   ```
3. Add patch application to the appropriate Makefile target

### For Replacement Files

Simply add the complete modified file to this directory and document the replacement instructions in the README.
