# Release Management

This directory contains scripts and documentation for managing dip-c releases.

## Release Workflow

### 1. Prepare a Release

Run the release preparation script:

```bash
./release/prepare_release.sh
```

This script will:
- Check for uncommitted changes
- Run the full test suite with coverage
- Regenerate the coverage badge
- Show current version/tag information
- Guide you through the release process

### 2. Version Numbering

Dip-C follows semantic versioning: `v0.MAJOR.MINOR`

- **MAJOR**: Increment for breaking API changes or major feature additions
- **MINOR**: Increment for backward-compatible improvements and bug fixes

Example versions:
- `v0.1.0` - Initial release with basic functionality
- `v0.2.0` - Added new features, maintains backward compatibility
- `v0.3.0` - Bug fixes and performance improvements
- `v1.0.0` - (Future) Stable API, ready for production use

### 3. Creating a Release

After running `prepare_release.sh`:

1. **Commit any changes** (updated badges, documentation):
   ```bash
   git add images/coverage-badge.svg README.md
   git commit -m "Update coverage badge for v0.X.Y release"
   ```

2. **Create an annotated tag**:
   ```bash
   git tag -a v0.X.Y -m "Release v0.X.Y: Brief description of changes"
   ```

3. **Push changes and tag**:
   ```bash
   git push origin main
   git push origin v0.X.Y
   ```

4. **Create a GitHub release** (optional):
   - Go to the GitHub repository
   - Click "Releases" → "Draft a new release"
   - Select your tag
   - Add release notes describing changes
   - Publish the release

## Future: Package Distribution

### BioConda (Planned)

Future script: `./release/build_conda.sh`

Will create a conda package recipe and submit to BioConda:
- Package metadata in `conda/meta.yaml`
- Build and test in conda environment
- Submit PR to bioconda-recipes repository

### PyPI (Planned)

Future script: `./release/build_pypi.sh`

Will build and upload to Python Package Index:
- Use `setup.py` or `pyproject.toml`
- Build source distribution and wheel
- Upload with `twine upload`

### GitHub Actions (Planned)

Automated CI/CD workflow:
- Run tests on each push
- Build packages on tag creation
- Automatically publish to BioConda and PyPI
- Update documentation

## Release Checklist

Before each release, ensure:

- [ ] All tests pass (`pytest`)
- [ ] Coverage badge is up to date
- [ ] README.md is current
- [ ] CHANGELOG.md is updated (if exists)
- [ ] Version number follows semantic versioning
- [ ] Git working directory is clean
- [ ] Tag message is descriptive

## Current Status

**Testing**: ✓ Comprehensive unit tests with pytest  
**Coverage**: ✓ Automated badge generation  
**Packaging**: ⚠️ Planned (BioConda, PyPI)  
**CI/CD**: ⚠️ Planned (GitHub Actions)  
**Documentation**: ⚠️ Planned (Read the Docs)
