# Releasing a new version of dip-c to PyPI

This guide explains how to publish a new version of dip-c so that users can get
it with `pip install -U run-dipc`.

## How it works (one-time setup, already done)

- The PyPI package is **`run-dipc`**. Its version string comes from the
  **`VERSION`** file in the repo (read automatically by `pyproject.toml`).
- Publishing is automated by **`.github/workflows/publish.yml`**, which runs on
  the **`pypi-packaging`** branch.
- The workflow fires when a tag matching **`v*`** (e.g. `v1.9.10`) is pushed to
  **this repository** (`tanlongzhi/dip-c`).
- Authentication uses **PyPI Trusted Publishing (OIDC)** — there is **no API
  token or password to manage**. PyPI trusts releases that come from
  `tanlongzhi/dip-c` running `publish.yml`.

So a release is just: **bump `VERSION`, then create a matching `vX.Y.Z` tag on
the `pypi-packaging` branch.** GitHub Actions does the rest (test → build →
upload).

Anyone with **Write** access to `tanlongzhi/dip-c` can do this.

## The three rules

1. **The version must be new.** PyPI permanently rejects re-uploading an existing
   version. Always bump `VERSION`.
2. **The tag must match `VERSION`.** `VERSION` = `1.9.10` ⟷ tag = `v1.9.10`.
   (A workflow check enforces this and fails the run if they differ.)
3. **Target the `pypi-packaging` branch**, on `tanlongzhi/dip-c` (not a fork).
   That branch holds the packaging and the workflow; the trusted publisher is
   bound to this repo.

## Option A — GitHub web UI (no local git needed)

**Step 1 — bump the version.**
1. Open the `VERSION` file on the `pypi-packaging` branch:
   `https://github.com/tanlongzhi/dip-c/blob/pypi-packaging/VERSION`
2. Click the pencil (Edit) icon.
3. Change the number, e.g. `1.9.9` → `1.9.10`.
   - Patch bump (`1.9.9` → `1.9.10`) for fixes.
   - Minor bump (`1.9.9` → `1.10.0`) for new features.
4. Commit directly to `pypi-packaging` (or open a PR into it and merge).

**Step 2 — create the release.**
1. Go to the repo's **Releases** page → **Draft a new release**.
2. **Choose a tag** → type the matching tag, e.g. `v1.9.10` →
   **Create new tag: `v1.9.10` on publish**.
3. **Target:** select **`pypi-packaging`** (this is important — not `master`).
4. Add a title and release notes if you like (optional, makes a nice changelog).
5. Click **Publish release**.

Publishing the release creates the tag, which triggers the workflow.

**Step 3 — verify.**
1. Open the **Actions** tab and watch the **Publish to PyPI** run.
   It runs the test suite (Python 3.10–3.13), then builds and uploads.
2. When it finishes green, confirm at
   `https://pypi.org/project/run-dipc/` and with:
   ```
   pip install -U run-dipc
   ```

## Option B — git command line

```bash
git checkout pypi-packaging
git pull

# 1. bump the version
echo 1.9.10 > VERSION
git commit -am "v1.9.10"

# 2. push the bump to this repo
git push origin pypi-packaging          # 'origin' must be tanlongzhi/dip-c

# 3. tag and push the tag — this triggers the publish
git tag v1.9.10
git push origin v1.9.10
```

> If you work from a fork, push the tag to the `tanlongzhi/dip-c` remote, not
> your fork — only tags on `tanlongzhi/dip-c` trigger the trusted publish.

Then verify as in Step 3 above.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| Workflow run fails at **"Verify tag matches VERSION"** | The tag (`vX.Y.Z`) and the `VERSION` file disagree. Bump `VERSION` to match, or delete and recreate the tag. |
| Upload fails: **"File already exists"** | That version is already on PyPI. Bump `VERSION` to a new number and release again. |
| **Nothing happened** after creating the release | Check the tag looks like `v*`, the release **targeted `pypi-packaging`**, and it was created on `tanlongzhi/dip-c` (not a fork). |
| **Tests failed** | The publish step only runs after tests pass. Fix the failing test, then release again. |
| **"trusted publisher" / OIDC error** during upload | The trusted publisher on PyPI must point at `tanlongzhi/dip-c`, workflow `publish.yml`. Check PyPI → Manage `run-dipc` → Publishing. |
