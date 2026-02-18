#!/bin/bash
# Release preparation script for dip-c
# This script runs tests, updates coverage badges, and prepares for a new release

set -e  # Exit on error

echo "========================================"
echo "Dip-C Release Preparation"
echo "========================================"
echo ""

# Get the repository root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$REPO_ROOT"

# Step 1: Check for uncommitted changes
echo "Step 1: Checking for uncommitted changes..."
if ! git diff-index --quiet HEAD --; then
    echo "⚠️  WARNING: You have uncommitted changes."
    echo "Please commit or stash your changes before preparing a release."
    echo ""
    git status --short
    echo ""
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 1
    fi
fi
echo "✓ Working directory is clean (or continuing with changes)"
echo ""

# Step 2: Run tests
echo "Step 2: Running test suite..."
cd tests
if pytest --cov=.. --cov-report=xml --cov-report=term; then
    echo "✓ All tests passed"
else
    echo "❌ Tests failed. Fix issues before releasing."
    exit 1
fi
cd "$REPO_ROOT"
echo ""

# Step 3: Regenerate coverage badge
echo "Step 3: Regenerating coverage badge..."
cd tests
if [ -f ../images/coverage-badge.svg ]; then
    rm ../images/coverage-badge.svg
    echo "  Removed old badge"
fi
if genbadge coverage -i ../coverage.xml -o ../images/coverage-badge.svg; then
    echo "✓ Coverage badge updated"
else
    echo "❌ Failed to generate coverage badge"
    exit 1
fi
cd "$REPO_ROOT"
echo ""

# Step 4: Show current version/tags and suggest next version
echo "Step 4: Version numbering..."
LATEST_TAG=$(git describe --tags --abbrev=0 2>/dev/null || echo "")
if [ -z "$LATEST_TAG" ]; then
    echo "  No existing tags found."
    SUGGESTED_VERSION="v0.1.0"
else
    echo "  Current version: $LATEST_TAG"
    COMMITS_SINCE=$(git rev-list ${LATEST_TAG}..HEAD --count 2>/dev/null || echo "0")
    echo "  Commits since tag: $COMMITS_SINCE"
    
    # Parse version and suggest next minor version
    if [[ $LATEST_TAG =~ ^v([0-9]+)\.([0-9]+)\.([0-9]+)$ ]]; then
        MAJOR="${BASH_REMATCH[1]}"
        MINOR="${BASH_REMATCH[2]}"
        PATCH="${BASH_REMATCH[3]}"
        SUGGESTED_VERSION="v${MAJOR}.${MINOR}.$((PATCH + 1))"
    elif [[ $LATEST_TAG =~ ^v([0-9]+)\.([0-9]+)$ ]]; then
        MAJOR="${BASH_REMATCH[1]}"
        MINOR="${BASH_REMATCH[2]}"
        SUGGESTED_VERSION="v${MAJOR}.$((MINOR + 1))"
    else
        SUGGESTED_VERSION="v0.1.0"
    fi
fi
echo "  Suggested next version: $SUGGESTED_VERSION"
echo ""

# Prompt for version number
echo "Enter new version number (include 'v' prefix, e.g., v0.2.4)"
echo "Or press Enter to use suggested version [$SUGGESTED_VERSION]:"
read -r NEW_VERSION

# Use suggested version if empty
if [ -z "$NEW_VERSION" ]; then
    NEW_VERSION="$SUGGESTED_VERSION"
    echo "Using suggested version: $NEW_VERSION"
fi

# Validate version format
if [[ ! $NEW_VERSION =~ ^v[0-9]+\.[0-9]+(\.[0-9]+)?$ ]]; then
    echo "❌ Invalid version format. Expected format: v0.1.0 or v0.1"
    echo "   (Must start with 'v' followed by numbers and dots)"
    exit 1
fi
echo "✓ Version validated: $NEW_VERSION"
echo ""

# Step 5: Create tag (optional)
echo "Step 5: Create git tag..."
read -p "Create tag '$NEW_VERSION' now? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Enter release message (brief description of changes):"
    read -r RELEASE_MSG
    if [ -z "$RELEASE_MSG" ]; then
        RELEASE_MSG="Release $NEW_VERSION"
    fi
    
    if git tag -a "$NEW_VERSION" -m "$RELEASE_MSG"; then
        echo "✓ Tag '$NEW_VERSION' created successfully"
        TAG_CREATED=true
    else
        echo "❌ Failed to create tag"
        exit 1
    fi
else
    echo "Skipping tag creation (you can create it manually later)"
    TAG_CREATED=false
fi
echo ""

# Step 6: Guide next steps
echo "========================================"
echo "Release Preparation Complete! ✓"
echo "========================================"
echo ""
if [ "$TAG_CREATED" = true ]; then
    echo "Next steps:"
    echo "  1. Review the tag you just created:"
    echo "     git show $NEW_VERSION"
    echo ""
    echo "  2. Commit any updated badges or documentation:"
    echo "     git add images/coverage-badge.svg"
    echo "     git commit -m 'Update coverage badge for $NEW_VERSION release'"
    echo ""
    echo "  3. Push changes and tag:"
    echo "     git push origin main"
    echo "     git push origin $NEW_VERSION"
    echo ""
    echo "  4. (Future) Build and upload packages:"
    echo "     - BioConda: ./release/build_conda.sh"
    echo "     - PyPI: ./release/build_pypi.sh"
else
    echo "Next steps:"
    echo "  1. Review the test coverage report (see above)"
    echo "  2. Commit any updated badges or documentation:"
    echo "     git add images/coverage-badge.svg"
    echo "     git commit -m 'Update coverage badge for release'"
    echo ""
    echo "  3. Create the release tag manually:"
    echo "     git tag -a $NEW_VERSION -m 'Release $NEW_VERSION'"
    echo ""
    echo "  4. Push changes and tag:"
    echo "     git push origin main"
    echo "     git push origin $NEW_VERSION"
    echo ""
    echo "  5. (Future) Build and upload packages:"
    echo "     - BioConda: ./release/build_conda.sh"
    echo "     - PyPI: ./release/build_pypi.sh"
fi
echo ""
