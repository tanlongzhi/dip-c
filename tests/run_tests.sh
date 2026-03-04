#!/bin/bash
# Run unit tests and generate coverage badge

set -e  # Exit on error

# Change to repository root
cd "$(dirname "$0")/.."

echo "Running unit tests with coverage..."
pytest --cov=. --cov-report=xml --cov-report=term tests/

echo "Generating coverage badge..."
# Remove old badge if it exists
rm -f images/coverage-badge.svg

# Generate new badge
genbadge coverage -i coverage.xml -o images/coverage-badge.svg

echo "Tests complete! Coverage badge saved to images/coverage-badge.svg"
echo "View detailed coverage report: open htmlcov/index.html"
