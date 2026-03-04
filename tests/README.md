# Running Tests for dip-c

## Installation

Install test dependencies:
```bash
pip install -r tests/requirements.txt
# or if using the vendored dependencies:
python -m pip install --target vendor pytest pytest-cov
```

## Running Tests

Run all tests:
```bash
pytest
```

Run with coverage report:
```bash
pytest --cov=. --cov-report=html --cov-report=term
```

Run specific test file:
```bash
pytest tests/test_classes.py
```

Run specific test class:
```bash
pytest tests/test_classes.py::TestHaplotypes
```

Run specific test:
```bash
pytest tests/test_classes.py::TestHaplotypes::test_haplotype_values
```

Run with verbose output:
```bash
pytest -v
```

## Test Structure

- `test_classes.py` - Unit tests for core data structures in `classes.py`
  - Haplotype utilities
  - Leg (genomic locus) operations
  - Con (contact) operations
  - G3dParticle (3D structure) operations
  - File parsing functions
  - Name conversion utilities

## Adding New Tests

1. Create test file in `tests/` directory with name `test_*.py`
2. Create test classes with name `Test*`
3. Create test functions with name `test_*`
4. Use pytest fixtures for reusable test data
5. Run `pytest` to verify

Example:
```python
def test_my_function():
    result = my_function(input_data)
    assert result == expected_output
```
