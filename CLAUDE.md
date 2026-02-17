# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

"Think Like A Physicist" — a terse undergraduate physics textbook by Casey Handmer and Jonathan Whitmore. Quarto book that renders to both a website and PDF. Source material includes Casey's "Think Like A Physicist" PDF and 169 handwritten IPhO (International Physics Olympiad) training note photos in `images/IPhO notes/`.

Scope: physics thinking and problem-solving only. No lab sections, no test-taking content.

## Workflow

Commit early and often, and push to GitHub after each logical chunk of work.

### TODO Hygiene

`TODO.md` tracks current and upcoming work. When a logical section (e.g., an entire Phase) is fully complete, move its items from `TODO.md` into `DONE.md` to keep the TODO file clean and focused on what's left.

## Build & Development Commands

```bash
# Preview site locally (live reload)
quarto preview

# Render full site (output goes to _site/)
quarto render

# Run derivation verification tests
uv run pytest tests/ -v
```

## Architecture

Quarto book project at repo root. Chapter files are `.qmd` at top level, organized into parts in `_quarto.yml`. Old attempts archived in `_archive/` (gitignored). Source images in `images/`. Output goes to `_book/` (gitignored).

## Content Conventions

### The 7 Ds Problem Template

Every worked problem MUST follow the 7 Ds (and the little s) structure. Each step must appear as a subsection, even if the content is brief. If a step is not applicable, write a short explanation (e.g., "Not needed for this problem" or "Diagram not needed — purely algebraic").

```markdown
## Problem Title {.unnumbered}

Problem statement goes here.

::: {.callout-warning}
## TODO: Diagram
Description of what the diagram should show. Include key elements, labels, and coordinate system.
:::

### Diagram {.unnumbered}
[Diagram or placeholder callout. Should be big — 2/3 of a page. Load the problem into your GPU.]

### Directions {.unnumbered}
[Define positive/negative directions. State coordinate system choice.]

### Definitions & Given Data {.unnumbered}
[Table or list of ALL variables with descriptions and units. Put everything on the page.]

### Diagnosis {.unnumbered}
[What type of problem is this? Identify the relevant physics: conservation laws, force laws, etc.]

### Derivation {.unnumbered}
[Write fundamental equations. Express diagnosis in symbols. Need as many equations as unknowns. Check dimensions as you go.]

### Determination (D'algebra) {.unnumbered}
[Algebraic manipulation to isolate the answer. Box the final result with $$\boxed{...}$$]

### Dimensions {.unnumbered}
[Verify dimensions of the answer. Check limiting cases / sanity check. Use smiley or checkmark for each case that makes sense.]

### Substitution {.unnumbered}
[Only if needed. Rough numerical calculation by hand. Check units. Include error estimate.]
```

### Diagram Placeholders

When a diagram is needed but not yet created, use a warning callout describing what the diagram should contain:

```markdown
::: {.callout-warning}
## TODO: Diagram
Description of diagram contents, key elements, labels, etc.
:::
```

### IPhO Note Photos

The handwritten IPhO note photos in `images/IPhO notes/` are reference material only — do NOT embed them in the book. Transcribe their content into proper Quarto markdown with LaTeX math.

### Derivation Verification Tests

Every boxed result (`$$\boxed{...}$$`) in the book MUST have a corresponding sympy test in `tests/test_<chapter>.py`. When transcribing a new chapter or adding a worked problem:

1. Add a test that derives the boxed result symbolically and asserts it matches.
2. Add tests for limiting cases mentioned in the Dimensions section.
3. Run `uv run pytest tests/ -v` to verify before committing.

Tests use `pytest` + `sympy`. One test file per chapter (e.g., `tests/test_mechanics.py`).

## Deployment

Push to `main` triggers `.github/workflows/publish.yml`:
1. Quarto renders the site to `_site/`
2. `actions/deploy-pages` deploys to GitHub Pages

GitHub Pages source must be set to "GitHub Actions" (not "Deploy from a branch").

## Theming

Dual light/dark mode using flatly (light) and darkly (dark) Bootstrap themes, with custom SCSS overrides.
# Python-Specific Agent Instructions

These instructions apply when working with Python code.

## Package Management with uv

Use **uv** exclusively for Python package management.

- All Python dependencies **must be installed, synchronized, and locked** using uv
- Never use pip, pip-tools, poetry, or conda directly for dependency management

### Common uv Commands

```bash
# Install dependencies
uv add <package>

# Remove dependencies
uv remove <package>

# Sync dependencies from lockfile
uv sync

# Run Python scripts
uv run script.py

# Run Python tools
uv run pytest
uv run ruff check .
uv run mypy .

# Launch Python REPL
uv run python
```

### PEP 723 Inline Script Metadata

For standalone scripts with dependencies, use PEP 723 inline metadata:

```python
# /// script
# dependencies = [
#   "requests",
#   "pandas>=2.0",
# ]
# ///

import requests
import pandas as pd

# Your script here
```

Run with: `uv run script.py`

Manage script dependencies:
```bash
uv add requests --script script.py
uv remove requests --script script.py
```

## Python Style Guidelines

- **Follow PEP 8** style guide for Python code
- Use **4 spaces** for indentation (not tabs)
- **Maximum line length of 88 characters** (Black formatter standard)
- Use **snake_case** for functions and variables
- Use **PascalCase** for classes
- Use **UPPER_CASE** for constants

### Formatting and Linting

```bash
# Format code with ruff
uv run ruff format .

# Lint code
uv run ruff check .

# Auto-fix issues
uv run ruff check --fix .

# Type checking with mypy
uv run mypy .
```

## Python Best Practices

### Type Hints

Always use type hints for function signatures and complex variables:

```python
from typing import Optional, List, Dict, Union
from pathlib import Path

def process_data(
    items: List[str],
    config: Dict[str, Any],
    output_path: Optional[Path] = None
) -> List[int]:
    """Process items according to config."""
    results: List[int] = []
    for item in items:
        results.append(len(item))
    return results
```

### String Formatting

Prefer f-strings for string formatting (Python 3.6+):

```python
# Good
name = "Alice"
age = 30
message = f"Hello, {name}! You are {age} years old."

# Avoid
message = "Hello, {}! You are {} years old.".format(name, age)
message = "Hello, %s! You are %d years old." % (name, age)
```

### File Operations with pathlib

Always use `pathlib` for file operations instead of `os.path`:

```python
from pathlib import Path

# Good
data_dir = Path("data")
file_path = data_dir / "input.txt"

if file_path.exists():
    content = file_path.read_text()
    lines = file_path.read_text().splitlines()

# Write to file
output_path = data_dir / "output.txt"
output_path.write_text("Hello, world!")

# Avoid
import os
file_path = os.path.join("data", "input.txt")
if os.path.exists(file_path):
    with open(file_path, 'r') as f:
        content = f.read()
```

### Context Managers

Use context managers (`with` statement) for resource management:

```python
# File handling
with open('file.txt') as f:
    content = f.read()

# Or with pathlib
from pathlib import Path
content = Path('file.txt').read_text()

# Database connections
with get_db_connection() as conn:
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM users")

# Custom context managers
from contextlib import contextmanager

@contextmanager
def timer(name: str):
    start = time.time()
    yield
    print(f"{name} took {time.time() - start:.2f}s")

with timer("data processing"):
    process_data()
```

### Comprehensions and Generators

Leverage list comprehensions and generator expressions where appropriate:

```python
# List comprehension
squares = [x**2 for x in range(10)]
evens = [x for x in numbers if x % 2 == 0]

# Dict comprehension
word_lengths = {word: len(word) for word in words}

# Set comprehension
unique_lengths = {len(word) for word in words}

# Generator expression (memory efficient)
sum_of_squares = sum(x**2 for x in range(1000000))

# Generator function
def read_large_file(file_path: Path):
    """Read file line by line without loading all into memory."""
    with file_path.open() as f:
        for line in f:
            yield line.strip()
```

### Common Python Patterns

```python
# Enumerate for indexed loops
for idx, item in enumerate(items):
    print(f"Item {idx}: {item}")

# Zip for parallel iteration
for x, y in zip(list1, list2):
    print(f"{x} -> {y}")

# Dict.get() with defaults
value = my_dict.get('key', default_value)

# Unpacking
first, *middle, last = items
head, *tail = items

# Multiple assignment
x, y = y, x  # Swap

# Any/all for collections
has_negatives = any(x < 0 for x in numbers)
all_positive = all(x > 0 for x in numbers)
```

## Error Handling

### Use Specific Exceptions

```python
# Good - specific exceptions
try:
    result = int(user_input)
except ValueError as e:
    logger.error(f"Invalid input: {user_input}")
    raise

except KeyError as e:
    logger.error(f"Missing key: {e}")
    return default_value

# Avoid - bare except
try:
    risky_operation()
except:  # Don't do this!
    pass
```

### Custom Exceptions

Create custom exceptions for domain-specific errors:

```python
class ValidationError(Exception):
    """Raised when data validation fails."""
    pass

class DataNotFoundError(Exception):
    """Raised when required data is not found."""
    pass

def validate_user(user_data: dict) -> None:
    if 'email' not in user_data:
        raise ValidationError("Email is required")
```

### Error Context

Include helpful context in error messages:

```python
# Good
if not file_path.exists():
    raise FileNotFoundError(
        f"Data file not found: {file_path}. "
        f"Expected location: {file_path.absolute()}"
    )

# Less helpful
if not file_path.exists():
    raise FileNotFoundError("File not found")
```

## Logging

Use the logging module, not print statements:

```python
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger(__name__)

# Use appropriate log levels
logger.debug("Detailed information for debugging")
logger.info("General information about program execution")
logger.warning("Something unexpected but not critical")
logger.error("Error occurred but program continues")
logger.exception("Error with full traceback")  # Use in except blocks
logger.critical("Critical error, program may not continue")

# Example usage
def process_file(file_path: Path) -> None:
    logger.info(f"Processing file: {file_path}")
    try:
        data = file_path.read_text()
        logger.debug(f"Read {len(data)} bytes")
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        raise
```

## Testing

### Use pytest

```bash
# Install pytest
uv add --dev pytest pytest-cov

# Run tests
uv run pytest

# Run with coverage
uv run pytest --cov=src tests/

# Run specific test
uv run pytest tests/test_module.py::test_function
```

### Test Structure

```python
import pytest
from pathlib import Path

def test_simple_case():
    """Test basic functionality."""
    result = add(2, 3)
    assert result == 5

def test_with_fixture(tmp_path: Path):
    """Test using pytest's tmp_path fixture."""
    test_file = tmp_path / "test.txt"
    test_file.write_text("Hello")
    assert test_file.read_text() == "Hello"

@pytest.mark.parametrize("input,expected", [
    (2, 4),
    (3, 9),
    (4, 16),
])
def test_square(input, expected):
    """Test square function with multiple inputs."""
    assert square(input) == expected

def test_exception():
    """Test that exception is raised."""
    with pytest.raises(ValueError, match="Invalid input"):
        process_invalid_data()
```

### Fixtures

```python
import pytest

@pytest.fixture
def sample_data():
    """Provide sample data for tests."""
    return {"name": "Alice", "age": 30}

@pytest.fixture
def database_connection():
    """Provide database connection with cleanup."""
    conn = create_connection()
    yield conn
    conn.close()

def test_with_fixture(sample_data):
    """Test using the fixture."""
    assert sample_data["name"] == "Alice"
```

## Documentation

### Docstrings

Write docstrings for all public modules, functions, classes, and methods:

```python
def calculate_total(items: List[float], tax_rate: float = 0.1) -> float:
    """
    Calculate total price including tax.

    Args:
        items: List of item prices
        tax_rate: Tax rate as decimal (default: 0.1 for 10%)

    Returns:
        Total price including tax

    Raises:
        ValueError: If tax_rate is negative

    Examples:
        >>> calculate_total([10.0, 20.0], tax_rate=0.1)
        33.0
    """
    if tax_rate < 0:
        raise ValueError("Tax rate cannot be negative")

    subtotal = sum(items)
    return subtotal * (1 + tax_rate)
```

## Virtual Environments

- Always use virtual environments (uv handles this automatically)
- Document required Python version in README or pyproject.toml
- Keep `pyproject.toml` and lock file in version control

## Performance Tips

- Use built-in functions (they're implemented in C and fast)
- Use `set` for membership testing instead of `list`
- Use `dict` for lookups instead of searching through lists
- Use generators for large datasets to save memory
- Profile before optimizing: `uv run python -m cProfile script.py`

```python
# Fast membership test
valid_ids = {1, 2, 3, 4, 5}  # set
if user_id in valid_ids:  # O(1)
    pass

# Slow membership test
valid_ids = [1, 2, 3, 4, 5]  # list
if user_id in valid_ids:  # O(n)
    pass
```
