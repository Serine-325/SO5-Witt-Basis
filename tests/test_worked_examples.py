
"""
Automated verification of the three worked examples.

1. Two equal quadratic extensions
2. Two distinct quadratic extensions
3. Biquadratic extension

Each example checks its Witt basis using exact SymPy arithmetic.
"""

from pathlib import Path
import runpy

import pytest


EXAMPLES_DIR = Path(__file__).resolve().parents[1] / "examples"


@pytest.mark.parametrize(
    "filename",
    [
        "sum_quadratics_equal_extensions.py",
        "sum_quadratics_distinct_extensions.py",
        "biquadratic_canonical.py",
    ],
)
def test_worked_example(filename):
    """Run the example and verify that all its assertions pass."""

    path = EXAMPLES_DIR / filename

    assert path.is_file(), f"Missing worked example: {path}"

    # Execute the example, including its Gram-matrix assertions.
    runpy.run_path(str(path), run_name="__main__")
