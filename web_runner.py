"""Web adapter for the verified command-line Witt basis constructor.

This adapter deliberately does not reimplement any Witt basis formulas.
It validates square-class choices, supplies CLI answers, and parses output.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import subprocess
import sys

CLASSES = ("epsilon", "pi", "epsilon*pi")
MU_OPTIONS = {
    "epsilon": ("1", "pi"),
    "pi": ("1", "epsilon"),
    "epsilon*pi": ("1", "epsilon"),
}
SCRIPT_NAME = "Witt_Basis_Sum_Quadratics.py"


@dataclass(frozen=True)
class WittResult:
    valid: bool
    message: str
    alpha: str
    beta: str
    mu1: str
    mu2: str
    minus_one_square: bool
    lam: str
    case: str
    basis: dict[str, str]
    witnesses: dict[str, str]
    raw_output: str


def hilbert(a: str, b: str, minus_one_square: bool) -> int:
    """Mirror the square-class Hilbert table used by the CLI program."""
    if a == "1" or b == "1":
        return 1
    if a == b == "epsilon":
        return 1
    if minus_one_square:
        if a == b:
            return 1
        if {a, b} in (
            {"epsilon", "pi"},
            {"epsilon", "epsilon*pi"},
            {"pi", "epsilon*pi"},
        ):
            return -1
    else:
        table = {
            frozenset(("epsilon", "pi")): -1,
            frozenset(("epsilon", "epsilon*pi")): -1,
            frozenset(("pi",)): -1,
            frozenset(("pi", "epsilon*pi")): 1,
            frozenset(("epsilon*pi",)): -1,
        }
        if frozenset((a, b)) in table:
            return table[frozenset((a, b))]
    raise ValueError(f"Hilbert symbol not covered: {a}, {b}")


def lambda_class(a: str, b: str) -> str:
    if a == b:
        return "1"
    return ({"epsilon", "pi", "epsilon*pi"} - {a, b}).pop()


def valid_parameters(minus_one_square: bool, alpha: str, beta: str, mu1: str, mu2: str) -> bool:
    if alpha not in CLASSES or beta not in CLASSES:
        raise ValueError("Alpha and beta must be nontrivial square classes.")
    if mu1 not in MU_OPTIONS[alpha] or mu2 not in MU_OPTIONS[beta]:
        raise ValueError("Mu parameters must be among the choices for their quadratic extensions.")
    return (
        hilbert(mu1, alpha, minus_one_square) * hilbert(mu2, beta, minus_one_square)
        == hilbert(alpha, beta, minus_one_square)
    )


def _lines(output: str, names: tuple[str, ...]) -> dict[str, str]:
    found = {}
    pattern = re.compile(r"^\s*(" + "|".join(map(re.escape, names)) + r")\s*=\s*(.+)$")
    for line in output.splitlines():
        match = pattern.match(line)
        if match:
            found[match.group(1)] = match.group(2).strip()
    return found


def run_constructor(
    minus_one_square: bool,
    alpha: str,
    beta: str,
    mu1: str,
    mu2: str,
    script_path: Path | None = None,
) -> WittResult:
    lam = lambda_class(alpha, beta)
    if not valid_parameters(minus_one_square, alpha, beta, mu1, mu2):
        return WittResult(
            False,
            "These parameters do not satisfy the split SO(5) Hilbert-symbol criterion: "
            "(mu1, alpha)(mu2, beta) = (alpha, beta). Please change mu1 or mu2.",
            alpha, beta, mu1, mu2, minus_one_square, lam, "", {}, {}, "",
        )

    path = Path(script_path) if script_path else Path(__file__).resolve().parent / SCRIPT_NAME
    if not path.is_file():
        raise FileNotFoundError(
            f"The verified constructor '{SCRIPT_NAME}' is missing. "
            "Upload it beside app.py and web_runner.py."
        )

    responses = [
        "yes" if minus_one_square else "no",
        str(CLASSES.index(alpha) + 1),
        str(CLASSES.index(beta) + 1),
        str(MU_OPTIONS[alpha].index(mu1) + 1),
        str(MU_OPTIONS[beta].index(mu2) + 1),
        "2",  # exit after the first construction
    ]
    completed = subprocess.run(
        [sys.executable, "-u", str(path)],
        input="\n".join(responses) + "\n",
        text=True,
        capture_output=True,
        cwd=str(path.parent),
        timeout=15,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "The original constructor exited with an error:\n"
            + (completed.stderr or completed.stdout)[-2000:]
        )

    output = completed.stdout
    basis = _lines(output, ("e1", "e2", "f1", "f2", "g1"))
    if set(basis) != {"e1", "e2", "f1", "f2", "g1"}:
        raise RuntimeError("The original constructor did not return all five Witt basis vectors.")

    case_match = re.search(r"^Case (4\.[123]):", output, flags=re.MULTILINE)
    case = case_match.group(1) if case_match else "Not identified"
    witnesses = _lines(output, ("r", "s", "u"))
    return WittResult(
        True, "Witt basis constructed by the verified Stage 1 program.",
        alpha, beta, mu1, mu2, minus_one_square, lam,
        case, basis, witnesses, output,
    )
