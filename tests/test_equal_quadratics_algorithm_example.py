import sympy as sp


def test_equal_quadratics_norm_witness(monkeypatch):
    import Witt_Basis_Sum_Quadratics as w

    # Over Q_3, -1 is not a square.
    monkeypatch.setattr(w, "minus_one_square", False, raising=False)

    u, condition = w.case_41_u("epsilon", "1", "1")

    assert u == "s^2 + i*t^2"
    assert "s^4 + t^4 = -1" in condition

    # Interpret i^2 = -1 symbolically.
    s, t = sp.symbols("s t", real=True)
    expression = sp.sympify(
        u.replace("^", "**"),
        locals={"s": s, "t": t, "i": sp.I},
    )

    norm = sp.expand(expression * sp.conjugate(expression))

    assert sp.simplify(norm - (s**4 + t**4)) == 0
    assert sp.simplify(norm.subs(t**4, -1 - s**4) + 1) == 0
