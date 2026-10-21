def test_distinct_quadratics_witnesses_match_example(monkeypatch):
    import Witt_Basis_Sum_Quadratics as w

    # Over Q_3, -1 is not a square.
    monkeypatch.setattr(w, "minus_one_square", False, raising=False)

    witnesses = w.case_42_plus_witnesses("pi", "epsilon*pi")

    assert witnesses == {
        "r": "1",
        "s": "1",
        "c": "sqrt(pi)",
        "a": "0",
        "b": "1/pi",
    }
