def mu4_options_for_case(quartic_case, subfield):
    """Return the two mu4 choices: admissible first, rejected second."""

    if quartic_case == "Exceptional twisted: d = eta*pi":
        return ("eta", "1")

    if subfield in ("F(sqrt(epsilon))", "F(sqrt(-1))"):
        return ("1", "varpi")

    return ("1", "epsilon")


def is_admissible_mu4(quartic_case, subfield, mu4_choice):
    """Check the mu4 admissibility rules used by the app."""

    admissible, rejected = mu4_options_for_case(
        quartic_case, subfield
    )

    return mu4_choice == admissible
