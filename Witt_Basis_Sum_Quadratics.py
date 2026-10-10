# Square classes in F*/F*2

square_classes = ["1", "epsilon", "pi", "epsilon*pi"]



print("Square classes:")

print(square_classes)



print()

print("==============================================")

print(" WITT BASIS CONSTRUCTOR FOR ELLIPTIC TORI")

print("                  OF SO(5)")

print("==============================================")









def hilbert(a, b):

    # Anything paired with 1 has Hilbert symbol 1

    if a == "1" or b == "1":

        return 1



    # Hilbert symbol is symmetric

    pair = {a, b}



    # (epsilon, epsilon) is always 1

    if a == "epsilon" and b == "epsilon":

        return 1



    if minus_one_square:

        # Table 2.4 when -1 IS a square

        if a == b:

            return 1



        if pair in [

            {"epsilon", "pi"},

            {"epsilon", "epsilon*pi"},

            {"pi", "epsilon*pi"}

        ]:

            return -1



    else:

        # Table 2.4 when -1 is NOT a square

        table = {

            ("epsilon", "pi"): -1,

            ("epsilon", "epsilon*pi"): -1,

            ("pi", "pi"): -1,

            ("pi", "epsilon*pi"): 1,

            ("epsilon*pi", "epsilon*pi"): -1

        }



        if (a, b) in table:

            return table[(a, b)]



        if (b, a) in table:

            return table[(b, a)]



    return None





def multiply_classes(a, b):

    if a == "1":

        return b



    if b == "1":

        return a



    if a == b:

        return "1"



    pair = {a, b}



    if pair == {"epsilon", "pi"}:

        return "epsilon*pi"



    if pair == {"epsilon", "epsilon*pi"}:

        return "pi"



    if pair == {"pi", "epsilon*pi"}:

        return "epsilon"

    return None



def square_factor(alpha, beta):

    if alpha == "1":

        return "1"



    if beta == "1":

        return "1"



    if alpha == beta:

        return f"1/({alpha})"



    pair = {alpha, beta}



    if pair == {"epsilon", "epsilon*pi"}:

        return "1/epsilon"



    if pair == {"pi", "epsilon*pi"}:

        return "1/pi"



    return "1"



print()

print("Please enter the parameters defining the torus.")




def inverse_square_factor(alpha, beta):
    """Return t=1/c for lambda=alpha*beta*c^2, as a readable string."""
    if alpha == "1" or beta == "1":
        return "1"
    if alpha == beta:
        return f"({alpha})"
    pair = {alpha, beta}
    if pair == {"epsilon", "pi"}:
        return "1"
    if pair == {"epsilon", "epsilon*pi"}:
        return "epsilon"
    if pair == {"pi", "epsilon*pi"}:
        return "pi"
    raise ValueError("Unknown square-class pair")

def choose_square_class(name):

    print()

    print("Choose", name)



    # alpha and beta cannot be 1

    if name == "alpha" or name == "beta":



        while True:

            print("1. epsilon")

            print("2. pi")

            print("3. epsilon*pi")



            choice = input("Your choice: ")



            if choice == "1":

                return "epsilon"

            elif choice == "2":

                return "pi"

            elif choice == "3":

                return "epsilon*pi"

            else:

                print("Invalid choice. Please enter 1, 2, or 3.")

                print()



    # mu1, mu2, and lambda can be any square class

    else:



        while True:

            print("1. 1")

            print("2. epsilon")

            print("3. pi")

            print("4. epsilon*pi")



            choice = input("Your choice: ")



            if choice == "1":

                return "1"

            elif choice == "2":

                return "epsilon"

            elif choice == "3":

                return "pi"

            elif choice == "4":

                return "epsilon*pi"

            else:

                print("Invalid choice. Please enter 1, 2, 3, or 4.")

                print()  

  # 2. PUT THE NEW FUNCTION HERE

def choose_mu(name, gamma):

    print()

    print("Choose", name)



    if gamma == "epsilon":

        options = ["1", "pi"]



    elif gamma == "pi":

        options = ["1", "epsilon"]



    elif gamma == "epsilon*pi":

        options = ["1", "epsilon"]



    while True:

        print("1.", options[0])

        print("2.", options[1])



        choice = input("Your choice: ")



        if choice == "1":

            return options[0]

        elif choice == "2":

            return options[1]

        else:

            print("Invalid choice. Please enter 1 or 2.")

            print()

def case_41_u(alpha, mu1, mu2):

    if alpha == "epsilon":

        if minus_one_square:

            return "i", "where i^2 = -1 and i is in F"

        else:

            return (

                "s^2 + i*t^2",

                "where i = sqrt(epsilon) and s^4 + t^4 = -1"

            )



    elif alpha == "pi" or alpha == "epsilon*pi":

        if minus_one_square:

            if mu1 == mu2:

                return "i", "where i^2 = -1 and i is in F"

        else:

            if (

                (mu1 == "epsilon" and mu2 == "1")

                or

                (mu1 == "1" and mu2 == "epsilon")

            ):

                return "1", "since the required norm is 1"



    return None, None

def case_42_plus_witnesses(alpha, beta):

    if minus_one_square:

        return None



    if (

        (alpha == "pi" and beta == "epsilon*pi")

        or

        (alpha == "epsilon*pi" and beta == "pi")

    ):

        return {

            "r": "1",

            "s": "1",

            "c": f"sqrt({alpha})",

            "a": "0",

            "b": f"1/{alpha}"

        }



    return None            

def norm_witness(gamma, target):

    # We want x in F(sqrt(gamma)) such that

    # N_gamma(x) has square class target.



    if target == "1":

        return "1"



    if minus_one_square:

        if target == gamma:

            return f"i*sqrt({gamma})"



    else:

        # Here epsilon = -1.



        if gamma == "epsilon":

            if target == "epsilon":

                return "s0^2 + i*t0^2"



        elif gamma == "pi":

            if target == "epsilon*pi":

                return "sqrt(pi)"



        elif gamma == "epsilon*pi":

            if target == "pi":

                return "sqrt(epsilon*pi)"



    return None

def construct_witt_basis():

    global minus_one_square



    print()

    answer = input("Is -1 a square in F? (yes/no): ").lower()

    minus_one_square = answer == "yes"

    alpha = choose_square_class("alpha")

    beta = choose_square_class("beta")

    while True:

       mu1 = choose_mu("mu1", alpha)

       mu2 = choose_mu("mu2", beta)



       left_side = hilbert(mu1, alpha) * hilbert(mu2, beta)

       right_side = hilbert(alpha, beta)



       if left_side == right_side:

          break



       print()

       print("The selected mu1 and mu2 do NOT satisfy the criterion.")

       print("(mu1, alpha)(mu2, beta) =", left_side)

       print("(alpha, beta) =", right_side)

       print("Please choose mu1 and mu2 again.")



    lam = multiply_classes(alpha, beta)

    d = square_factor(alpha, beta)

    print("Square factor d =", d)

    lambda_required = lam



    print()

    print("Checking the criterion...")

    print("lambda is determined automatically as alpha*beta:", lam)

    print("(mu1, alpha)(mu2, beta) =", left_side)

    print("(alpha, beta) =", right_side)



    print("The criterion is satisfied.")







    print()

    print("Determining the Witt basis case...")

    print("DEBUG alpha =", alpha)

    print("DEBUG beta =", beta)

    print("DEBUG alpha == beta =", alpha == beta)



    if alpha == beta:

        print("Case 4.1: alpha = beta")

        print("Use Proposition 4.1")

        u, u_condition = case_41_u(alpha, mu1, mu2)



        if u is not None:

            print()

            print("Explicit norm witness:")

            print(f"u = {u}")

            print(u_condition)

            print()

            print(f"Choose u in F(sqrt({alpha})) satisfying:")

            print(f"N_{alpha}(u) = -{mu1}/{mu2}")



            print()

            print("The Witt basis is:")

            print("e1 = (1, u, 0)")

            print(f"f1 = (1/(2*{mu1})) * (1, -u, 0)")

            print(f"e2 = (sqrt({alpha}), sqrt({alpha})*u, 0)")

            print(

                f"f2 = (1/(2*{alpha}*{mu1})) * "

                f"(-sqrt({alpha}), sqrt({alpha})*u, 0)"

            )

            print("g1 = (0, 0, 1)")





    else:

        ab = hilbert(alpha, beta)

        sign1 = hilbert(mu1, alpha)

        sign2 = hilbert(mu2, beta)



        print("alpha and beta are distinct.")

        print("(alpha, beta) =", ab)

        print("(mu1, alpha) =", sign1)

        print("(mu2, beta) =", sign2)



        if ab == 1:

            print("Case 4.2: (alpha, beta) = 1")



            if sign1 == 1 and sign2 == 1:



                print("Subcase: (+1, +1)")

                print("Use Proposition 4.3")

                witnesses = case_42_plus_witnesses(alpha, beta)



                if witnesses is not None:

                    print()

                    print("Explicit witnesses:")

                    print(f"r = {witnesses['r']}")

                    print(f"s = {witnesses['s']}")

                    print(f"c = {witnesses['c']}")

                    print(f"a = {witnesses['a']}")

                    print(f"b = {witnesses['b']}")



                    print()

                    print("The explicit Witt basis is:")



                    if alpha == "pi":

                        print("e1 = (1, 0, 1)")

                        print("f1 = (1/2, 0, -1/2)")

                        print("e2 = (sqrt(pi), sqrt(epsilon*pi), 0)")

                        print("f2 = (-1/(2*sqrt(pi)), sqrt(epsilon*pi)/(2*pi), 0)")

                        print("g1 = (0, 1, 0)")



                    elif alpha == "epsilon*pi":

                        print("e1 = (1, 0, 1)")

                        print("f1 = (1/2, 0, -1/2)")

                        print("e2 = (sqrt(epsilon*pi), sqrt(pi), 0)")

                        print("f2 = (-1/(2*sqrt(epsilon*pi)), sqrt(pi)/(2*(epsilon*pi)), 0)")

                        print("g1 = (0, 1, 0)")



            elif sign1 == -1 and sign2 == -1:

                print("Subcase: (-1, -1)")

                print("Use Proposition 4.7")



                print()

                print("This is the exceptional case with -1 not a square in F.")



                print()

                print("Choose s, t in F satisfying:")

                print("s^4 + t^4 = -1")



                print()



                if alpha == "pi" and beta == "epsilon*pi":

                    print("Canonical ordering: alpha = pi, beta = epsilon*pi")

                    print()

                    print("The Witt basis is:")

                    print("e1 = (sqrt(pi), i*sqrt(pi), 0)")

                    print("f1 = (1/(2*sqrt(pi)), -i/(2*sqrt(pi)), 0)")

                    print("e2 = (t^2, -s^2, 1)")

                    print("f2 = (t^2/2, -s^2/2, -1/2)")

                    print("g1 = (s^2, t^2, 0)")



                elif alpha == "epsilon*pi" and beta == "pi":

                    print("Reversed ordering: alpha = epsilon*pi, beta = pi")

                    print()

                    print("The Witt basis is:")

                    print("e1 = (i*sqrt(pi), sqrt(pi), 0)")

                    print("f1 = (-i/(2*sqrt(pi)), 1/(2*sqrt(pi)), 0)")

                    print("e2 = (-s^2, t^2, 1)")

                    print("f2 = (-s^2/2, t^2/2, -1/2)")

                    print("g1 = (t^2, s^2, 0)")

        elif ab == -1:
            print("Case 4.3: (alpha, beta) = -1")

            if (sign1, sign2) in ((1, -1), (-1, 1)):
                swapped = (sign1 == -1)
                print("Subcase:", "(-1, +1)" if swapped else "(+1, -1)")
                print("Use Proposition 4.8")

                # Proposition 4.8 always uses (+1, -1) after any swap.
                A, B = (beta, alpha) if swapped else (alpha, beta)
                M1, M2 = (mu2, mu1) if swapped else (mu1, mu2)

                # The chosen lambda is a square-class representative.
                # lambda = A*B*c^2; t = 1/c is used in the z-coordinate.
                c = square_factor(A, B)
                t = inverse_square_factor(A, B)

                # An exact norm equation is stronger than a square-class
                # equation. Never pass a square class to norm_witness() and
                # claim that the returned element solves the exact equation.
                r = "1" if M1 == "1" else "r"
                # N_B(1) = 1, so s = 1 whenever M2 = A exactly.
                s = "1" if M2 == A else "s"
                exact_special = (
                    minus_one_square
                    and A == "pi" and B == "epsilon*pi"
                    and M1 == "1" and M2 == "epsilon"
                )
                reverse_special = (
                    minus_one_square
                    and A == "epsilon*pi" and B == "pi"
                    and M1 == "1" and M2 == "epsilon"
                )
                epsilon_special = (
                    minus_one_square
                    and A == "epsilon*pi" and B == "epsilon"
                    and M1 == "1" and M2 == "pi"
                )
                unit_witness_special = (
                    A == "epsilon" and B == "epsilon*pi"
                    and M1 == "1" and M2 == "epsilon"
                )
                quartic_witness_special = (
                    not minus_one_square
                    and A == "epsilon*pi" and B == "epsilon"
                    and M1 == "1" and M2 == "pi"
                )
                if quartic_witness_special:
                    # Lemma 2.1.5: epsilon = -1, i^2 = -1,
                    # s0^4 + t0^4 = -1. Hence N_epsilon(s0^2+i*t0^2)
                    # = s0^4+t0^4 = epsilon, so pi*N_epsilon(s) = A.
                    s = "s0**2 + i*t0**2"
                elif unit_witness_special:
                    # M1*N_A(1) = 1 and M2*N_B(1) = epsilon = A.
                    s = "1"
                elif exact_special:
                    s = "i*sqrt(epsilon*pi)/epsilon"
                elif reverse_special:
                    s = "i*sqrt(pi)"
                elif epsilon_special:
                    # N_epsilon(i*sqrt(epsilon)) = epsilon, since i^2 = -1.
                    # Therefore pi*N_epsilon(s) = epsilon*pi = A exactly.
                    s = "i*sqrt(epsilon)"

                if swapped:
                    print("Interchanging the two quadratic summands.")
                    print(f"After interchange: alpha = {A}, beta = {B}, "
                          f"mu1 = {M1}, mu2 = {M2}")

                print("Required exact norm equations:")
                print(f"{M1} * N_{A}(r) = 1")
                print(f"{M2} * N_{B}(s) = {A}")
                print(f"lambda = {lam} = ({A})*({B})*({c})^2")
                if r == "r" or s == "s":
                    print("Symbolic witnesses: choose r and s satisfying "
                          "the displayed exact norm equations.")
                    print("Their existence follows from the Hilbert-symbol criterion;")
                    print("their explicit values depend on the field F.")
                else:
                    print("Exact witnesses:")
                print(f"r = {r}")
                print(f"s = {s}")
                if quartic_witness_special:
                    print("Lemma 2.1.5: epsilon = -1, i^2 = -1, "
                          "s0^4 + t0^4 = -1")
                    print("N_epsilon(s) = s0^4 + t0^4 = epsilon")

                print("The Witt basis in the original coordinate order is:")
                if not swapped:
                    print(f"e1 = (({r})*sqrt({A}), ({s}), 0)")
                    print(f"f1 = (-({r})*sqrt({A})/(2*({A})), "
                          f"({s})/(2*({A})), 0)")
                    if exact_special or reverse_special:
                        # Rescale the second hyperbolic pair to remove pi.
                        print("e2 = (0, i, 1)")
                        print("f2 = (0, -i/(2*epsilon), 1/(2*epsilon))")
                    else:
                        print(f"e2 = (0, ({s})*sqrt({B}), ({t}))")
                        print(f"f2 = (0, -({s})*sqrt({B})/(2*({A})*({B})), "
                              f"({t})/(2*({A})*({B})))")
                    print(f"g1 = ({r}, 0, 0)")
                else:
                    print(f"e1 = (({s}), ({r})*sqrt({A}), 0)")
                    print(f"f1 = (({s})/(2*({A})), "
                          f"-({r})*sqrt({A})/(2*({A})), 0)")
                    if exact_special or reverse_special:
                        print("e2 = (i, 0, 1)")
                        print("f2 = (-i/(2*epsilon), 0, 1/(2*epsilon))")
                    else:
                        print(f"e2 = (({s})*sqrt({B}), 0, ({t}))")
                        print(f"f2 = (-({s})*sqrt({B})/(2*({A})*({B})), "
                              f"0, ({t})/(2*({A})*({B})))")
                    print(f"g1 = (0, {r}, 0)")
            else:
                print("Unsupported sign combination for Proposition 4.8")

while True:

    construct_witt_basis()



    print()

    print("----------------------------------------------")

    print("Would you like to construct another Witt basis?")

    print("1. Yes")

    print("2. Exit")



    choice = input("Your choice: ")



    if choice == "2":

        print()

        print("Program ended.")

        break
