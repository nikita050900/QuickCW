"""

This will house small calculators and other useful utility functions.

"""

import argparse
import math

def calc_mass_ratio(mc,mt):

"""

Code to save you from ever again having to do the quadratic equation
while converting a pair of Mc and Mtot to a mass ratio.

Go in peace.

"""
    term1 = math.sqrt(
        ((mt**4) * ((mc/mt)**(2/3)))
        -
        (4* (mc**2) * (mt**2) * (mc/mt)**(1/3))
    )

    term2 = 2* (mc**2)

    term3 = (mt**2) * ((mc/mt)**(1/3))

    term4 = 2 * (mc**2)

    qplus = (term1 - term2 + term3)/term4

    qminus = (-1*term1 - term2 + term3)/term4

    if (mc/mt > 0.435275282):
        print("\n\t------- ERROR -------\nThe inputs you provided imply a mass ratio of q > 1\n(in other words, Mc is too big compared to Mtot).\n\nRemember q is defined 0 < q < 1. Please check your values and try again.\n\nIf you are inputing the results from a CW upper-limit run, this result implies your results are not constraining on q.\n\t---------------------\n")
        exit()

    if (qplus>qminus):
        q = qminus
        multiplier = qplus
    else:
        q = qplus
        multiplier = qminus
        

    return q,multiplier




if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Calculate the q value for a Mc, Mtot pair. q is defined 0 < q < 1. Mc and Mtot must be given in log(solar mass), i.e. 5 x 10^9 Msun would be given as 9.69897 .")
    parser.add_argument('--mc', '-c', dest='mc', type=float, default=None,
                        help="log(Chirp mass in solar mass)")
    parser.add_argument('--mtot', '-t', dest='mt', type=float, default=None,
                        help="log(Total mass in solar mass)")
    args = parser.parse_args()

    if args.mc == None or args.mt == None:
        print("\n\t------- ERROR -------\nYou must supply both Mc and Mtot values as -c and -t, respectively.\n\t---------------------\n")
        exit()

    mc = 10**(args.mc)
    mt = 10**(args.mt)

    q,multiplier = calc_mass_ratio(mc,mt)

    print("Your q = "+str(q))
    print("This implies one mass "+str(multiplier)+" times the mass of the other.")

