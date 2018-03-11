# -*- coding: utf-8 -*-
"""
entropy.py - functions for boltzmann probabilities and shannon entropies
"""

import bigfloat

def compute_entropy(probabilities):
    # Compute Shannon Entropy of the given distribution
    # https://en.wikipedia.org/wiki/Entropy_(information_theory)
    # SE = - ∑ P_i * log_b( P_i )
    # where
    #  P_i is the probability of sample i
    #  log_b is log-base-b
    # It is not clear in the Pereyra paper what units were used, so we
    # try the two most common bases: 2 and e
    entropy = 0.0
    for probability in probabilities:
        # base e: entropy measured in nats
        entropy -= probability * bigfloat.log(probability)
    return entropy

def compute_boltzmann(energies, debug = False):
    # create Boltzmann distribution
    # https://en.wikipedia.org/wiki/Boltzmann_distribution
    # p_i = exp( -E_i / kT ) / ∑ exp( -E_j / kT )
    # where:
    #    p_i is the probability of state i occuring
    #    E_i is the energy of state i
    #    E_j is the energy of state j, where the ∑ in the denominator
    #    iterates over all states j
    # first we calculate exp( -E_i / kT ) for all states
    if len(energies) == 0:
        return []
    if debug:
        print("Boltzmann Distribution")
    divisor = bigfloat.BigFloat.exact(0.0)
    for energy in energies:
        ep = bigfloat.exp(-energy)
        if debug:
            print("energy, ep = %g, %g" % (energy, ep))
        # divisor = ∑ exp( -E_j / kT )
        divisor += ep
    probabilities = []
    for energy in energies:
        # p_i = exp( -E_i / kT ) / divisor
        numerator = bigfloat.exp(-energy)
        probability = numerator / divisor
        # save probability to dictionary
        probabilities.append(float(probability))
        if debug:
            print("%s / %s = %s" % (numerator, divisor, probability))

    return probabilities
