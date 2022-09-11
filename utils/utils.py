import numpy as np


def roulette_wheel_pop(pop, p, number):
    chosen=np.random.choice(len(pop),number,False,p)
    chosen = [pop[a] for a in chosen]
    return chosen
