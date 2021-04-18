from math import acos, pi

def get_kF (mu):
    return acos(-0.5*mu)

mu = -1.9
kF = get_kF (mu)
print (kF / pi)
