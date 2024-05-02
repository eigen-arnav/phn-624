def oscillator_strength(energy):
    # energy in MeV 
    return 197.33/(940*energy)**.5

A = 36 
b1 = oscillator_strength(41*A**(-1/3))
b2 = oscillator_strength(45*A**(-1/3) - 25*A**(-2/3))
print(b1) ; print(b2) ; print(b1-b2)
