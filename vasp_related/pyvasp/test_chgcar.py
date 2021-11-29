import chgcar

f = chgcar.chgcar().read('Al/CHGCAR')
print f.get_nions()
f.get_charges()
