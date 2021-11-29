import poscar

pos = poscar.poscar().read('Al/POSCAR')
print pos.get_spacegroup()
