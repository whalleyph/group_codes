import utils, poscar

pos = poscar.poscar('test_files/POSCAR', types='C Fe1 Fe2 M'.split())

centered_poscar = utils.center_around_333 (pos, point=pos.coord[15])

for i in centered_poscar[0]:
  print i
