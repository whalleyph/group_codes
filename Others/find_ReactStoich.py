import sympy
from sympy import Symbol
from molmass import Formula

x = Symbol('x')

Reactant1 = 'LiVO2'
Reactant2 = 'Li3NbO4'
Product = 'Li4NbVO6'

r1 = Formula(Reactant1)
r2 = Formula(Reactant2)
p = Formula(Product)

my_list = []
comp_list = []
for n in range(1,10,1):
    cal = sympy.solve(n*r1.mass*x+n*(1.0-x)*r2.mass-p.mass,"x")
    z = float(format(cal[0], '.2f'))
    if (cal[0] >= 0 and z >= 0.0 and z <= 1.0):
        val = n*r1.mass*z+n*(1.0-z)*r2.mass-p.mass
        my_list.append(val)
        comp_list.append(z)

if len(my_list) > 1:
    min_value = min([i for i in my_list if i >= 0])
    min_index = my_list.index(min_value)
    comp = comp_list[min_index]
    print ('The value of x in xLiVO2 + (1-x)Li3NbO4 = {:s}'.format(Product), 'is: ' "{:.2f}".format(comp))
    
else:
    print ('The value of x in xLiVO2 + (1-x)Li3NbO4 = {:s}'.format(Product), 'is: ' "{:.2f}".format(min(comp_list)))
