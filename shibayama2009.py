import sympy
import Hamilton
import Birkhoff

x = sympy.IndexedBase("x")
y = sympy.IndexedBase("y")
p = sympy.IndexedBase("p")
q = sympy.IndexedBase("q")
t = sympy.IndexedBase("tau")
o = sympy.Symbol("omega", positive=True)
a = sympy.Symbol("alpha", positive=True)
b = sympy.Symbol("beta")
l = sympy.Symbol("lambda", positive=True)

H = p[1]*p[1]/4 + o*o/4/q[1]/q[1] - 1/a/2/q[1] - 2/sympy.sqrt(q[1]*q[1]+q[2]*q[2]) + (a+2)*p[2]*p[2]/a/4
ham = Hamilton.Hamiltonian(H, [p[1],p[2]], [q[1],q[2]])
equilibrium_points = [0, 0, o*o*a/(4*a+1), 0]
ham.expand_around_equilibrium(equilibrium_points, max_degree=4)
ham.rescale()
ham.coeff_subs([(a, (8-l*l)/(4*l*l-4))])
ham.rotate45()

birkhoff = Birkhoff.LieTransform.fromHamiltonian(ham)
birkhoff.exec()
print(birkhoff.normalform())
