import sympy

class LieTransform:
    ps = []
    qs = []
    variables = []
    new_variables = []
    dim = 0
    dims = 0 # dims = 2*dim
    hamiltonian = 0
    mat_hamiltonian = [[]]
    frequency = []
    generator_list= [0]
    generator_function = 0
    max_degree = 0
    normalform_flag = False
    normalform_calculated = False
    normalform_function = 0 

    def __const_sum_list__(self, num, max_sum, min_sum = 0):
        ans = {0: [tuple([0 for _ in range(num)])]}
        for i in range(1, max_sum+1):
            prev_list = ans[i-1]
            new_set = set()
            for p in prev_list:
                for j in range(num):
                    tmp_p = [(p[k] if k != j else p[k]+1) for k in range(num)]
                    new_set.add(tuple(tmp_p))
            new_list = list(new_set)
            ans[i] = new_list
        if min_sum > 0:
            for i in range(min_sum):
                ans.pop(i)
        return ans

    def __const_sum__(self, num, s):
        return self.__const_sum_list__(num, s, s)[s]

    def __deg2monomial__(self, deg_list):
        monomial = 1
        for i in range(self.dims):
            monomial *= self.variables[i]**deg_list[i]
        return monomial
    
    def __decomposeH__(self, sympoly, deg_lists):
        self.mat_hamiltonian = [[0 for _ in range(self.max_degree-d-1)] for d in range(self.max_degree-1)]
        for deg_list in deg_lists:
            d = 0
            for i in range(self.dims):
                d += deg_list[i]
            if 2 <= d and d <= self.max_degree:
                self.mat_hamiltonian[0][d-2] += sympoly.coeff_monomial(deg_list)*self.__deg2monomial__(deg_list)*sympy.factorial(d-2)
        self.frequency = [0 for _ in range(self.dim)]
        for i in range(self.dim):
            self.frequency[i] = self.mat_hamiltonian[0][0].coeff(self.variables[i],1).coeff(self.variables[i+self.dim],1)
        return 0

    def __PoissonBracket__(self, F, G):
        eq = 0
        for i in range(self.dim):
            eq += sympy.diff(F, self.ps[i])*sympy.diff(G, self.qs[i]) - sympy.diff(F, self.qs[i])*sympy.diff(G, self.ps[i])
        return eq

    def __PoissonBracketWithH00__(self, deg_list):
        ans = 0
        for i in range(self.dim):
            ans += self.frequency[i] * (deg_list[i+self.dim] - deg_list[i])
        return ans

    def __solveW__(self, n):
        eq = self.mat_hamiltonian[0][n]
        if n>1:
            for i in range(1,n):
                for k in range(n-i+1):
                    eq += sympy.binomial(n-k-1, i-1)*self.__PoissonBracket__(self.mat_hamiltonian[k][n-i-k], self.generator_list[i])
        eq = sympy.Poly(eq, *self.variables)
        deg_lists = eq.monoms()
        ans = 0
        for deg_list in deg_lists:
            # check if degredd of ps and qs is equal
            # if equal, continue
            isSame = True
            for i in range(self.dim):
                if deg_list[i] != deg_list[i+self.dim]:
                    isSame = False
            if not isSame:
                tmp_coeff = - eq.coeff_monomial(deg_list)
                tmp_coeff /= self.__PoissonBracketWithH00__(deg_list)
                ans += tmp_coeff*self.__deg2monomial__(deg_list)
        return ans

    def __calcH__(self, i, j):
        eq = self.mat_hamiltonian[i-1][j+1]
        for k in range(j+1):
            eq += sympy.binomial(j, k)*self.__PoissonBracket__(self.mat_hamiltonian[i-1][j-k], self.generator_list[k+1])
        return eq

    def __init__(self, H, ps, qs):
        self.ps = ps.copy()
        self.qs = qs.copy()
        self.variables = [*ps, *qs]
        self.dim = len(ps)
        self.dims = 2*self.dim
        t = sympy.IndexedBase("tau")
        self.new_variables = [t[i+1] for i in range(self.dim)]
        self.hamiltonian = H
        sympolyH = sympy.Poly(H, *self.variables)
        self.max_degree = sympolyH.degree()
        deg_lists = sympolyH.monoms()
        self.__decomposeH__(sympolyH, deg_lists)
        print("Initialized")
        print("do exec() and calculate normal form")

    @classmethod
    def fromHamiltonian(cls, hamiltonian):
        ham = hamiltonian.approx_function
        ps = hamiltonian.ps
        qs = hamiltonian.qs
        return cls(ham, ps, qs)


    def __del__(self):
        print("DEL!")

    def exec(self):
        for d in range(1, self.max_degree-1):
            self.generator_list.append(self.__solveW__(d))
            for i in range(1, d+1):
                self.mat_hamiltonian[i][d-i] = sympy.simplify(self.__calcH__(i, d-i))
        self.normalform_flag = True
        print("Lie transform completed!!")
        print("do normalform() to print normal form of input Hamiltonian")
    
    def normalform(self, new_var = None):
        if new_var is None:
            new_var = self.new_variables
        if len(new_var) != self.dim:
            print("new variable length not correct!!")
            return 0
        if not self.normalform_flag:
            print("Lie transform is not yet done!!")
            print("do exec() to execute Lie transform")
            return 0
        if not self.normalform_calculated:
            self.normalform_function = 0
            for d in range(self.max_degree-1):
                eq = sympy.Poly(self.mat_hamiltonian[d][0], *self.variables)
                deg_lists = eq.monoms()
                for deg_list in deg_lists:
                    monomial = 1
                    for i in range(self.dim):
                        monomial *= (sympy.I*new_var[i])**deg_list[i]
                    self.normalform_function += sympy.factor(eq.coeff_monomial(deg_list))*monomial/sympy.factorial(d)
            self.normalform_calculated = True    
        return self.normalform_function

    def __getitem__(self,ind):
        return self.mat_hamiltonian[ind]