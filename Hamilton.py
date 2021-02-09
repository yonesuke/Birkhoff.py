import sympy

class Hamiltonian:
    hamiltonian = 0
    ps = []
    qs = []
    variables = []
    dim = 0
    dims = 0 # dims = 2*dim
    approx_function = 0
    approx_function_set = {}
    isExpanded = False
    equilibrium_points = []
    rescaling_factor = []
    diff_set = {}

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

    def __deg2monomial__(self, deg_list):
        monomial = 1
        for i in range(self.dims):
            monomial *= self.variables[i]**deg_list[i]
        return monomial

    def __init__(self, hamiltonian, ps, qs):
        self.hamiltonian = hamiltonian
        self.ps = ps.copy()
        self.qs = qs.copy()
        self.variables = [*ps, *qs]
        self.dim = len(ps)
        self.dims = 2*self.dim
        self.equilibrium_points = [0 for _ in range(self.dims)]
        self.rescaling_factor = [1 for _ in range(self.dims)]
        self.diff_set[0] = [(tuple([0 for _ in range(self.dims)]), hamiltonian)]

    def expand_around_equilibrium(self, equilibrium_points, max_degree = 5):
        self.equilibrium_points = equilibrium_points.copy()
        degree_set = self.__const_sum_list__(self.dims, max_degree)
        for d in range(1, max_degree+1):
            degree_lists = degree_set[d]
            prev_diff_set = dict(self.diff_set[d-1])
            diff_list = []
            for degree_list in degree_lists:
                ind = 0
                while degree_list[ind]==0:
                    ind += 1
                prev_degree = tuple([(degree_list[k]-1 if k==ind else degree_list[k]) for k in range(self.dims)])
                prev_diff = prev_diff_set[prev_degree]
                diff_func = sympy.simplify(sympy.diff(prev_diff, self.variables[ind]))
                diff_list.append((degree_list, diff_func))
            self.diff_set[d] = diff_list
        subs_point = [(self.variables[i], self.equilibrium_points[i]) for i in range(self.dims)]
        self.approx_function = 0
        self.approx_function_set = {}
        for d in range(max_degree+1):
            for diffs in self.diff_set[d]:
                degree_list, diff_func = diffs
                tmp_coeff = sympy.simplify(diff_func.subs(subs_point))
                for i in range(self.dims):
                    tmp_coeff /= sympy.factorial(degree_list[i])
                tmp_term = tmp_coeff*self.__deg2monomial__(degree_list)
                self.approx_function += tmp_term
                if tmp_coeff != 0:
                    self.approx_function_set[degree_list] = tmp_coeff
        self.isExpanded = True
        return self.approx_function

    def rescale(self):
        if not self.isExpanded:
            print("Expand around an equilibrium point, then rescale!!")
            return 0
        for i in range(self.dim):
            # get coeffs of p1^2 and q1^2 and determine rescaling factor
            pi_coeff = self.approx_function_set[tuple([(2 if j==i else 0) for j in range(self.dims)])]
            qi_coeff = self.approx_function_set[tuple([(2 if j==i+self.dim else 0) for j in range(self.dims)])]
            r4 = pi_coeff/qi_coeff
            r = sympy.sqrt(sympy.sqrt(r4))
            self.rescaling_factor[i] = 1/r
            self.rescaling_factor[i+self.dim] = r
        self.approx_function = 0
        for degree_list in self.approx_function_set.keys():
            tmp_coeff = self.approx_function_set[degree_list]
            for i in range(self.dims):
                tmp_coeff *= self.rescaling_factor[i]**degree_list[i]
            self.approx_function_set[degree_list] = tmp_coeff
            self.approx_function += tmp_coeff*self.__deg2monomial__(degree_list)
        return self.approx_function

    def rotate45(self):
        rotate = [
            *[(self.ps[i]-sympy.I*self.qs[i])/sympy.sqrt(2) for i in range(self.dim)],
            *[(-sympy.I*self.ps[i]+self.qs[i])/sympy.sqrt(2) for i in range(self.dim)]
        ]
        eq = 0
        for degree_list in self.approx_function_set.keys():
            tmp_term = self.approx_function_set[degree_list]
            for i in range(self.dims):
                tmp_term *= sympy.expand(rotate[i]**degree_list[i])
            eq += tmp_term
        new_func_set = {}
        self.approx_function = 0
        eq = sympy.Poly(eq, *self.variables)
        degree_lists = eq.monoms()
        for degree_list in degree_lists:
            tmp_coeff = eq.coeff_monomial(degree_list)
            tmp_coeff = sympy.factor(sympy.expand(tmp_coeff))
            new_func_set[degree_list] = tmp_coeff
            self.approx_function += tmp_coeff * self.__deg2monomial__(degree_list)
        self.approx_function_set = new_func_set.copy()
        return self.approx_function

    def coeff_subs(self, subs_list):
        self.approx_function = 0
        for degree_list in self.approx_function_set.keys():
            tmp_coeff = self.approx_function_set[degree_list]
            tmp_coeff = sympy.factor(tmp_coeff.subs(subs_list))
            self.approx_function_set[degree_list] = tmp_coeff
            self.approx_function += tmp_coeff * self.__deg2monomial__(degree_list)
        return self.approx_function

