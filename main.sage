class OurCurves(object):
    def __init__(self, elt):
        self._normalized = {}
        self._normalization = {}
        self._reduction_type = {}
        if type(elt) == list:
            A,B,C,a,b,c = elt
            cm = sage.structure.element.get_coercion_model()
            K = cm.common_parent(A,B,C,a,b,c)
            K = K.fraction_field()
            R.<x,y,z> = PolynomialRing(K)
            self._equation = A*x^4 + B*y^4 + C*z^4 + a*y^2*z^2 + b*x^2*z^2 + c*x^2*y^2
            self._coeffs = {'A':A, 'B':B, 'C':C, 'a':a, 'b':b, 'c':c}
            self._K = K
        else: #given elt ins an equation
            R = elt.parent()
            K = R.base_ring()
            K = K.fraction_field()
            x,y,z = R.gens()
            A = elt.coefficient({x:4})
            B = elt.coefficient({y:4})
            C = elt.coefficient({z:4})
            a = elt.coefficient({y:2, z:2})
            b = elt.coefficient({x:2, z:2})
            c = elt.coefficient({x:2, y:2})
            self._equation = elt
            self._coeffs = {'A':A, 'B':B, 'C':C, 'a':a, 'b':b, 'c':c}
            self._K = K

    def __repr__(self):
        return str(self._equation)

    def conic(self):
        try:
            return self._conic
        except AttributeError:
            K = self.K
            R.<u,v,w> = PolynomialRing(K)
            self._conic = A*u^2 + B*v^2 + C*w^2 + a*v*w + b*u*w + c*u*v
            return self._conic

    def invariants(self):
        ## Compute the invariants of the given model
        try:
            return self._invariants
        except AttributeError:
            A = self._coeffs['A']
            B = self._coeffs['B']
            C = self._coeffs['C']
            a = self._coeffs['a']
            b = self._coeffs['b']
            c = self._coeffs['c']
            Da = a^2 - 4*B*C
            Db = b^2 - 4*A*C
            Dc = c^2 - 4*A*B
            I3 = A*B*C
            II3 = A*Da + B*Db + C*Dc
            III3 = - 4*A*B*C + A*a^2 + B*b^2 + C*c^2 - a*b*c #Delta X #a*b*c == - III3 + II3 + 8*I3
            I6 = Da*Db*Dc
            I = A*B*Da*Db + A*C*Da*Dc + B*C*Db*Dc
            self._invariants = [I3, II3, III3, I6, I]
            return self._invariants

    def is_normalized(self, p):
        ##return if the given model is normalized wrt the valuation associated to p as defined in Prop 3.1
        try:
            return self._normalized[p]
        except KeyError:
            A = self._coeffs['A']
            B = self._coeffs['B']
            C = self._coeffs['C']
            a = self._coeffs['a']
            b = self._coeffs['b']
            c = self._coeffs['c']
            K = self._K
            v = K.valuation(p)
            for L in [[A, B, c], [A, b, C], [a, B, C],[A, b, c], [a, B, c], [a, b, C]]:
                if all(map(lambda i : v(i) > 0, L)):
                    self._normalized[p] = False
                    return False
            self._normalized[p] = True
            return True
    
    def normalized_model(self, p):
        ##given a model of the curve, returns an isomorphic model that is normalized wrt the valuation
        # associated to p as defined in Prop 3.1
        if self.is_normalized(p):
            return self
        else:
            try:
                return self._normalization[p]
            except KeyError:
                D = self._coeffs
                K = self._K
                v = K.valuation(p)
                F = self._equation
                for L in [['A', 'b', 'c', 0], ['a', 'B', 'c', 1], ['a', 'b', 'C', 2]]:
                    if all(map(lambda i : v(D[i]) > 0, L[:-1])):
                        var = list(F.parent().gens())
                        val = list(map(lambda i : v(D[i])/2, L[:-1]))
                        val[L[-1]] /= 2
                        r = min(val)
                        try:
                            elt = v.element_with_valuation(r)
                        except ValueError:
                            P.<X> = PolynomialRing(p.parent())
                            den = r.denominator()
                            K.<m> = K.extension(X^den - p, 'aux').absolute_field()
                            v = K.valuation(p)
                            elt = v.element_with_valuation(r)
                        var[L[-1]] /= elt
                        F = F(var)
                        D = coeffs(F)
                        assert not all(map(lambda i : v(D[i])  > 0, L[:-1]))
                for L in [['A', 'B', 'c', 2], ['A', 'b', 'C', 1], ['a', 'B', 'C', 0]]:
                    if all(map(lambda i : v(D[i]) > 0, L[:-1])):
                        var = list(F.parent().gens())
                        val = list(map(lambda i : v(D[i])/4, L[:-1]))
                        val[L[-1]] *= 2
                        r = min(val)
                        try:
                            elt = v.element_with_valuation(r)
                        except ValueError:
                            den = r.denominator()
                            K.<m> = K.extension(X^den - p, 'aux').absolute_field()
                            v = K.valuation(p)
                            elt = v.element_with_valuation(r)
                        for i in range(3):
                            if i == L[-1]:
                                var[i] *= elt
                            else:
                                var[i] /= elt
                        F = F(var)
                        [A,B,C,a,b,c] = coeffs(F)
                        print(L)
                        assert not all(map(lambda i : v(D[i])  > 0, L[:-1]))
                self._normalization[p] = OurCurves(F)
                return self._normalization[p]

    def reduction_type(self, p):
        ##Using the decisional tree, compute the reduction type of the curve, and return the corresponding
        try:
            return self._reduction_type[p]
        except KeyError:
            [I3, II3, III3, I6, I] = self.invariants()
            K = self._K
            v = K.valuation(p)
            L = [v(I3/II3), v(II3/III3), v(I6/III3^2)]
            if min(L)>=0:
                if L[2] == 0:
                    if L[0] == 0:
                        return "Good reduction"
                    else:
                        if L[1] == 0:
                            if v(I/III3^2) == 0:
                                return "3.8 (d)"
                            else:
                                return "3.8 (g)"
                        else:
                            return "3.8 (h)"
                else:
                    if L[0] == 0:
                        if v(I/III3^2) == 0:
                            return "3.8 (a)"
                        else:
                            if L[1] == 0:
                                return "3.8 (b)"
                            else:
                                return "3.8 (c)"
                    else:
                        if v(I/III3^2) == 0:
                            return "3.8 (e)"
                        else:
                            if v(I^2/(I6*I3*III3)) > 0:
                                x = v(I6/(I3*II3))
                                if x == 0:
                                    return "3.8 (f.iii)"
                                elif x > 0:
                                    return "3.8 (f.i)"
                                else:
                                    return "3.8 (f.ii)"
                            else:
                                x = v(I6/I)
                                if x == 0:
                                    if v(I/(I3*III3)) == 0:
                                        return "3.8 (f.iii)"
                                    else: #TODO check that this is only negative
                                        return "3.8 (f.vi)"
                                elif x > 0:
                                    y = v(I/(I3*III3))
                                    if y == 0:
                                        return "3.8 (f.v)"
                                    elif y > 0:
                                        return "3.8 (f.i)"
                                    else:
                                        return "3.8 (f.iv)"
                                else:
                                    return "3.8 (f.ii)"
            else:
                M = [v(II3/I3), v(III3/I3), v(I6/I3^2), v(I/I3^2)]
                x = min(M)
                if x == 0:
                    if v(I6/I)>0:
                        y = v(I6/(III3*(I3+II3)))
                        if y == 0:
                            return "3.9 (b.iii)"
                        elif y > 0:
                            return "3.9 (b.i)"
                        else:
                            return "3.0 (b.ii)"
                    else:
                        return "3.9 (a)"
                elif x > 0:
                    N = [v(II3^2/(I3*III3)), v(I6^2/(I3*III3^3))]
                    if min(N) >= 0:
                        if N[1] == 0:
                            return "3.9 (c.i)"
                        else:
                            if v((II3^2 - 16*I2*III3)/(I3*III3)) == 0:
                                return "3.9 (c.iii)"
                            else:
                                return "3.9 (c.ii)"
                    else:
                        if v((I6*I3)/II3^3) > 0:
                            y = v((I6*I3)/(I3*II3*III3))
                            if y == 0:
                                return "3.9 (c.vii)"
                            elif y > 0:
                                return "3.9 (c.v)"
                            else:
                                return "3.9 (c.vi)"
                        else:
                            return "3.9 (c.iv)"
                else:
                    z = v(I/II3^2)
                    if z == 0:
                        if v(I6/I) == 0:
                            return "3.9 (d)"
                        else:
                            y = v(I6/(III3*(I3+II3)))
                            if y == 0:
                                return "3.9 (b.iii)"
                            elif y > 0:
                                return "3.9 (b.i)"
                            else:
                                return "3.9 (b.ii)"
                    elif z > 0:
                        return "3.9 (e)"
                    else:
                        return "3.9 (d)"
            return "Error"
    
def coeffs(elt):
    R = elt.parent()
    x,y,z = R.gens()
    K = R.base_ring()
    A = elt.coefficient({x:4})
    B = elt.coefficient({y:4})
    C = elt.coefficient({z:4})
    a = elt.coefficient({y:2, z:2})
    b = elt.coefficient({x:2, z:2})
    c = elt.coefficient({x:2, y:2})
    return {'A':K(A), 'B':K(B), 'C':K(C), 'a':K(a), 'b':K(b), 'c':K(c)}
