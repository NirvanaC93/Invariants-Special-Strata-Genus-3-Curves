class OurCurves(object):
    """
    Creates a curve of the family

    .. MATH::
    
        A*x^4 + B*y^4 + C*z^4 + a*y^2*z^2 + b*x^2*z^2 + c*x^2*y^2

    INPUT:

    -  ``elt`` - A list of length 6 with the coeficients [A,B,C,a,b,c] or
            a polynomial A*x^4 + B*y^4 + C*z^4 + a*y^2*z^2 + b*x^2*z^2 + c*x^2*y^2

    EXAMPLES::

        sage: C = OurCurves([7,7,7,7,7,1]); C
        Curve over Rational Field defined by 7*x^4 + x^2*y^2 + 7*y^4 + 7*x^2*z^2 + 7*y^2*z^2 + 7*z^4

        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: F = 7*x^4 + 7*y^4 + 7*z^4 + 7*y^2*z^2 + 7*x^2*z^2 + x^2*y^2
        sage: D = OurCurves(F); D
        Curve over Rational Field defined by 7*x^4 + x^2*y^2 + 7*y^4 + 7*x^2*z^2 + 7*y^2*z^2 + 7*z^4
    """
    def __init__(self, elt):
        from sage.rings.polynomial.multi_polynomial_element import is_MPolynomial
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
        elif is_MPolynomial(elt): #given elt is an equation
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
        else:
            raise TypeError("Argument (=%s) must be a list or a polynomial."%T)

    def defining_polynomial(self):
        """
        Return the defining polynomial of the curve.
        """
        return self._equation

    def __repr__(self):
        """
        Return a string representation of this curve.
        """
        return "Curve over {} defined by {}".format(self._K,
                 self._equation)

    def conic(self):
        """
        Return the equation of the conic obtained by quotiening the curve
        by the Klein group.

        EXAMPLES::

            sage: C = OurCurves([7,7,7,7,7,1]); C
            sage: C.conic()
            Projective Conic Curve over Rational Field defined by 7*u^2 + u*v + 7*v^2 + 7*u*w + 7*v*w + 7*w^2
        """
        try:
            return self._conic
        except AttributeError:
            K = self._K
            R.<u,v,w> = PolynomialRing(K)
            A = self._coeffs['A']
            B = self._coeffs['B']
            C = self._coeffs['C']
            a = self._coeffs['a']
            b = self._coeffs['b']
            c = self._coeffs['c']
            self._conic = Conic(A*u^2 + B*v^2 + C*w^2 + a*v*w + b*u*w + c*u*v)
            return self._conic

    def invariants(self):
        """
        Return the projective invariants I_3, I'_3 I''_3, I_6 and I of the curve
        as defined in our paper. The invariants have weight 3, 3, 3, 6 and 6
        respectively.

        EXAMPLES::

            sage: C = OurCurves([7,7,7,7,7,1]); C
            sage: C.invariants()
            [343, -3423, -728, -4213755, 3868011]
        """
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

    def discriminant(self):
        """
        Return the discriminant of the curve.

        EXAMPLES::

            sage: C = OurCurves([7,7,7,7,7,1]); C
            sage: C.disctiminant()
            -417636311076891767890575/256
        """
        [I3, II3, III3, I6, I] = self.invariants()
        return -2^(-20)*I3*III3^4*I6^2
    
    def is_normalized(self, p):
        """
        Return whether the curve is normalized with respect to ``p``, as defined
        in Proposition 3.1 of the paper.

        INPUT:

        -  ``p`` - the prime (different from 2) determining the valuation.

        EXAMPLES::

            sage: C = OurCurves([7,7,7,7,7,1]); C
            sage: C.is_normalized(7)
            False
        """
        if p == 2:
            raise NotImplementedError
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
        """
        Return a normalized model of the curve that is normalized with respect
        to ``p``, as defined in Proposition 3.1 of the paper.

        INPUT:

        -  ``p`` - the prime (different from 2) determining the valuation.

        EXAMPLES::

            sage: C = OurCurves([7,7,7,7,7,1]); C
            sage: C.normalized_model(7)
            Curve over Number Field in m with defining polynomial X^4 - 7
              defined by 7*x^4 + x^2*y^2 + 7*y^4 + (m^2)*x^2*z^2 + (m^2)*y^2*z^2 + z^4
        """
        if p == 2:
            raise NotImplementedError
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
                            P.<X> = PolynomialRing(K)
                            den = r.denominator()
                            poly = X^den - p
                            poly2, _ = poly.factor()[0]
                            K.<m> = K.extension(poly2, 'aux').absolute_field()
                            v = K.valuation(p)
                            elt = v.element_with_valuation(r)
                        R = F.parent().change_ring(K)
                        var[L[-1]] = R(var[L[-1]])/elt
                        F = R(F)(var)
                        D = coeffs(F)
                        assert not all(map(lambda i : v(D[i])  > 0, L[:-1]))
                for L in [['A', 'B', 'c', 2], ['A', 'b', 'C', 1], ['a', 'B', 'C', 0]]:
                    if all(map(lambda i : v(D[i]) > 0, L[:-1])):
                        var = list(F.parent().gens())
                        val = list(map(lambda i : v(D[i])/4, L[:-1]))
                        r = min(val)
                        try:
                            elt = v.element_with_valuation(r)
                        except ValueError:
                            P.<X> = PolynomialRing(K)
                            den = r.denominator()
                            poly = X^den - p
                            poly2, _ = poly.factor()[0]
                            K.<m> = K.extension(poly2, 'aux').absolute_field()
                            v = K.valuation(p)
                            elt = v.element_with_valuation(r)
                        R = F.parent().change_ring(K)
                        for i in range(3):
                            if i == L[-1]:
                                var[i] = elt*R(var[i])
                            else:
                                var[i] = R(var[i])/elt
                        F = R(F)(var)
                        D = coeffs(F)
                        print(F)
                        assert not all(map(lambda i : v(D[i])  > 0, L[:-1]))
                self._normalization[p] = OurCurves(F)
                return self._normalization[p]

    def reduction_type(self, p):
        """
        Return the reduction type of the curve with respect to ``p``,
        using the decision tree in the repository.

        INPUT:

        -  ``p`` - the prime (different from 2) determining the valuation.

        EXAMPLES::

            sage: C = OurCurves([7,7,7,7,7,1]); C
            sage: C.reduction_type(7)
            '3.8 (f.iii)'
        """
        if p == 2:
            raise NotImplementedError
        if self.discriminant() == 0:
            raise ValueError("The curve is singular.")
        try:
            return self._reduction_type[p]
        except KeyError:
            [I3, II3, III3, I6, I] = self.invariants()
            K = self._K
            v = K.valuation(p)
            L = [v(I3/III3), v(II3/III3), v(I6/III3^2)]
            if min(L)>=0:
                if L[2] == 0:
                    if L[0] == 0:
                        return "Good reduction"
                    else:
                        if v(I/III3^2) == 0:
                            return "3.8 (d)"
                        else:
                            if L[1] == 0:
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
                                x = v((I3*II3)/I6)
                                if x == 0:
                                    return "3.8 (f.iii)"
                                elif x < 0:
                                    return "3.8 (f.i)"
                                else:
                                    return "3.8 (f.ii)"
                            else:
                                x = v(I/I6)
                                if x == 0:
                                    if v(I/(I3*III3)) == 0:
                                        return "3.8 (f.iii)"
                                    else:
                                        return "3.8 (f.vi)"
                                elif x < 0:
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
                    if v(I/I6) < 0:
                        y = min(v((III3*I3)/I6), v((III3*II3)/I6))
                        if y == 0:
                            return "3.9 (b.iii)"
                        elif y < 0:
                            return "3.9 (b.i)"
                        else:
                            return "3.9 (b.ii)"
                    else:
                        return "3.9 (a)"
                elif x > 0:
                    N = [v(II3^2/(I3*III3)), v(I6^2/(I3*III3^3))]
                    if min(N) >= 0:
                        if N[1] == 0:
                            return "3.9 (c.i)"
                        else:
                            if v((II3^2 - 16*I3*III3)/(I3*III3)) == 0:
                                return "3.9 (c.ii)"
                            else:
                                return "3.9 (c.iii)"
                    else:
                        if v(II3^3/(I6*I3)) < 0:
                            y = v((I3*II3*III3)/(I6*I3))
                            if y == 0:
                                return "3.9 (c.vii)"
                            elif y < 0:
                                return "3.9 (c.v)"
                            else:
                                return "3.9 (c.vi)"
                        else:
                            return "3.9 (c.iv)"
                else:
                    z = v(I) - 2*v(II3)
                    if z == 0:
                        if v(I/I6) == 0:
                            return "3.9 (d)"
                        else:
                            y = min(v((III3*I3)/I6), v((III3*II3)/I6))
                            if y == 0:
                                return "3.9 (b.iii)"
                            elif y < 0:
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
