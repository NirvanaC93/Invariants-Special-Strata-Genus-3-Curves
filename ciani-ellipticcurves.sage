""" 
Code accompanying the paper
"Reduction types of genus-3 curves in a special stratum of their moduli space"


Here, we compute the elliptic curve E_a in the cover Y -> E_a -> X
and in particular derive the formula for the j-invariant provided in Section 2.1
We use the same notation as in that section.
"""

K.<A,B,C,a,b,c> = QQ[]
L.<alpha,beta,gamma> = Frac(K)[]
R.<u,v,w> = L[]

relations = [alpha^2 - 2*a*alpha + 4*B*C, beta^2 - 2*b*beta + 4*A*C, gamma^2 - 2*c*gamma + 4*A*B]
I = L.ideal(relations)

# branch points of the cover 
Pa1 = (0,alpha,-2*B)
Pa2 = (0,-2*C,alpha)
Pb1 = (-2*C,0,beta)
Pb2 = (beta,0,-2*A)
Pc1 = (gamma,-2*A,0)
Pc2 = (-2*B,gamma,0)

# xi is a coordinate that sends (Pb1,Pb2,Pc1) -> (0,1,infty)
xi1 = (-2*A)/(gamma*(beta^2-4*A*C)) * (gamma*beta*u + 2*B*beta*v + 2*C*gamma*w)/w
xi2 = (-2*A)/(gamma*(beta^2-4*A*C)) * (-2*gamma*beta*(C*w+a*v+b*u) + 2*C*gamma*(2*A*u+gamma*v))/(2*A*u+gamma*v)

#check conditions on xi
assert xi1(Pb1).numerator().reduce(I) == 0
assert (xi1(Pb2) - 1).numerator().reduce(I) == 0
assert xi2(Pb1).numerator().reduce(I) == 0
assert (xi2(Pb2) - 1).numerator().reduce(I) == 0
#xi1(Pc1) gives division by zero error

#check that lambda is correct
lam = xi2(Pc2)
lam_check = 1/2 + (2*A*a-b*c)*beta*gamma/(2*(b*beta-4*A*C)*(c*gamma-4*A*B))
assert (lam - lam_check).numerator().reduce(I) == 0

#we compute the j-invariant in terms of lambda = xi(Pc2)
j_inv = 256*(lam^2-lam+1)^3/(lam^2*(lam-1)^2)
Da = a^2 - 4*B*C
Db = b^2 - 4*A*C
Dc = c^2 - 4*A*B
DX = A*a^2 + B*b^2 + C*c^2 - 4*A*B*C - a*b*c
j_inv_paper = -256 * (3*A*DX - (2*A*a-b*c)^2)^3/(A^2*Db*Dc*DX^2)

assert (j_inv-j_inv_paper).numerator().reduce(I) == 0