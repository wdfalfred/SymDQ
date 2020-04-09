# Reference : 
# DOI: 10.1177/02783649922066213
from __future__ import print_function

from sympy import S
from sympy import simplify, trigsimp
from sympy import conjugate, sin, cos, sqrt
from sympy.algebras.quaternion import Quaternion
from sympy.core.expr import Expr


class DualQuaternion(Expr):
    _op_priority = 11.0

    is_commutative = False

    def __new__(cls, p=Quaternion(0), q=Quaternion(0)):
        if not isinstance(p, Quaternion):
            p = Quaternion(p)
        if not isinstance(q, Quaternion):
            q = Quaternion(q)

        obj = Expr.__new__(cls, p, q)
        obj._p = p
        obj._q = q
        return obj

    @property
    def real(self):
        return self._p

    @property
    def dual(self):
        return self._q

    def __add__(self, other):
        return self.add(other)

    def __radd__(self, other):
        return self.add(other)

    def __sub__(self, other):
        return self.add(other*-1)

    def __mul__(self, other):
        return self._generic_mul(self, other)

    def __rmul__(self, other):
        return self._generic_mul(other, self)

    def __neg__(self):
        return DualQuaternion(-self._p, -self._q)
    
    def add(self, other):
        dq1 = self
        dq2 = simplify(other)

        if not isinstance(dq2, DualQuaternion):
            if isinstance(dq2, Quaternion) or dq2.is_commutative:
                return DualQuaternion(dq1.real + dq2, dq1.dual)
            else:
                raise ValueError("Only DualQuaternions or Quaternions can be added with a DualQuaternion.")
        return DualQuaternion(dq1.real + dq2.real, dq1.dual + dq2.dual)
    
    def mul(self, other):
        return self._generic_mul(self, other)

    @staticmethod
    def _generic_mul(dq1, dq2):
        dq1 = simplify(dq1)
        dq2 = simplify(dq2)

        if not isinstance(dq1, DualQuaternion) and not isinstance(dq2, DualQuaternion):
            return dq1 * dq2

        if not isinstance(dq1, DualQuaternion):
            if dq1.is_commutative:
                return DualQuaternion(dq1 * dq2._p, dq1 * dq2._q)
            else:
                raise ValueError("Only DualQuaternions or commutative expressions can be multiplied with a DualQuaternion.")

        if not isinstance(dq2, DualQuaternion):
            if dq2.is_commutative:
                return DualQuaternion(dq1._p * dq2, dq1._q * dq2)
            else:
                raise ValueError("Only DualQuaternions or commutative expressions can be multiplied with a DualQuaternion.")

        return DualQuaternion(dq1._p * dq2._p, dq1._p * dq2._q + dq1._q * dq2._p)
    
    def quaternion_conjugate(self):
        return DualQuaternion(conjugate(self._p), conjugate(self._q))

    def dual_number_conjugate(self):
        return DualQuaternion(self._p, -self._q)

    def combined_conjugate(self):
        return DualQuaternion(conjugate(self._p), -conjugate(self._q))

    def norm(self):
        """Returns the norm of the dual quaternion."""
        return self * self.quaternion_conjugate()

    def is_unit(self):
        """Returns True if is a unit dual quaternion and False otherwise."""
        p = self._p
        q = self._q
        return trigsimp(p.a**2 + p.b**2 + p.c**2 + p.d**2) == 1 and trigsimp(p.a*q.a + p.b*q.b + p.c*q.c + p.d*q.d) == 0

    def rotation_matrix(self):
        return self._p.to_rotation_matrix()

    @staticmethod
    def transform_point(pin, t):
        """Returns the coordinates of the point pin(a 3 tuple) after transformation.
        Parameters
        ==========
        pin : tuple
            A 3-element tuple of coordinates of a point which needs to be
            transformed.
        r : DualQuaternion
            Screw axis and dual angle of rotation.

        Returns
        =======
        tuple
            The coordinates of the point after transformation.
        """
        pout =  (t * DualQuaternion(1, Quaternion(0, *pin)) * t.combined_conjugate()).dual
        return (pout.b, pout.c, pout.d)

    @classmethod
    def from_screw(cls, l, m, theta, d):
        """Returns the unit dual quaternion corresponding to a screw.

        Parameters
        ==========
        l : tuple
            unit vector parallel to the screw axis
        m : tuple
            moment of l 
        theta : number
        d : number

        Returns
        =======
        DualQuaternion

        """
        (x, y, z) = l
        (a, b, c) = m
        if trigsimp(x**2 + y**2 + z**2) != 1 or trigsimp(x*a + y*b + z*c) != 0:
            raise ValueError("Expected l to be a unit vector and m perpendicular to l!")

        q_r = Quaternion(cos(theta * S.Half), *(sin(theta * S.Half) * i for i in l))
        q_d = Quaternion(-d * S.Half * sin(theta * S.Half), 
        *(d * S.Half * cos(theta * S.Half) * i + sin(theta * S.Half) * j for i, j in zip(l, m)))
        return DualQuaternion(q_r, q_d)
