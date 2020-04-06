from __future__ import print_function

from sympy import simplify, conjugate, sin, cos
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
        p = simplify(p)
        q = simplify(q)
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
        pass

    @classmethod
    def from_screw(cls, l, m, theta, d):
        pass