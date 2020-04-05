import pytest
import sympy
from sympy.algebras.quaternion import Quaternion
from .context import DualQuaternion


class TestDualQuaternion():
    a, b, c, d = sympy.symbols('a b c d', real=True)
    x, y, z, w = sympy.symbols('x y z w', real=True)
    q1 = Quaternion(a, b, c, d)
    q2 = Quaternion(x, y, z, w)

    def test_should_create_dual_quaternion(self):
        dq = DualQuaternion(self.q1, self.q2)
        assert isinstance(dq, DualQuaternion)
        assert dq.real == self.q1
        assert dq.dual == self.q2
    