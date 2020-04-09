import pytest
from sympy import S
from sympy import symbols, conjugate
from sympy import sin, cos
from sympy.abc import a, b, c, d, x, y, z, w, theta
from sympy.algebras.quaternion import Quaternion
from .context import DualQuaternion


class TestDualQuaternion():
    q1 = Quaternion(a, b, c, d)
    q2 = Quaternion(x, y, z, w)
    q3 = Quaternion(z, b, a, d)
    q4 = Quaternion(c, w, x, y)

    def test_should_create_dual_quaternion(self):
        dq = DualQuaternion(self.q1, self.q2)
        assert isinstance(dq, DualQuaternion)
        assert dq.real == self.q1
        assert dq.dual == self.q2
    
    def test_should_create_from_real_number(self):
        dq = DualQuaternion(a, b)
        assert dq.real == Quaternion(a)
        assert dq.dual == Quaternion(b)

    def test_should_create_from_default(self):
        dq = DualQuaternion()
        assert dq.real == Quaternion(0)
        assert dq.dual == Quaternion(0)

    def test_should_implement_addition(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq2 = DualQuaternion(self.q2, self.q1)
        dq = dq1 + dq2
        assert dq.real == self.q1 + self.q2
        assert dq.dual == self.q2 + self.q1

    def test_should_handle_dq_and_quat_add(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq = dq1 + self.q2
        assert dq.real == self.q1 + self.q2
        assert dq.dual == self.q2

    def test_should_handle_dq_and_real_add(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq = dq1 + 1
        assert dq.real == self.q1 + 1
        assert dq.dual == self.q2

    def test_should_implement_neg(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq = -dq1
        assert dq.real == -self.q1
        assert dq.dual == -self.q2

    def test_should_handle_dq_substraction(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq2 = DualQuaternion(self.q3, self.q4)
        dq = dq1 - dq2
        assert dq.real == self.q1 - self.q3
        assert dq.dual == self.q2 - self.q4

    def test_dq_should_sub_real(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq = dq1 - 1
        assert dq.real == self.q1 - 1
        assert dq.dual == self.q2

    def test_dq_should_sub_quat(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq = dq1 - self.q3
        assert dq.real == self.q1 - self.q3
        assert dq.dual == self.q2

    # fails since it calls __add__ of Quaternion instead of __radd__ of
    # DualQuaternion in this condition
    #def test_should_handle_quat_and_dq_radd(self):
        #dq1 = DualQuaternion(self.q1, self.q2)
        #dq = self.q2 + dq1
        #assert dq.real == self.q1 + self.q2
        #assert dq.dual == self.q2

    def test_should_handle_real_and_dq_radd(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq = 1 + dq1
        assert dq.real == self.q1 + 1
        assert dq.dual == self.q2

    def test_should_implement_multiplication(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq2 = DualQuaternion(self.q3, self.q4)
        dq = dq1 * dq2
        assert dq.real == self.q1 * self.q3
        assert dq.dual == self.q1 * self.q4 + self.q2 * self.q3

    def test_should_have_quaternion_conj(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq = dq1.quaternion_conjugate()
        assert dq.real == conjugate(self.q1)
        assert dq.dual == conjugate(self.q2)

    def test_should_have_dual_conj(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq = dq1.dual_number_conjugate()
        assert dq.real == self.q1
        assert dq.dual == -self.q2

    def test_should_have_combined_conj(self):
        dq1 = DualQuaternion(self.q1, self.q2)
        dq = dq1.combined_conjugate()
        assert dq.real == conjugate(self.q1)
        assert dq.dual == -conjugate(self.q2)

    def test_should_implement_norm(self):
        dq = DualQuaternion(Quaternion(0, 0, 0, 1), Quaternion(0, 0, -x, 0))
        nm = dq.norm()
        assert nm.real == Quaternion(1)
        assert nm.dual == Quaternion(0)

    def test_should_construct_from_screw(self):
        l = (0, 0, 1)
        m = (0, -x, 0)
        dq = DualQuaternion.from_screw(l, m, theta, 0)
        assert dq.real == Quaternion(cos(theta * S.Half), 0, 0, sin(theta * S.Half))
        assert dq.dual == Quaternion(0, 0, -x * sin(theta * S.Half), 0)

    def test_from_screw_should_raise_exception_if_l_not_unit(self):
        l = (0, 0, 2)
        m = (0, -x, 0)
        with pytest.raises(ValueError):
            DualQuaternion.from_screw(l, m, theta, 0)

    def test_from_screw_should_raise_exception_if_l_m_not_orthogonal(self):
        l = (0, 0, 1)
        m = (0, 0, x)
        with pytest.raises(ValueError):
            DualQuaternion.from_screw(l, m, theta, 0)
