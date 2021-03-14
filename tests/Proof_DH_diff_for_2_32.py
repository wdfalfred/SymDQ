import pytest
from sympy import *
from sympy import symbols, conjugate
from sympy import sin, cos
from sympy.abc import a, b, c, d, x, y, z, w, theta
from sympy.algebras.quaternion import Quaternion
from context import DualQuaternion
from sympy import simplify
from sympy.abc import theta, alpha

print("test")
d, a = symbols('d  a')
q1 = Quaternion(1, 0, 0, 0)
q0 = Quaternion(0, 0, 0, 0)
rdi = Quaternion(0, 0, 0, d)
qti = Quaternion(cos(theta * 0.5), 0, 0, sin(theta * 0.5))
rai = Quaternion(0, a, 0, 0)
qai = Quaternion(cos(alpha * 0.5), 0, 0, sin(alpha * 0.5))

# dqti = diff(qti, theta)
# print(dqti)

dq1 = DualQuaternion(q1, 0.5 * rdi)
dq2 = DualQuaternion(qti, q0)
dq3 = DualQuaternion(q1, 0.5 * rai)
dq4 = DualQuaternion(qai, q0)

dq12 = dq1 * dq2
dq123 = dq12 * dq3
dq1234 = dq123 * dq4

ddqr = diff(dq1234.real, theta)

ddqd = diff(dq1234.dual, theta)

ddq1234 = DualQuaternion(ddqr, ddqd)

pprint(simplify(2 * ddq1234 * dq1234.quaternion_conjugate()))
 

pprint(diff(sin(2 * theta) * cos(alpha), theta))
