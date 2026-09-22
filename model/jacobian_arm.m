function Jm = jacobian_arm(phi, L)

q1 = phi(1);
q2 = phi(2);
q3 = phi(3);

L1 = L(1);
L2 = L(2);
L3 = L(3);

Jm = [
    -L1*sin(q1) - L2*sin(q1+q2) - L3*sin(q1+q2+q3), ...
    -L2*sin(q1+q2) - L3*sin(q1+q2+q3), ...
    -L3*sin(q1+q2+q3);

    L1*cos(q1) + L2*cos(q1+q2) + L3*cos(q1+q2+q3), ...
    L2*cos(q1+q2) + L3*cos(q1+q2+q3), ...
    L3*cos(q1+q2+q3)
    ];

end