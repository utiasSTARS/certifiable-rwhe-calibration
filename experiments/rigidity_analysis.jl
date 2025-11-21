using LinearAlgebra

p1 = [0; 0; 0]
p2 = [1; 0; 0]
p3 = [0; 1; 0]
p4 = [1; 1; 0]
P = [p1 p2 p3 p4]

z1 = p1 - p2
z2 = p1 - p3
z3 = p1 - p4 
z4 = p2 - p3 
z5 = p2 - p4 
z6 = p3 - p4

Z = [z1 z2 z3 z4 z5 z6]

dZ = [z1 zeros(3, 5);
      zeros(3) z2 zeros(3, 4);
      zeros(3, 2) z3 zeros(3, 3);
      zeros(3, 3) z4 zeros(3, 2);
      zeros(3, 4) z5 zeros(3)
      zeros(3, 5) z6]

incidence = [1 -1 0 0;
             1 0 -1 0;
             1 0 0 -1;
             0 1 -1 0;
             0 1 0 -1;
             0 0 1 -1]

rigidity_test = dZ'*kron(incidence, I(3))

rigidity = [z1 z2 z3 zeros(3, 3);
           -z1 zeros(3, 2) z4 z5 zeros(3);
           zeros(3) -z2 zeros(3) -z4 zeros(3) z6;
           zeros(3, 2) -z3 zeros(3) -z5 -z6]'
