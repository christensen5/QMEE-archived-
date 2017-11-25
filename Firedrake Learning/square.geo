Point(1) = { -2., -2., 0, 0.1}; 
Point(2) = { -2., 2., 0, 0.1}; 
Point(3) = { 2., 2., 0, 0.1}; 
Point(4) = { 2., -2. , 0, 0.1}; 
Line(1) = { 1,2}; 
Line(2) = { 2,3};
Line(3) = { 3,4};
Line(4) = { 4,1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(6) = {1};

Physical Line(1) = { 1, 2, 3, 4 };

 
Physical Surface(11) = {6};
Mesh.Algorithm = 6; // frontal=6, delannay=5, meshadapt=1
