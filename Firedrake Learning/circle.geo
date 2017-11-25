d = 2;
Point(1) = { 0, 0, 0, 0.1}; 
Point(2) = { -d/2, 0.0, 0, 0.1}; 
Point(3) = { 0 , d/2, 0, 0.1}; 
Point(4) = { d/2, 0. , 0, 0.1}; 
Point(5) = { 0, -d/2, 0, 0.1}; 
Circle(1) = { 2, 1, 3}; 
Circle(2) = { 3, 1, 4};
Circle(3) = { 4, 1, 5};
Circle(4) = { 5, 1, 2};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(6) = {1};

Physical Line(1) = { 1, 2, 3, 4 };

 
Physical Surface(11) = {6};
Mesh.Algorithm = 6; // frontal=6, delannay=5, meshadapt=1
