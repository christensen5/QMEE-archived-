cl__1 = 0.3;
cl__2 = 0.1;
Point(1) = {0, 0, 0, 0.3};
Point(2) = {2.8, 0, 0, 0.3};
Point(3) = {10.8, 0, 0, 0.3};
Point(4) = {13.6, 0, 0, 0.3};
Point(5) = {10.8, 4.8, 0, 0.3};
Point(6) = {10.8, 8.800000000000001, 0, 0.3};
Point(7) = {10.8, 13.6, 0, 0.3};
Point(8) = {13.6, 13.6, 0, 0.3};
Point(9) = {0, 13.6, 0, 0.3};
Point(10) = {2.8, 13.6, 0, 0.3};
Point(11) = {2.8, 4.8, 0, 0.3};
Point(12) = {2.8, 8.800000000000001, 0, 0.3};
Line(1) = {9, 1};
Line(2) = {1, 2};
Line(3) = {2, 11};
Line(4) = {11, 5};
Line(5) = {5, 3};
Line(6) = {3, 4};
Line(7) = {4, 8};
Line(8) = {8, 7};
Line(9) = {7, 6};
Line(10) = {6, 12};
Line(11) = {12, 10};
Line(12) = {10, 9};
Line Loop(14) = {10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9};
Plane Surface(14) = {14};
Physical Line(15) = {1, 3, 4, 5, 7, 9, 10, 11};
Physical Line(16) = {2, 6};
Physical Line(17) = {8, 12};
Physical Surface(18) = {14};
