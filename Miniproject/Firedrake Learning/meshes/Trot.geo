cl__1 = 0.05;
Point(1) = {0, 0, 0, 0.1};
Point(2) = {1, 0, 0, 0.1};
Point(3) = {1, 2, 0, 0.1};
Point(4) = {2, 2, 0, 0.1};
Point(5) = {2, 3, 0, 0.1};
Point(6) = {1, 3, 0, 0.1};
Point(7) = {1, 5, 0, 0.1};
Point(8) = {0, 5, 0, 0.1};
Line(10) = {1, 2};
Line(11) = {2, 3};
Line(12) = {3, 4};
Line(13) = {4, 5};
Line(14) = {5, 6};
Line(15) = {6, 7};
Line(16) = {7, 8};
Line(17) = {8, 1};
Line Loop(19) = {17, 10, 11, 12, 13, 14, 15, 16};
Plane Surface(19) = {19};