cl1 = 1.0;
cl2 = 10.0;
r = 50;
Point(1) = {0, 0, 0, cl1};
Point(2) = {r, 0, 0, cl2};
Point(3) = {-r, 0, 0, cl2};
Point(4) = {0, r, 0, cl2};
Point(5) = {0, -r, 0, cl2};
Point(6) = {-1, 1, 0, cl1};
Point(7) = {1, 1, 0, cl1};
Point(8) = {1, -1, 0, cl1};
Point(9) = {-1, -1, 0, cl1};
Circle(1) = {4, 1, 2};
Circle(2) = {2, 1, 5};
Circle(3) = {5, 1, 3};
Circle(4) = {3, 1, 4};
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 6};
Line Loop(14) = {5, 6, 7, 8};
Plane Surface(14) = {14};
Line Loop(15) = {1, 2, 3, 4, -8, -7, -6, -5};
Plane Surface(15) = {15};
Physical Surface(16) = {14};
Physical Surface(17) = {15};
