function [alok1, alok2, alok3, alok4, alok6] = calc_aloks_Tang(p1, p2, p3, Q)

[alok1, alok2] = exercise24_Tang(p1, p2, p3, Q);
[alok3, alok4] = exercise43_Tang(p1, p2, p3, Q);
alok6 = project_work11_alok6_Tang(p1, p2, p3, Q);