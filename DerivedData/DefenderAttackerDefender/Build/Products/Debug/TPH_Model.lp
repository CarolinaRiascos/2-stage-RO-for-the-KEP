\ENCODING=ISO-8859-1
\Problem name: IloCplex

Maximize
 obj1: 4 r.1 + 4 r.2 + 2 r.6 + 3 s.3 + 3 r.4 + 2 r.5 + 3 s.2 + 2 s.4 + s.5
       + 3 r.3 + 3 s.1
Subject To
 c1: r.1 + r.2 + r.6 + s.3 <= 1
 c2: r.1 + r.2 + r.6 + s.3 + r.4 + r.5 + s.2 + s.4 + s.5 <= 1
 c3: r.1 + r.2 + r.4 + s.4 + r.3 <= 1
 c4: r.2 + r.4 + r.5 + s.5 + r.3 <= 1
 c5: s.3 + r.4 + r.5 + s.2 + s.4 + s.5 + s.1 <= 1
 c6: s.3 + s.2 + s.4 + s.5 + s.1 <= 1
Bounds
 0 <= r.1 <= 1
 0 <= r.2 <= 1
 0 <= r.6 <= 1
 0 <= s.3 <= 1
 0 <= r.4 <= 1
 0 <= r.5 <= 1
 0 <= s.2 <= 1
 0 <= s.4 <= 1
 0 <= s.5 <= 1
 0 <= r.3 <= 1
 0 <= s.1 <= 1
Binaries
 r.1  r.2  r.6  s.3  r.4  r.5  s.2  s.4  s.5  r.3  s.1 
End
