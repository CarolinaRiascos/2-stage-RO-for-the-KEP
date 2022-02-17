\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 Obj_GrandSubP: _Excess.0+ Beta.0
Subject To
 _ExcessConst       _Excess.0+ selv.3+ selv.4+ selv.6+ selv.5 + selv.2+ selv.1 = 6
 Selected.3         selv.3- x.1- x.2- x.3- y.2 - x.4 <= 0
 Selected.4         selv.4- x.1- x.2- y.3 - x.4 <= 0
 Selected.6         selv.6- y.2 - y.3 - y.1 - y.4 <= 0
 Selected.5         selv.5- x.2- y.2 - y.3 - y.1 - y.4 <= 0
 Selected.2         selv.2- x.1- x.2- x.3- y.2 - y.3 - y.1 - x.5 <= 0
 Selected.1         selv.1- x.1- x.3- y.1 - x.5 <= 0
 BetaLB             - 4 x.1- 4 x.2- 3 x.3- 3 y.2 - 3 y.3 - 3 y.1 + Beta.0 - x.4 - x.5 - y.4 >= 0
 Active_cyvar2nd.0  x.1+ x.2+ x.3+ y.2 + y.3 + y.1
                       + vertex.3 + arc.3.1 + vertex.1 + arc.1.2
                       + vertex.2 + arc.2.4 + vertex.4 + arc.4.3
                       - x.4 - x.5 >= 1
 Active_cyvar2nd.1  x.1+ x.2+ x.3+ y.2 + y.3 + y.1
                       + vertex.3 + vertex.2 + vertex.4 + arc.3.4
                       + arc.4.5 + vertex.5 + arc.5.2 + arc.2.3
                       - x.4 - x.5 >= 1
 Active_cyvar2nd.2: x.1+ x.2+ x.3+ y.2 + y.3 + y.1
                       + vertex.3 + arc.3.1 + vertex.1 + arc.1.2
                       + vertex.2 + arc.2.3 - x.4 - x.5 >= 1
 Active_chvar2nd.0: x.1+ x.2+ x.3+ y.2 + y.3 + y.1
                       + vertex.1 + vertex.2 + vertex.5 + arc.5.2
                       + vertex.6 + arc.6.5 + arc.2.1 - y.4 >= 1
 Active_chvar2nd.1: x.1+ x.2+ x.3+ y.2 + y.3 + y.1
                       + vertex.3 + vertex.2 + vertex.5 + arc.5.2
                       + arc.2.3 + vertex.6 + arc.6.5 - y.4 >= 1
 Active_chvar2nd.2: x.1+ x.2+ x.3+ y.2 + y.3 + y.1
                       + vertex.2 + arc.2.4 + vertex.4 + vertex.5
                       + arc.5.2 + vertex.6 + arc.6.5 - y.4 >= 1
 VtxSum:            vertex.3 + vertex.1 + vertex.2 + vertex.4
                       + vertex.5 + vertex.6 <= 0
 ArcSum:            arc.3.1 + arc.1.2 + arc.2.4 + arc.4.3
                       + arc.3.4 + arc.4.5 + arc.5.2 + arc.2.3
                       + arc.6.5 + arc.2.1 <= 0
 Active_cyvar2nd.3: x.1+ x.2+ x.3+ y.2 + y.3 + vertex.3
                       + vertex.4 + arc.4.3 + arc.3.4 + x.4 >= 1
 Active_cyvar2nd.4: x.1+ x.2+ x.3+ y.2 + y.3 + y.1
                       + vertex.1 + arc.1.2 + vertex.2 + arc.2.1
                       + x.5 >= 1
 Active_chvar2nd.3: x.2+ y.2 + y.3 + y.1 + vertex.5
                       + vertex.6 + arc.6.5 + y.4 >= 1
 NewBound.1:        Beta.0 + vertex.3 + arc.3.1 + vertex.1
                       + arc.1.2 + vertex.2 + arc.2.4 + vertex.4
                       + arc.4.3 + arc.3.4 + arc.4.5 + vertex.5
                       + arc.5.2 + arc.2.3 + vertex.6 + arc.6.5
                       + arc.2.1 >= 5
Bounds
 0 <= selv.3<= 1
 0 <= selv.4<= 1
 0 <= selv.6<= 1
 0 <= selv.5<= 1
 0 <= selv.2<= 1
 0 <= selv.1<= 1
 0 <= x.1<= 1
 0 <= x.2<= 1
 0 <= x.3<= 1
 0 <= y.2 <= 1
 0 <= y.3 <= 1
 0 <= y.1 <= 1
 0 <= vertex.3 <= 1
 0 <= arc.3.1 <= 1
 0 <= vertex.1 <= 1
 0 <= arc.1.2 <= 1
 0 <= vertex.2 <= 1
 0 <= arc.2.4 <= 1
 0 <= vertex.4 <= 1
 0 <= arc.4.3 <= 1
 0 <= arc.3.4 <= 1
 0 <= arc.4.5 <= 1
 0 <= vertex.5 <= 1
 0 <= arc.5.2 <= 1
 0 <= arc.2.3 <= 1
 0 <= vertex.6 <= 1
 0 <= arc.6.5 <= 1
 0 <= arc.2.1 <= 1
 0 <= x.4 <= 1
 0 <= x.5 <= 1
 0 <= y.4 <= 1
Binaries
 selv.3 selv.4 selv.6 selv.5 selv.2 selv.1 x.1 x.2
 x.3 y.2  y.3  y.1  vertex.3  arc.3.1  vertex.1 
 arc.1.2  vertex.2  arc.2.4  vertex.4  arc.4.3  arc.3.4 
 arc.4.5  vertex.5  arc.5.2  arc.2.3  vertex.6  arc.6.5 
 arc.2.1  x.4  x.5  y.4 
End
