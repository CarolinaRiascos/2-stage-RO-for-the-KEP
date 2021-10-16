\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj: Beta.0
Subject To
 BoundingC.1:   Beta.0 - 3 r.0 - 3 r.1 - 3 r.2 - 3 r.3 - 3 r.4 - 2 r.5 - 2 r.6
                - 3 r.7 - 2 r.8 - 3 r.9 >= 0
 ActiveCycle.0: r.0 + vertex.1 + arc.1.37 + vertex.37 + arc.37.35 + vertex.35
                + arc.35.1 >= 1
 ActiveCycle.1: r.1 + vertex.6 + arc.6.20 + vertex.20 + arc.20.45 + vertex.45
                + arc.45.6 >= 1
 ActiveCycle.2: r.2 + vertex.10 + arc.10.32 + vertex.32 + arc.32.14 + vertex.14
                + arc.14.10 >= 1
 ActiveCycle.3: r.3 + vertex.12 + arc.12.51 + vertex.51 + arc.51.18 + vertex.18
                + arc.18.12 >= 1
 ActiveCycle.4: r.4 + vertex.16 + arc.16.21 + vertex.21 + arc.21.36 + vertex.36
                + arc.36.16 >= 1
 ActiveCycle.5: r.5 + vertex.30 + arc.30.52 + vertex.52 + arc.52.30 >= 1
 ActiveCycle.6: r.6 + vertex.34 + arc.34.46 + vertex.46 + arc.46.34 >= 1
 ActiveCycle.7: r.7 + vertex.38 + arc.38.48 + vertex.48 + arc.48.57 + vertex.57
                + arc.57.38 >= 1
 ActiveCycle.8: r.8 + vertex.39 + arc.39.58 + vertex.58 + arc.58.39 >= 1
 ActiveCycle.9: r.9 + vertex.42 + arc.42.62 + vertex.62 + arc.62.50 + vertex.50
                + arc.50.42 >= 1
 MakeOneFail.1: vertex.1 + arc.1.37 + vertex.37 + arc.37.35 + vertex.35
                + arc.35.1 + vertex.6 + arc.6.20 + vertex.20 + arc.20.45
                + vertex.45 + arc.45.6 + vertex.10 + arc.10.32 + vertex.32
                + arc.32.14 + vertex.14 + arc.14.10 + vertex.12 + arc.12.51
                + vertex.51 + arc.51.18 + vertex.18 + arc.18.12 + vertex.16
                + arc.16.21 + vertex.21 + arc.21.36 + vertex.36 + arc.36.16
                + vertex.30 + arc.30.52 + vertex.52 + arc.52.30 + vertex.34
                + arc.34.46 + vertex.46 + arc.46.34 + vertex.38 + arc.38.48
                + vertex.48 + arc.48.57 + vertex.57 + arc.57.38 + vertex.39
                + arc.39.58 + vertex.58 + arc.58.39 + vertex.42 + arc.42.62
                + vertex.62 + arc.62.50 + vertex.50 + arc.50.42 >= 1
 c13:           vertex.1 + vertex.37 + vertex.35 + vertex.6 + vertex.20
                + vertex.45 + vertex.10 + vertex.32 + vertex.14 + vertex.12
                + vertex.51 + vertex.18 + vertex.16 + vertex.21 + vertex.36
                + vertex.30 + vertex.52 + vertex.34 + vertex.46 + vertex.38
                + vertex.48 + vertex.57 + vertex.39 + vertex.58 + vertex.42
                + vertex.62 + vertex.50 + vertex.0 + vertex.2 + vertex.3
                + vertex.4 + vertex.5 + vertex.7 + vertex.8 + vertex.9
                + vertex.11 + vertex.13 + vertex.15 + vertex.17 + vertex.19
                + vertex.22 + vertex.23 + vertex.24 + vertex.25 + vertex.26
                + vertex.27 + vertex.28 + vertex.29 + vertex.31 + vertex.33
                + vertex.40 + vertex.41 + vertex.43 + vertex.44 + vertex.47
                + vertex.49 + vertex.53 + vertex.54 + vertex.55 + vertex.56
                + vertex.59 + vertex.60 + vertex.61 + vertex.63 <= 1
 c14:           arc.1.37 + arc.37.35 + arc.35.1 + arc.6.20 + arc.20.45
                + arc.45.6 + arc.10.32 + arc.32.14 + arc.14.10 + arc.12.51
                + arc.51.18 + arc.18.12 + arc.16.21 + arc.21.36 + arc.36.16
                + arc.30.52 + arc.52.30 + arc.34.46 + arc.46.34 + arc.38.48
                + arc.48.57 + arc.57.38 + arc.39.58 + arc.58.39 + arc.42.62
                + arc.62.50 + arc.50.42 + arc.0.36 + arc.0.18 + arc.0.50
                + arc.0.48 + arc.0.58 + arc.0.35 + arc.1.63 + arc.1.19
                + arc.1.57 + arc.1.39 + arc.3.5 + arc.3.24 + arc.3.48
                + arc.3.62 + arc.3.11 + arc.3.25 + arc.3.35 + arc.3.0 + arc.3.8
                + arc.3.18 + arc.3.61 + arc.3.56 + arc.3.37 + arc.3.57
                + arc.3.4 + arc.3.53 + arc.3.43 + arc.3.49 + arc.3.55 + arc.3.6
                + arc.3.63 + arc.3.23 + arc.3.51 + arc.3.59 + arc.3.47
                + arc.3.2 + arc.3.32 + arc.3.50 + arc.3.29 + arc.3.36
                + arc.3.19 + arc.3.30 + arc.3.44 + arc.3.7 + arc.3.9 + arc.3.38
                + arc.4.35 + arc.4.18 + arc.4.58 + arc.4.48 + arc.4.38
                + arc.4.50 + arc.4.36 + arc.5.35 + arc.5.18 + arc.5.48
                + arc.5.58 + arc.5.38 + arc.5.36 + arc.6.63 + arc.6.19
                + arc.6.39 + arc.6.21 + arc.6.14 + arc.6.32 + arc.6.37
                + arc.7.20 + arc.7.19 + arc.7.42 + arc.7.37 + arc.7.21
                + arc.7.10 + arc.9.37 + arc.9.19 + arc.9.63 + arc.9.39
                + arc.9.14 + arc.9.32 + arc.10.14 + arc.10.37 + arc.10.63
                + arc.10.57 + arc.10.19 + arc.11.1 + arc.11.37 + arc.11.19
                + arc.11.16 + arc.11.32 + arc.11.14 + arc.12.10 + arc.12.21
                + arc.12.19 + arc.12.37 + arc.12.32 + arc.12.63 + arc.12.14
                + arc.13.26 + arc.13.44 + arc.13.6 + arc.13.50 + arc.13.48
                + arc.13.8 + arc.13.32 + arc.13.18 + arc.13.14 + arc.13.36
                + arc.13.0 + arc.13.54 + arc.13.4 + arc.13.62 + arc.13.60
                + arc.13.34 + arc.13.56 + arc.13.2 + arc.13.11 + arc.13.39
                + arc.13.38 + arc.13.29 + arc.13.35 + arc.13.21 + arc.13.33
                + arc.13.37 + arc.13.43 + arc.13.25 + arc.13.9 + arc.13.41
                + arc.13.61 + arc.13.1 + arc.13.59 + arc.13.5 + arc.13.53
                + arc.13.19 + arc.13.63 + arc.13.7 + arc.14.19 + arc.14.51
                + arc.14.63 + arc.14.57 + arc.14.37 + arc.15.36 + arc.15.4
                + arc.15.32 + arc.15.0 + arc.15.58 + arc.15.51 + arc.15.49
                + arc.15.2 + arc.15.55 + arc.15.10 + arc.15.14 + arc.15.57
                + arc.15.19 + arc.15.43 + arc.15.11 + arc.15.21 + arc.15.37
                + arc.15.61 + arc.15.29 + arc.15.35 + arc.15.24 + arc.15.60
                + arc.15.7 + arc.15.48 + arc.15.62 + arc.15.47 + arc.15.6
                + arc.15.5 + arc.15.9 + arc.15.44 + arc.15.56 + arc.15.34
                + arc.15.63 + arc.15.25 + arc.15.54 + arc.16.32 + arc.16.10
                + arc.16.20 + arc.16.19 + arc.16.37 + arc.17.48 + arc.17.40
                + arc.17.36 + arc.17.18 + arc.17.35 + arc.18.19 + arc.18.37
                + arc.18.63 + arc.18.57 + arc.20.61 + arc.20.28 + arc.20.0
                + arc.20.24 + arc.20.34 + arc.20.4 + arc.20.60 + arc.20.56
                + arc.20.18 + arc.20.2 + arc.20.6 + arc.20.42 + arc.20.8
                + arc.20.44 + arc.20.37 + arc.20.9 + arc.20.50 + arc.20.35
                + arc.20.19 + arc.20.29 + arc.20.14 + arc.20.5 + arc.20.36
                + arc.20.53 + arc.20.11 + arc.20.47 + arc.20.43 + arc.20.63
                + arc.20.25 + arc.21.18 + arc.21.35 + arc.21.48 + arc.21.38
                + arc.22.37 + arc.22.55 + arc.22.19 + arc.22.61 + arc.22.50
                + arc.22.59 + arc.22.43 + arc.22.11 + arc.22.6 + arc.22.49
                + arc.22.63 + arc.22.4 + arc.22.8 + arc.22.0 + arc.22.58
                + arc.22.18 + arc.22.62 + arc.22.48 + arc.22.44 + arc.22.54
                + arc.22.26 + arc.22.14 + arc.22.34 + arc.22.56 + arc.22.36
                + arc.22.60 + arc.22.30 + arc.22.38 + arc.22.2 + arc.22.10
                + arc.22.5 + arc.22.9 + arc.22.21 + arc.22.25 + arc.22.53
                + arc.22.35 + arc.22.17 + arc.23.35 + arc.23.48 + arc.23.36
                + arc.23.18 + arc.23.58 + arc.24.38 + arc.24.48 + arc.24.36
                + arc.24.40 + arc.24.58 + arc.24.18 + arc.24.35 + arc.25.37
                + arc.25.14 + arc.25.19 + arc.25.63 + arc.26.57 + arc.26.19
                + arc.26.37 + arc.26.63 + arc.26.10 + arc.26.42 + arc.26.46
                + arc.26.32 + arc.27.9 + arc.27.37 + arc.27.59 + arc.27.61
                + arc.27.43 + arc.27.53 + arc.27.57 + arc.27.33 + arc.27.36
                + arc.27.22 + arc.27.55 + arc.27.32 + arc.27.19 + arc.27.34
                + arc.27.20 + arc.27.63 + arc.27.0 + arc.27.11 + arc.27.38
                + arc.27.24 + arc.27.8 + arc.27.6 + arc.27.44 + arc.27.10
                + arc.27.2 + arc.27.7 + arc.27.62 + arc.27.12 + arc.27.25
                + arc.27.56 + arc.27.18 + arc.27.29 + arc.27.48 + arc.27.5
                + arc.27.35 + arc.27.60 + arc.27.14 + arc.27.42 + arc.29.35
                + arc.29.48 + arc.29.36 + arc.29.18 + arc.29.50 + arc.30.54
                + arc.30.48 + arc.30.50 + arc.30.2 + arc.30.44 + arc.30.56
                + arc.30.18 + arc.30.60 + arc.30.26 + arc.30.34 + arc.30.61
                + arc.30.9 + arc.30.62 + arc.30.14 + arc.30.43 + arc.30.63
                + arc.30.5 + arc.30.53 + arc.30.13 + arc.30.35 + arc.30.21
                + arc.30.47 + arc.30.11 + arc.30.19 + arc.30.23 + arc.30.36
                + arc.30.55 + arc.30.8 + arc.30.24 + arc.30.58 + arc.30.4
                + arc.30.37 + arc.30.0 + arc.30.25 + arc.31.46 + arc.31.10
                + arc.31.32 + arc.31.20 + arc.31.1 + arc.31.63 + arc.31.37
                + arc.31.19 + arc.32.19 + arc.32.20 + arc.32.57 + arc.33.63
                + arc.33.37 + arc.33.19 + arc.33.42 + arc.34.20 + arc.34.21
                + arc.34.19 + arc.34.37 + arc.34.63 + arc.34.32 + arc.35.14
                + arc.35.32 + arc.35.37 + arc.35.42 + arc.35.57 + arc.35.19
                + arc.35.63 + arc.36.20 + arc.36.39 + arc.36.21 + arc.36.14
                + arc.36.42 + arc.36.19 + arc.36.41 + arc.36.37 + arc.36.63
                + arc.36.32 + arc.37.36 + arc.37.18 + arc.37.48 + arc.38.35
                + arc.38.36 + arc.38.18 + arc.38.58 + arc.39.36 + arc.39.35
                + arc.39.18 + arc.41.32 + arc.41.63 + arc.41.57 + arc.41.19
                + arc.41.14 + arc.41.37 + arc.42.5 + arc.42.56 + arc.42.44
                + arc.42.30 + arc.42.24 + arc.42.35 + arc.42.50 + arc.42.58
                + arc.42.48 + arc.42.54 + arc.42.36 + arc.42.18 + arc.42.8
                + arc.42.33 + arc.42.25 + arc.42.37 + arc.42.4 + arc.42.0
                + arc.42.31 + arc.42.39 + arc.42.23 + arc.42.53 + arc.42.63
                + arc.42.2 + arc.42.6 + arc.42.11 + arc.42.29 + arc.42.61
                + arc.42.19 + arc.42.47 + arc.42.59 + arc.42.43 + arc.42.60
                + arc.42.34 + arc.43.37 + arc.43.63 + arc.43.21 + arc.43.57
                + arc.43.19 + arc.43.39 + arc.43.14 + arc.44.63 + arc.44.20
                + arc.44.14 + arc.44.32 + arc.44.37 + arc.44.19 + arc.44.21
                + arc.44.39 + arc.45.19 + arc.45.44 + arc.45.29 + arc.45.0
                + arc.45.35 + arc.45.48 + arc.45.8 + arc.45.62 + arc.45.2
                + arc.45.4 + arc.45.37 + arc.45.57 + arc.45.11 + arc.45.42
                + arc.45.60 + arc.45.59 + arc.45.54 + arc.45.14 + arc.45.36
                + arc.45.18 + arc.45.50 + arc.45.34 + arc.45.58 + arc.45.43
                + arc.45.61 + arc.45.5 + arc.45.7 + arc.45.56 + arc.45.41
                + arc.45.52 + arc.45.38 + arc.45.20 + arc.45.9 + arc.45.24
                + arc.45.63 + arc.45.21 + arc.45.25 + arc.45.53 + arc.46.19
                + arc.46.6 + arc.46.56 + arc.46.18 + arc.46.36 + arc.46.30
                + arc.46.60 + arc.46.5 + arc.46.44 + arc.46.48 + arc.46.9
                + arc.46.0 + arc.46.4 + arc.46.62 + arc.46.58 + arc.46.54
                + arc.46.10 + arc.46.14 + arc.46.2 + arc.46.23 + arc.46.8
                + arc.46.37 + arc.46.29 + arc.46.25 + arc.46.59 + arc.46.55
                + arc.46.11 + arc.46.7 + arc.46.35 + arc.46.63 + arc.46.53
                + arc.46.61 + arc.46.20 + arc.46.38 + arc.46.43 + arc.47.18
                + arc.47.35 + arc.47.58 + arc.47.48 + arc.47.36 + arc.48.20
                + arc.48.63 + arc.48.19 + arc.48.37 + arc.48.14 + arc.48.10
                + arc.49.48 + arc.49.56 + arc.49.20 + arc.49.44 + arc.49.8
                + arc.49.6 + arc.49.24 + arc.49.36 + arc.49.30 + arc.49.60
                + arc.49.38 + arc.49.35 + arc.49.5 + arc.49.53 + arc.49.61
                + arc.49.50 + arc.49.29 + arc.49.63 + arc.49.58 + arc.49.7
                + arc.49.43 + arc.49.25 + arc.49.57 + arc.49.59 + arc.49.37
                + arc.49.0 + arc.49.19 + arc.49.9 + arc.49.34 + arc.49.11
                + arc.49.18 + arc.49.4 + arc.49.54 + arc.49.2 + arc.50.10
                + arc.50.63 + arc.50.32 + arc.50.37 + arc.50.21 + arc.50.19
                + arc.50.12 + arc.51.40 + arc.51.35 + arc.51.48 + arc.51.58
                + arc.51.38 + arc.51.36 + arc.52.4 + arc.52.21 + arc.52.0
                + arc.52.60 + arc.52.43 + arc.52.53 + arc.52.55 + arc.52.2
                + arc.52.63 + arc.52.19 + arc.52.61 + arc.52.20 + arc.52.36
                + arc.52.34 + arc.52.58 + arc.52.6 + arc.52.50 + arc.52.5
                + arc.52.14 + arc.52.11 + arc.52.44 + arc.52.56 + arc.52.54
                + arc.52.48 + arc.52.46 + arc.52.25 + arc.52.35 + arc.52.42
                + arc.52.18 + arc.52.9 + arc.52.29 + arc.52.10 + arc.52.33
                + arc.52.59 + arc.52.23 + arc.52.37 + arc.53.35 + arc.53.58
                + arc.53.36 + arc.53.18 + arc.53.50 + arc.53.38 + arc.53.48
                + arc.54.32 + arc.54.42 + arc.54.21 + arc.54.19 + arc.54.37
                + arc.54.10 + arc.54.57 + arc.54.63 + arc.54.20 + arc.55.1
                + arc.55.20 + arc.55.32 + arc.55.19 + arc.55.37 + arc.55.21
                + arc.55.63 + arc.55.57 + arc.56.19 + arc.56.63 + arc.56.21
                + arc.56.57 + arc.56.10 + arc.56.32 + arc.56.37 + arc.57.36
                + arc.57.48 + arc.57.35 + arc.58.20 + arc.58.14 + arc.58.10
                + arc.58.32 + arc.58.19 + arc.58.63 + arc.58.21 + arc.58.37
                + arc.58.42 + arc.59.20 + arc.59.16 + arc.59.19 + arc.59.57
                + arc.59.10 + arc.59.37 + arc.59.63 + arc.59.42 + arc.60.18
                + arc.60.35 + arc.60.58 + arc.60.36 + arc.60.48 + arc.60.38
                + arc.62.48 + arc.62.58 + arc.62.36 + arc.62.35 + arc.62.18
                + arc.63.38 + arc.63.48 + arc.63.35 + arc.63.18 + arc.63.36
                + arc.64.2 + arc.64.43 + arc.64.9 + arc.64.61 + arc.64.58
                + arc.64.53 + arc.64.5 + arc.64.11 + arc.64.25 + arc.64.48
                + arc.64.60 + arc.64.62 + arc.64.49 + arc.64.34 + arc.64.44
                + arc.64.0 + arc.64.8 + arc.64.36 + arc.64.18 + arc.64.54
                + arc.64.6 + arc.64.35 + arc.64.23 + arc.64.4 + arc.64.30
                + arc.64.29 + arc.64.7 + arc.64.55 + arc.65.9 + arc.65.61
                + arc.65.43 + arc.65.7 + arc.65.5 + arc.65.24 + arc.65.59
                + arc.65.29 + arc.65.30 + arc.65.56 + arc.65.2 + arc.65.62
                + arc.65.54 + arc.65.53 + arc.65.25 + arc.65.28 + arc.65.60
                + arc.65.44 + arc.65.4 + arc.65.34 + arc.65.8 + arc.65.0
                + arc.65.11 + arc.65.47 + arc.66.8 + arc.66.24 + arc.66.44
                + arc.66.7 + arc.66.47 + arc.66.4 + arc.66.34 + arc.66.60
                + arc.66.62 + arc.66.56 + arc.66.2 + arc.66.0 + arc.66.11
                + arc.66.23 + arc.66.53 + arc.66.5 + arc.66.29 + arc.66.61
                + arc.66.9 + arc.66.6 + arc.66.43 <= 2
 c15:           21 vertex.0 + arc.0.36 + arc.0.18 + arc.0.50 + arc.0.48
                + arc.0.58 + arc.0.35 + arc.3.0 + arc.13.0 + arc.15.0
                + arc.20.0 + arc.22.0 + arc.27.0 + arc.30.0 + arc.42.0
                + arc.45.0 + arc.46.0 + arc.49.0 + arc.52.0 + arc.64.0
                + arc.65.0 + arc.66.0 <= 21
 c16:           10 vertex.1 + arc.1.37 + arc.35.1 + arc.1.63 + arc.1.19
                + arc.1.57 + arc.1.39 + arc.11.1 + arc.13.1 + arc.31.1
                + arc.55.1 <= 10
 c17:           15 vertex.2 + arc.3.2 + arc.13.2 + arc.15.2 + arc.20.2
                + arc.22.2 + arc.27.2 + arc.30.2 + arc.42.2 + arc.45.2
                + arc.46.2 + arc.49.2 + arc.52.2 + arc.64.2 + arc.65.2
                + arc.66.2 <= 15
 c18:           36 vertex.3 + arc.3.5 + arc.3.24 + arc.3.48 + arc.3.62
                + arc.3.11 + arc.3.25 + arc.3.35 + arc.3.0 + arc.3.8 + arc.3.18
                + arc.3.61 + arc.3.56 + arc.3.37 + arc.3.57 + arc.3.4
                + arc.3.53 + arc.3.43 + arc.3.49 + arc.3.55 + arc.3.6
                + arc.3.63 + arc.3.23 + arc.3.51 + arc.3.59 + arc.3.47
                + arc.3.2 + arc.3.32 + arc.3.50 + arc.3.29 + arc.3.36
                + arc.3.19 + arc.3.30 + arc.3.44 + arc.3.7 + arc.3.9 + arc.3.38
                <= 36
 c19:           21 vertex.4 + arc.3.4 + arc.4.35 + arc.4.18 + arc.4.58
                + arc.4.48 + arc.4.38 + arc.4.50 + arc.4.36 + arc.13.4
                + arc.15.4 + arc.20.4 + arc.22.4 + arc.30.4 + arc.42.4
                + arc.45.4 + arc.46.4 + arc.49.4 + arc.52.4 + arc.64.4
                + arc.65.4 + arc.66.4 <= 21
 c20:           21 vertex.5 + arc.3.5 + arc.5.35 + arc.5.18 + arc.5.48
                + arc.5.58 + arc.5.38 + arc.5.36 + arc.13.5 + arc.15.5
                + arc.20.5 + arc.22.5 + arc.27.5 + arc.30.5 + arc.42.5
                + arc.45.5 + arc.46.5 + arc.49.5 + arc.52.5 + arc.64.5
                + arc.65.5 + arc.66.5 <= 21
 c21:           21 vertex.6 + arc.6.20 + arc.45.6 + arc.3.6 + arc.6.63
                + arc.6.19 + arc.6.39 + arc.6.21 + arc.6.14 + arc.6.32
                + arc.6.37 + arc.13.6 + arc.15.6 + arc.20.6 + arc.22.6
                + arc.27.6 + arc.42.6 + arc.46.6 + arc.49.6 + arc.52.6
                + arc.64.6 + arc.66.6 <= 21
 c22:           16 vertex.7 + arc.3.7 + arc.7.20 + arc.7.19 + arc.7.42
                + arc.7.37 + arc.7.21 + arc.7.10 + arc.13.7 + arc.15.7
                + arc.27.7 + arc.45.7 + arc.46.7 + arc.49.7 + arc.64.7
                + arc.65.7 + arc.66.7 <= 16
 c23:           13 vertex.8 + arc.3.8 + arc.13.8 + arc.20.8 + arc.22.8
                + arc.27.8 + arc.30.8 + arc.42.8 + arc.45.8 + arc.46.8
                + arc.49.8 + arc.64.8 + arc.65.8 + arc.66.8 <= 13
 c24:           20 vertex.9 + arc.3.9 + arc.9.37 + arc.9.19 + arc.9.63
                + arc.9.39 + arc.9.14 + arc.9.32 + arc.13.9 + arc.15.9
                + arc.20.9 + arc.22.9 + arc.27.9 + arc.30.9 + arc.45.9
                + arc.46.9 + arc.49.9 + arc.52.9 + arc.64.9 + arc.65.9
                + arc.66.9 <= 20
 c25:           23 vertex.10 + arc.10.32 + arc.14.10 + arc.7.10 + arc.10.14
                + arc.10.37 + arc.10.63 + arc.10.57 + arc.10.19 + arc.12.10
                + arc.15.10 + arc.16.10 + arc.22.10 + arc.26.10 + arc.27.10
                + arc.31.10 + arc.46.10 + arc.48.10 + arc.50.10 + arc.52.10
                + arc.54.10 + arc.56.10 + arc.58.10 + arc.59.10 <= 23
 c26:           21 vertex.11 + arc.3.11 + arc.11.1 + arc.11.37 + arc.11.19
                + arc.11.16 + arc.11.32 + arc.11.14 + arc.13.11 + arc.15.11
                + arc.20.11 + arc.22.11 + arc.27.11 + arc.30.11 + arc.42.11
                + arc.45.11 + arc.46.11 + arc.49.11 + arc.52.11 + arc.64.11
                + arc.65.11 + arc.66.11 <= 21
 c27:           11 vertex.12 + arc.12.51 + arc.18.12 + arc.12.10 + arc.12.21
                + arc.12.19 + arc.12.37 + arc.12.32 + arc.12.63 + arc.12.14
                + arc.27.12 + arc.50.12 <= 11
 c28:           39 vertex.13 + arc.13.26 + arc.13.44 + arc.13.6 + arc.13.50
                + arc.13.48 + arc.13.8 + arc.13.32 + arc.13.18 + arc.13.14
                + arc.13.36 + arc.13.0 + arc.13.54 + arc.13.4 + arc.13.62
                + arc.13.60 + arc.13.34 + arc.13.56 + arc.13.2 + arc.13.11
                + arc.13.39 + arc.13.38 + arc.13.29 + arc.13.35 + arc.13.21
                + arc.13.33 + arc.13.37 + arc.13.43 + arc.13.25 + arc.13.9
                + arc.13.41 + arc.13.61 + arc.13.1 + arc.13.59 + arc.13.5
                + arc.13.53 + arc.13.19 + arc.13.63 + arc.13.7 + arc.30.13
                <= 39
 c29:           arc.32.14 + 29 vertex.14 + arc.14.10 + arc.6.14 + arc.9.14
                + arc.10.14 + arc.11.14 + arc.12.14 + arc.13.14 + arc.14.19
                + arc.14.51 + arc.14.63 + arc.14.57 + arc.14.37 + arc.15.14
                + arc.20.14 + arc.22.14 + arc.25.14 + arc.27.14 + arc.30.14
                + arc.35.14 + arc.36.14 + arc.41.14 + arc.43.14 + arc.44.14
                + arc.45.14 + arc.46.14 + arc.48.14 + arc.52.14 + arc.58.14
                <= 29
 c30:           35 vertex.15 + arc.15.36 + arc.15.4 + arc.15.32 + arc.15.0
                + arc.15.58 + arc.15.51 + arc.15.49 + arc.15.2 + arc.15.55
                + arc.15.10 + arc.15.14 + arc.15.57 + arc.15.19 + arc.15.43
                + arc.15.11 + arc.15.21 + arc.15.37 + arc.15.61 + arc.15.29
                + arc.15.35 + arc.15.24 + arc.15.60 + arc.15.7 + arc.15.48
                + arc.15.62 + arc.15.47 + arc.15.6 + arc.15.5 + arc.15.9
                + arc.15.44 + arc.15.56 + arc.15.34 + arc.15.63 + arc.15.25
                + arc.15.54 <= 35
 c31:           9 vertex.16 + arc.16.21 + arc.36.16 + arc.11.16 + arc.16.32
                + arc.16.10 + arc.16.20 + arc.16.19 + arc.16.37 + arc.59.16
                <= 9
 c32:           6 vertex.17 + arc.17.48 + arc.17.40 + arc.17.36 + arc.17.18
                + arc.17.35 + arc.22.17 <= 6
 c33:           arc.51.18 + 34 vertex.18 + arc.18.12 + arc.0.18 + arc.3.18
                + arc.4.18 + arc.5.18 + arc.13.18 + arc.17.18 + arc.18.19
                + arc.18.37 + arc.18.63 + arc.18.57 + arc.20.18 + arc.21.18
                + arc.22.18 + arc.23.18 + arc.24.18 + arc.27.18 + arc.29.18
                + arc.30.18 + arc.37.18 + arc.38.18 + arc.39.18 + arc.42.18
                + arc.45.18 + arc.46.18 + arc.47.18 + arc.49.18 + arc.52.18
                + arc.53.18 + arc.60.18 + arc.62.18 + arc.63.18 + arc.64.18
                <= 34
 c34:           40 vertex.19 + arc.1.19 + arc.3.19 + arc.6.19 + arc.7.19
                + arc.9.19 + arc.10.19 + arc.11.19 + arc.12.19 + arc.13.19
                + arc.14.19 + arc.15.19 + arc.16.19 + arc.18.19 + arc.20.19
                + arc.22.19 + arc.25.19 + arc.26.19 + arc.27.19 + arc.30.19
                + arc.31.19 + arc.32.19 + arc.33.19 + arc.34.19 + arc.35.19
                + arc.36.19 + arc.41.19 + arc.42.19 + arc.43.19 + arc.44.19
                + arc.45.19 + arc.46.19 + arc.48.19 + arc.49.19 + arc.50.19
                + arc.52.19 + arc.54.19 + arc.55.19 + arc.56.19 + arc.58.19
                + arc.59.19 <= 40
 c35:           arc.6.20 + 48 vertex.20 + arc.20.45 + arc.7.20 + arc.16.20
                + arc.20.61 + arc.20.28 + arc.20.0 + arc.20.24 + arc.20.34
                + arc.20.4 + arc.20.60 + arc.20.56 + arc.20.18 + arc.20.2
                + arc.20.6 + arc.20.42 + arc.20.8 + arc.20.44 + arc.20.37
                + arc.20.9 + arc.20.50 + arc.20.35 + arc.20.19 + arc.20.29
                + arc.20.14 + arc.20.5 + arc.20.36 + arc.20.53 + arc.20.11
                + arc.20.47 + arc.20.43 + arc.20.63 + arc.20.25 + arc.27.20
                + arc.31.20 + arc.32.20 + arc.34.20 + arc.36.20 + arc.44.20
                + arc.45.20 + arc.46.20 + arc.48.20 + arc.49.20 + arc.52.20
                + arc.54.20 + arc.55.20 + arc.58.20 + arc.59.20 <= 48
 c36:           arc.16.21 + 24 vertex.21 + arc.21.36 + arc.6.21 + arc.7.21
                + arc.12.21 + arc.13.21 + arc.15.21 + arc.21.18 + arc.21.35
                + arc.21.48 + arc.21.38 + arc.22.21 + arc.30.21 + arc.34.21
                + arc.36.21 + arc.43.21 + arc.44.21 + arc.45.21 + arc.50.21
                + arc.52.21 + arc.54.21 + arc.55.21 + arc.56.21 + arc.58.21
                <= 24
 c37:           38 vertex.22 + arc.22.37 + arc.22.55 + arc.22.19 + arc.22.61
                + arc.22.50 + arc.22.59 + arc.22.43 + arc.22.11 + arc.22.6
                + arc.22.49 + arc.22.63 + arc.22.4 + arc.22.8 + arc.22.0
                + arc.22.58 + arc.22.18 + arc.22.62 + arc.22.48 + arc.22.44
                + arc.22.54 + arc.22.26 + arc.22.14 + arc.22.34 + arc.22.56
                + arc.22.36 + arc.22.60 + arc.22.30 + arc.22.38 + arc.22.2
                + arc.22.10 + arc.22.5 + arc.22.9 + arc.22.21 + arc.22.25
                + arc.22.53 + arc.22.35 + arc.22.17 + arc.27.22 <= 38
 c38:           12 vertex.23 + arc.3.23 + arc.23.35 + arc.23.48 + arc.23.36
                + arc.23.18 + arc.23.58 + arc.30.23 + arc.42.23 + arc.46.23
                + arc.52.23 + arc.64.23 + arc.66.23 <= 12
 c39:           17 vertex.24 + arc.3.24 + arc.15.24 + arc.20.24 + arc.24.38
                + arc.24.48 + arc.24.36 + arc.24.40 + arc.24.58 + arc.24.18
                + arc.24.35 + arc.27.24 + arc.30.24 + arc.42.24 + arc.45.24
                + arc.49.24 + arc.65.24 + arc.66.24 <= 17
 c40:           18 vertex.25 + arc.3.25 + arc.13.25 + arc.15.25 + arc.20.25
                + arc.22.25 + arc.25.37 + arc.25.14 + arc.25.19 + arc.25.63
                + arc.27.25 + arc.30.25 + arc.42.25 + arc.45.25 + arc.46.25
                + arc.49.25 + arc.52.25 + arc.64.25 + arc.65.25 <= 18
 c41:           11 vertex.26 + arc.13.26 + arc.22.26 + arc.26.57 + arc.26.19
                + arc.26.37 + arc.26.63 + arc.26.10 + arc.26.42 + arc.26.46
                + arc.26.32 + arc.30.26 <= 11
 c42:           38 vertex.27 + arc.27.9 + arc.27.37 + arc.27.59 + arc.27.61
                + arc.27.43 + arc.27.53 + arc.27.57 + arc.27.33 + arc.27.36
                + arc.27.22 + arc.27.55 + arc.27.32 + arc.27.19 + arc.27.34
                + arc.27.20 + arc.27.63 + arc.27.0 + arc.27.11 + arc.27.38
                + arc.27.24 + arc.27.8 + arc.27.6 + arc.27.44 + arc.27.10
                + arc.27.2 + arc.27.7 + arc.27.62 + arc.27.12 + arc.27.25
                + arc.27.56 + arc.27.18 + arc.27.29 + arc.27.48 + arc.27.5
                + arc.27.35 + arc.27.60 + arc.27.14 + arc.27.42 <= 38
 c43:           2 vertex.28 + arc.20.28 + arc.65.28 <= 2
 c44:           18 vertex.29 + arc.3.29 + arc.13.29 + arc.15.29 + arc.20.29
                + arc.27.29 + arc.29.35 + arc.29.48 + arc.29.36 + arc.29.18
                + arc.29.50 + arc.42.29 + arc.45.29 + arc.46.29 + arc.49.29
                + arc.52.29 + arc.64.29 + arc.65.29 + arc.66.29 <= 18
 c45:           43 vertex.30 + arc.30.52 + arc.52.30 + arc.3.30 + arc.22.30
                + arc.30.54 + arc.30.48 + arc.30.50 + arc.30.2 + arc.30.44
                + arc.30.56 + arc.30.18 + arc.30.60 + arc.30.26 + arc.30.34
                + arc.30.61 + arc.30.9 + arc.30.62 + arc.30.14 + arc.30.43
                + arc.30.63 + arc.30.5 + arc.30.53 + arc.30.13 + arc.30.35
                + arc.30.21 + arc.30.47 + arc.30.11 + arc.30.19 + arc.30.23
                + arc.30.36 + arc.30.55 + arc.30.8 + arc.30.24 + arc.30.58
                + arc.30.4 + arc.30.37 + arc.30.0 + arc.30.25 + arc.42.30
                + arc.46.30 + arc.49.30 + arc.64.30 + arc.65.30 <= 43
 c46:           9 vertex.31 + arc.31.46 + arc.31.10 + arc.31.32 + arc.31.20
                + arc.31.1 + arc.31.63 + arc.31.37 + arc.31.19 + arc.42.31 <= 9
 c47:           arc.10.32 + 26 vertex.32 + arc.32.14 + arc.3.32 + arc.6.32
                + arc.9.32 + arc.11.32 + arc.12.32 + arc.13.32 + arc.15.32
                + arc.16.32 + arc.26.32 + arc.27.32 + arc.31.32 + arc.32.19
                + arc.32.20 + arc.32.57 + arc.34.32 + arc.35.32 + arc.36.32
                + arc.41.32 + arc.44.32 + arc.50.32 + arc.54.32 + arc.55.32
                + arc.56.32 + arc.58.32 <= 26
 c48:           8 vertex.33 + arc.13.33 + arc.27.33 + arc.33.63 + arc.33.37
                + arc.33.19 + arc.33.42 + arc.42.33 + arc.52.33 <= 8
 c49:           21 vertex.34 + arc.34.46 + arc.46.34 + arc.13.34 + arc.15.34
                + arc.20.34 + arc.22.34 + arc.27.34 + arc.30.34 + arc.34.20
                + arc.34.21 + arc.34.19 + arc.34.37 + arc.34.63 + arc.34.32
                + arc.42.34 + arc.45.34 + arc.49.34 + arc.52.34 + arc.64.34
                + arc.65.34 + arc.66.34 <= 21
 c50:           arc.37.35 + 39 vertex.35 + arc.35.1 + arc.0.35 + arc.3.35
                + arc.4.35 + arc.5.35 + arc.13.35 + arc.15.35 + arc.17.35
                + arc.20.35 + arc.21.35 + arc.22.35 + arc.23.35 + arc.24.35
                + arc.27.35 + arc.29.35 + arc.30.35 + arc.35.14 + arc.35.32
                + arc.35.37 + arc.35.42 + arc.35.57 + arc.35.19 + arc.35.63
                + arc.38.35 + arc.39.35 + arc.42.35 + arc.45.35 + arc.46.35
                + arc.47.35 + arc.49.35 + arc.51.35 + arc.52.35 + arc.53.35
                + arc.57.35 + arc.60.35 + arc.62.35 + arc.63.35 + arc.64.35
                <= 39
 c51:           arc.21.36 + 42 vertex.36 + arc.36.16 + arc.0.36 + arc.3.36
                + arc.4.36 + arc.5.36 + arc.13.36 + arc.15.36 + arc.17.36
                + arc.20.36 + arc.22.36 + arc.23.36 + arc.24.36 + arc.27.36
                + arc.29.36 + arc.30.36 + arc.36.20 + arc.36.39 + arc.36.21
                + arc.36.14 + arc.36.42 + arc.36.19 + arc.36.41 + arc.36.37
                + arc.36.63 + arc.36.32 + arc.37.36 + arc.38.36 + arc.39.36
                + arc.42.36 + arc.45.36 + arc.46.36 + arc.47.36 + arc.49.36
                + arc.51.36 + arc.52.36 + arc.53.36 + arc.57.36 + arc.60.36
                + arc.62.36 + arc.63.36 + arc.64.36 <= 42
 c52:           arc.1.37 + 43 vertex.37 + arc.37.35 + arc.3.37 + arc.6.37
                + arc.7.37 + arc.9.37 + arc.10.37 + arc.11.37 + arc.12.37
                + arc.13.37 + arc.14.37 + arc.15.37 + arc.16.37 + arc.18.37
                + arc.20.37 + arc.22.37 + arc.25.37 + arc.26.37 + arc.27.37
                + arc.30.37 + arc.31.37 + arc.33.37 + arc.34.37 + arc.35.37
                + arc.36.37 + arc.37.36 + arc.37.18 + arc.37.48 + arc.41.37
                + arc.42.37 + arc.43.37 + arc.44.37 + arc.45.37 + arc.46.37
                + arc.48.37 + arc.49.37 + arc.50.37 + arc.52.37 + arc.54.37
                + arc.55.37 + arc.56.37 + arc.58.37 + arc.59.37 <= 43
 c53:           21 vertex.38 + arc.38.48 + arc.57.38 + arc.3.38 + arc.4.38
                + arc.5.38 + arc.13.38 + arc.21.38 + arc.22.38 + arc.24.38
                + arc.27.38 + arc.38.35 + arc.38.36 + arc.38.18 + arc.38.58
                + arc.45.38 + arc.46.38 + arc.49.38 + arc.51.38 + arc.53.38
                + arc.60.38 + arc.63.38 <= 21
 c54:           13 vertex.39 + arc.39.58 + arc.58.39 + arc.1.39 + arc.6.39
                + arc.9.39 + arc.13.39 + arc.36.39 + arc.39.36 + arc.39.35
                + arc.39.18 + arc.42.39 + arc.43.39 + arc.44.39 <= 13
 c55:           3 vertex.40 + arc.17.40 + arc.24.40 + arc.51.40 <= 3
 c56:           9 vertex.41 + arc.13.41 + arc.36.41 + arc.41.32 + arc.41.63
                + arc.41.57 + arc.41.19 + arc.41.14 + arc.41.37 + arc.45.41
                <= 9
 c57:           48 vertex.42 + arc.42.62 + arc.50.42 + arc.7.42 + arc.20.42
                + arc.26.42 + arc.27.42 + arc.33.42 + arc.35.42 + arc.36.42
                + arc.42.5 + arc.42.56 + arc.42.44 + arc.42.30 + arc.42.24
                + arc.42.35 + arc.42.50 + arc.42.58 + arc.42.48 + arc.42.54
                + arc.42.36 + arc.42.18 + arc.42.8 + arc.42.33 + arc.42.25
                + arc.42.37 + arc.42.4 + arc.42.0 + arc.42.31 + arc.42.39
                + arc.42.23 + arc.42.53 + arc.42.63 + arc.42.2 + arc.42.6
                + arc.42.11 + arc.42.29 + arc.42.61 + arc.42.19 + arc.42.47
                + arc.42.59 + arc.42.43 + arc.42.60 + arc.42.34 + arc.45.42
                + arc.52.42 + arc.54.42 + arc.58.42 + arc.59.42 <= 48
 c58:           22 vertex.43 + arc.3.43 + arc.13.43 + arc.15.43 + arc.20.43
                + arc.22.43 + arc.27.43 + arc.30.43 + arc.42.43 + arc.43.37
                + arc.43.63 + arc.43.21 + arc.43.57 + arc.43.19 + arc.43.39
                + arc.43.14 + arc.45.43 + arc.46.43 + arc.49.43 + arc.52.43
                + arc.64.43 + arc.65.43 + arc.66.43 <= 22
 c59:           23 vertex.44 + arc.3.44 + arc.13.44 + arc.15.44 + arc.20.44
                + arc.22.44 + arc.27.44 + arc.30.44 + arc.42.44 + arc.44.63
                + arc.44.20 + arc.44.14 + arc.44.32 + arc.44.37 + arc.44.19
                + arc.44.21 + arc.44.39 + arc.45.44 + arc.46.44 + arc.49.44
                + arc.52.44 + arc.64.44 + arc.65.44 + arc.66.44 <= 23
 c60:           arc.20.45 + 40 vertex.45 + arc.45.6 + arc.45.19 + arc.45.44
                + arc.45.29 + arc.45.0 + arc.45.35 + arc.45.48 + arc.45.8
                + arc.45.62 + arc.45.2 + arc.45.4 + arc.45.37 + arc.45.57
                + arc.45.11 + arc.45.42 + arc.45.60 + arc.45.59 + arc.45.54
                + arc.45.14 + arc.45.36 + arc.45.18 + arc.45.50 + arc.45.34
                + arc.45.58 + arc.45.43 + arc.45.61 + arc.45.5 + arc.45.7
                + arc.45.56 + arc.45.41 + arc.45.52 + arc.45.38 + arc.45.20
                + arc.45.9 + arc.45.24 + arc.45.63 + arc.45.21 + arc.45.25
                + arc.45.53 <= 40
 c61:           arc.34.46 + 40 vertex.46 + arc.46.34 + arc.26.46 + arc.31.46
                + arc.46.19 + arc.46.6 + arc.46.56 + arc.46.18 + arc.46.36
                + arc.46.30 + arc.46.60 + arc.46.5 + arc.46.44 + arc.46.48
                + arc.46.9 + arc.46.0 + arc.46.4 + arc.46.62 + arc.46.58
                + arc.46.54 + arc.46.10 + arc.46.14 + arc.46.2 + arc.46.23
                + arc.46.8 + arc.46.37 + arc.46.29 + arc.46.25 + arc.46.59
                + arc.46.55 + arc.46.11 + arc.46.7 + arc.46.35 + arc.46.63
                + arc.46.53 + arc.46.61 + arc.46.20 + arc.46.38 + arc.46.43
                + arc.52.46 <= 40
 c62:           12 vertex.47 + arc.3.47 + arc.15.47 + arc.20.47 + arc.30.47
                + arc.42.47 + arc.47.18 + arc.47.35 + arc.47.58 + arc.47.48
                + arc.47.36 + arc.65.47 + arc.66.47 <= 12
 c63:           arc.38.48 + 36 vertex.48 + arc.48.57 + arc.0.48 + arc.3.48
                + arc.4.48 + arc.5.48 + arc.13.48 + arc.15.48 + arc.17.48
                + arc.21.48 + arc.22.48 + arc.23.48 + arc.24.48 + arc.27.48
                + arc.29.48 + arc.30.48 + arc.37.48 + arc.42.48 + arc.45.48
                + arc.46.48 + arc.47.48 + arc.48.20 + arc.48.63 + arc.48.19
                + arc.48.37 + arc.48.14 + arc.48.10 + arc.49.48 + arc.51.48
                + arc.52.48 + arc.53.48 + arc.57.48 + arc.60.48 + arc.62.48
                + arc.63.48 + arc.64.48 <= 36
 c64:           38 vertex.49 + arc.3.49 + arc.15.49 + arc.22.49 + arc.49.48
                + arc.49.56 + arc.49.20 + arc.49.44 + arc.49.8 + arc.49.6
                + arc.49.24 + arc.49.36 + arc.49.30 + arc.49.60 + arc.49.38
                + arc.49.35 + arc.49.5 + arc.49.53 + arc.49.61 + arc.49.50
                + arc.49.29 + arc.49.63 + arc.49.58 + arc.49.7 + arc.49.43
                + arc.49.25 + arc.49.57 + arc.49.59 + arc.49.37 + arc.49.0
                + arc.49.19 + arc.49.9 + arc.49.34 + arc.49.11 + arc.49.18
                + arc.49.4 + arc.49.54 + arc.49.2 + arc.64.49 <= 38
 c65:           arc.62.50 + 22 vertex.50 + arc.50.42 + arc.0.50 + arc.3.50
                + arc.4.50 + arc.13.50 + arc.20.50 + arc.22.50 + arc.29.50
                + arc.30.50 + arc.42.50 + arc.45.50 + arc.49.50 + arc.50.10
                + arc.50.63 + arc.50.32 + arc.50.37 + arc.50.21 + arc.50.19
                + arc.50.12 + arc.52.50 + arc.53.50 <= 22
 c66:           arc.12.51 + 11 vertex.51 + arc.51.18 + arc.3.51 + arc.14.51
                + arc.15.51 + arc.51.40 + arc.51.35 + arc.51.48 + arc.51.58
                + arc.51.38 + arc.51.36 <= 11
 c67:           arc.30.52 + 39 vertex.52 + arc.52.30 + arc.45.52 + arc.52.4
                + arc.52.21 + arc.52.0 + arc.52.60 + arc.52.43 + arc.52.53
                + arc.52.55 + arc.52.2 + arc.52.63 + arc.52.19 + arc.52.61
                + arc.52.20 + arc.52.36 + arc.52.34 + arc.52.58 + arc.52.6
                + arc.52.50 + arc.52.5 + arc.52.14 + arc.52.11 + arc.52.44
                + arc.52.56 + arc.52.54 + arc.52.48 + arc.52.46 + arc.52.25
                + arc.52.35 + arc.52.42 + arc.52.18 + arc.52.9 + arc.52.29
                + arc.52.10 + arc.52.33 + arc.52.59 + arc.52.23 + arc.52.37
                <= 39
 c68:           21 vertex.53 + arc.3.53 + arc.13.53 + arc.20.53 + arc.22.53
                + arc.27.53 + arc.30.53 + arc.42.53 + arc.45.53 + arc.46.53
                + arc.49.53 + arc.52.53 + arc.53.35 + arc.53.58 + arc.53.36
                + arc.53.18 + arc.53.50 + arc.53.38 + arc.53.48 + arc.64.53
                + arc.65.53 + arc.66.53 <= 21
 c69:           20 vertex.54 + arc.13.54 + arc.15.54 + arc.22.54 + arc.30.54
                + arc.42.54 + arc.45.54 + arc.46.54 + arc.49.54 + arc.52.54
                + arc.54.32 + arc.54.42 + arc.54.21 + arc.54.19 + arc.54.37
                + arc.54.10 + arc.54.57 + arc.54.63 + arc.54.20 + arc.64.54
                + arc.65.54 <= 20
 c70:           16 vertex.55 + arc.3.55 + arc.15.55 + arc.22.55 + arc.27.55
                + arc.30.55 + arc.46.55 + arc.52.55 + arc.55.1 + arc.55.20
                + arc.55.32 + arc.55.19 + arc.55.37 + arc.55.21 + arc.55.63
                + arc.55.57 + arc.64.55 <= 16
 c71:           21 vertex.56 + arc.3.56 + arc.13.56 + arc.15.56 + arc.20.56
                + arc.22.56 + arc.27.56 + arc.30.56 + arc.42.56 + arc.45.56
                + arc.46.56 + arc.49.56 + arc.52.56 + arc.56.19 + arc.56.63
                + arc.56.21 + arc.56.57 + arc.56.10 + arc.56.32 + arc.56.37
                + arc.65.56 + arc.66.56 <= 21
 c72:           arc.48.57 + 23 vertex.57 + arc.57.38 + arc.1.57 + arc.3.57
                + arc.10.57 + arc.14.57 + arc.15.57 + arc.18.57 + arc.26.57
                + arc.27.57 + arc.32.57 + arc.35.57 + arc.41.57 + arc.43.57
                + arc.45.57 + arc.49.57 + arc.54.57 + arc.55.57 + arc.56.57
                + arc.57.36 + arc.57.48 + arc.57.35 + arc.59.57 <= 23
 c73:           arc.39.58 + 31 vertex.58 + arc.58.39 + arc.0.58 + arc.4.58
                + arc.5.58 + arc.15.58 + arc.22.58 + arc.23.58 + arc.24.58
                + arc.30.58 + arc.38.58 + arc.42.58 + arc.45.58 + arc.46.58
                + arc.47.58 + arc.49.58 + arc.51.58 + arc.52.58 + arc.53.58
                + arc.58.20 + arc.58.14 + arc.58.10 + arc.58.32 + arc.58.19
                + arc.58.63 + arc.58.21 + arc.58.37 + arc.58.42 + arc.60.58
                + arc.62.58 + arc.64.58 <= 31
 c74:           18 vertex.59 + arc.3.59 + arc.13.59 + arc.22.59 + arc.27.59
                + arc.42.59 + arc.45.59 + arc.46.59 + arc.49.59 + arc.52.59
                + arc.59.20 + arc.59.16 + arc.59.19 + arc.59.57 + arc.59.10
                + arc.59.37 + arc.59.63 + arc.59.42 + arc.65.59 <= 18
 c75:           20 vertex.60 + arc.13.60 + arc.15.60 + arc.20.60 + arc.22.60
                + arc.27.60 + arc.30.60 + arc.42.60 + arc.45.60 + arc.46.60
                + arc.49.60 + arc.52.60 + arc.60.18 + arc.60.35 + arc.60.58
                + arc.60.36 + arc.60.48 + arc.60.38 + arc.64.60 + arc.65.60
                + arc.66.60 <= 20
 c76:           15 vertex.61 + arc.3.61 + arc.13.61 + arc.15.61 + arc.20.61
                + arc.22.61 + arc.27.61 + arc.30.61 + arc.42.61 + arc.45.61
                + arc.46.61 + arc.49.61 + arc.52.61 + arc.64.61 + arc.65.61
                + arc.66.61 <= 15
 c77:           arc.42.62 + 18 vertex.62 + arc.62.50 + arc.3.62 + arc.13.62
                + arc.15.62 + arc.22.62 + arc.27.62 + arc.30.62 + arc.45.62
                + arc.46.62 + arc.62.48 + arc.62.58 + arc.62.36 + arc.62.35
                + arc.62.18 + arc.64.62 + arc.65.62 + arc.66.62 <= 18
 c78:           41 vertex.63 + arc.1.63 + arc.3.63 + arc.6.63 + arc.9.63
                + arc.10.63 + arc.12.63 + arc.13.63 + arc.14.63 + arc.15.63
                + arc.18.63 + arc.20.63 + arc.22.63 + arc.25.63 + arc.26.63
                + arc.27.63 + arc.30.63 + arc.31.63 + arc.33.63 + arc.34.63
                + arc.35.63 + arc.36.63 + arc.41.63 + arc.42.63 + arc.43.63
                + arc.44.63 + arc.45.63 + arc.46.63 + arc.48.63 + arc.49.63
                + arc.50.63 + arc.52.63 + arc.54.63 + arc.55.63 + arc.56.63
                + arc.58.63 + arc.59.63 + arc.63.38 + arc.63.48 + arc.63.35
                + arc.63.18 + arc.63.36 <= 41
Bounds
 0 <= r.0 <= 1
 0 <= r.1 <= 1
 0 <= r.2 <= 1
 0 <= r.3 <= 1
 0 <= r.4 <= 1
 0 <= r.5 <= 1
 0 <= r.6 <= 1
 0 <= r.7 <= 1
 0 <= r.8 <= 1
 0 <= r.9 <= 1
 0 <= vertex.1 <= 1
 0 <= arc.1.37 <= 1
 0 <= vertex.37 <= 1
 0 <= arc.37.35 <= 1
 0 <= vertex.35 <= 1
 0 <= arc.35.1 <= 1
 0 <= vertex.6 <= 1
 0 <= arc.6.20 <= 1
 0 <= vertex.20 <= 1
 0 <= arc.20.45 <= 1
 0 <= vertex.45 <= 1
 0 <= arc.45.6 <= 1
 0 <= vertex.10 <= 1
 0 <= arc.10.32 <= 1
 0 <= vertex.32 <= 1
 0 <= arc.32.14 <= 1
 0 <= vertex.14 <= 1
 0 <= arc.14.10 <= 1
 0 <= vertex.12 <= 1
 0 <= arc.12.51 <= 1
 0 <= vertex.51 <= 1
 0 <= arc.51.18 <= 1
 0 <= vertex.18 <= 1
 0 <= arc.18.12 <= 1
 0 <= vertex.16 <= 1
 0 <= arc.16.21 <= 1
 0 <= vertex.21 <= 1
 0 <= arc.21.36 <= 1
 0 <= vertex.36 <= 1
 0 <= arc.36.16 <= 1
 0 <= vertex.30 <= 1
 0 <= arc.30.52 <= 1
 0 <= vertex.52 <= 1
 0 <= arc.52.30 <= 1
 0 <= vertex.34 <= 1
 0 <= arc.34.46 <= 1
 0 <= vertex.46 <= 1
 0 <= arc.46.34 <= 1
 0 <= vertex.38 <= 1
 0 <= arc.38.48 <= 1
 0 <= vertex.48 <= 1
 0 <= arc.48.57 <= 1
 0 <= vertex.57 <= 1
 0 <= arc.57.38 <= 1
 0 <= vertex.39 <= 1
 0 <= arc.39.58 <= 1
 0 <= vertex.58 <= 1
 0 <= arc.58.39 <= 1
 0 <= vertex.42 <= 1
 0 <= arc.42.62 <= 1
 0 <= vertex.62 <= 1
 0 <= arc.62.50 <= 1
 0 <= vertex.50 <= 1
 0 <= arc.50.42 <= 1
 0 <= vertex.0 <= 1
 0 <= vertex.2 <= 1
 0 <= vertex.3 <= 1
 0 <= vertex.4 <= 1
 0 <= vertex.5 <= 1
 0 <= vertex.7 <= 1
 0 <= vertex.8 <= 1
 0 <= vertex.9 <= 1
 0 <= vertex.11 <= 1
 0 <= vertex.13 <= 1
 0 <= vertex.15 <= 1
 0 <= vertex.17 <= 1
 0 <= vertex.19 <= 1
 0 <= vertex.22 <= 1
 0 <= vertex.23 <= 1
 0 <= vertex.24 <= 1
 0 <= vertex.25 <= 1
 0 <= vertex.26 <= 1
 0 <= vertex.27 <= 1
 0 <= vertex.28 <= 1
 0 <= vertex.29 <= 1
 0 <= vertex.31 <= 1
 0 <= vertex.33 <= 1
 0 <= vertex.40 <= 1
 0 <= vertex.41 <= 1
 0 <= vertex.43 <= 1
 0 <= vertex.44 <= 1
 0 <= vertex.47 <= 1
 0 <= vertex.49 <= 1
 0 <= vertex.53 <= 1
 0 <= vertex.54 <= 1
 0 <= vertex.55 <= 1
 0 <= vertex.56 <= 1
 0 <= vertex.59 <= 1
 0 <= vertex.60 <= 1
 0 <= vertex.61 <= 1
 0 <= vertex.63 <= 1
 0 <= arc.0.36 <= 1
 0 <= arc.0.18 <= 1
 0 <= arc.0.50 <= 1
 0 <= arc.0.48 <= 1
 0 <= arc.0.58 <= 1
 0 <= arc.0.35 <= 1
 0 <= arc.1.63 <= 1
 0 <= arc.1.19 <= 1
 0 <= arc.1.57 <= 1
 0 <= arc.1.39 <= 1
 0 <= arc.3.5 <= 1
 0 <= arc.3.24 <= 1
 0 <= arc.3.48 <= 1
 0 <= arc.3.62 <= 1
 0 <= arc.3.11 <= 1
 0 <= arc.3.25 <= 1
 0 <= arc.3.35 <= 1
 0 <= arc.3.0 <= 1
 0 <= arc.3.8 <= 1
 0 <= arc.3.18 <= 1
 0 <= arc.3.61 <= 1
 0 <= arc.3.56 <= 1
 0 <= arc.3.37 <= 1
 0 <= arc.3.57 <= 1
 0 <= arc.3.4 <= 1
 0 <= arc.3.53 <= 1
 0 <= arc.3.43 <= 1
 0 <= arc.3.49 <= 1
 0 <= arc.3.55 <= 1
 0 <= arc.3.6 <= 1
 0 <= arc.3.63 <= 1
 0 <= arc.3.23 <= 1
 0 <= arc.3.51 <= 1
 0 <= arc.3.59 <= 1
 0 <= arc.3.47 <= 1
 0 <= arc.3.2 <= 1
 0 <= arc.3.32 <= 1
 0 <= arc.3.50 <= 1
 0 <= arc.3.29 <= 1
 0 <= arc.3.36 <= 1
 0 <= arc.3.19 <= 1
 0 <= arc.3.30 <= 1
 0 <= arc.3.44 <= 1
 0 <= arc.3.7 <= 1
 0 <= arc.3.9 <= 1
 0 <= arc.3.38 <= 1
 0 <= arc.4.35 <= 1
 0 <= arc.4.18 <= 1
 0 <= arc.4.58 <= 1
 0 <= arc.4.48 <= 1
 0 <= arc.4.38 <= 1
 0 <= arc.4.50 <= 1
 0 <= arc.4.36 <= 1
 0 <= arc.5.35 <= 1
 0 <= arc.5.18 <= 1
 0 <= arc.5.48 <= 1
 0 <= arc.5.58 <= 1
 0 <= arc.5.38 <= 1
 0 <= arc.5.36 <= 1
 0 <= arc.6.63 <= 1
 0 <= arc.6.19 <= 1
 0 <= arc.6.39 <= 1
 0 <= arc.6.21 <= 1
 0 <= arc.6.14 <= 1
 0 <= arc.6.32 <= 1
 0 <= arc.6.37 <= 1
 0 <= arc.7.20 <= 1
 0 <= arc.7.19 <= 1
 0 <= arc.7.42 <= 1
 0 <= arc.7.37 <= 1
 0 <= arc.7.21 <= 1
 0 <= arc.7.10 <= 1
 0 <= arc.9.37 <= 1
 0 <= arc.9.19 <= 1
 0 <= arc.9.63 <= 1
 0 <= arc.9.39 <= 1
 0 <= arc.9.14 <= 1
 0 <= arc.9.32 <= 1
 0 <= arc.10.14 <= 1
 0 <= arc.10.37 <= 1
 0 <= arc.10.63 <= 1
 0 <= arc.10.57 <= 1
 0 <= arc.10.19 <= 1
 0 <= arc.11.1 <= 1
 0 <= arc.11.37 <= 1
 0 <= arc.11.19 <= 1
 0 <= arc.11.16 <= 1
 0 <= arc.11.32 <= 1
 0 <= arc.11.14 <= 1
 0 <= arc.12.10 <= 1
 0 <= arc.12.21 <= 1
 0 <= arc.12.19 <= 1
 0 <= arc.12.37 <= 1
 0 <= arc.12.32 <= 1
 0 <= arc.12.63 <= 1
 0 <= arc.12.14 <= 1
 0 <= arc.13.26 <= 1
 0 <= arc.13.44 <= 1
 0 <= arc.13.6 <= 1
 0 <= arc.13.50 <= 1
 0 <= arc.13.48 <= 1
 0 <= arc.13.8 <= 1
 0 <= arc.13.32 <= 1
 0 <= arc.13.18 <= 1
 0 <= arc.13.14 <= 1
 0 <= arc.13.36 <= 1
 0 <= arc.13.0 <= 1
 0 <= arc.13.54 <= 1
 0 <= arc.13.4 <= 1
 0 <= arc.13.62 <= 1
 0 <= arc.13.60 <= 1
 0 <= arc.13.34 <= 1
 0 <= arc.13.56 <= 1
 0 <= arc.13.2 <= 1
 0 <= arc.13.11 <= 1
 0 <= arc.13.39 <= 1
 0 <= arc.13.38 <= 1
 0 <= arc.13.29 <= 1
 0 <= arc.13.35 <= 1
 0 <= arc.13.21 <= 1
 0 <= arc.13.33 <= 1
 0 <= arc.13.37 <= 1
 0 <= arc.13.43 <= 1
 0 <= arc.13.25 <= 1
 0 <= arc.13.9 <= 1
 0 <= arc.13.41 <= 1
 0 <= arc.13.61 <= 1
 0 <= arc.13.1 <= 1
 0 <= arc.13.59 <= 1
 0 <= arc.13.5 <= 1
 0 <= arc.13.53 <= 1
 0 <= arc.13.19 <= 1
 0 <= arc.13.63 <= 1
 0 <= arc.13.7 <= 1
 0 <= arc.14.19 <= 1
 0 <= arc.14.51 <= 1
 0 <= arc.14.63 <= 1
 0 <= arc.14.57 <= 1
 0 <= arc.14.37 <= 1
 0 <= arc.15.36 <= 1
 0 <= arc.15.4 <= 1
 0 <= arc.15.32 <= 1
 0 <= arc.15.0 <= 1
 0 <= arc.15.58 <= 1
 0 <= arc.15.51 <= 1
 0 <= arc.15.49 <= 1
 0 <= arc.15.2 <= 1
 0 <= arc.15.55 <= 1
 0 <= arc.15.10 <= 1
 0 <= arc.15.14 <= 1
 0 <= arc.15.57 <= 1
 0 <= arc.15.19 <= 1
 0 <= arc.15.43 <= 1
 0 <= arc.15.11 <= 1
 0 <= arc.15.21 <= 1
 0 <= arc.15.37 <= 1
 0 <= arc.15.61 <= 1
 0 <= arc.15.29 <= 1
 0 <= arc.15.35 <= 1
 0 <= arc.15.24 <= 1
 0 <= arc.15.60 <= 1
 0 <= arc.15.7 <= 1
 0 <= arc.15.48 <= 1
 0 <= arc.15.62 <= 1
 0 <= arc.15.47 <= 1
 0 <= arc.15.6 <= 1
 0 <= arc.15.5 <= 1
 0 <= arc.15.9 <= 1
 0 <= arc.15.44 <= 1
 0 <= arc.15.56 <= 1
 0 <= arc.15.34 <= 1
 0 <= arc.15.63 <= 1
 0 <= arc.15.25 <= 1
 0 <= arc.15.54 <= 1
 0 <= arc.16.32 <= 1
 0 <= arc.16.10 <= 1
 0 <= arc.16.20 <= 1
 0 <= arc.16.19 <= 1
 0 <= arc.16.37 <= 1
 0 <= arc.17.48 <= 1
 0 <= arc.17.40 <= 1
 0 <= arc.17.36 <= 1
 0 <= arc.17.18 <= 1
 0 <= arc.17.35 <= 1
 0 <= arc.18.19 <= 1
 0 <= arc.18.37 <= 1
 0 <= arc.18.63 <= 1
 0 <= arc.18.57 <= 1
 0 <= arc.20.61 <= 1
 0 <= arc.20.28 <= 1
 0 <= arc.20.0 <= 1
 0 <= arc.20.24 <= 1
 0 <= arc.20.34 <= 1
 0 <= arc.20.4 <= 1
 0 <= arc.20.60 <= 1
 0 <= arc.20.56 <= 1
 0 <= arc.20.18 <= 1
 0 <= arc.20.2 <= 1
 0 <= arc.20.6 <= 1
 0 <= arc.20.42 <= 1
 0 <= arc.20.8 <= 1
 0 <= arc.20.44 <= 1
 0 <= arc.20.37 <= 1
 0 <= arc.20.9 <= 1
 0 <= arc.20.50 <= 1
 0 <= arc.20.35 <= 1
 0 <= arc.20.19 <= 1
 0 <= arc.20.29 <= 1
 0 <= arc.20.14 <= 1
 0 <= arc.20.5 <= 1
 0 <= arc.20.36 <= 1
 0 <= arc.20.53 <= 1
 0 <= arc.20.11 <= 1
 0 <= arc.20.47 <= 1
 0 <= arc.20.43 <= 1
 0 <= arc.20.63 <= 1
 0 <= arc.20.25 <= 1
 0 <= arc.21.18 <= 1
 0 <= arc.21.35 <= 1
 0 <= arc.21.48 <= 1
 0 <= arc.21.38 <= 1
 0 <= arc.22.37 <= 1
 0 <= arc.22.55 <= 1
 0 <= arc.22.19 <= 1
 0 <= arc.22.61 <= 1
 0 <= arc.22.50 <= 1
 0 <= arc.22.59 <= 1
 0 <= arc.22.43 <= 1
 0 <= arc.22.11 <= 1
 0 <= arc.22.6 <= 1
 0 <= arc.22.49 <= 1
 0 <= arc.22.63 <= 1
 0 <= arc.22.4 <= 1
 0 <= arc.22.8 <= 1
 0 <= arc.22.0 <= 1
 0 <= arc.22.58 <= 1
 0 <= arc.22.18 <= 1
 0 <= arc.22.62 <= 1
 0 <= arc.22.48 <= 1
 0 <= arc.22.44 <= 1
 0 <= arc.22.54 <= 1
 0 <= arc.22.26 <= 1
 0 <= arc.22.14 <= 1
 0 <= arc.22.34 <= 1
 0 <= arc.22.56 <= 1
 0 <= arc.22.36 <= 1
 0 <= arc.22.60 <= 1
 0 <= arc.22.30 <= 1
 0 <= arc.22.38 <= 1
 0 <= arc.22.2 <= 1
 0 <= arc.22.10 <= 1
 0 <= arc.22.5 <= 1
 0 <= arc.22.9 <= 1
 0 <= arc.22.21 <= 1
 0 <= arc.22.25 <= 1
 0 <= arc.22.53 <= 1
 0 <= arc.22.35 <= 1
 0 <= arc.22.17 <= 1
 0 <= arc.23.35 <= 1
 0 <= arc.23.48 <= 1
 0 <= arc.23.36 <= 1
 0 <= arc.23.18 <= 1
 0 <= arc.23.58 <= 1
 0 <= arc.24.38 <= 1
 0 <= arc.24.48 <= 1
 0 <= arc.24.36 <= 1
 0 <= arc.24.40 <= 1
 0 <= arc.24.58 <= 1
 0 <= arc.24.18 <= 1
 0 <= arc.24.35 <= 1
 0 <= arc.25.37 <= 1
 0 <= arc.25.14 <= 1
 0 <= arc.25.19 <= 1
 0 <= arc.25.63 <= 1
 0 <= arc.26.57 <= 1
 0 <= arc.26.19 <= 1
 0 <= arc.26.37 <= 1
 0 <= arc.26.63 <= 1
 0 <= arc.26.10 <= 1
 0 <= arc.26.42 <= 1
 0 <= arc.26.46 <= 1
 0 <= arc.26.32 <= 1
 0 <= arc.27.9 <= 1
 0 <= arc.27.37 <= 1
 0 <= arc.27.59 <= 1
 0 <= arc.27.61 <= 1
 0 <= arc.27.43 <= 1
 0 <= arc.27.53 <= 1
 0 <= arc.27.57 <= 1
 0 <= arc.27.33 <= 1
 0 <= arc.27.36 <= 1
 0 <= arc.27.22 <= 1
 0 <= arc.27.55 <= 1
 0 <= arc.27.32 <= 1
 0 <= arc.27.19 <= 1
 0 <= arc.27.34 <= 1
 0 <= arc.27.20 <= 1
 0 <= arc.27.63 <= 1
 0 <= arc.27.0 <= 1
 0 <= arc.27.11 <= 1
 0 <= arc.27.38 <= 1
 0 <= arc.27.24 <= 1
 0 <= arc.27.8 <= 1
 0 <= arc.27.6 <= 1
 0 <= arc.27.44 <= 1
 0 <= arc.27.10 <= 1
 0 <= arc.27.2 <= 1
 0 <= arc.27.7 <= 1
 0 <= arc.27.62 <= 1
 0 <= arc.27.12 <= 1
 0 <= arc.27.25 <= 1
 0 <= arc.27.56 <= 1
 0 <= arc.27.18 <= 1
 0 <= arc.27.29 <= 1
 0 <= arc.27.48 <= 1
 0 <= arc.27.5 <= 1
 0 <= arc.27.35 <= 1
 0 <= arc.27.60 <= 1
 0 <= arc.27.14 <= 1
 0 <= arc.27.42 <= 1
 0 <= arc.29.35 <= 1
 0 <= arc.29.48 <= 1
 0 <= arc.29.36 <= 1
 0 <= arc.29.18 <= 1
 0 <= arc.29.50 <= 1
 0 <= arc.30.54 <= 1
 0 <= arc.30.48 <= 1
 0 <= arc.30.50 <= 1
 0 <= arc.30.2 <= 1
 0 <= arc.30.44 <= 1
 0 <= arc.30.56 <= 1
 0 <= arc.30.18 <= 1
 0 <= arc.30.60 <= 1
 0 <= arc.30.26 <= 1
 0 <= arc.30.34 <= 1
 0 <= arc.30.61 <= 1
 0 <= arc.30.9 <= 1
 0 <= arc.30.62 <= 1
 0 <= arc.30.14 <= 1
 0 <= arc.30.43 <= 1
 0 <= arc.30.63 <= 1
 0 <= arc.30.5 <= 1
 0 <= arc.30.53 <= 1
 0 <= arc.30.13 <= 1
 0 <= arc.30.35 <= 1
 0 <= arc.30.21 <= 1
 0 <= arc.30.47 <= 1
 0 <= arc.30.11 <= 1
 0 <= arc.30.19 <= 1
 0 <= arc.30.23 <= 1
 0 <= arc.30.36 <= 1
 0 <= arc.30.55 <= 1
 0 <= arc.30.8 <= 1
 0 <= arc.30.24 <= 1
 0 <= arc.30.58 <= 1
 0 <= arc.30.4 <= 1
 0 <= arc.30.37 <= 1
 0 <= arc.30.0 <= 1
 0 <= arc.30.25 <= 1
 0 <= arc.31.46 <= 1
 0 <= arc.31.10 <= 1
 0 <= arc.31.32 <= 1
 0 <= arc.31.20 <= 1
 0 <= arc.31.1 <= 1
 0 <= arc.31.63 <= 1
 0 <= arc.31.37 <= 1
 0 <= arc.31.19 <= 1
 0 <= arc.32.19 <= 1
 0 <= arc.32.20 <= 1
 0 <= arc.32.57 <= 1
 0 <= arc.33.63 <= 1
 0 <= arc.33.37 <= 1
 0 <= arc.33.19 <= 1
 0 <= arc.33.42 <= 1
 0 <= arc.34.20 <= 1
 0 <= arc.34.21 <= 1
 0 <= arc.34.19 <= 1
 0 <= arc.34.37 <= 1
 0 <= arc.34.63 <= 1
 0 <= arc.34.32 <= 1
 0 <= arc.35.14 <= 1
 0 <= arc.35.32 <= 1
 0 <= arc.35.37 <= 1
 0 <= arc.35.42 <= 1
 0 <= arc.35.57 <= 1
 0 <= arc.35.19 <= 1
 0 <= arc.35.63 <= 1
 0 <= arc.36.20 <= 1
 0 <= arc.36.39 <= 1
 0 <= arc.36.21 <= 1
 0 <= arc.36.14 <= 1
 0 <= arc.36.42 <= 1
 0 <= arc.36.19 <= 1
 0 <= arc.36.41 <= 1
 0 <= arc.36.37 <= 1
 0 <= arc.36.63 <= 1
 0 <= arc.36.32 <= 1
 0 <= arc.37.36 <= 1
 0 <= arc.37.18 <= 1
 0 <= arc.37.48 <= 1
 0 <= arc.38.35 <= 1
 0 <= arc.38.36 <= 1
 0 <= arc.38.18 <= 1
 0 <= arc.38.58 <= 1
 0 <= arc.39.36 <= 1
 0 <= arc.39.35 <= 1
 0 <= arc.39.18 <= 1
 0 <= arc.41.32 <= 1
 0 <= arc.41.63 <= 1
 0 <= arc.41.57 <= 1
 0 <= arc.41.19 <= 1
 0 <= arc.41.14 <= 1
 0 <= arc.41.37 <= 1
 0 <= arc.42.5 <= 1
 0 <= arc.42.56 <= 1
 0 <= arc.42.44 <= 1
 0 <= arc.42.30 <= 1
 0 <= arc.42.24 <= 1
 0 <= arc.42.35 <= 1
 0 <= arc.42.50 <= 1
 0 <= arc.42.58 <= 1
 0 <= arc.42.48 <= 1
 0 <= arc.42.54 <= 1
 0 <= arc.42.36 <= 1
 0 <= arc.42.18 <= 1
 0 <= arc.42.8 <= 1
 0 <= arc.42.33 <= 1
 0 <= arc.42.25 <= 1
 0 <= arc.42.37 <= 1
 0 <= arc.42.4 <= 1
 0 <= arc.42.0 <= 1
 0 <= arc.42.31 <= 1
 0 <= arc.42.39 <= 1
 0 <= arc.42.23 <= 1
 0 <= arc.42.53 <= 1
 0 <= arc.42.63 <= 1
 0 <= arc.42.2 <= 1
 0 <= arc.42.6 <= 1
 0 <= arc.42.11 <= 1
 0 <= arc.42.29 <= 1
 0 <= arc.42.61 <= 1
 0 <= arc.42.19 <= 1
 0 <= arc.42.47 <= 1
 0 <= arc.42.59 <= 1
 0 <= arc.42.43 <= 1
 0 <= arc.42.60 <= 1
 0 <= arc.42.34 <= 1
 0 <= arc.43.37 <= 1
 0 <= arc.43.63 <= 1
 0 <= arc.43.21 <= 1
 0 <= arc.43.57 <= 1
 0 <= arc.43.19 <= 1
 0 <= arc.43.39 <= 1
 0 <= arc.43.14 <= 1
 0 <= arc.44.63 <= 1
 0 <= arc.44.20 <= 1
 0 <= arc.44.14 <= 1
 0 <= arc.44.32 <= 1
 0 <= arc.44.37 <= 1
 0 <= arc.44.19 <= 1
 0 <= arc.44.21 <= 1
 0 <= arc.44.39 <= 1
 0 <= arc.45.19 <= 1
 0 <= arc.45.44 <= 1
 0 <= arc.45.29 <= 1
 0 <= arc.45.0 <= 1
 0 <= arc.45.35 <= 1
 0 <= arc.45.48 <= 1
 0 <= arc.45.8 <= 1
 0 <= arc.45.62 <= 1
 0 <= arc.45.2 <= 1
 0 <= arc.45.4 <= 1
 0 <= arc.45.37 <= 1
 0 <= arc.45.57 <= 1
 0 <= arc.45.11 <= 1
 0 <= arc.45.42 <= 1
 0 <= arc.45.60 <= 1
 0 <= arc.45.59 <= 1
 0 <= arc.45.54 <= 1
 0 <= arc.45.14 <= 1
 0 <= arc.45.36 <= 1
 0 <= arc.45.18 <= 1
 0 <= arc.45.50 <= 1
 0 <= arc.45.34 <= 1
 0 <= arc.45.58 <= 1
 0 <= arc.45.43 <= 1
 0 <= arc.45.61 <= 1
 0 <= arc.45.5 <= 1
 0 <= arc.45.7 <= 1
 0 <= arc.45.56 <= 1
 0 <= arc.45.41 <= 1
 0 <= arc.45.52 <= 1
 0 <= arc.45.38 <= 1
 0 <= arc.45.20 <= 1
 0 <= arc.45.9 <= 1
 0 <= arc.45.24 <= 1
 0 <= arc.45.63 <= 1
 0 <= arc.45.21 <= 1
 0 <= arc.45.25 <= 1
 0 <= arc.45.53 <= 1
 0 <= arc.46.19 <= 1
 0 <= arc.46.6 <= 1
 0 <= arc.46.56 <= 1
 0 <= arc.46.18 <= 1
 0 <= arc.46.36 <= 1
 0 <= arc.46.30 <= 1
 0 <= arc.46.60 <= 1
 0 <= arc.46.5 <= 1
 0 <= arc.46.44 <= 1
 0 <= arc.46.48 <= 1
 0 <= arc.46.9 <= 1
 0 <= arc.46.0 <= 1
 0 <= arc.46.4 <= 1
 0 <= arc.46.62 <= 1
 0 <= arc.46.58 <= 1
 0 <= arc.46.54 <= 1
 0 <= arc.46.10 <= 1
 0 <= arc.46.14 <= 1
 0 <= arc.46.2 <= 1
 0 <= arc.46.23 <= 1
 0 <= arc.46.8 <= 1
 0 <= arc.46.37 <= 1
 0 <= arc.46.29 <= 1
 0 <= arc.46.25 <= 1
 0 <= arc.46.59 <= 1
 0 <= arc.46.55 <= 1
 0 <= arc.46.11 <= 1
 0 <= arc.46.7 <= 1
 0 <= arc.46.35 <= 1
 0 <= arc.46.63 <= 1
 0 <= arc.46.53 <= 1
 0 <= arc.46.61 <= 1
 0 <= arc.46.20 <= 1
 0 <= arc.46.38 <= 1
 0 <= arc.46.43 <= 1
 0 <= arc.47.18 <= 1
 0 <= arc.47.35 <= 1
 0 <= arc.47.58 <= 1
 0 <= arc.47.48 <= 1
 0 <= arc.47.36 <= 1
 0 <= arc.48.20 <= 1
 0 <= arc.48.63 <= 1
 0 <= arc.48.19 <= 1
 0 <= arc.48.37 <= 1
 0 <= arc.48.14 <= 1
 0 <= arc.48.10 <= 1
 0 <= arc.49.48 <= 1
 0 <= arc.49.56 <= 1
 0 <= arc.49.20 <= 1
 0 <= arc.49.44 <= 1
 0 <= arc.49.8 <= 1
 0 <= arc.49.6 <= 1
 0 <= arc.49.24 <= 1
 0 <= arc.49.36 <= 1
 0 <= arc.49.30 <= 1
 0 <= arc.49.60 <= 1
 0 <= arc.49.38 <= 1
 0 <= arc.49.35 <= 1
 0 <= arc.49.5 <= 1
 0 <= arc.49.53 <= 1
 0 <= arc.49.61 <= 1
 0 <= arc.49.50 <= 1
 0 <= arc.49.29 <= 1
 0 <= arc.49.63 <= 1
 0 <= arc.49.58 <= 1
 0 <= arc.49.7 <= 1
 0 <= arc.49.43 <= 1
 0 <= arc.49.25 <= 1
 0 <= arc.49.57 <= 1
 0 <= arc.49.59 <= 1
 0 <= arc.49.37 <= 1
 0 <= arc.49.0 <= 1
 0 <= arc.49.19 <= 1
 0 <= arc.49.9 <= 1
 0 <= arc.49.34 <= 1
 0 <= arc.49.11 <= 1
 0 <= arc.49.18 <= 1
 0 <= arc.49.4 <= 1
 0 <= arc.49.54 <= 1
 0 <= arc.49.2 <= 1
 0 <= arc.50.10 <= 1
 0 <= arc.50.63 <= 1
 0 <= arc.50.32 <= 1
 0 <= arc.50.37 <= 1
 0 <= arc.50.21 <= 1
 0 <= arc.50.19 <= 1
 0 <= arc.50.12 <= 1
 0 <= arc.51.40 <= 1
 0 <= arc.51.35 <= 1
 0 <= arc.51.48 <= 1
 0 <= arc.51.58 <= 1
 0 <= arc.51.38 <= 1
 0 <= arc.51.36 <= 1
 0 <= arc.52.4 <= 1
 0 <= arc.52.21 <= 1
 0 <= arc.52.0 <= 1
 0 <= arc.52.60 <= 1
 0 <= arc.52.43 <= 1
 0 <= arc.52.53 <= 1
 0 <= arc.52.55 <= 1
 0 <= arc.52.2 <= 1
 0 <= arc.52.63 <= 1
 0 <= arc.52.19 <= 1
 0 <= arc.52.61 <= 1
 0 <= arc.52.20 <= 1
 0 <= arc.52.36 <= 1
 0 <= arc.52.34 <= 1
 0 <= arc.52.58 <= 1
 0 <= arc.52.6 <= 1
 0 <= arc.52.50 <= 1
 0 <= arc.52.5 <= 1
 0 <= arc.52.14 <= 1
 0 <= arc.52.11 <= 1
 0 <= arc.52.44 <= 1
 0 <= arc.52.56 <= 1
 0 <= arc.52.54 <= 1
 0 <= arc.52.48 <= 1
 0 <= arc.52.46 <= 1
 0 <= arc.52.25 <= 1
 0 <= arc.52.35 <= 1
 0 <= arc.52.42 <= 1
 0 <= arc.52.18 <= 1
 0 <= arc.52.9 <= 1
 0 <= arc.52.29 <= 1
 0 <= arc.52.10 <= 1
 0 <= arc.52.33 <= 1
 0 <= arc.52.59 <= 1
 0 <= arc.52.23 <= 1
 0 <= arc.52.37 <= 1
 0 <= arc.53.35 <= 1
 0 <= arc.53.58 <= 1
 0 <= arc.53.36 <= 1
 0 <= arc.53.18 <= 1
 0 <= arc.53.50 <= 1
 0 <= arc.53.38 <= 1
 0 <= arc.53.48 <= 1
 0 <= arc.54.32 <= 1
 0 <= arc.54.42 <= 1
 0 <= arc.54.21 <= 1
 0 <= arc.54.19 <= 1
 0 <= arc.54.37 <= 1
 0 <= arc.54.10 <= 1
 0 <= arc.54.57 <= 1
 0 <= arc.54.63 <= 1
 0 <= arc.54.20 <= 1
 0 <= arc.55.1 <= 1
 0 <= arc.55.20 <= 1
 0 <= arc.55.32 <= 1
 0 <= arc.55.19 <= 1
 0 <= arc.55.37 <= 1
 0 <= arc.55.21 <= 1
 0 <= arc.55.63 <= 1
 0 <= arc.55.57 <= 1
 0 <= arc.56.19 <= 1
 0 <= arc.56.63 <= 1
 0 <= arc.56.21 <= 1
 0 <= arc.56.57 <= 1
 0 <= arc.56.10 <= 1
 0 <= arc.56.32 <= 1
 0 <= arc.56.37 <= 1
 0 <= arc.57.36 <= 1
 0 <= arc.57.48 <= 1
 0 <= arc.57.35 <= 1
 0 <= arc.58.20 <= 1
 0 <= arc.58.14 <= 1
 0 <= arc.58.10 <= 1
 0 <= arc.58.32 <= 1
 0 <= arc.58.19 <= 1
 0 <= arc.58.63 <= 1
 0 <= arc.58.21 <= 1
 0 <= arc.58.37 <= 1
 0 <= arc.58.42 <= 1
 0 <= arc.59.20 <= 1
 0 <= arc.59.16 <= 1
 0 <= arc.59.19 <= 1
 0 <= arc.59.57 <= 1
 0 <= arc.59.10 <= 1
 0 <= arc.59.37 <= 1
 0 <= arc.59.63 <= 1
 0 <= arc.59.42 <= 1
 0 <= arc.60.18 <= 1
 0 <= arc.60.35 <= 1
 0 <= arc.60.58 <= 1
 0 <= arc.60.36 <= 1
 0 <= arc.60.48 <= 1
 0 <= arc.60.38 <= 1
 0 <= arc.62.48 <= 1
 0 <= arc.62.58 <= 1
 0 <= arc.62.36 <= 1
 0 <= arc.62.35 <= 1
 0 <= arc.62.18 <= 1
 0 <= arc.63.38 <= 1
 0 <= arc.63.48 <= 1
 0 <= arc.63.35 <= 1
 0 <= arc.63.18 <= 1
 0 <= arc.63.36 <= 1
 0 <= arc.64.2 <= 1
 0 <= arc.64.43 <= 1
 0 <= arc.64.9 <= 1
 0 <= arc.64.61 <= 1
 0 <= arc.64.58 <= 1
 0 <= arc.64.53 <= 1
 0 <= arc.64.5 <= 1
 0 <= arc.64.11 <= 1
 0 <= arc.64.25 <= 1
 0 <= arc.64.48 <= 1
 0 <= arc.64.60 <= 1
 0 <= arc.64.62 <= 1
 0 <= arc.64.49 <= 1
 0 <= arc.64.34 <= 1
 0 <= arc.64.44 <= 1
 0 <= arc.64.0 <= 1
 0 <= arc.64.8 <= 1
 0 <= arc.64.36 <= 1
 0 <= arc.64.18 <= 1
 0 <= arc.64.54 <= 1
 0 <= arc.64.6 <= 1
 0 <= arc.64.35 <= 1
 0 <= arc.64.23 <= 1
 0 <= arc.64.4 <= 1
 0 <= arc.64.30 <= 1
 0 <= arc.64.29 <= 1
 0 <= arc.64.7 <= 1
 0 <= arc.64.55 <= 1
 0 <= arc.65.9 <= 1
 0 <= arc.65.61 <= 1
 0 <= arc.65.43 <= 1
 0 <= arc.65.7 <= 1
 0 <= arc.65.5 <= 1
 0 <= arc.65.24 <= 1
 0 <= arc.65.59 <= 1
 0 <= arc.65.29 <= 1
 0 <= arc.65.30 <= 1
 0 <= arc.65.56 <= 1
 0 <= arc.65.2 <= 1
 0 <= arc.65.62 <= 1
 0 <= arc.65.54 <= 1
 0 <= arc.65.53 <= 1
 0 <= arc.65.25 <= 1
 0 <= arc.65.28 <= 1
 0 <= arc.65.60 <= 1
 0 <= arc.65.44 <= 1
 0 <= arc.65.4 <= 1
 0 <= arc.65.34 <= 1
 0 <= arc.65.8 <= 1
 0 <= arc.65.0 <= 1
 0 <= arc.65.11 <= 1
 0 <= arc.65.47 <= 1
 0 <= arc.66.8 <= 1
 0 <= arc.66.24 <= 1
 0 <= arc.66.44 <= 1
 0 <= arc.66.7 <= 1
 0 <= arc.66.47 <= 1
 0 <= arc.66.4 <= 1
 0 <= arc.66.34 <= 1
 0 <= arc.66.60 <= 1
 0 <= arc.66.62 <= 1
 0 <= arc.66.56 <= 1
 0 <= arc.66.2 <= 1
 0 <= arc.66.0 <= 1
 0 <= arc.66.11 <= 1
 0 <= arc.66.23 <= 1
 0 <= arc.66.53 <= 1
 0 <= arc.66.5 <= 1
 0 <= arc.66.29 <= 1
 0 <= arc.66.61 <= 1
 0 <= arc.66.9 <= 1
 0 <= arc.66.6 <= 1
 0 <= arc.66.43 <= 1
Binaries
 r.0  r.1  r.2  r.3  r.4  r.5  r.6  r.7  r.8  r.9  vertex.1  arc.1.37 
 vertex.37  arc.37.35  vertex.35  arc.35.1  vertex.6  arc.6.20  vertex.20 
 arc.20.45  vertex.45  arc.45.6  vertex.10  arc.10.32  vertex.32  arc.32.14 
 vertex.14  arc.14.10  vertex.12  arc.12.51  vertex.51  arc.51.18  vertex.18 
 arc.18.12  vertex.16  arc.16.21  vertex.21  arc.21.36  vertex.36  arc.36.16 
 vertex.30  arc.30.52  vertex.52  arc.52.30  vertex.34  arc.34.46  vertex.46 
 arc.46.34  vertex.38  arc.38.48  vertex.48  arc.48.57  vertex.57  arc.57.38 
 vertex.39  arc.39.58  vertex.58  arc.58.39  vertex.42  arc.42.62  vertex.62 
 arc.62.50  vertex.50  arc.50.42  vertex.0  vertex.2  vertex.3  vertex.4 
 vertex.5  vertex.7  vertex.8  vertex.9  vertex.11  vertex.13  vertex.15 
 vertex.17  vertex.19  vertex.22  vertex.23  vertex.24  vertex.25  vertex.26 
 vertex.27  vertex.28  vertex.29  vertex.31  vertex.33  vertex.40  vertex.41 
 vertex.43  vertex.44  vertex.47  vertex.49  vertex.53  vertex.54  vertex.55 
 vertex.56  vertex.59  vertex.60  vertex.61  vertex.63  arc.0.36  arc.0.18 
 arc.0.50  arc.0.48  arc.0.58  arc.0.35  arc.1.63  arc.1.19  arc.1.57 
 arc.1.39  arc.3.5  arc.3.24  arc.3.48  arc.3.62  arc.3.11  arc.3.25  arc.3.35 
 arc.3.0  arc.3.8  arc.3.18  arc.3.61  arc.3.56  arc.3.37  arc.3.57  arc.3.4 
 arc.3.53  arc.3.43  arc.3.49  arc.3.55  arc.3.6  arc.3.63  arc.3.23  arc.3.51 
 arc.3.59  arc.3.47  arc.3.2  arc.3.32  arc.3.50  arc.3.29  arc.3.36  arc.3.19 
 arc.3.30  arc.3.44  arc.3.7  arc.3.9  arc.3.38  arc.4.35  arc.4.18  arc.4.58 
 arc.4.48  arc.4.38  arc.4.50  arc.4.36  arc.5.35  arc.5.18  arc.5.48 
 arc.5.58  arc.5.38  arc.5.36  arc.6.63  arc.6.19  arc.6.39  arc.6.21 
 arc.6.14  arc.6.32  arc.6.37  arc.7.20  arc.7.19  arc.7.42  arc.7.37 
 arc.7.21  arc.7.10  arc.9.37  arc.9.19  arc.9.63  arc.9.39  arc.9.14 
 arc.9.32  arc.10.14  arc.10.37  arc.10.63  arc.10.57  arc.10.19  arc.11.1 
 arc.11.37  arc.11.19  arc.11.16  arc.11.32  arc.11.14  arc.12.10  arc.12.21 
 arc.12.19  arc.12.37  arc.12.32  arc.12.63  arc.12.14  arc.13.26  arc.13.44 
 arc.13.6  arc.13.50  arc.13.48  arc.13.8  arc.13.32  arc.13.18  arc.13.14 
 arc.13.36  arc.13.0  arc.13.54  arc.13.4  arc.13.62  arc.13.60  arc.13.34 
 arc.13.56  arc.13.2  arc.13.11  arc.13.39  arc.13.38  arc.13.29  arc.13.35 
 arc.13.21  arc.13.33  arc.13.37  arc.13.43  arc.13.25  arc.13.9  arc.13.41 
 arc.13.61  arc.13.1  arc.13.59  arc.13.5  arc.13.53  arc.13.19  arc.13.63 
 arc.13.7  arc.14.19  arc.14.51  arc.14.63  arc.14.57  arc.14.37  arc.15.36 
 arc.15.4  arc.15.32  arc.15.0  arc.15.58  arc.15.51  arc.15.49  arc.15.2 
 arc.15.55  arc.15.10  arc.15.14  arc.15.57  arc.15.19  arc.15.43  arc.15.11 
 arc.15.21  arc.15.37  arc.15.61  arc.15.29  arc.15.35  arc.15.24  arc.15.60 
 arc.15.7  arc.15.48  arc.15.62  arc.15.47  arc.15.6  arc.15.5  arc.15.9 
 arc.15.44  arc.15.56  arc.15.34  arc.15.63  arc.15.25  arc.15.54  arc.16.32 
 arc.16.10  arc.16.20  arc.16.19  arc.16.37  arc.17.48  arc.17.40  arc.17.36 
 arc.17.18  arc.17.35  arc.18.19  arc.18.37  arc.18.63  arc.18.57  arc.20.61 
 arc.20.28  arc.20.0  arc.20.24  arc.20.34  arc.20.4  arc.20.60  arc.20.56 
 arc.20.18  arc.20.2  arc.20.6  arc.20.42  arc.20.8  arc.20.44  arc.20.37 
 arc.20.9  arc.20.50  arc.20.35  arc.20.19  arc.20.29  arc.20.14  arc.20.5 
 arc.20.36  arc.20.53  arc.20.11  arc.20.47  arc.20.43  arc.20.63  arc.20.25 
 arc.21.18  arc.21.35  arc.21.48  arc.21.38  arc.22.37  arc.22.55  arc.22.19 
 arc.22.61  arc.22.50  arc.22.59  arc.22.43  arc.22.11  arc.22.6  arc.22.49 
 arc.22.63  arc.22.4  arc.22.8  arc.22.0  arc.22.58  arc.22.18  arc.22.62 
 arc.22.48  arc.22.44  arc.22.54  arc.22.26  arc.22.14  arc.22.34  arc.22.56 
 arc.22.36  arc.22.60  arc.22.30  arc.22.38  arc.22.2  arc.22.10  arc.22.5 
 arc.22.9  arc.22.21  arc.22.25  arc.22.53  arc.22.35  arc.22.17  arc.23.35 
 arc.23.48  arc.23.36  arc.23.18  arc.23.58  arc.24.38  arc.24.48  arc.24.36 
 arc.24.40  arc.24.58  arc.24.18  arc.24.35  arc.25.37  arc.25.14  arc.25.19 
 arc.25.63  arc.26.57  arc.26.19  arc.26.37  arc.26.63  arc.26.10  arc.26.42 
 arc.26.46  arc.26.32  arc.27.9  arc.27.37  arc.27.59  arc.27.61  arc.27.43 
 arc.27.53  arc.27.57  arc.27.33  arc.27.36  arc.27.22  arc.27.55  arc.27.32 
 arc.27.19  arc.27.34  arc.27.20  arc.27.63  arc.27.0  arc.27.11  arc.27.38 
 arc.27.24  arc.27.8  arc.27.6  arc.27.44  arc.27.10  arc.27.2  arc.27.7 
 arc.27.62  arc.27.12  arc.27.25  arc.27.56  arc.27.18  arc.27.29  arc.27.48 
 arc.27.5  arc.27.35  arc.27.60  arc.27.14  arc.27.42  arc.29.35  arc.29.48 
 arc.29.36  arc.29.18  arc.29.50  arc.30.54  arc.30.48  arc.30.50  arc.30.2 
 arc.30.44  arc.30.56  arc.30.18  arc.30.60  arc.30.26  arc.30.34  arc.30.61 
 arc.30.9  arc.30.62  arc.30.14  arc.30.43  arc.30.63  arc.30.5  arc.30.53 
 arc.30.13  arc.30.35  arc.30.21  arc.30.47  arc.30.11  arc.30.19  arc.30.23 
 arc.30.36  arc.30.55  arc.30.8  arc.30.24  arc.30.58  arc.30.4  arc.30.37 
 arc.30.0  arc.30.25  arc.31.46  arc.31.10  arc.31.32  arc.31.20  arc.31.1 
 arc.31.63  arc.31.37  arc.31.19  arc.32.19  arc.32.20  arc.32.57  arc.33.63 
 arc.33.37  arc.33.19  arc.33.42  arc.34.20  arc.34.21  arc.34.19  arc.34.37 
 arc.34.63  arc.34.32  arc.35.14  arc.35.32  arc.35.37  arc.35.42  arc.35.57 
 arc.35.19  arc.35.63  arc.36.20  arc.36.39  arc.36.21  arc.36.14  arc.36.42 
 arc.36.19  arc.36.41  arc.36.37  arc.36.63  arc.36.32  arc.37.36  arc.37.18 
 arc.37.48  arc.38.35  arc.38.36  arc.38.18  arc.38.58  arc.39.36  arc.39.35 
 arc.39.18  arc.41.32  arc.41.63  arc.41.57  arc.41.19  arc.41.14  arc.41.37 
 arc.42.5  arc.42.56  arc.42.44  arc.42.30  arc.42.24  arc.42.35  arc.42.50 
 arc.42.58  arc.42.48  arc.42.54  arc.42.36  arc.42.18  arc.42.8  arc.42.33 
 arc.42.25  arc.42.37  arc.42.4  arc.42.0  arc.42.31  arc.42.39  arc.42.23 
 arc.42.53  arc.42.63  arc.42.2  arc.42.6  arc.42.11  arc.42.29  arc.42.61 
 arc.42.19  arc.42.47  arc.42.59  arc.42.43  arc.42.60  arc.42.34  arc.43.37 
 arc.43.63  arc.43.21  arc.43.57  arc.43.19  arc.43.39  arc.43.14  arc.44.63 
 arc.44.20  arc.44.14  arc.44.32  arc.44.37  arc.44.19  arc.44.21  arc.44.39 
 arc.45.19  arc.45.44  arc.45.29  arc.45.0  arc.45.35  arc.45.48  arc.45.8 
 arc.45.62  arc.45.2  arc.45.4  arc.45.37  arc.45.57  arc.45.11  arc.45.42 
 arc.45.60  arc.45.59  arc.45.54  arc.45.14  arc.45.36  arc.45.18  arc.45.50 
 arc.45.34  arc.45.58  arc.45.43  arc.45.61  arc.45.5  arc.45.7  arc.45.56 
 arc.45.41  arc.45.52  arc.45.38  arc.45.20  arc.45.9  arc.45.24  arc.45.63 
 arc.45.21  arc.45.25  arc.45.53  arc.46.19  arc.46.6  arc.46.56  arc.46.18 
 arc.46.36  arc.46.30  arc.46.60  arc.46.5  arc.46.44  arc.46.48  arc.46.9 
 arc.46.0  arc.46.4  arc.46.62  arc.46.58  arc.46.54  arc.46.10  arc.46.14 
 arc.46.2  arc.46.23  arc.46.8  arc.46.37  arc.46.29  arc.46.25  arc.46.59 
 arc.46.55  arc.46.11  arc.46.7  arc.46.35  arc.46.63  arc.46.53  arc.46.61 
 arc.46.20  arc.46.38  arc.46.43  arc.47.18  arc.47.35  arc.47.58  arc.47.48 
 arc.47.36  arc.48.20  arc.48.63  arc.48.19  arc.48.37  arc.48.14  arc.48.10 
 arc.49.48  arc.49.56  arc.49.20  arc.49.44  arc.49.8  arc.49.6  arc.49.24 
 arc.49.36  arc.49.30  arc.49.60  arc.49.38  arc.49.35  arc.49.5  arc.49.53 
 arc.49.61  arc.49.50  arc.49.29  arc.49.63  arc.49.58  arc.49.7  arc.49.43 
 arc.49.25  arc.49.57  arc.49.59  arc.49.37  arc.49.0  arc.49.19  arc.49.9 
 arc.49.34  arc.49.11  arc.49.18  arc.49.4  arc.49.54  arc.49.2  arc.50.10 
 arc.50.63  arc.50.32  arc.50.37  arc.50.21  arc.50.19  arc.50.12  arc.51.40 
 arc.51.35  arc.51.48  arc.51.58  arc.51.38  arc.51.36  arc.52.4  arc.52.21 
 arc.52.0  arc.52.60  arc.52.43  arc.52.53  arc.52.55  arc.52.2  arc.52.63 
 arc.52.19  arc.52.61  arc.52.20  arc.52.36  arc.52.34  arc.52.58  arc.52.6 
 arc.52.50  arc.52.5  arc.52.14  arc.52.11  arc.52.44  arc.52.56  arc.52.54 
 arc.52.48  arc.52.46  arc.52.25  arc.52.35  arc.52.42  arc.52.18  arc.52.9 
 arc.52.29  arc.52.10  arc.52.33  arc.52.59  arc.52.23  arc.52.37  arc.53.35 
 arc.53.58  arc.53.36  arc.53.18  arc.53.50  arc.53.38  arc.53.48  arc.54.32 
 arc.54.42  arc.54.21  arc.54.19  arc.54.37  arc.54.10  arc.54.57  arc.54.63 
 arc.54.20  arc.55.1  arc.55.20  arc.55.32  arc.55.19  arc.55.37  arc.55.21 
 arc.55.63  arc.55.57  arc.56.19  arc.56.63  arc.56.21  arc.56.57  arc.56.10 
 arc.56.32  arc.56.37  arc.57.36  arc.57.48  arc.57.35  arc.58.20  arc.58.14 
 arc.58.10  arc.58.32  arc.58.19  arc.58.63  arc.58.21  arc.58.37  arc.58.42 
 arc.59.20  arc.59.16  arc.59.19  arc.59.57  arc.59.10  arc.59.37  arc.59.63 
 arc.59.42  arc.60.18  arc.60.35  arc.60.58  arc.60.36  arc.60.48  arc.60.38 
 arc.62.48  arc.62.58  arc.62.36  arc.62.35  arc.62.18  arc.63.38  arc.63.48 
 arc.63.35  arc.63.18  arc.63.36  arc.64.2  arc.64.43  arc.64.9  arc.64.61 
 arc.64.58  arc.64.53  arc.64.5  arc.64.11  arc.64.25  arc.64.48  arc.64.60 
 arc.64.62  arc.64.49  arc.64.34  arc.64.44  arc.64.0  arc.64.8  arc.64.36 
 arc.64.18  arc.64.54  arc.64.6  arc.64.35  arc.64.23  arc.64.4  arc.64.30 
 arc.64.29  arc.64.7  arc.64.55  arc.65.9  arc.65.61  arc.65.43  arc.65.7 
 arc.65.5  arc.65.24  arc.65.59  arc.65.29  arc.65.30  arc.65.56  arc.65.2 
 arc.65.62  arc.65.54  arc.65.53  arc.65.25  arc.65.28  arc.65.60  arc.65.44 
 arc.65.4  arc.65.34  arc.65.8  arc.65.0  arc.65.11  arc.65.47  arc.66.8 
 arc.66.24  arc.66.44  arc.66.7  arc.66.47  arc.66.4  arc.66.34  arc.66.60 
 arc.66.62  arc.66.56  arc.66.2  arc.66.0  arc.66.11  arc.66.23  arc.66.53 
 arc.66.5  arc.66.29  arc.66.61  arc.66.9  arc.66.6  arc.66.43 
End
