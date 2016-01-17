# Symbolic discretization of F_c, Q_c
from sympy import MatrixSymbol, symbols, ZeroMatrix, Identity, \
    BlockMatrix, block_collapse, BlockDiagMatrix, simplify

class Discretizer:
    def discretizeF(self):
        Z33 = ZeroMatrix(3, 3)
        Z31 = ZeroMatrix(3, 1)
        Z13 = ZeroMatrix(1, 3)
        Z11 = ZeroMatrix(1, 1)

        I33 = Identity(3)
        I11 = Identity(1)

        F11 = Z33
        F12 = MatrixSymbol('F12', 3, 3)
        F13 = MatrixSymbol('F13', 3, 3)
        F14 = Z33
        F15 = Z31
        F16 = Z31
        F17 = Z33

        F21 = Z33
        F22 = MatrixSymbol('F22', 3, 3)
        F23 = MatrixSymbol('F23', 3, 3)
        F24 = MatrixSymbol('F24', 3, 3)
        F25 = MatrixSymbol('F25', 3, 1)
        F26 = Z31
        F27 = Z33

        F31 = Z33
        F32 = Z33
        F33 = MatrixSymbol('F33', 3, 3)
        F34 = I33
        F35 = Z31
        F36 = Z31
        F37 = Z33

        F41 = Z33
        F42 = Z33
        F43 = Z33
        F44 = MatrixSymbol('F44', 3, 3)
        F45 = MatrixSymbol('F45', 3, 1)
        F46 = MatrixSymbol('F46', 3, 1)
        F47 = MatrixSymbol('F47', 3, 3)

        F51 = Z13
        F52 = Z13
        F53 = Z13
        F54 = Z13
        F55 = Z11
        F56 = Z11
        F57 = Z13

        F61 = Z13
        F62 = Z13
        F63 = Z13
        F64 = Z13
        F65 = Z11
        F66 = Z11
        F67 = Z13

        F71 = Z33
        F72 = Z33
        F73 = Z33
        F74 = Z33
        F75 = Z31
        F76 = Z31
        F77 = Z33

        F_c = BlockMatrix([[F11, F12, F13, F14, F15, F16, F17],
                           [F21, F22, F23, F24, F25, F26, F27],
                           [F31, F32, F33, F34, F35, F36, F37],
                           [F41, F42, F43, F44, F45, F46, F47],
                           [F51, F52, F53, F54, F55, F56, F57],
                           [F61, F62, F63, F64, F65, F66, F67],
                           [F71, F72, F73, F74, F75, F76, F77]])
        print(F_c)

        I_d = Identity(17)

        Delta_t = symbols('Delta_t')
        F_k = block_collapse(I_d + F_c * Delta_t)

        print(F_k)

        # k=6
        G11 = Z33
        G12 = Z33
        G13 = Z33
        G14 = Z33
        G15 = Z33
        G16 = Z33
        G17 = Z33
        G18 = Z33

        GZ3 = BlockMatrix([[Z33, Z33, Z33, Z33, Z33, Z33, Z33, Z33]])

        G21 = MatrixSymbol('G21', 3, 3)
        G22 = MatrixSymbol('G22', 3, 3)
        G23 = MatrixSymbol('G23', 3, 3)
        G24 = MatrixSymbol('G24', 3, 3)
        G25 = MatrixSymbol('G25', 3, 3)
        G26 = MatrixSymbol('G26', 3, 3)
        G27 = I33
        G28 = Z33

        G41 = MatrixSymbol('G41', 3, 3)
        G42 = MatrixSymbol('G42', 3, 3)
        G43 = MatrixSymbol('G43', 3, 3)
        G44 = MatrixSymbol('G44', 3, 3)
        G45 = MatrixSymbol('G45', 3, 3)
        G46 = MatrixSymbol('G46', 3, 3)
        G47 = Z33
        G48 = I33

        G51 = Z13
        G52 = Z13
        G53 = Z13
        G54 = Z13
        G55 = Z13
        G56 = Z13
        G57 = Z13
        G58 = Z13

        G_c = BlockMatrix([[G11, G12, G13, G14, G15, G16, G17, G18],
                           [G21, G22, G23, G24, G25, G26, G27, G28],
                           [G11, G12, G13, G14, G15, G16, G17, G18],
                           [G41, G42, G43, G44, G45, G46, G47, G48],
                           [G51, G52, G53, G54, G55, G56, G57, G58],
                           [G51, G52, G53, G54, G55, G56, G57, G58],
                           [G11, G12, G13, G14, G15, G16, G17, G18]])

        print(G_c)

        Q11 = MatrixSymbol('sigma2_T', 3, 3)
        Q22 = Q11
        Q33 = Q11
        Q44 = Q11
        Q55 = Q11
        Q66 = Q11
        Q77 = MatrixSymbol('sigma2_A', 3, 3)
        Q88 = MatrixSymbol('sigma2_M', 3, 3)
        Q_c = BlockDiagMatrix(Q11, Q22, Q33, Q44, Q55, Q66, Q77, Q88)

        Q_k = block_collapse(G_c * Q_c * G_c.T)
        Delta_t_2 = symbols('Delta_t_2')
        Q_k = block_collapse(Q_k * Delta_t_2)
        print(simplify(Q_k))

        # Discussion
        P_p = MatrixSymbol('P_p', 3, 3)
        P_v = MatrixSymbol('P_v', 3, 3)
        P_q = MatrixSymbol('P_q', 3, 3)
        P_a = MatrixSymbol('P_a', 3, 3)
        P_t = MatrixSymbol('P_t', 1, 1)
        P_m = MatrixSymbol('P_m', 1, 1)
        P_j = MatrixSymbol('P_j', 3, 3)
        P0 = BlockDiagMatrix(P_p, P_v, P_q, P_a, P_t, P_m, P_j)

        F_kd = block_collapse(F_k - F_c*Delta_t + F_c)
        Q_kd = block_collapse(G_c * Q_c * G_c.T)

        P1 = block_collapse(F_kd * P0 * F_kd.T + Q_kd)
        #print(P1)
        P2 = block_collapse(F_kd * P1 * F_kd.T + Q_kd)
        print(P2)
        H = BlockMatrix([[I33, Z33, Z33, Z33, Z31, Z31, Z33],
                         [Z33, Z33, I33, Z33, Z31, Z31, Z33]])
        S2 = block_collapse(H * P2 * H.T + Identity(6))
        #print(S2)
        K2 = block_collapse(P2 * H.T * S2.I)
        #print(K2)
        P2_upd = block_collapse((Identity(17) - K2*H)*P2)
        print(P2_upd)


Discretizer().discretizeF()