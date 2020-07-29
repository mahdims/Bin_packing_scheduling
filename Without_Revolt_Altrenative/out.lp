\ Model packing
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  60 gamma[0,0] + 60 gamma[1,1] + 60 gamma[2,2] + 60 gamma[3,3]
   + 60 gamma[4,4] + 60 gamma[5,5] + 60 gamma[6,6] + 60 gamma[7,7]
   + 60 gamma[8,8] + 60 gamma[9,9] + 0.0662904 Q[0] + 0.0662904 Q[1]
   + 0.0662904 Q[2] + 0.0662904 Q[3] + 0.0662904 Q[4] + 0.0662904 Q[5]
   + 0.0662904 Q[6] + 0.0662904 Q[7] + 0.0662904 Q[8] + 0.0662904 Q[9]
Subject To
 R0: alpha[0,0] = 1
 R1: alpha[1,1] = 1
 R2: alpha[0,2] + alpha[2,2] = 1
 R3: alpha[1,3] + alpha[3,3] = 1
 R4: alpha[4,4] = 1
 R5: alpha[5,5] = 1
 R6: alpha[0,6] + alpha[2,6] + alpha[6,6] = 1
 R7: alpha[1,7] + alpha[3,7] + alpha[7,7] = 1
 R8: alpha[0,8] + alpha[2,8] + alpha[6,8] + alpha[8,8] = 1
 R9: alpha[0,9] + alpha[2,9] + alpha[6,9] + alpha[8,9] + alpha[9,9] = 1
 R10: - 8 alpha[0,0] + alpha[0,2] + alpha[0,6] + alpha[0,8] + alpha[0,9]
   <= 0
 R11: - 7 alpha[1,1] + alpha[1,3] + alpha[1,7] <= 0
 R12: - 6 alpha[2,2] + alpha[2,6] + alpha[2,8] + alpha[2,9] <= 0
 R13: - 5 alpha[3,3] + alpha[3,7] <= 0
 R14: - 4 alpha[4,4] <= 0
 R15: - 3 alpha[5,5] <= 0
 R16: - 2 alpha[6,6] + alpha[6,8] + alpha[6,9] <= 0
 R17: - alpha[7,7] <= 0
 R18: alpha[8,9] <= 0
 R19: - alpha[0,0] + beta[0,0] + beta[1,0] + beta[2,0] + beta[3,0]
   + beta[4,0] + beta[5,0] + beta[6,0] + beta[7,0] + beta[8,0] + beta[9,0]
   = 0
 R20: - alpha[1,1] + beta[0,1] + beta[1,1] + beta[2,1] + beta[3,1]
   + beta[4,1] + beta[5,1] + beta[6,1] + beta[7,1] + beta[8,1] + beta[9,1]
   = 0
 R21: - alpha[2,2] + beta[0,2] + beta[1,2] + beta[2,2] + beta[3,2]
   + beta[4,2] + beta[5,2] + beta[6,2] + beta[7,2] + beta[8,2] + beta[9,2]
   = 0
 R22: - alpha[3,3] + beta[0,3] + beta[1,3] + beta[2,3] + beta[3,3]
   + beta[4,3] + beta[5,3] + beta[6,3] + beta[7,3] + beta[8,3] + beta[9,3]
   = 0
 R23: - alpha[4,4] + beta[0,4] + beta[1,4] + beta[2,4] + beta[3,4]
   + beta[4,4] + beta[5,4] + beta[6,4] + beta[7,4] + beta[8,4] + beta[9,4]
   = 0
 R24: - alpha[5,5] + beta[0,5] + beta[1,5] + beta[2,5] + beta[3,5]
   + beta[4,5] + beta[5,5] + beta[6,5] + beta[7,5] + beta[8,5] + beta[9,5]
   = 0
 R25: - alpha[6,6] + beta[0,6] + beta[1,6] + beta[2,6] + beta[3,6]
   + beta[4,6] + beta[5,6] + beta[6,6] + beta[7,6] + beta[8,6] + beta[9,6]
   = 0
 R26: - alpha[7,7] + beta[0,7] + beta[1,7] + beta[2,7] + beta[3,7]
   + beta[4,7] + beta[5,7] + beta[6,7] + beta[7,7] + beta[8,7] + beta[9,7]
   = 0
 R27: - alpha[8,8] + beta[0,8] + beta[1,8] + beta[2,8] + beta[3,8]
   + beta[4,8] + beta[5,8] + beta[6,8] + beta[7,8] + beta[8,8] + beta[9,8]
   = 0
 R28: - alpha[9,9] + beta[0,9] + beta[1,9] + beta[2,9] + beta[3,9]
   + beta[4,9] + beta[5,9] + beta[6,9] + beta[7,9] + beta[8,9] + beta[9,9]
   = 0
 R29: 21 alpha[0,0] + 21 alpha[0,2] + 21 alpha[0,6] + 21 alpha[0,8]
   + 21 alpha[0,9] - 21 alpha[1,1] - 21 alpha[1,3] - 21 alpha[1,7]
   + 131 beta[1,0] <= 130.997
 R30: 21 alpha[0,0] + 21 alpha[0,2] + 21 alpha[0,6] + 21 alpha[0,8]
   + 21 alpha[0,9] - 21 alpha[2,2] - 21 alpha[2,6] - 21 alpha[2,8]
   - 21 alpha[2,9] + 131 beta[2,0] <= 130.997
 R31: 21 alpha[1,1] + 21 alpha[1,3] + 21 alpha[1,7] - 21 alpha[2,2]
   - 21 alpha[2,6] - 21 alpha[2,8] - 21 alpha[2,9] + 131 beta[2,1]
   <= 130.997
 R32: 21 alpha[0,0] + 21 alpha[0,2] + 21 alpha[0,6] + 21 alpha[0,8]
   + 21 alpha[0,9] - 21 alpha[3,3] - 21 alpha[3,7] + 131 beta[3,0]
   <= 130.997
 R33: 21 alpha[1,1] + 21 alpha[1,3] + 21 alpha[1,7] - 21 alpha[3,3]
   - 21 alpha[3,7] + 131 beta[3,1] <= 130.997
 R34: 21 alpha[2,2] + 21 alpha[2,6] + 21 alpha[2,8] + 21 alpha[2,9]
   - 21 alpha[3,3] - 21 alpha[3,7] + 131 beta[3,2] <= 130.997
 R35: 21 alpha[0,0] + 21 alpha[0,2] + 21 alpha[0,6] + 21 alpha[0,8]
   + 21 alpha[0,9] - 21 alpha[4,4] + 131 beta[4,0] <= 130.997
 R36: 21 alpha[1,1] + 21 alpha[1,3] + 21 alpha[1,7] - 21 alpha[4,4]
   + 131 beta[4,1] <= 130.997
 R37: 21 alpha[2,2] + 21 alpha[2,6] + 21 alpha[2,8] + 21 alpha[2,9]
   - 21 alpha[4,4] + 131 beta[4,2] <= 130.997
 R38: 21 alpha[3,3] + 21 alpha[3,7] - 21 alpha[4,4] + 131 beta[4,3]
   <= 130.997
 R39: 21 alpha[0,0] + 21 alpha[0,2] + 21 alpha[0,6] + 21 alpha[0,8]
   + 21 alpha[0,9] - 10 alpha[5,5] + 131 beta[5,0] <= 130.997
 R40: 21 alpha[1,1] + 21 alpha[1,3] + 21 alpha[1,7] - 10 alpha[5,5]
   + 131 beta[5,1] <= 130.997
 R41: 21 alpha[2,2] + 21 alpha[2,6] + 21 alpha[2,8] + 21 alpha[2,9]
   - 10 alpha[5,5] + 131 beta[5,2] <= 130.997
 R42: 21 alpha[3,3] + 21 alpha[3,7] - 10 alpha[5,5] + 131 beta[5,3]
   <= 130.997
 R43: 21 alpha[4,4] - 10 alpha[5,5] + 131 beta[5,4] <= 130.997
 R44: 21 alpha[0,0] + 21 alpha[0,2] + 21 alpha[0,6] + 21 alpha[0,8]
   + 21 alpha[0,9] - 21 alpha[6,6] - 21 alpha[6,8] - 21 alpha[6,9]
   + 131 beta[6,0] <= 130.997
 R45: 21 alpha[1,1] + 21 alpha[1,3] + 21 alpha[1,7] - 21 alpha[6,6]
   - 21 alpha[6,8] - 21 alpha[6,9] + 131 beta[6,1] <= 130.997
 R46: 21 alpha[2,2] + 21 alpha[2,6] + 21 alpha[2,8] + 21 alpha[2,9]
   - 21 alpha[6,6] - 21 alpha[6,8] - 21 alpha[6,9] + 131 beta[6,2]
   <= 130.997
 R47: 21 alpha[3,3] + 21 alpha[3,7] - 21 alpha[6,6] - 21 alpha[6,8]
   - 21 alpha[6,9] + 131 beta[6,3] <= 130.997
 R48: 21 alpha[4,4] - 21 alpha[6,6] - 21 alpha[6,8] - 21 alpha[6,9]
   + 131 beta[6,4] <= 130.997
 R49: 10 alpha[5,5] - 21 alpha[6,6] - 21 alpha[6,8] - 21 alpha[6,9]
   + 131 beta[6,5] <= 130.997
 R50: 21 alpha[0,0] + 21 alpha[0,2] + 21 alpha[0,6] + 21 alpha[0,8]
   + 21 alpha[0,9] - 21 alpha[7,7] + 131 beta[7,0] <= 130.997
 R51: 21 alpha[1,1] + 21 alpha[1,3] + 21 alpha[1,7] - 21 alpha[7,7]
   + 131 beta[7,1] <= 130.997
 R52: 21 alpha[2,2] + 21 alpha[2,6] + 21 alpha[2,8] + 21 alpha[2,9]
   - 21 alpha[7,7] + 131 beta[7,2] <= 130.997
 R53: 21 alpha[3,3] + 21 alpha[3,7] - 21 alpha[7,7] + 131 beta[7,3]
   <= 130.997
 R54: 21 alpha[4,4] - 21 alpha[7,7] + 131 beta[7,4] <= 130.997
 R55: 10 alpha[5,5] - 21 alpha[7,7] + 131 beta[7,5] <= 130.997
 R56: 21 alpha[6,6] + 21 alpha[6,8] + 21 alpha[6,9] - 21 alpha[7,7]
   + 131 beta[7,6] <= 130.997
 R57: 21 alpha[0,0] + 21 alpha[0,2] + 21 alpha[0,6] + 21 alpha[0,8]
   + 21 alpha[0,9] - 21 alpha[8,8] - 21 alpha[8,9] + 131 beta[8,0]
   <= 130.997
 R58: 21 alpha[1,1] + 21 alpha[1,3] + 21 alpha[1,7] - 21 alpha[8,8]
   - 21 alpha[8,9] + 131 beta[8,1] <= 130.997
 R59: 21 alpha[2,2] + 21 alpha[2,6] + 21 alpha[2,8] + 21 alpha[2,9]
   - 21 alpha[8,8] - 21 alpha[8,9] + 131 beta[8,2] <= 130.997
 R60: 21 alpha[3,3] + 21 alpha[3,7] - 21 alpha[8,8] - 21 alpha[8,9]
   + 131 beta[8,3] <= 130.997
 R61: 21 alpha[4,4] - 21 alpha[8,8] - 21 alpha[8,9] + 131 beta[8,4]
   <= 130.997
 R62: 10 alpha[5,5] - 21 alpha[8,8] - 21 alpha[8,9] + 131 beta[8,5]
   <= 130.997
 R63: 21 alpha[6,6] + 21 alpha[6,8] + 21 alpha[6,9] - 21 alpha[8,8]
   - 21 alpha[8,9] + 131 beta[8,6] <= 130.997
 R64: 21 alpha[7,7] - 21 alpha[8,8] - 21 alpha[8,9] + 131 beta[8,7]
   <= 130.997
 R65: 21 alpha[0,0] + 21 alpha[0,2] + 21 alpha[0,6] + 21 alpha[0,8]
   + 21 alpha[0,9] - 21 alpha[9,9] + 131 beta[9,0] <= 130.997
 R66: 21 alpha[1,1] + 21 alpha[1,3] + 21 alpha[1,7] - 21 alpha[9,9]
   + 131 beta[9,1] <= 130.997
 R67: 21 alpha[2,2] + 21 alpha[2,6] + 21 alpha[2,8] + 21 alpha[2,9]
   - 21 alpha[9,9] + 131 beta[9,2] <= 130.997
 R68: 21 alpha[3,3] + 21 alpha[3,7] - 21 alpha[9,9] + 131 beta[9,3]
   <= 130.997
 R69: 21 alpha[4,4] - 21 alpha[9,9] + 131 beta[9,4] <= 130.997
 R70: 10 alpha[5,5] - 21 alpha[9,9] + 131 beta[9,5] <= 130.997
 R71: 21 alpha[6,6] + 21 alpha[6,8] + 21 alpha[6,9] - 21 alpha[9,9]
   + 131 beta[9,6] <= 130.997
 R72: 21 alpha[7,7] - 21 alpha[9,9] + 131 beta[9,7] <= 130.997
 R73: 21 alpha[8,8] + 21 alpha[8,9] - 21 alpha[9,9] + 131 beta[9,8]
   <= 130.997
 R74: 130 beta[0,0] <= 130
 R75: - 21 alpha[0,0] - 21 alpha[0,2] - 21 alpha[0,6] - 21 alpha[0,8]
   - 21 alpha[0,9] + 21 alpha[1,1] + 21 alpha[1,3] + 21 alpha[1,7]
   + 130 beta[0,1] <= 130
 R76: - 21 alpha[0,0] - 21 alpha[0,2] - 21 alpha[0,6] - 21 alpha[0,8]
   - 21 alpha[0,9] + 21 alpha[2,2] + 21 alpha[2,6] + 21 alpha[2,8]
   + 21 alpha[2,9] + 130 beta[0,2] <= 130
 R77: - 21 alpha[0,0] - 21 alpha[0,2] - 21 alpha[0,6] - 21 alpha[0,8]
   - 21 alpha[0,9] + 21 alpha[3,3] + 21 alpha[3,7] + 130 beta[0,3] <= 130
 R78: - 21 alpha[0,0] - 21 alpha[0,2] - 21 alpha[0,6] - 21 alpha[0,8]
   - 21 alpha[0,9] + 21 alpha[4,4] + 130 beta[0,4] <= 130
 R79: - 21 alpha[0,0] - 21 alpha[0,2] - 21 alpha[0,6] - 21 alpha[0,8]
   - 21 alpha[0,9] + 10 alpha[5,5] + 130 beta[0,5] <= 130
 R80: - 21 alpha[0,0] - 21 alpha[0,2] - 21 alpha[0,6] - 21 alpha[0,8]
   - 21 alpha[0,9] + 21 alpha[6,6] + 21 alpha[6,8] + 21 alpha[6,9]
   + 130 beta[0,6] <= 130
 R81: - 21 alpha[0,0] - 21 alpha[0,2] - 21 alpha[0,6] - 21 alpha[0,8]
   - 21 alpha[0,9] + 21 alpha[7,7] + 130 beta[0,7] <= 130
 R82: - 21 alpha[0,0] - 21 alpha[0,2] - 21 alpha[0,6] - 21 alpha[0,8]
   - 21 alpha[0,9] + 21 alpha[8,8] + 21 alpha[8,9] + 130 beta[0,8] <= 130
 R83: - 21 alpha[0,0] - 21 alpha[0,2] - 21 alpha[0,6] - 21 alpha[0,8]
   - 21 alpha[0,9] + 21 alpha[9,9] + 130 beta[0,9] <= 130
 R84: 130 beta[1,1] <= 130
 R85: - 21 alpha[1,1] - 21 alpha[1,3] - 21 alpha[1,7] + 21 alpha[2,2]
   + 21 alpha[2,6] + 21 alpha[2,8] + 21 alpha[2,9] + 130 beta[1,2] <= 130
 R86: - 21 alpha[1,1] - 21 alpha[1,3] - 21 alpha[1,7] + 21 alpha[3,3]
   + 21 alpha[3,7] + 130 beta[1,3] <= 130
 R87: - 21 alpha[1,1] - 21 alpha[1,3] - 21 alpha[1,7] + 21 alpha[4,4]
   + 130 beta[1,4] <= 130
 R88: - 21 alpha[1,1] - 21 alpha[1,3] - 21 alpha[1,7] + 10 alpha[5,5]
   + 130 beta[1,5] <= 130
 R89: - 21 alpha[1,1] - 21 alpha[1,3] - 21 alpha[1,7] + 21 alpha[6,6]
   + 21 alpha[6,8] + 21 alpha[6,9] + 130 beta[1,6] <= 130
 R90: - 21 alpha[1,1] - 21 alpha[1,3] - 21 alpha[1,7] + 21 alpha[7,7]
   + 130 beta[1,7] <= 130
 R91: - 21 alpha[1,1] - 21 alpha[1,3] - 21 alpha[1,7] + 21 alpha[8,8]
   + 21 alpha[8,9] + 130 beta[1,8] <= 130
 R92: - 21 alpha[1,1] - 21 alpha[1,3] - 21 alpha[1,7] + 21 alpha[9,9]
   + 130 beta[1,9] <= 130
 R93: 130 beta[2,2] <= 130
 R94: - 21 alpha[2,2] - 21 alpha[2,6] - 21 alpha[2,8] - 21 alpha[2,9]
   + 21 alpha[3,3] + 21 alpha[3,7] + 130 beta[2,3] <= 130
 R95: - 21 alpha[2,2] - 21 alpha[2,6] - 21 alpha[2,8] - 21 alpha[2,9]
   + 21 alpha[4,4] + 130 beta[2,4] <= 130
 R96: - 21 alpha[2,2] - 21 alpha[2,6] - 21 alpha[2,8] - 21 alpha[2,9]
   + 10 alpha[5,5] + 130 beta[2,5] <= 130
 R97: - 21 alpha[2,2] - 21 alpha[2,6] - 21 alpha[2,8] - 21 alpha[2,9]
   + 21 alpha[6,6] + 21 alpha[6,8] + 21 alpha[6,9] + 130 beta[2,6] <= 130
 R98: - 21 alpha[2,2] - 21 alpha[2,6] - 21 alpha[2,8] - 21 alpha[2,9]
   + 21 alpha[7,7] + 130 beta[2,7] <= 130
 R99: - 21 alpha[2,2] - 21 alpha[2,6] - 21 alpha[2,8] - 21 alpha[2,9]
   + 21 alpha[8,8] + 21 alpha[8,9] + 130 beta[2,8] <= 130
 R100: - 21 alpha[2,2] - 21 alpha[2,6] - 21 alpha[2,8] - 21 alpha[2,9]
   + 21 alpha[9,9] + 130 beta[2,9] <= 130
 R101: 130 beta[3,3] <= 130
 R102: - 21 alpha[3,3] - 21 alpha[3,7] + 21 alpha[4,4] + 130 beta[3,4]
   <= 130
 R103: - 21 alpha[3,3] - 21 alpha[3,7] + 10 alpha[5,5] + 130 beta[3,5]
   <= 130
 R104: - 21 alpha[3,3] - 21 alpha[3,7] + 21 alpha[6,6] + 21 alpha[6,8]
   + 21 alpha[6,9] + 130 beta[3,6] <= 130
 R105: - 21 alpha[3,3] - 21 alpha[3,7] + 21 alpha[7,7] + 130 beta[3,7]
   <= 130
 R106: - 21 alpha[3,3] - 21 alpha[3,7] + 21 alpha[8,8] + 21 alpha[8,9]
   + 130 beta[3,8] <= 130
 R107: - 21 alpha[3,3] - 21 alpha[3,7] + 21 alpha[9,9] + 130 beta[3,9]
   <= 130
 R108: 130 beta[4,4] <= 130
 R109: - 21 alpha[4,4] + 10 alpha[5,5] + 130 beta[4,5] <= 130
 R110: - 21 alpha[4,4] + 21 alpha[6,6] + 21 alpha[6,8] + 21 alpha[6,9]
   + 130 beta[4,6] <= 130
 R111: - 21 alpha[4,4] + 21 alpha[7,7] + 130 beta[4,7] <= 130
 R112: - 21 alpha[4,4] + 21 alpha[8,8] + 21 alpha[8,9] + 130 beta[4,8]
   <= 130
 R113: - 21 alpha[4,4] + 21 alpha[9,9] + 130 beta[4,9] <= 130
 R114: 130 beta[5,5] <= 130
 R115: - 10 alpha[5,5] + 21 alpha[6,6] + 21 alpha[6,8] + 21 alpha[6,9]
   + 130 beta[5,6] <= 130
 R116: - 10 alpha[5,5] + 21 alpha[7,7] + 130 beta[5,7] <= 130
 R117: - 10 alpha[5,5] + 21 alpha[8,8] + 21 alpha[8,9] + 130 beta[5,8]
   <= 130
 R118: - 10 alpha[5,5] + 21 alpha[9,9] + 130 beta[5,9] <= 130
 R119: 130 beta[6,6] <= 130
 R120: - 21 alpha[6,6] - 21 alpha[6,8] - 21 alpha[6,9] + 21 alpha[7,7]
   + 130 beta[6,7] <= 130
 R121: - 21 alpha[6,6] - 21 alpha[6,8] - 21 alpha[6,9] + 21 alpha[8,8]
   + 21 alpha[8,9] + 130 beta[6,8] <= 130
 R122: - 21 alpha[6,6] - 21 alpha[6,8] - 21 alpha[6,9] + 21 alpha[9,9]
   + 130 beta[6,9] <= 130
 R123: 130 beta[7,7] <= 130
 R124: - 21 alpha[7,7] + 21 alpha[8,8] + 21 alpha[8,9] + 130 beta[7,8]
   <= 130
 R125: - 21 alpha[7,7] + 21 alpha[9,9] + 130 beta[7,9] <= 130
 R126: 130 beta[8,8] <= 130
 R127: - 21 alpha[8,8] - 21 alpha[8,9] + 21 alpha[9,9] + 130 beta[8,9]
   <= 130
 R128: - 30.2 beta[0,0] + 29.7 beta[0,1] + 14.8 beta[0,2] + 29.7 beta[0,3]
   + 10 beta[0,4] + 5 beta[0,5] + 14.8 beta[0,6] + 29.7 beta[0,7]
   + 14.8 beta[0,8] + 14.8 beta[0,9] <= 0
 R129: 14.8 beta[1,0] - 15.3 beta[1,1] + 14.8 beta[1,2] + 29.7 beta[1,3]
   + 10 beta[1,4] + 5 beta[1,5] + 14.8 beta[1,6] + 29.7 beta[1,7]
   + 14.8 beta[1,8] + 14.8 beta[1,9] <= 0
 R130: 14.8 beta[2,0] + 29.7 beta[2,1] - 30.2 beta[2,2] + 29.7 beta[2,3]
   + 10 beta[2,4] + 5 beta[2,5] + 14.8 beta[2,6] + 29.7 beta[2,7]
   + 14.8 beta[2,8] + 14.8 beta[2,9] <= 0
 R131: 14.8 beta[3,0] + 29.7 beta[3,1] + 14.8 beta[3,2] - 15.3 beta[3,3]
   + 10 beta[3,4] + 5 beta[3,5] + 14.8 beta[3,6] + 29.7 beta[3,7]
   + 14.8 beta[3,8] + 14.8 beta[3,9] <= 0
 R132: 14.8 beta[4,0] + 29.7 beta[4,1] + 14.8 beta[4,2] + 29.7 beta[4,3]
   - 35 beta[4,4] + 5 beta[4,5] + 14.8 beta[4,6] + 29.7 beta[4,7]
   + 14.8 beta[4,8] + 14.8 beta[4,9] <= 0
 R133: 14.8 beta[5,0] + 29.7 beta[5,1] + 14.8 beta[5,2] + 29.7 beta[5,3]
   + 10 beta[5,4] - 40 beta[5,5] + 14.8 beta[5,6] + 29.7 beta[5,7]
   + 14.8 beta[5,8] + 14.8 beta[5,9] <= 0
 R134: 14.8 beta[6,0] + 29.7 beta[6,1] + 14.8 beta[6,2] + 29.7 beta[6,3]
   + 10 beta[6,4] + 5 beta[6,5] - 30.2 beta[6,6] + 29.7 beta[6,7]
   + 14.8 beta[6,8] + 14.8 beta[6,9] <= 0
 R135: 14.8 beta[7,0] + 29.7 beta[7,1] + 14.8 beta[7,2] + 29.7 beta[7,3]
   + 10 beta[7,4] + 5 beta[7,5] + 14.8 beta[7,6] - 15.3 beta[7,7]
   + 14.8 beta[7,8] + 14.8 beta[7,9] <= 0
 R136: 14.8 beta[8,0] + 29.7 beta[8,1] + 14.8 beta[8,2] + 29.7 beta[8,3]
   + 10 beta[8,4] + 5 beta[8,5] + 14.8 beta[8,6] + 29.7 beta[8,7]
   - 30.2 beta[8,8] + 14.8 beta[8,9] <= 0
 R137: 14.8 beta[9,0] + 29.7 beta[9,1] + 14.8 beta[9,2] + 29.7 beta[9,3]
   + 10 beta[9,4] + 5 beta[9,5] + 14.8 beta[9,6] + 29.7 beta[9,7]
   + 14.8 beta[9,8] - 30.2 beta[9,9] <= 0
 R138: - beta[0,0] + gamma[0,0] = 0
 R139: - beta[1,1] + gamma[0,1] + gamma[1,1] = 0
 R140: - beta[2,2] + gamma[0,2] + gamma[1,2] + gamma[2,2] = 0
 R141: - beta[3,3] + gamma[0,3] + gamma[1,3] + gamma[2,3] + gamma[3,3] = 0
 R142: - beta[4,4] + gamma[0,4] + gamma[1,4] + gamma[2,4] + gamma[3,4]
   + gamma[4,4] = 0
 R143: - beta[5,5] + gamma[0,5] + gamma[1,5] + gamma[2,5] + gamma[3,5]
   + gamma[4,5] + gamma[5,5] = 0
 R144: - beta[6,6] + gamma[0,6] + gamma[1,6] + gamma[2,6] + gamma[3,6]
   + gamma[4,6] + gamma[5,6] + gamma[6,6] = 0
 R145: - beta[7,7] + gamma[0,7] + gamma[1,7] + gamma[2,7] + gamma[3,7]
   + gamma[4,7] + gamma[5,7] + gamma[6,7] + gamma[7,7] = 0
 R146: - beta[8,8] + gamma[0,8] + gamma[1,8] + gamma[2,8] + gamma[3,8]
   + gamma[4,8] + gamma[5,8] + gamma[6,8] + gamma[7,8] + gamma[8,8] = 0
 R147: - beta[9,9] + gamma[0,9] + gamma[1,9] + gamma[2,9] + gamma[3,9]
   + gamma[4,9] + gamma[5,9] + gamma[6,9] + gamma[7,9] + gamma[8,9]
   + gamma[9,9] = 0
 R148: - 109 gamma[0,0] + 21 gamma[0,1] + 21 gamma[0,2] + 21 gamma[0,3]
   + 21 gamma[0,4] + 10 gamma[0,5] + 21 gamma[0,6] + 21 gamma[0,7]
   + 21 gamma[0,8] + 21 gamma[0,9] + 21 delta[0,1,0] + 21 delta[0,2,0]
   + 21 delta[0,2,1] + 21 delta[0,3,0] + 21 delta[0,3,1] + 21 delta[0,3,2]
   + 21 delta[0,4,0] + 21 delta[0,4,1] + 21 delta[0,4,2] + 21 delta[0,4,3]
   + 10 delta[0,5,0] + 10 delta[0,5,1] + 10 delta[0,5,2] + 10 delta[0,5,3]
   + 10 delta[0,5,4] + 21 delta[0,6,0] + 21 delta[0,6,1] + 21 delta[0,6,2]
   + 21 delta[0,6,3] + 21 delta[0,6,4] + 21 delta[0,6,5] + 21 delta[0,7,0]
   + 21 delta[0,7,1] + 21 delta[0,7,2] + 21 delta[0,7,3] + 21 delta[0,7,4]
   + 21 delta[0,7,5] + 21 delta[0,7,6] + 21 delta[0,8,0] + 21 delta[0,8,1]
   + 21 delta[0,8,2] + 21 delta[0,8,3] + 21 delta[0,8,4] + 21 delta[0,8,5]
   + 21 delta[0,8,6] + 21 delta[0,8,7] + 21 delta[0,9,0] + 21 delta[0,9,1]
   + 21 delta[0,9,2] + 21 delta[0,9,3] + 21 delta[0,9,4] + 21 delta[0,9,5]
   + 21 delta[0,9,6] + 21 delta[0,9,7] + 21 delta[0,9,8] <= 0
 R149: - 109 gamma[1,1] + 21 gamma[1,2] + 21 gamma[1,3] + 21 gamma[1,4]
   + 10 gamma[1,5] + 21 gamma[1,6] + 21 gamma[1,7] + 21 gamma[1,8]
   + 21 gamma[1,9] + 21 delta[1,2,1] + 21 delta[1,3,1] + 21 delta[1,3,2]
   + 21 delta[1,4,1] + 21 delta[1,4,2] + 21 delta[1,4,3] + 10 delta[1,5,1]
   + 10 delta[1,5,2] + 10 delta[1,5,3] + 10 delta[1,5,4] + 21 delta[1,6,1]
   + 21 delta[1,6,2] + 21 delta[1,6,3] + 21 delta[1,6,4] + 21 delta[1,6,5]
   + 21 delta[1,7,1] + 21 delta[1,7,2] + 21 delta[1,7,3] + 21 delta[1,7,4]
   + 21 delta[1,7,5] + 21 delta[1,7,6] + 21 delta[1,8,1] + 21 delta[1,8,2]
   + 21 delta[1,8,3] + 21 delta[1,8,4] + 21 delta[1,8,5] + 21 delta[1,8,6]
   + 21 delta[1,8,7] + 21 delta[1,9,1] + 21 delta[1,9,2] + 21 delta[1,9,3]
   + 21 delta[1,9,4] + 21 delta[1,9,5] + 21 delta[1,9,6] + 21 delta[1,9,7]
   + 21 delta[1,9,8] <= 0
 R150: - 109 gamma[2,2] + 21 gamma[2,3] + 21 gamma[2,4] + 10 gamma[2,5]
   + 21 gamma[2,6] + 21 gamma[2,7] + 21 gamma[2,8] + 21 gamma[2,9]
   + 21 delta[2,3,2] + 21 delta[2,4,2] + 21 delta[2,4,3] + 10 delta[2,5,2]
   + 10 delta[2,5,3] + 10 delta[2,5,4] + 21 delta[2,6,2] + 21 delta[2,6,3]
   + 21 delta[2,6,4] + 21 delta[2,6,5] + 21 delta[2,7,2] + 21 delta[2,7,3]
   + 21 delta[2,7,4] + 21 delta[2,7,5] + 21 delta[2,7,6] + 21 delta[2,8,2]
   + 21 delta[2,8,3] + 21 delta[2,8,4] + 21 delta[2,8,5] + 21 delta[2,8,6]
   + 21 delta[2,8,7] + 21 delta[2,9,2] + 21 delta[2,9,3] + 21 delta[2,9,4]
   + 21 delta[2,9,5] + 21 delta[2,9,6] + 21 delta[2,9,7] + 21 delta[2,9,8]
   <= 0
 R151: - 109 gamma[3,3] + 21 gamma[3,4] + 10 gamma[3,5] + 21 gamma[3,6]
   + 21 gamma[3,7] + 21 gamma[3,8] + 21 gamma[3,9] + 21 delta[3,4,3]
   + 10 delta[3,5,3] + 10 delta[3,5,4] + 21 delta[3,6,3] + 21 delta[3,6,4]
   + 21 delta[3,6,5] + 21 delta[3,7,3] + 21 delta[3,7,4] + 21 delta[3,7,5]
   + 21 delta[3,7,6] + 21 delta[3,8,3] + 21 delta[3,8,4] + 21 delta[3,8,5]
   + 21 delta[3,8,6] + 21 delta[3,8,7] + 21 delta[3,9,3] + 21 delta[3,9,4]
   + 21 delta[3,9,5] + 21 delta[3,9,6] + 21 delta[3,9,7] + 21 delta[3,9,8]
   <= 0
 R152: - 109 gamma[4,4] + 10 gamma[4,5] + 21 gamma[4,6] + 21 gamma[4,7]
   + 21 gamma[4,8] + 21 gamma[4,9] + 10 delta[4,5,4] + 21 delta[4,6,4]
   + 21 delta[4,6,5] + 21 delta[4,7,4] + 21 delta[4,7,5] + 21 delta[4,7,6]
   + 21 delta[4,8,4] + 21 delta[4,8,5] + 21 delta[4,8,6] + 21 delta[4,8,7]
   + 21 delta[4,9,4] + 21 delta[4,9,5] + 21 delta[4,9,6] + 21 delta[4,9,7]
   + 21 delta[4,9,8] <= 0
 R153: - 120 gamma[5,5] + 21 gamma[5,6] + 21 gamma[5,7] + 21 gamma[5,8]
   + 21 gamma[5,9] + 21 delta[5,6,5] + 21 delta[5,7,5] + 21 delta[5,7,6]
   + 21 delta[5,8,5] + 21 delta[5,8,6] + 21 delta[5,8,7] + 21 delta[5,9,5]
   + 21 delta[5,9,6] + 21 delta[5,9,7] + 21 delta[5,9,8] <= 0
 R154: - 109 gamma[6,6] + 21 gamma[6,7] + 21 gamma[6,8] + 21 gamma[6,9]
   + 21 delta[6,7,6] + 21 delta[6,8,6] + 21 delta[6,8,7] + 21 delta[6,9,6]
   + 21 delta[6,9,7] + 21 delta[6,9,8] <= 0
 R155: - 109 gamma[7,7] + 21 gamma[7,8] + 21 gamma[7,9] + 21 delta[7,8,7]
   + 21 delta[7,9,7] + 21 delta[7,9,8] <= 0
 R156: - 109 gamma[8,8] + 21 gamma[8,9] + 21 delta[8,9,8] <= 0
 R157: alpha[0,2] + gamma[0,0] - delta[0,2,0] <= 1
 R158: alpha[1,3] + gamma[0,1] - delta[0,3,1] <= 1
 R159: alpha[0,6] + gamma[0,0] - delta[0,6,0] <= 1
 R160: alpha[2,6] + gamma[0,2] - delta[0,6,2] <= 1
 R161: alpha[1,7] + gamma[0,1] - delta[0,7,1] <= 1
 R162: alpha[3,7] + gamma[0,3] - delta[0,7,3] <= 1
 R163: alpha[0,8] + gamma[0,0] - delta[0,8,0] <= 1
 R164: alpha[2,8] + gamma[0,2] - delta[0,8,2] <= 1
 R165: alpha[6,8] + gamma[0,6] - delta[0,8,6] <= 1
 R166: alpha[0,9] + gamma[0,0] - delta[0,9,0] <= 1
 R167: alpha[2,9] + gamma[0,2] - delta[0,9,2] <= 1
 R168: alpha[6,9] + gamma[0,6] - delta[0,9,6] <= 1
 R169: alpha[8,9] + gamma[0,8] - delta[0,9,8] <= 1
 R170: alpha[1,3] + gamma[1,1] - delta[1,3,1] <= 1
 R171: alpha[2,6] + gamma[1,2] - delta[1,6,2] <= 1
 R172: alpha[1,7] + gamma[1,1] - delta[1,7,1] <= 1
 R173: alpha[3,7] + gamma[1,3] - delta[1,7,3] <= 1
 R174: alpha[2,8] + gamma[1,2] - delta[1,8,2] <= 1
 R175: alpha[6,8] + gamma[1,6] - delta[1,8,6] <= 1
 R176: alpha[2,9] + gamma[1,2] - delta[1,9,2] <= 1
 R177: alpha[6,9] + gamma[1,6] - delta[1,9,6] <= 1
 R178: alpha[8,9] + gamma[1,8] - delta[1,9,8] <= 1
 R179: alpha[2,6] + gamma[2,2] - delta[2,6,2] <= 1
 R180: alpha[3,7] + gamma[2,3] - delta[2,7,3] <= 1
 R181: alpha[2,8] + gamma[2,2] - delta[2,8,2] <= 1
 R182: alpha[6,8] + gamma[2,6] - delta[2,8,6] <= 1
 R183: alpha[2,9] + gamma[2,2] - delta[2,9,2] <= 1
 R184: alpha[6,9] + gamma[2,6] - delta[2,9,6] <= 1
 R185: alpha[8,9] + gamma[2,8] - delta[2,9,8] <= 1
 R186: alpha[3,7] + gamma[3,3] - delta[3,7,3] <= 1
 R187: alpha[6,8] + gamma[3,6] - delta[3,8,6] <= 1
 R188: alpha[6,9] + gamma[3,6] - delta[3,9,6] <= 1
 R189: alpha[8,9] + gamma[3,8] - delta[3,9,8] <= 1
 R190: alpha[6,8] + gamma[4,6] - delta[4,8,6] <= 1
 R191: alpha[6,9] + gamma[4,6] - delta[4,9,6] <= 1
 R192: alpha[8,9] + gamma[4,8] - delta[4,9,8] <= 1
 R193: alpha[6,8] + gamma[5,6] - delta[5,8,6] <= 1
 R194: alpha[6,9] + gamma[5,6] - delta[5,9,6] <= 1
 R195: alpha[8,9] + gamma[5,8] - delta[5,9,8] <= 1
 R196: alpha[6,8] + gamma[6,6] - delta[6,8,6] <= 1
 R197: alpha[6,9] + gamma[6,6] - delta[6,9,6] <= 1
 R198: alpha[8,9] + gamma[6,8] - delta[6,9,8] <= 1
 R199: alpha[8,9] + gamma[7,8] - delta[7,9,8] <= 1
 R200: alpha[8,9] + gamma[8,8] - delta[8,9,8] <= 1
 R201: 0.5 alpha[0,2] + 0.5 gamma[0,0] - delta[0,2,0] >= 0
 R202: 0.5 alpha[1,3] + 0.5 gamma[0,1] - delta[0,3,1] >= 0
 R203: 0.5 alpha[0,6] + 0.5 gamma[0,0] - delta[0,6,0] >= 0
 R204: 0.5 alpha[2,6] + 0.5 gamma[0,2] - delta[0,6,2] >= 0
 R205: 0.5 alpha[1,7] + 0.5 gamma[0,1] - delta[0,7,1] >= 0
 R206: 0.5 alpha[3,7] + 0.5 gamma[0,3] - delta[0,7,3] >= 0
 R207: 0.5 alpha[0,8] + 0.5 gamma[0,0] - delta[0,8,0] >= 0
 R208: 0.5 alpha[2,8] + 0.5 gamma[0,2] - delta[0,8,2] >= 0
 R209: 0.5 alpha[6,8] + 0.5 gamma[0,6] - delta[0,8,6] >= 0
 R210: 0.5 alpha[0,9] + 0.5 gamma[0,0] - delta[0,9,0] >= 0
 R211: 0.5 alpha[2,9] + 0.5 gamma[0,2] - delta[0,9,2] >= 0
 R212: 0.5 alpha[6,9] + 0.5 gamma[0,6] - delta[0,9,6] >= 0
 R213: 0.5 alpha[8,9] + 0.5 gamma[0,8] - delta[0,9,8] >= 0
 R214: 0.5 alpha[1,3] + 0.5 gamma[1,1] - delta[1,3,1] >= 0
 R215: 0.5 alpha[2,6] + 0.5 gamma[1,2] - delta[1,6,2] >= 0
 R216: 0.5 alpha[1,7] + 0.5 gamma[1,1] - delta[1,7,1] >= 0
 R217: 0.5 alpha[3,7] + 0.5 gamma[1,3] - delta[1,7,3] >= 0
 R218: 0.5 alpha[2,8] + 0.5 gamma[1,2] - delta[1,8,2] >= 0
 R219: 0.5 alpha[6,8] + 0.5 gamma[1,6] - delta[1,8,6] >= 0
 R220: 0.5 alpha[2,9] + 0.5 gamma[1,2] - delta[1,9,2] >= 0
 R221: 0.5 alpha[6,9] + 0.5 gamma[1,6] - delta[1,9,6] >= 0
 R222: 0.5 alpha[8,9] + 0.5 gamma[1,8] - delta[1,9,8] >= 0
 R223: 0.5 alpha[2,6] + 0.5 gamma[2,2] - delta[2,6,2] >= 0
 R224: 0.5 alpha[3,7] + 0.5 gamma[2,3] - delta[2,7,3] >= 0
 R225: 0.5 alpha[2,8] + 0.5 gamma[2,2] - delta[2,8,2] >= 0
 R226: 0.5 alpha[6,8] + 0.5 gamma[2,6] - delta[2,8,6] >= 0
 R227: 0.5 alpha[2,9] + 0.5 gamma[2,2] - delta[2,9,2] >= 0
 R228: 0.5 alpha[6,9] + 0.5 gamma[2,6] - delta[2,9,6] >= 0
 R229: 0.5 alpha[8,9] + 0.5 gamma[2,8] - delta[2,9,8] >= 0
 R230: 0.5 alpha[3,7] + 0.5 gamma[3,3] - delta[3,7,3] >= 0
 R231: 0.5 alpha[6,8] + 0.5 gamma[3,6] - delta[3,8,6] >= 0
 R232: 0.5 alpha[6,9] + 0.5 gamma[3,6] - delta[3,9,6] >= 0
 R233: 0.5 alpha[8,9] + 0.5 gamma[3,8] - delta[3,9,8] >= 0
 R234: 0.5 alpha[6,8] + 0.5 gamma[4,6] - delta[4,8,6] >= 0
 R235: 0.5 alpha[6,9] + 0.5 gamma[4,6] - delta[4,9,6] >= 0
 R236: 0.5 alpha[8,9] + 0.5 gamma[4,8] - delta[4,9,8] >= 0
 R237: 0.5 alpha[6,8] + 0.5 gamma[5,6] - delta[5,8,6] >= 0
 R238: 0.5 alpha[6,9] + 0.5 gamma[5,6] - delta[5,9,6] >= 0
 R239: 0.5 alpha[8,9] + 0.5 gamma[5,8] - delta[5,9,8] >= 0
 R240: 0.5 alpha[6,8] + 0.5 gamma[6,6] - delta[6,8,6] >= 0
 R241: 0.5 alpha[6,9] + 0.5 gamma[6,6] - delta[6,9,6] >= 0
 R242: 0.5 alpha[8,9] + 0.5 gamma[6,8] - delta[6,9,8] >= 0
 R243: 0.5 alpha[8,9] + 0.5 gamma[7,8] - delta[7,9,8] >= 0
 R244: 0.5 alpha[8,9] + 0.5 gamma[8,8] - delta[8,9,8] >= 0
 R245: - 10 gamma[0,0] + gamma[0,1] + gamma[0,2] + gamma[0,3] + gamma[0,4]
   + gamma[0,5] + gamma[0,6] + gamma[0,7] + gamma[0,8] + gamma[0,9] <= 0
 R246: - 9 gamma[1,1] + gamma[1,2] + gamma[1,3] + gamma[1,4] + gamma[1,5]
   + gamma[1,6] + gamma[1,7] + gamma[1,8] + gamma[1,9] <= 0
 R247: - 8 gamma[2,2] + gamma[2,3] + gamma[2,4] + gamma[2,5] + gamma[2,6]
   + gamma[2,7] + gamma[2,8] + gamma[2,9] <= 0
 R248: - 7 gamma[3,3] + gamma[3,4] + gamma[3,5] + gamma[3,6] + gamma[3,7]
   + gamma[3,8] + gamma[3,9] <= 0
 R249: - 6 gamma[4,4] + gamma[4,5] + gamma[4,6] + gamma[4,7] + gamma[4,8]
   + gamma[4,9] <= 0
 R250: - 5 gamma[5,5] + gamma[5,6] + gamma[5,7] + gamma[5,8] + gamma[5,9]
   <= 0
 R251: - 4 gamma[6,6] + gamma[6,7] + gamma[6,8] + gamma[6,9] <= 0
 R252: - 3 gamma[7,7] + gamma[7,8] + gamma[7,9] <= 0
 R253: - 2 gamma[8,8] + gamma[8,9] <= 0
 R254: - gamma[0,0] + x[0,0] + x[0,1] + x[0,2] = 0
 R255: - gamma[1,1] + x[1,0] + x[1,1] + x[1,2] = 0
 R256: - gamma[2,2] + x[2,0] + x[2,1] + x[2,2] = 0
 R257: - gamma[3,3] + x[3,0] + x[3,1] + x[3,2] = 0
 R258: - gamma[4,4] + x[4,0] + x[4,1] + x[4,2] = 0
 R259: - gamma[5,5] + x[5,0] + x[5,1] + x[5,2] = 0
 R260: - gamma[6,6] + x[6,0] + x[6,1] + x[6,2] = 0
 R261: - gamma[7,7] + x[7,0] + x[7,1] + x[7,2] = 0
 R262: - gamma[8,8] + x[8,0] + x[8,1] + x[8,2] = 0
 R263: - gamma[9,9] + x[9,0] + x[9,1] + x[9,2] = 0
 R264: alpha[0,0] + beta[0,0] + gamma[0,0] - a[0,0,0,0] <= 2
 R265: alpha[0,0] + beta[1,0] + gamma[0,1] - a[0,1,0,0] <= 2
 R266: alpha[0,0] + beta[1,0] + gamma[1,1] - a[1,1,0,0] <= 2
 R267: alpha[0,0] + beta[2,0] + gamma[0,2] - a[0,2,0,0] <= 2
 R268: alpha[0,0] + beta[2,0] + gamma[1,2] - a[1,2,0,0] <= 2
 R269: alpha[0,0] + beta[2,0] + gamma[2,2] - a[2,2,0,0] <= 2
 R270: alpha[0,0] + beta[3,0] + gamma[0,3] - a[0,3,0,0] <= 2
 R271: alpha[0,0] + beta[3,0] + gamma[1,3] - a[1,3,0,0] <= 2
 R272: alpha[0,0] + beta[3,0] + gamma[2,3] - a[2,3,0,0] <= 2
 R273: alpha[0,0] + beta[3,0] + gamma[3,3] - a[3,3,0,0] <= 2
 R274: alpha[0,0] + beta[4,0] + gamma[0,4] - a[0,4,0,0] <= 2
 R275: alpha[0,0] + beta[4,0] + gamma[1,4] - a[1,4,0,0] <= 2
 R276: alpha[0,0] + beta[4,0] + gamma[2,4] - a[2,4,0,0] <= 2
 R277: alpha[0,0] + beta[4,0] + gamma[3,4] - a[3,4,0,0] <= 2
 R278: alpha[0,0] + beta[4,0] + gamma[4,4] - a[4,4,0,0] <= 2
 R279: alpha[0,0] + beta[5,0] + gamma[0,5] - a[0,5,0,0] <= 2
 R280: alpha[0,0] + beta[5,0] + gamma[1,5] - a[1,5,0,0] <= 2
 R281: alpha[0,0] + beta[5,0] + gamma[2,5] - a[2,5,0,0] <= 2
 R282: alpha[0,0] + beta[5,0] + gamma[3,5] - a[3,5,0,0] <= 2
 R283: alpha[0,0] + beta[5,0] + gamma[4,5] - a[4,5,0,0] <= 2
 R284: alpha[0,0] + beta[5,0] + gamma[5,5] - a[5,5,0,0] <= 2
 R285: alpha[0,0] + beta[6,0] + gamma[0,6] - a[0,6,0,0] <= 2
 R286: alpha[0,0] + beta[6,0] + gamma[1,6] - a[1,6,0,0] <= 2
 R287: alpha[0,0] + beta[6,0] + gamma[2,6] - a[2,6,0,0] <= 2
 R288: alpha[0,0] + beta[6,0] + gamma[3,6] - a[3,6,0,0] <= 2
 R289: alpha[0,0] + beta[6,0] + gamma[4,6] - a[4,6,0,0] <= 2
 R290: alpha[0,0] + beta[6,0] + gamma[5,6] - a[5,6,0,0] <= 2
 R291: alpha[0,0] + beta[6,0] + gamma[6,6] - a[6,6,0,0] <= 2
 R292: alpha[0,0] + beta[7,0] + gamma[0,7] - a[0,7,0,0] <= 2
 R293: alpha[0,0] + beta[7,0] + gamma[1,7] - a[1,7,0,0] <= 2
 R294: alpha[0,0] + beta[7,0] + gamma[2,7] - a[2,7,0,0] <= 2
 R295: alpha[0,0] + beta[7,0] + gamma[3,7] - a[3,7,0,0] <= 2
 R296: alpha[0,0] + beta[7,0] + gamma[4,7] - a[4,7,0,0] <= 2
 R297: alpha[0,0] + beta[7,0] + gamma[5,7] - a[5,7,0,0] <= 2
 R298: alpha[0,0] + beta[7,0] + gamma[6,7] - a[6,7,0,0] <= 2
 R299: alpha[0,0] + beta[7,0] + gamma[7,7] - a[7,7,0,0] <= 2
 R300: alpha[0,0] + beta[8,0] + gamma[0,8] - a[0,8,0,0] <= 2
 R301: alpha[0,0] + beta[8,0] + gamma[1,8] - a[1,8,0,0] <= 2
 R302: alpha[0,0] + beta[8,0] + gamma[2,8] - a[2,8,0,0] <= 2
 R303: alpha[0,0] + beta[8,0] + gamma[3,8] - a[3,8,0,0] <= 2
 R304: alpha[0,0] + beta[8,0] + gamma[4,8] - a[4,8,0,0] <= 2
 R305: alpha[0,0] + beta[8,0] + gamma[5,8] - a[5,8,0,0] <= 2
 R306: alpha[0,0] + beta[8,0] + gamma[6,8] - a[6,8,0,0] <= 2
 R307: alpha[0,0] + beta[8,0] + gamma[7,8] - a[7,8,0,0] <= 2
 R308: alpha[0,0] + beta[8,0] + gamma[8,8] - a[8,8,0,0] <= 2
 R309: alpha[0,0] + beta[9,0] + gamma[0,9] - a[0,9,0,0] <= 2
 R310: alpha[0,0] + beta[9,0] + gamma[1,9] - a[1,9,0,0] <= 2
 R311: alpha[0,0] + beta[9,0] + gamma[2,9] - a[2,9,0,0] <= 2
 R312: alpha[0,0] + beta[9,0] + gamma[3,9] - a[3,9,0,0] <= 2
 R313: alpha[0,0] + beta[9,0] + gamma[4,9] - a[4,9,0,0] <= 2
 R314: alpha[0,0] + beta[9,0] + gamma[5,9] - a[5,9,0,0] <= 2
 R315: alpha[0,0] + beta[9,0] + gamma[6,9] - a[6,9,0,0] <= 2
 R316: alpha[0,0] + beta[9,0] + gamma[7,9] - a[7,9,0,0] <= 2
 R317: alpha[0,0] + beta[9,0] + gamma[8,9] - a[8,9,0,0] <= 2
 R318: alpha[0,0] + beta[9,0] + gamma[9,9] - a[9,9,0,0] <= 2
 R319: alpha[0,2] + beta[0,0] + gamma[0,0] - a[0,0,0,2] <= 2
 R320: alpha[0,2] + beta[1,0] + gamma[0,1] - a[0,1,0,2] <= 2
 R321: alpha[0,2] + beta[1,0] + gamma[1,1] - a[1,1,0,2] <= 2
 R322: alpha[0,2] + beta[2,0] + gamma[0,2] - a[0,2,0,2] <= 2
 R323: alpha[0,2] + beta[2,0] + gamma[1,2] - a[1,2,0,2] <= 2
 R324: alpha[0,2] + beta[2,0] + gamma[2,2] - a[2,2,0,2] <= 2
 R325: alpha[0,2] + beta[3,0] + gamma[0,3] - a[0,3,0,2] <= 2
 R326: alpha[0,2] + beta[3,0] + gamma[1,3] - a[1,3,0,2] <= 2
 R327: alpha[0,2] + beta[3,0] + gamma[2,3] - a[2,3,0,2] <= 2
 R328: alpha[0,2] + beta[3,0] + gamma[3,3] - a[3,3,0,2] <= 2
 R329: alpha[0,2] + beta[4,0] + gamma[0,4] - a[0,4,0,2] <= 2
 R330: alpha[0,2] + beta[4,0] + gamma[1,4] - a[1,4,0,2] <= 2
 R331: alpha[0,2] + beta[4,0] + gamma[2,4] - a[2,4,0,2] <= 2
 R332: alpha[0,2] + beta[4,0] + gamma[3,4] - a[3,4,0,2] <= 2
 R333: alpha[0,2] + beta[4,0] + gamma[4,4] - a[4,4,0,2] <= 2
 R334: alpha[0,2] + beta[5,0] + gamma[0,5] - a[0,5,0,2] <= 2
 R335: alpha[0,2] + beta[5,0] + gamma[1,5] - a[1,5,0,2] <= 2
 R336: alpha[0,2] + beta[5,0] + gamma[2,5] - a[2,5,0,2] <= 2
 R337: alpha[0,2] + beta[5,0] + gamma[3,5] - a[3,5,0,2] <= 2
 R338: alpha[0,2] + beta[5,0] + gamma[4,5] - a[4,5,0,2] <= 2
 R339: alpha[0,2] + beta[5,0] + gamma[5,5] - a[5,5,0,2] <= 2
 R340: alpha[0,2] + beta[6,0] + gamma[0,6] - a[0,6,0,2] <= 2
 R341: alpha[0,2] + beta[6,0] + gamma[1,6] - a[1,6,0,2] <= 2
 R342: alpha[0,2] + beta[6,0] + gamma[2,6] - a[2,6,0,2] <= 2
 R343: alpha[0,2] + beta[6,0] + gamma[3,6] - a[3,6,0,2] <= 2
 R344: alpha[0,2] + beta[6,0] + gamma[4,6] - a[4,6,0,2] <= 2
 R345: alpha[0,2] + beta[6,0] + gamma[5,6] - a[5,6,0,2] <= 2
 R346: alpha[0,2] + beta[6,0] + gamma[6,6] - a[6,6,0,2] <= 2
 R347: alpha[0,2] + beta[7,0] + gamma[0,7] - a[0,7,0,2] <= 2
 R348: alpha[0,2] + beta[7,0] + gamma[1,7] - a[1,7,0,2] <= 2
 R349: alpha[0,2] + beta[7,0] + gamma[2,7] - a[2,7,0,2] <= 2
 R350: alpha[0,2] + beta[7,0] + gamma[3,7] - a[3,7,0,2] <= 2
 R351: alpha[0,2] + beta[7,0] + gamma[4,7] - a[4,7,0,2] <= 2
 R352: alpha[0,2] + beta[7,0] + gamma[5,7] - a[5,7,0,2] <= 2
 R353: alpha[0,2] + beta[7,0] + gamma[6,7] - a[6,7,0,2] <= 2
 R354: alpha[0,2] + beta[7,0] + gamma[7,7] - a[7,7,0,2] <= 2
 R355: alpha[0,2] + beta[8,0] + gamma[0,8] - a[0,8,0,2] <= 2
 R356: alpha[0,2] + beta[8,0] + gamma[1,8] - a[1,8,0,2] <= 2
 R357: alpha[0,2] + beta[8,0] + gamma[2,8] - a[2,8,0,2] <= 2
 R358: alpha[0,2] + beta[8,0] + gamma[3,8] - a[3,8,0,2] <= 2
 R359: alpha[0,2] + beta[8,0] + gamma[4,8] - a[4,8,0,2] <= 2
 R360: alpha[0,2] + beta[8,0] + gamma[5,8] - a[5,8,0,2] <= 2
 R361: alpha[0,2] + beta[8,0] + gamma[6,8] - a[6,8,0,2] <= 2
 R362: alpha[0,2] + beta[8,0] + gamma[7,8] - a[7,8,0,2] <= 2
 R363: alpha[0,2] + beta[8,0] + gamma[8,8] - a[8,8,0,2] <= 2
 R364: alpha[0,2] + beta[9,0] + gamma[0,9] - a[0,9,0,2] <= 2
 R365: alpha[0,2] + beta[9,0] + gamma[1,9] - a[1,9,0,2] <= 2
 R366: alpha[0,2] + beta[9,0] + gamma[2,9] - a[2,9,0,2] <= 2
 R367: alpha[0,2] + beta[9,0] + gamma[3,9] - a[3,9,0,2] <= 2
 R368: alpha[0,2] + beta[9,0] + gamma[4,9] - a[4,9,0,2] <= 2
 R369: alpha[0,2] + beta[9,0] + gamma[5,9] - a[5,9,0,2] <= 2
 R370: alpha[0,2] + beta[9,0] + gamma[6,9] - a[6,9,0,2] <= 2
 R371: alpha[0,2] + beta[9,0] + gamma[7,9] - a[7,9,0,2] <= 2
 R372: alpha[0,2] + beta[9,0] + gamma[8,9] - a[8,9,0,2] <= 2
 R373: alpha[0,2] + beta[9,0] + gamma[9,9] - a[9,9,0,2] <= 2
 R374: alpha[0,6] + beta[0,0] + gamma[0,0] - a[0,0,0,6] <= 2
 R375: alpha[0,6] + beta[1,0] + gamma[0,1] - a[0,1,0,6] <= 2
 R376: alpha[0,6] + beta[1,0] + gamma[1,1] - a[1,1,0,6] <= 2
 R377: alpha[0,6] + beta[2,0] + gamma[0,2] - a[0,2,0,6] <= 2
 R378: alpha[0,6] + beta[2,0] + gamma[1,2] - a[1,2,0,6] <= 2
 R379: alpha[0,6] + beta[2,0] + gamma[2,2] - a[2,2,0,6] <= 2
 R380: alpha[0,6] + beta[3,0] + gamma[0,3] - a[0,3,0,6] <= 2
 R381: alpha[0,6] + beta[3,0] + gamma[1,3] - a[1,3,0,6] <= 2
 R382: alpha[0,6] + beta[3,0] + gamma[2,3] - a[2,3,0,6] <= 2
 R383: alpha[0,6] + beta[3,0] + gamma[3,3] - a[3,3,0,6] <= 2
 R384: alpha[0,6] + beta[4,0] + gamma[0,4] - a[0,4,0,6] <= 2
 R385: alpha[0,6] + beta[4,0] + gamma[1,4] - a[1,4,0,6] <= 2
 R386: alpha[0,6] + beta[4,0] + gamma[2,4] - a[2,4,0,6] <= 2
 R387: alpha[0,6] + beta[4,0] + gamma[3,4] - a[3,4,0,6] <= 2
 R388: alpha[0,6] + beta[4,0] + gamma[4,4] - a[4,4,0,6] <= 2
 R389: alpha[0,6] + beta[5,0] + gamma[0,5] - a[0,5,0,6] <= 2
 R390: alpha[0,6] + beta[5,0] + gamma[1,5] - a[1,5,0,6] <= 2
 R391: alpha[0,6] + beta[5,0] + gamma[2,5] - a[2,5,0,6] <= 2
 R392: alpha[0,6] + beta[5,0] + gamma[3,5] - a[3,5,0,6] <= 2
 R393: alpha[0,6] + beta[5,0] + gamma[4,5] - a[4,5,0,6] <= 2
 R394: alpha[0,6] + beta[5,0] + gamma[5,5] - a[5,5,0,6] <= 2
 R395: alpha[0,6] + beta[6,0] + gamma[0,6] - a[0,6,0,6] <= 2
 R396: alpha[0,6] + beta[6,0] + gamma[1,6] - a[1,6,0,6] <= 2
 R397: alpha[0,6] + beta[6,0] + gamma[2,6] - a[2,6,0,6] <= 2
 R398: alpha[0,6] + beta[6,0] + gamma[3,6] - a[3,6,0,6] <= 2
 R399: alpha[0,6] + beta[6,0] + gamma[4,6] - a[4,6,0,6] <= 2
 R400: alpha[0,6] + beta[6,0] + gamma[5,6] - a[5,6,0,6] <= 2
 R401: alpha[0,6] + beta[6,0] + gamma[6,6] - a[6,6,0,6] <= 2
 R402: alpha[0,6] + beta[7,0] + gamma[0,7] - a[0,7,0,6] <= 2
 R403: alpha[0,6] + beta[7,0] + gamma[1,7] - a[1,7,0,6] <= 2
 R404: alpha[0,6] + beta[7,0] + gamma[2,7] - a[2,7,0,6] <= 2
 R405: alpha[0,6] + beta[7,0] + gamma[3,7] - a[3,7,0,6] <= 2
 R406: alpha[0,6] + beta[7,0] + gamma[4,7] - a[4,7,0,6] <= 2
 R407: alpha[0,6] + beta[7,0] + gamma[5,7] - a[5,7,0,6] <= 2
 R408: alpha[0,6] + beta[7,0] + gamma[6,7] - a[6,7,0,6] <= 2
 R409: alpha[0,6] + beta[7,0] + gamma[7,7] - a[7,7,0,6] <= 2
 R410: alpha[0,6] + beta[8,0] + gamma[0,8] - a[0,8,0,6] <= 2
 R411: alpha[0,6] + beta[8,0] + gamma[1,8] - a[1,8,0,6] <= 2
 R412: alpha[0,6] + beta[8,0] + gamma[2,8] - a[2,8,0,6] <= 2
 R413: alpha[0,6] + beta[8,0] + gamma[3,8] - a[3,8,0,6] <= 2
 R414: alpha[0,6] + beta[8,0] + gamma[4,8] - a[4,8,0,6] <= 2
 R415: alpha[0,6] + beta[8,0] + gamma[5,8] - a[5,8,0,6] <= 2
 R416: alpha[0,6] + beta[8,0] + gamma[6,8] - a[6,8,0,6] <= 2
 R417: alpha[0,6] + beta[8,0] + gamma[7,8] - a[7,8,0,6] <= 2
 R418: alpha[0,6] + beta[8,0] + gamma[8,8] - a[8,8,0,6] <= 2
 R419: alpha[0,6] + beta[9,0] + gamma[0,9] - a[0,9,0,6] <= 2
 R420: alpha[0,6] + beta[9,0] + gamma[1,9] - a[1,9,0,6] <= 2
 R421: alpha[0,6] + beta[9,0] + gamma[2,9] - a[2,9,0,6] <= 2
 R422: alpha[0,6] + beta[9,0] + gamma[3,9] - a[3,9,0,6] <= 2
 R423: alpha[0,6] + beta[9,0] + gamma[4,9] - a[4,9,0,6] <= 2
 R424: alpha[0,6] + beta[9,0] + gamma[5,9] - a[5,9,0,6] <= 2
 R425: alpha[0,6] + beta[9,0] + gamma[6,9] - a[6,9,0,6] <= 2
 R426: alpha[0,6] + beta[9,0] + gamma[7,9] - a[7,9,0,6] <= 2
 R427: alpha[0,6] + beta[9,0] + gamma[8,9] - a[8,9,0,6] <= 2
 R428: alpha[0,6] + beta[9,0] + gamma[9,9] - a[9,9,0,6] <= 2
 R429: alpha[0,8] + beta[0,0] + gamma[0,0] - a[0,0,0,8] <= 2
 R430: alpha[0,8] + beta[1,0] + gamma[0,1] - a[0,1,0,8] <= 2
 R431: alpha[0,8] + beta[1,0] + gamma[1,1] - a[1,1,0,8] <= 2
 R432: alpha[0,8] + beta[2,0] + gamma[0,2] - a[0,2,0,8] <= 2
 R433: alpha[0,8] + beta[2,0] + gamma[1,2] - a[1,2,0,8] <= 2
 R434: alpha[0,8] + beta[2,0] + gamma[2,2] - a[2,2,0,8] <= 2
 R435: alpha[0,8] + beta[3,0] + gamma[0,3] - a[0,3,0,8] <= 2
 R436: alpha[0,8] + beta[3,0] + gamma[1,3] - a[1,3,0,8] <= 2
 R437: alpha[0,8] + beta[3,0] + gamma[2,3] - a[2,3,0,8] <= 2
 R438: alpha[0,8] + beta[3,0] + gamma[3,3] - a[3,3,0,8] <= 2
 R439: alpha[0,8] + beta[4,0] + gamma[0,4] - a[0,4,0,8] <= 2
 R440: alpha[0,8] + beta[4,0] + gamma[1,4] - a[1,4,0,8] <= 2
 R441: alpha[0,8] + beta[4,0] + gamma[2,4] - a[2,4,0,8] <= 2
 R442: alpha[0,8] + beta[4,0] + gamma[3,4] - a[3,4,0,8] <= 2
 R443: alpha[0,8] + beta[4,0] + gamma[4,4] - a[4,4,0,8] <= 2
 R444: alpha[0,8] + beta[5,0] + gamma[0,5] - a[0,5,0,8] <= 2
 R445: alpha[0,8] + beta[5,0] + gamma[1,5] - a[1,5,0,8] <= 2
 R446: alpha[0,8] + beta[5,0] + gamma[2,5] - a[2,5,0,8] <= 2
 R447: alpha[0,8] + beta[5,0] + gamma[3,5] - a[3,5,0,8] <= 2
 R448: alpha[0,8] + beta[5,0] + gamma[4,5] - a[4,5,0,8] <= 2
 R449: alpha[0,8] + beta[5,0] + gamma[5,5] - a[5,5,0,8] <= 2
 R450: alpha[0,8] + beta[6,0] + gamma[0,6] - a[0,6,0,8] <= 2
 R451: alpha[0,8] + beta[6,0] + gamma[1,6] - a[1,6,0,8] <= 2
 R452: alpha[0,8] + beta[6,0] + gamma[2,6] - a[2,6,0,8] <= 2
 R453: alpha[0,8] + beta[6,0] + gamma[3,6] - a[3,6,0,8] <= 2
 R454: alpha[0,8] + beta[6,0] + gamma[4,6] - a[4,6,0,8] <= 2
 R455: alpha[0,8] + beta[6,0] + gamma[5,6] - a[5,6,0,8] <= 2
 R456: alpha[0,8] + beta[6,0] + gamma[6,6] - a[6,6,0,8] <= 2
 R457: alpha[0,8] + beta[7,0] + gamma[0,7] - a[0,7,0,8] <= 2
 R458: alpha[0,8] + beta[7,0] + gamma[1,7] - a[1,7,0,8] <= 2
 R459: alpha[0,8] + beta[7,0] + gamma[2,7] - a[2,7,0,8] <= 2
 R460: alpha[0,8] + beta[7,0] + gamma[3,7] - a[3,7,0,8] <= 2
 R461: alpha[0,8] + beta[7,0] + gamma[4,7] - a[4,7,0,8] <= 2
 R462: alpha[0,8] + beta[7,0] + gamma[5,7] - a[5,7,0,8] <= 2
 R463: alpha[0,8] + beta[7,0] + gamma[6,7] - a[6,7,0,8] <= 2
 R464: alpha[0,8] + beta[7,0] + gamma[7,7] - a[7,7,0,8] <= 2
 R465: alpha[0,8] + beta[8,0] + gamma[0,8] - a[0,8,0,8] <= 2
 R466: alpha[0,8] + beta[8,0] + gamma[1,8] - a[1,8,0,8] <= 2
 R467: alpha[0,8] + beta[8,0] + gamma[2,8] - a[2,8,0,8] <= 2
 R468: alpha[0,8] + beta[8,0] + gamma[3,8] - a[3,8,0,8] <= 2
 R469: alpha[0,8] + beta[8,0] + gamma[4,8] - a[4,8,0,8] <= 2
 R470: alpha[0,8] + beta[8,0] + gamma[5,8] - a[5,8,0,8] <= 2
 R471: alpha[0,8] + beta[8,0] + gamma[6,8] - a[6,8,0,8] <= 2
 R472: alpha[0,8] + beta[8,0] + gamma[7,8] - a[7,8,0,8] <= 2
 R473: alpha[0,8] + beta[8,0] + gamma[8,8] - a[8,8,0,8] <= 2
 R474: alpha[0,8] + beta[9,0] + gamma[0,9] - a[0,9,0,8] <= 2
 R475: alpha[0,8] + beta[9,0] + gamma[1,9] - a[1,9,0,8] <= 2
 R476: alpha[0,8] + beta[9,0] + gamma[2,9] - a[2,9,0,8] <= 2
 R477: alpha[0,8] + beta[9,0] + gamma[3,9] - a[3,9,0,8] <= 2
 R478: alpha[0,8] + beta[9,0] + gamma[4,9] - a[4,9,0,8] <= 2
 R479: alpha[0,8] + beta[9,0] + gamma[5,9] - a[5,9,0,8] <= 2
 R480: alpha[0,8] + beta[9,0] + gamma[6,9] - a[6,9,0,8] <= 2
 R481: alpha[0,8] + beta[9,0] + gamma[7,9] - a[7,9,0,8] <= 2
 R482: alpha[0,8] + beta[9,0] + gamma[8,9] - a[8,9,0,8] <= 2
 R483: alpha[0,8] + beta[9,0] + gamma[9,9] - a[9,9,0,8] <= 2
 R484: alpha[0,9] + beta[0,0] + gamma[0,0] - a[0,0,0,9] <= 2
 R485: alpha[0,9] + beta[1,0] + gamma[0,1] - a[0,1,0,9] <= 2
 R486: alpha[0,9] + beta[1,0] + gamma[1,1] - a[1,1,0,9] <= 2
 R487: alpha[0,9] + beta[2,0] + gamma[0,2] - a[0,2,0,9] <= 2
 R488: alpha[0,9] + beta[2,0] + gamma[1,2] - a[1,2,0,9] <= 2
 R489: alpha[0,9] + beta[2,0] + gamma[2,2] - a[2,2,0,9] <= 2
 R490: alpha[0,9] + beta[3,0] + gamma[0,3] - a[0,3,0,9] <= 2
 R491: alpha[0,9] + beta[3,0] + gamma[1,3] - a[1,3,0,9] <= 2
 R492: alpha[0,9] + beta[3,0] + gamma[2,3] - a[2,3,0,9] <= 2
 R493: alpha[0,9] + beta[3,0] + gamma[3,3] - a[3,3,0,9] <= 2
 R494: alpha[0,9] + beta[4,0] + gamma[0,4] - a[0,4,0,9] <= 2
 R495: alpha[0,9] + beta[4,0] + gamma[1,4] - a[1,4,0,9] <= 2
 R496: alpha[0,9] + beta[4,0] + gamma[2,4] - a[2,4,0,9] <= 2
 R497: alpha[0,9] + beta[4,0] + gamma[3,4] - a[3,4,0,9] <= 2
 R498: alpha[0,9] + beta[4,0] + gamma[4,4] - a[4,4,0,9] <= 2
 R499: alpha[0,9] + beta[5,0] + gamma[0,5] - a[0,5,0,9] <= 2
 R500: alpha[0,9] + beta[5,0] + gamma[1,5] - a[1,5,0,9] <= 2
 R501: alpha[0,9] + beta[5,0] + gamma[2,5] - a[2,5,0,9] <= 2
 R502: alpha[0,9] + beta[5,0] + gamma[3,5] - a[3,5,0,9] <= 2
 R503: alpha[0,9] + beta[5,0] + gamma[4,5] - a[4,5,0,9] <= 2
 R504: alpha[0,9] + beta[5,0] + gamma[5,5] - a[5,5,0,9] <= 2
 R505: alpha[0,9] + beta[6,0] + gamma[0,6] - a[0,6,0,9] <= 2
 R506: alpha[0,9] + beta[6,0] + gamma[1,6] - a[1,6,0,9] <= 2
 R507: alpha[0,9] + beta[6,0] + gamma[2,6] - a[2,6,0,9] <= 2
 R508: alpha[0,9] + beta[6,0] + gamma[3,6] - a[3,6,0,9] <= 2
 R509: alpha[0,9] + beta[6,0] + gamma[4,6] - a[4,6,0,9] <= 2
 R510: alpha[0,9] + beta[6,0] + gamma[5,6] - a[5,6,0,9] <= 2
 R511: alpha[0,9] + beta[6,0] + gamma[6,6] - a[6,6,0,9] <= 2
 R512: alpha[0,9] + beta[7,0] + gamma[0,7] - a[0,7,0,9] <= 2
 R513: alpha[0,9] + beta[7,0] + gamma[1,7] - a[1,7,0,9] <= 2
 R514: alpha[0,9] + beta[7,0] + gamma[2,7] - a[2,7,0,9] <= 2
 R515: alpha[0,9] + beta[7,0] + gamma[3,7] - a[3,7,0,9] <= 2
 R516: alpha[0,9] + beta[7,0] + gamma[4,7] - a[4,7,0,9] <= 2
 R517: alpha[0,9] + beta[7,0] + gamma[5,7] - a[5,7,0,9] <= 2
 R518: alpha[0,9] + beta[7,0] + gamma[6,7] - a[6,7,0,9] <= 2
 R519: alpha[0,9] + beta[7,0] + gamma[7,7] - a[7,7,0,9] <= 2
 R520: alpha[0,9] + beta[8,0] + gamma[0,8] - a[0,8,0,9] <= 2
 R521: alpha[0,9] + beta[8,0] + gamma[1,8] - a[1,8,0,9] <= 2
 R522: alpha[0,9] + beta[8,0] + gamma[2,8] - a[2,8,0,9] <= 2
 R523: alpha[0,9] + beta[8,0] + gamma[3,8] - a[3,8,0,9] <= 2
 R524: alpha[0,9] + beta[8,0] + gamma[4,8] - a[4,8,0,9] <= 2
 R525: alpha[0,9] + beta[8,0] + gamma[5,8] - a[5,8,0,9] <= 2
 R526: alpha[0,9] + beta[8,0] + gamma[6,8] - a[6,8,0,9] <= 2
 R527: alpha[0,9] + beta[8,0] + gamma[7,8] - a[7,8,0,9] <= 2
 R528: alpha[0,9] + beta[8,0] + gamma[8,8] - a[8,8,0,9] <= 2
 R529: alpha[0,9] + beta[9,0] + gamma[0,9] - a[0,9,0,9] <= 2
 R530: alpha[0,9] + beta[9,0] + gamma[1,9] - a[1,9,0,9] <= 2
 R531: alpha[0,9] + beta[9,0] + gamma[2,9] - a[2,9,0,9] <= 2
 R532: alpha[0,9] + beta[9,0] + gamma[3,9] - a[3,9,0,9] <= 2
 R533: alpha[0,9] + beta[9,0] + gamma[4,9] - a[4,9,0,9] <= 2
 R534: alpha[0,9] + beta[9,0] + gamma[5,9] - a[5,9,0,9] <= 2
 R535: alpha[0,9] + beta[9,0] + gamma[6,9] - a[6,9,0,9] <= 2
 R536: alpha[0,9] + beta[9,0] + gamma[7,9] - a[7,9,0,9] <= 2
 R537: alpha[0,9] + beta[9,0] + gamma[8,9] - a[8,9,0,9] <= 2
 R538: alpha[0,9] + beta[9,0] + gamma[9,9] - a[9,9,0,9] <= 2
 R539: alpha[1,1] + beta[0,1] + gamma[0,0] - a[0,0,1,1] <= 2
 R540: alpha[1,1] + beta[1,1] + gamma[0,1] - a[0,1,1,1] <= 2
 R541: alpha[1,1] + beta[1,1] + gamma[1,1] - a[1,1,1,1] <= 2
 R542: alpha[1,1] + beta[2,1] + gamma[0,2] - a[0,2,1,1] <= 2
 R543: alpha[1,1] + beta[2,1] + gamma[1,2] - a[1,2,1,1] <= 2
 R544: alpha[1,1] + beta[2,1] + gamma[2,2] - a[2,2,1,1] <= 2
 R545: alpha[1,1] + beta[3,1] + gamma[0,3] - a[0,3,1,1] <= 2
 R546: alpha[1,1] + beta[3,1] + gamma[1,3] - a[1,3,1,1] <= 2
 R547: alpha[1,1] + beta[3,1] + gamma[2,3] - a[2,3,1,1] <= 2
 R548: alpha[1,1] + beta[3,1] + gamma[3,3] - a[3,3,1,1] <= 2
 R549: alpha[1,1] + beta[4,1] + gamma[0,4] - a[0,4,1,1] <= 2
 R550: alpha[1,1] + beta[4,1] + gamma[1,4] - a[1,4,1,1] <= 2
 R551: alpha[1,1] + beta[4,1] + gamma[2,4] - a[2,4,1,1] <= 2
 R552: alpha[1,1] + beta[4,1] + gamma[3,4] - a[3,4,1,1] <= 2
 R553: alpha[1,1] + beta[4,1] + gamma[4,4] - a[4,4,1,1] <= 2
 R554: alpha[1,1] + beta[5,1] + gamma[0,5] - a[0,5,1,1] <= 2
 R555: alpha[1,1] + beta[5,1] + gamma[1,5] - a[1,5,1,1] <= 2
 R556: alpha[1,1] + beta[5,1] + gamma[2,5] - a[2,5,1,1] <= 2
 R557: alpha[1,1] + beta[5,1] + gamma[3,5] - a[3,5,1,1] <= 2
 R558: alpha[1,1] + beta[5,1] + gamma[4,5] - a[4,5,1,1] <= 2
 R559: alpha[1,1] + beta[5,1] + gamma[5,5] - a[5,5,1,1] <= 2
 R560: alpha[1,1] + beta[6,1] + gamma[0,6] - a[0,6,1,1] <= 2
 R561: alpha[1,1] + beta[6,1] + gamma[1,6] - a[1,6,1,1] <= 2
 R562: alpha[1,1] + beta[6,1] + gamma[2,6] - a[2,6,1,1] <= 2
 R563: alpha[1,1] + beta[6,1] + gamma[3,6] - a[3,6,1,1] <= 2
 R564: alpha[1,1] + beta[6,1] + gamma[4,6] - a[4,6,1,1] <= 2
 R565: alpha[1,1] + beta[6,1] + gamma[5,6] - a[5,6,1,1] <= 2
 R566: alpha[1,1] + beta[6,1] + gamma[6,6] - a[6,6,1,1] <= 2
 R567: alpha[1,1] + beta[7,1] + gamma[0,7] - a[0,7,1,1] <= 2
 R568: alpha[1,1] + beta[7,1] + gamma[1,7] - a[1,7,1,1] <= 2
 R569: alpha[1,1] + beta[7,1] + gamma[2,7] - a[2,7,1,1] <= 2
 R570: alpha[1,1] + beta[7,1] + gamma[3,7] - a[3,7,1,1] <= 2
 R571: alpha[1,1] + beta[7,1] + gamma[4,7] - a[4,7,1,1] <= 2
 R572: alpha[1,1] + beta[7,1] + gamma[5,7] - a[5,7,1,1] <= 2
 R573: alpha[1,1] + beta[7,1] + gamma[6,7] - a[6,7,1,1] <= 2
 R574: alpha[1,1] + beta[7,1] + gamma[7,7] - a[7,7,1,1] <= 2
 R575: alpha[1,1] + beta[8,1] + gamma[0,8] - a[0,8,1,1] <= 2
 R576: alpha[1,1] + beta[8,1] + gamma[1,8] - a[1,8,1,1] <= 2
 R577: alpha[1,1] + beta[8,1] + gamma[2,8] - a[2,8,1,1] <= 2
 R578: alpha[1,1] + beta[8,1] + gamma[3,8] - a[3,8,1,1] <= 2
 R579: alpha[1,1] + beta[8,1] + gamma[4,8] - a[4,8,1,1] <= 2
 R580: alpha[1,1] + beta[8,1] + gamma[5,8] - a[5,8,1,1] <= 2
 R581: alpha[1,1] + beta[8,1] + gamma[6,8] - a[6,8,1,1] <= 2
 R582: alpha[1,1] + beta[8,1] + gamma[7,8] - a[7,8,1,1] <= 2
 R583: alpha[1,1] + beta[8,1] + gamma[8,8] - a[8,8,1,1] <= 2
 R584: alpha[1,1] + beta[9,1] + gamma[0,9] - a[0,9,1,1] <= 2
 R585: alpha[1,1] + beta[9,1] + gamma[1,9] - a[1,9,1,1] <= 2
 R586: alpha[1,1] + beta[9,1] + gamma[2,9] - a[2,9,1,1] <= 2
 R587: alpha[1,1] + beta[9,1] + gamma[3,9] - a[3,9,1,1] <= 2
 R588: alpha[1,1] + beta[9,1] + gamma[4,9] - a[4,9,1,1] <= 2
 R589: alpha[1,1] + beta[9,1] + gamma[5,9] - a[5,9,1,1] <= 2
 R590: alpha[1,1] + beta[9,1] + gamma[6,9] - a[6,9,1,1] <= 2
 R591: alpha[1,1] + beta[9,1] + gamma[7,9] - a[7,9,1,1] <= 2
 R592: alpha[1,1] + beta[9,1] + gamma[8,9] - a[8,9,1,1] <= 2
 R593: alpha[1,1] + beta[9,1] + gamma[9,9] - a[9,9,1,1] <= 2
 R594: alpha[1,3] + beta[0,1] + gamma[0,0] - a[0,0,1,3] <= 2
 R595: alpha[1,3] + beta[1,1] + gamma[0,1] - a[0,1,1,3] <= 2
 R596: alpha[1,3] + beta[1,1] + gamma[1,1] - a[1,1,1,3] <= 2
 R597: alpha[1,3] + beta[2,1] + gamma[0,2] - a[0,2,1,3] <= 2
 R598: alpha[1,3] + beta[2,1] + gamma[1,2] - a[1,2,1,3] <= 2
 R599: alpha[1,3] + beta[2,1] + gamma[2,2] - a[2,2,1,3] <= 2
 R600: alpha[1,3] + beta[3,1] + gamma[0,3] - a[0,3,1,3] <= 2
 R601: alpha[1,3] + beta[3,1] + gamma[1,3] - a[1,3,1,3] <= 2
 R602: alpha[1,3] + beta[3,1] + gamma[2,3] - a[2,3,1,3] <= 2
 R603: alpha[1,3] + beta[3,1] + gamma[3,3] - a[3,3,1,3] <= 2
 R604: alpha[1,3] + beta[4,1] + gamma[0,4] - a[0,4,1,3] <= 2
 R605: alpha[1,3] + beta[4,1] + gamma[1,4] - a[1,4,1,3] <= 2
 R606: alpha[1,3] + beta[4,1] + gamma[2,4] - a[2,4,1,3] <= 2
 R607: alpha[1,3] + beta[4,1] + gamma[3,4] - a[3,4,1,3] <= 2
 R608: alpha[1,3] + beta[4,1] + gamma[4,4] - a[4,4,1,3] <= 2
 R609: alpha[1,3] + beta[5,1] + gamma[0,5] - a[0,5,1,3] <= 2
 R610: alpha[1,3] + beta[5,1] + gamma[1,5] - a[1,5,1,3] <= 2
 R611: alpha[1,3] + beta[5,1] + gamma[2,5] - a[2,5,1,3] <= 2
 R612: alpha[1,3] + beta[5,1] + gamma[3,5] - a[3,5,1,3] <= 2
 R613: alpha[1,3] + beta[5,1] + gamma[4,5] - a[4,5,1,3] <= 2
 R614: alpha[1,3] + beta[5,1] + gamma[5,5] - a[5,5,1,3] <= 2
 R615: alpha[1,3] + beta[6,1] + gamma[0,6] - a[0,6,1,3] <= 2
 R616: alpha[1,3] + beta[6,1] + gamma[1,6] - a[1,6,1,3] <= 2
 R617: alpha[1,3] + beta[6,1] + gamma[2,6] - a[2,6,1,3] <= 2
 R618: alpha[1,3] + beta[6,1] + gamma[3,6] - a[3,6,1,3] <= 2
 R619: alpha[1,3] + beta[6,1] + gamma[4,6] - a[4,6,1,3] <= 2
 R620: alpha[1,3] + beta[6,1] + gamma[5,6] - a[5,6,1,3] <= 2
 R621: alpha[1,3] + beta[6,1] + gamma[6,6] - a[6,6,1,3] <= 2
 R622: alpha[1,3] + beta[7,1] + gamma[0,7] - a[0,7,1,3] <= 2
 R623: alpha[1,3] + beta[7,1] + gamma[1,7] - a[1,7,1,3] <= 2
 R624: alpha[1,3] + beta[7,1] + gamma[2,7] - a[2,7,1,3] <= 2
 R625: alpha[1,3] + beta[7,1] + gamma[3,7] - a[3,7,1,3] <= 2
 R626: alpha[1,3] + beta[7,1] + gamma[4,7] - a[4,7,1,3] <= 2
 R627: alpha[1,3] + beta[7,1] + gamma[5,7] - a[5,7,1,3] <= 2
 R628: alpha[1,3] + beta[7,1] + gamma[6,7] - a[6,7,1,3] <= 2
 R629: alpha[1,3] + beta[7,1] + gamma[7,7] - a[7,7,1,3] <= 2
 R630: alpha[1,3] + beta[8,1] + gamma[0,8] - a[0,8,1,3] <= 2
 R631: alpha[1,3] + beta[8,1] + gamma[1,8] - a[1,8,1,3] <= 2
 R632: alpha[1,3] + beta[8,1] + gamma[2,8] - a[2,8,1,3] <= 2
 R633: alpha[1,3] + beta[8,1] + gamma[3,8] - a[3,8,1,3] <= 2
 R634: alpha[1,3] + beta[8,1] + gamma[4,8] - a[4,8,1,3] <= 2
 R635: alpha[1,3] + beta[8,1] + gamma[5,8] - a[5,8,1,3] <= 2
 R636: alpha[1,3] + beta[8,1] + gamma[6,8] - a[6,8,1,3] <= 2
 R637: alpha[1,3] + beta[8,1] + gamma[7,8] - a[7,8,1,3] <= 2
 R638: alpha[1,3] + beta[8,1] + gamma[8,8] - a[8,8,1,3] <= 2
 R639: alpha[1,3] + beta[9,1] + gamma[0,9] - a[0,9,1,3] <= 2
 R640: alpha[1,3] + beta[9,1] + gamma[1,9] - a[1,9,1,3] <= 2
 R641: alpha[1,3] + beta[9,1] + gamma[2,9] - a[2,9,1,3] <= 2
 R642: alpha[1,3] + beta[9,1] + gamma[3,9] - a[3,9,1,3] <= 2
 R643: alpha[1,3] + beta[9,1] + gamma[4,9] - a[4,9,1,3] <= 2
 R644: alpha[1,3] + beta[9,1] + gamma[5,9] - a[5,9,1,3] <= 2
 R645: alpha[1,3] + beta[9,1] + gamma[6,9] - a[6,9,1,3] <= 2
 R646: alpha[1,3] + beta[9,1] + gamma[7,9] - a[7,9,1,3] <= 2
 R647: alpha[1,3] + beta[9,1] + gamma[8,9] - a[8,9,1,3] <= 2
 R648: alpha[1,3] + beta[9,1] + gamma[9,9] - a[9,9,1,3] <= 2
 R649: alpha[1,7] + beta[0,1] + gamma[0,0] - a[0,0,1,7] <= 2
 R650: alpha[1,7] + beta[1,1] + gamma[0,1] - a[0,1,1,7] <= 2
 R651: alpha[1,7] + beta[1,1] + gamma[1,1] - a[1,1,1,7] <= 2
 R652: alpha[1,7] + beta[2,1] + gamma[0,2] - a[0,2,1,7] <= 2
 R653: alpha[1,7] + beta[2,1] + gamma[1,2] - a[1,2,1,7] <= 2
 R654: alpha[1,7] + beta[2,1] + gamma[2,2] - a[2,2,1,7] <= 2
 R655: alpha[1,7] + beta[3,1] + gamma[0,3] - a[0,3,1,7] <= 2
 R656: alpha[1,7] + beta[3,1] + gamma[1,3] - a[1,3,1,7] <= 2
 R657: alpha[1,7] + beta[3,1] + gamma[2,3] - a[2,3,1,7] <= 2
 R658: alpha[1,7] + beta[3,1] + gamma[3,3] - a[3,3,1,7] <= 2
 R659: alpha[1,7] + beta[4,1] + gamma[0,4] - a[0,4,1,7] <= 2
 R660: alpha[1,7] + beta[4,1] + gamma[1,4] - a[1,4,1,7] <= 2
 R661: alpha[1,7] + beta[4,1] + gamma[2,4] - a[2,4,1,7] <= 2
 R662: alpha[1,7] + beta[4,1] + gamma[3,4] - a[3,4,1,7] <= 2
 R663: alpha[1,7] + beta[4,1] + gamma[4,4] - a[4,4,1,7] <= 2
 R664: alpha[1,7] + beta[5,1] + gamma[0,5] - a[0,5,1,7] <= 2
 R665: alpha[1,7] + beta[5,1] + gamma[1,5] - a[1,5,1,7] <= 2
 R666: alpha[1,7] + beta[5,1] + gamma[2,5] - a[2,5,1,7] <= 2
 R667: alpha[1,7] + beta[5,1] + gamma[3,5] - a[3,5,1,7] <= 2
 R668: alpha[1,7] + beta[5,1] + gamma[4,5] - a[4,5,1,7] <= 2
 R669: alpha[1,7] + beta[5,1] + gamma[5,5] - a[5,5,1,7] <= 2
 R670: alpha[1,7] + beta[6,1] + gamma[0,6] - a[0,6,1,7] <= 2
 R671: alpha[1,7] + beta[6,1] + gamma[1,6] - a[1,6,1,7] <= 2
 R672: alpha[1,7] + beta[6,1] + gamma[2,6] - a[2,6,1,7] <= 2
 R673: alpha[1,7] + beta[6,1] + gamma[3,6] - a[3,6,1,7] <= 2
 R674: alpha[1,7] + beta[6,1] + gamma[4,6] - a[4,6,1,7] <= 2
 R675: alpha[1,7] + beta[6,1] + gamma[5,6] - a[5,6,1,7] <= 2
 R676: alpha[1,7] + beta[6,1] + gamma[6,6] - a[6,6,1,7] <= 2
 R677: alpha[1,7] + beta[7,1] + gamma[0,7] - a[0,7,1,7] <= 2
 R678: alpha[1,7] + beta[7,1] + gamma[1,7] - a[1,7,1,7] <= 2
 R679: alpha[1,7] + beta[7,1] + gamma[2,7] - a[2,7,1,7] <= 2
 R680: alpha[1,7] + beta[7,1] + gamma[3,7] - a[3,7,1,7] <= 2
 R681: alpha[1,7] + beta[7,1] + gamma[4,7] - a[4,7,1,7] <= 2
 R682: alpha[1,7] + beta[7,1] + gamma[5,7] - a[5,7,1,7] <= 2
 R683: alpha[1,7] + beta[7,1] + gamma[6,7] - a[6,7,1,7] <= 2
 R684: alpha[1,7] + beta[7,1] + gamma[7,7] - a[7,7,1,7] <= 2
 R685: alpha[1,7] + beta[8,1] + gamma[0,8] - a[0,8,1,7] <= 2
 R686: alpha[1,7] + beta[8,1] + gamma[1,8] - a[1,8,1,7] <= 2
 R687: alpha[1,7] + beta[8,1] + gamma[2,8] - a[2,8,1,7] <= 2
 R688: alpha[1,7] + beta[8,1] + gamma[3,8] - a[3,8,1,7] <= 2
 R689: alpha[1,7] + beta[8,1] + gamma[4,8] - a[4,8,1,7] <= 2
 R690: alpha[1,7] + beta[8,1] + gamma[5,8] - a[5,8,1,7] <= 2
 R691: alpha[1,7] + beta[8,1] + gamma[6,8] - a[6,8,1,7] <= 2
 R692: alpha[1,7] + beta[8,1] + gamma[7,8] - a[7,8,1,7] <= 2
 R693: alpha[1,7] + beta[8,1] + gamma[8,8] - a[8,8,1,7] <= 2
 R694: alpha[1,7] + beta[9,1] + gamma[0,9] - a[0,9,1,7] <= 2
 R695: alpha[1,7] + beta[9,1] + gamma[1,9] - a[1,9,1,7] <= 2
 R696: alpha[1,7] + beta[9,1] + gamma[2,9] - a[2,9,1,7] <= 2
 R697: alpha[1,7] + beta[9,1] + gamma[3,9] - a[3,9,1,7] <= 2
 R698: alpha[1,7] + beta[9,1] + gamma[4,9] - a[4,9,1,7] <= 2
 R699: alpha[1,7] + beta[9,1] + gamma[5,9] - a[5,9,1,7] <= 2
 R700: alpha[1,7] + beta[9,1] + gamma[6,9] - a[6,9,1,7] <= 2
 R701: alpha[1,7] + beta[9,1] + gamma[7,9] - a[7,9,1,7] <= 2
 R702: alpha[1,7] + beta[9,1] + gamma[8,9] - a[8,9,1,7] <= 2
 R703: alpha[1,7] + beta[9,1] + gamma[9,9] - a[9,9,1,7] <= 2
 R704: alpha[2,2] + beta[0,2] + gamma[0,0] - a[0,0,2,2] <= 2
 R705: alpha[2,2] + beta[1,2] + gamma[0,1] - a[0,1,2,2] <= 2
 R706: alpha[2,2] + beta[1,2] + gamma[1,1] - a[1,1,2,2] <= 2
 R707: alpha[2,2] + beta[2,2] + gamma[0,2] - a[0,2,2,2] <= 2
 R708: alpha[2,2] + beta[2,2] + gamma[1,2] - a[1,2,2,2] <= 2
 R709: alpha[2,2] + beta[2,2] + gamma[2,2] - a[2,2,2,2] <= 2
 R710: alpha[2,2] + beta[3,2] + gamma[0,3] - a[0,3,2,2] <= 2
 R711: alpha[2,2] + beta[3,2] + gamma[1,3] - a[1,3,2,2] <= 2
 R712: alpha[2,2] + beta[3,2] + gamma[2,3] - a[2,3,2,2] <= 2
 R713: alpha[2,2] + beta[3,2] + gamma[3,3] - a[3,3,2,2] <= 2
 R714: alpha[2,2] + beta[4,2] + gamma[0,4] - a[0,4,2,2] <= 2
 R715: alpha[2,2] + beta[4,2] + gamma[1,4] - a[1,4,2,2] <= 2
 R716: alpha[2,2] + beta[4,2] + gamma[2,4] - a[2,4,2,2] <= 2
 R717: alpha[2,2] + beta[4,2] + gamma[3,4] - a[3,4,2,2] <= 2
 R718: alpha[2,2] + beta[4,2] + gamma[4,4] - a[4,4,2,2] <= 2
 R719: alpha[2,2] + beta[5,2] + gamma[0,5] - a[0,5,2,2] <= 2
 R720: alpha[2,2] + beta[5,2] + gamma[1,5] - a[1,5,2,2] <= 2
 R721: alpha[2,2] + beta[5,2] + gamma[2,5] - a[2,5,2,2] <= 2
 R722: alpha[2,2] + beta[5,2] + gamma[3,5] - a[3,5,2,2] <= 2
 R723: alpha[2,2] + beta[5,2] + gamma[4,5] - a[4,5,2,2] <= 2
 R724: alpha[2,2] + beta[5,2] + gamma[5,5] - a[5,5,2,2] <= 2
 R725: alpha[2,2] + beta[6,2] + gamma[0,6] - a[0,6,2,2] <= 2
 R726: alpha[2,2] + beta[6,2] + gamma[1,6] - a[1,6,2,2] <= 2
 R727: alpha[2,2] + beta[6,2] + gamma[2,6] - a[2,6,2,2] <= 2
 R728: alpha[2,2] + beta[6,2] + gamma[3,6] - a[3,6,2,2] <= 2
 R729: alpha[2,2] + beta[6,2] + gamma[4,6] - a[4,6,2,2] <= 2
 R730: alpha[2,2] + beta[6,2] + gamma[5,6] - a[5,6,2,2] <= 2
 R731: alpha[2,2] + beta[6,2] + gamma[6,6] - a[6,6,2,2] <= 2
 R732: alpha[2,2] + beta[7,2] + gamma[0,7] - a[0,7,2,2] <= 2
 R733: alpha[2,2] + beta[7,2] + gamma[1,7] - a[1,7,2,2] <= 2
 R734: alpha[2,2] + beta[7,2] + gamma[2,7] - a[2,7,2,2] <= 2
 R735: alpha[2,2] + beta[7,2] + gamma[3,7] - a[3,7,2,2] <= 2
 R736: alpha[2,2] + beta[7,2] + gamma[4,7] - a[4,7,2,2] <= 2
 R737: alpha[2,2] + beta[7,2] + gamma[5,7] - a[5,7,2,2] <= 2
 R738: alpha[2,2] + beta[7,2] + gamma[6,7] - a[6,7,2,2] <= 2
 R739: alpha[2,2] + beta[7,2] + gamma[7,7] - a[7,7,2,2] <= 2
 R740: alpha[2,2] + beta[8,2] + gamma[0,8] - a[0,8,2,2] <= 2
 R741: alpha[2,2] + beta[8,2] + gamma[1,8] - a[1,8,2,2] <= 2
 R742: alpha[2,2] + beta[8,2] + gamma[2,8] - a[2,8,2,2] <= 2
 R743: alpha[2,2] + beta[8,2] + gamma[3,8] - a[3,8,2,2] <= 2
 R744: alpha[2,2] + beta[8,2] + gamma[4,8] - a[4,8,2,2] <= 2
 R745: alpha[2,2] + beta[8,2] + gamma[5,8] - a[5,8,2,2] <= 2
 R746: alpha[2,2] + beta[8,2] + gamma[6,8] - a[6,8,2,2] <= 2
 R747: alpha[2,2] + beta[8,2] + gamma[7,8] - a[7,8,2,2] <= 2
 R748: alpha[2,2] + beta[8,2] + gamma[8,8] - a[8,8,2,2] <= 2
 R749: alpha[2,2] + beta[9,2] + gamma[0,9] - a[0,9,2,2] <= 2
 R750: alpha[2,2] + beta[9,2] + gamma[1,9] - a[1,9,2,2] <= 2
 R751: alpha[2,2] + beta[9,2] + gamma[2,9] - a[2,9,2,2] <= 2
 R752: alpha[2,2] + beta[9,2] + gamma[3,9] - a[3,9,2,2] <= 2
 R753: alpha[2,2] + beta[9,2] + gamma[4,9] - a[4,9,2,2] <= 2
 R754: alpha[2,2] + beta[9,2] + gamma[5,9] - a[5,9,2,2] <= 2
 R755: alpha[2,2] + beta[9,2] + gamma[6,9] - a[6,9,2,2] <= 2
 R756: alpha[2,2] + beta[9,2] + gamma[7,9] - a[7,9,2,2] <= 2
 R757: alpha[2,2] + beta[9,2] + gamma[8,9] - a[8,9,2,2] <= 2
 R758: alpha[2,2] + beta[9,2] + gamma[9,9] - a[9,9,2,2] <= 2
 R759: alpha[2,6] + beta[0,2] + gamma[0,0] - a[0,0,2,6] <= 2
 R760: alpha[2,6] + beta[1,2] + gamma[0,1] - a[0,1,2,6] <= 2
 R761: alpha[2,6] + beta[1,2] + gamma[1,1] - a[1,1,2,6] <= 2
 R762: alpha[2,6] + beta[2,2] + gamma[0,2] - a[0,2,2,6] <= 2
 R763: alpha[2,6] + beta[2,2] + gamma[1,2] - a[1,2,2,6] <= 2
 R764: alpha[2,6] + beta[2,2] + gamma[2,2] - a[2,2,2,6] <= 2
 R765: alpha[2,6] + beta[3,2] + gamma[0,3] - a[0,3,2,6] <= 2
 R766: alpha[2,6] + beta[3,2] + gamma[1,3] - a[1,3,2,6] <= 2
 R767: alpha[2,6] + beta[3,2] + gamma[2,3] - a[2,3,2,6] <= 2
 R768: alpha[2,6] + beta[3,2] + gamma[3,3] - a[3,3,2,6] <= 2
 R769: alpha[2,6] + beta[4,2] + gamma[0,4] - a[0,4,2,6] <= 2
 R770: alpha[2,6] + beta[4,2] + gamma[1,4] - a[1,4,2,6] <= 2
 R771: alpha[2,6] + beta[4,2] + gamma[2,4] - a[2,4,2,6] <= 2
 R772: alpha[2,6] + beta[4,2] + gamma[3,4] - a[3,4,2,6] <= 2
 R773: alpha[2,6] + beta[4,2] + gamma[4,4] - a[4,4,2,6] <= 2
 R774: alpha[2,6] + beta[5,2] + gamma[0,5] - a[0,5,2,6] <= 2
 R775: alpha[2,6] + beta[5,2] + gamma[1,5] - a[1,5,2,6] <= 2
 R776: alpha[2,6] + beta[5,2] + gamma[2,5] - a[2,5,2,6] <= 2
 R777: alpha[2,6] + beta[5,2] + gamma[3,5] - a[3,5,2,6] <= 2
 R778: alpha[2,6] + beta[5,2] + gamma[4,5] - a[4,5,2,6] <= 2
 R779: alpha[2,6] + beta[5,2] + gamma[5,5] - a[5,5,2,6] <= 2
 R780: alpha[2,6] + beta[6,2] + gamma[0,6] - a[0,6,2,6] <= 2
 R781: alpha[2,6] + beta[6,2] + gamma[1,6] - a[1,6,2,6] <= 2
 R782: alpha[2,6] + beta[6,2] + gamma[2,6] - a[2,6,2,6] <= 2
 R783: alpha[2,6] + beta[6,2] + gamma[3,6] - a[3,6,2,6] <= 2
 R784: alpha[2,6] + beta[6,2] + gamma[4,6] - a[4,6,2,6] <= 2
 R785: alpha[2,6] + beta[6,2] + gamma[5,6] - a[5,6,2,6] <= 2
 R786: alpha[2,6] + beta[6,2] + gamma[6,6] - a[6,6,2,6] <= 2
 R787: alpha[2,6] + beta[7,2] + gamma[0,7] - a[0,7,2,6] <= 2
 R788: alpha[2,6] + beta[7,2] + gamma[1,7] - a[1,7,2,6] <= 2
 R789: alpha[2,6] + beta[7,2] + gamma[2,7] - a[2,7,2,6] <= 2
 R790: alpha[2,6] + beta[7,2] + gamma[3,7] - a[3,7,2,6] <= 2
 R791: alpha[2,6] + beta[7,2] + gamma[4,7] - a[4,7,2,6] <= 2
 R792: alpha[2,6] + beta[7,2] + gamma[5,7] - a[5,7,2,6] <= 2
 R793: alpha[2,6] + beta[7,2] + gamma[6,7] - a[6,7,2,6] <= 2
 R794: alpha[2,6] + beta[7,2] + gamma[7,7] - a[7,7,2,6] <= 2
 R795: alpha[2,6] + beta[8,2] + gamma[0,8] - a[0,8,2,6] <= 2
 R796: alpha[2,6] + beta[8,2] + gamma[1,8] - a[1,8,2,6] <= 2
 R797: alpha[2,6] + beta[8,2] + gamma[2,8] - a[2,8,2,6] <= 2
 R798: alpha[2,6] + beta[8,2] + gamma[3,8] - a[3,8,2,6] <= 2
 R799: alpha[2,6] + beta[8,2] + gamma[4,8] - a[4,8,2,6] <= 2
 R800: alpha[2,6] + beta[8,2] + gamma[5,8] - a[5,8,2,6] <= 2
 R801: alpha[2,6] + beta[8,2] + gamma[6,8] - a[6,8,2,6] <= 2
 R802: alpha[2,6] + beta[8,2] + gamma[7,8] - a[7,8,2,6] <= 2
 R803: alpha[2,6] + beta[8,2] + gamma[8,8] - a[8,8,2,6] <= 2
 R804: alpha[2,6] + beta[9,2] + gamma[0,9] - a[0,9,2,6] <= 2
 R805: alpha[2,6] + beta[9,2] + gamma[1,9] - a[1,9,2,6] <= 2
 R806: alpha[2,6] + beta[9,2] + gamma[2,9] - a[2,9,2,6] <= 2
 R807: alpha[2,6] + beta[9,2] + gamma[3,9] - a[3,9,2,6] <= 2
 R808: alpha[2,6] + beta[9,2] + gamma[4,9] - a[4,9,2,6] <= 2
 R809: alpha[2,6] + beta[9,2] + gamma[5,9] - a[5,9,2,6] <= 2
 R810: alpha[2,6] + beta[9,2] + gamma[6,9] - a[6,9,2,6] <= 2
 R811: alpha[2,6] + beta[9,2] + gamma[7,9] - a[7,9,2,6] <= 2
 R812: alpha[2,6] + beta[9,2] + gamma[8,9] - a[8,9,2,6] <= 2
 R813: alpha[2,6] + beta[9,2] + gamma[9,9] - a[9,9,2,6] <= 2
 R814: alpha[2,8] + beta[0,2] + gamma[0,0] - a[0,0,2,8] <= 2
 R815: alpha[2,8] + beta[1,2] + gamma[0,1] - a[0,1,2,8] <= 2
 R816: alpha[2,8] + beta[1,2] + gamma[1,1] - a[1,1,2,8] <= 2
 R817: alpha[2,8] + beta[2,2] + gamma[0,2] - a[0,2,2,8] <= 2
 R818: alpha[2,8] + beta[2,2] + gamma[1,2] - a[1,2,2,8] <= 2
 R819: alpha[2,8] + beta[2,2] + gamma[2,2] - a[2,2,2,8] <= 2
 R820: alpha[2,8] + beta[3,2] + gamma[0,3] - a[0,3,2,8] <= 2
 R821: alpha[2,8] + beta[3,2] + gamma[1,3] - a[1,3,2,8] <= 2
 R822: alpha[2,8] + beta[3,2] + gamma[2,3] - a[2,3,2,8] <= 2
 R823: alpha[2,8] + beta[3,2] + gamma[3,3] - a[3,3,2,8] <= 2
 R824: alpha[2,8] + beta[4,2] + gamma[0,4] - a[0,4,2,8] <= 2
 R825: alpha[2,8] + beta[4,2] + gamma[1,4] - a[1,4,2,8] <= 2
 R826: alpha[2,8] + beta[4,2] + gamma[2,4] - a[2,4,2,8] <= 2
 R827: alpha[2,8] + beta[4,2] + gamma[3,4] - a[3,4,2,8] <= 2
 R828: alpha[2,8] + beta[4,2] + gamma[4,4] - a[4,4,2,8] <= 2
 R829: alpha[2,8] + beta[5,2] + gamma[0,5] - a[0,5,2,8] <= 2
 R830: alpha[2,8] + beta[5,2] + gamma[1,5] - a[1,5,2,8] <= 2
 R831: alpha[2,8] + beta[5,2] + gamma[2,5] - a[2,5,2,8] <= 2
 R832: alpha[2,8] + beta[5,2] + gamma[3,5] - a[3,5,2,8] <= 2
 R833: alpha[2,8] + beta[5,2] + gamma[4,5] - a[4,5,2,8] <= 2
 R834: alpha[2,8] + beta[5,2] + gamma[5,5] - a[5,5,2,8] <= 2
 R835: alpha[2,8] + beta[6,2] + gamma[0,6] - a[0,6,2,8] <= 2
 R836: alpha[2,8] + beta[6,2] + gamma[1,6] - a[1,6,2,8] <= 2
 R837: alpha[2,8] + beta[6,2] + gamma[2,6] - a[2,6,2,8] <= 2
 R838: alpha[2,8] + beta[6,2] + gamma[3,6] - a[3,6,2,8] <= 2
 R839: alpha[2,8] + beta[6,2] + gamma[4,6] - a[4,6,2,8] <= 2
 R840: alpha[2,8] + beta[6,2] + gamma[5,6] - a[5,6,2,8] <= 2
 R841: alpha[2,8] + beta[6,2] + gamma[6,6] - a[6,6,2,8] <= 2
 R842: alpha[2,8] + beta[7,2] + gamma[0,7] - a[0,7,2,8] <= 2
 R843: alpha[2,8] + beta[7,2] + gamma[1,7] - a[1,7,2,8] <= 2
 R844: alpha[2,8] + beta[7,2] + gamma[2,7] - a[2,7,2,8] <= 2
 R845: alpha[2,8] + beta[7,2] + gamma[3,7] - a[3,7,2,8] <= 2
 R846: alpha[2,8] + beta[7,2] + gamma[4,7] - a[4,7,2,8] <= 2
 R847: alpha[2,8] + beta[7,2] + gamma[5,7] - a[5,7,2,8] <= 2
 R848: alpha[2,8] + beta[7,2] + gamma[6,7] - a[6,7,2,8] <= 2
 R849: alpha[2,8] + beta[7,2] + gamma[7,7] - a[7,7,2,8] <= 2
 R850: alpha[2,8] + beta[8,2] + gamma[0,8] - a[0,8,2,8] <= 2
 R851: alpha[2,8] + beta[8,2] + gamma[1,8] - a[1,8,2,8] <= 2
 R852: alpha[2,8] + beta[8,2] + gamma[2,8] - a[2,8,2,8] <= 2
 R853: alpha[2,8] + beta[8,2] + gamma[3,8] - a[3,8,2,8] <= 2
 R854: alpha[2,8] + beta[8,2] + gamma[4,8] - a[4,8,2,8] <= 2
 R855: alpha[2,8] + beta[8,2] + gamma[5,8] - a[5,8,2,8] <= 2
 R856: alpha[2,8] + beta[8,2] + gamma[6,8] - a[6,8,2,8] <= 2
 R857: alpha[2,8] + beta[8,2] + gamma[7,8] - a[7,8,2,8] <= 2
 R858: alpha[2,8] + beta[8,2] + gamma[8,8] - a[8,8,2,8] <= 2
 R859: alpha[2,8] + beta[9,2] + gamma[0,9] - a[0,9,2,8] <= 2
 R860: alpha[2,8] + beta[9,2] + gamma[1,9] - a[1,9,2,8] <= 2
 R861: alpha[2,8] + beta[9,2] + gamma[2,9] - a[2,9,2,8] <= 2
 R862: alpha[2,8] + beta[9,2] + gamma[3,9] - a[3,9,2,8] <= 2
 R863: alpha[2,8] + beta[9,2] + gamma[4,9] - a[4,9,2,8] <= 2
 R864: alpha[2,8] + beta[9,2] + gamma[5,9] - a[5,9,2,8] <= 2
 R865: alpha[2,8] + beta[9,2] + gamma[6,9] - a[6,9,2,8] <= 2
 R866: alpha[2,8] + beta[9,2] + gamma[7,9] - a[7,9,2,8] <= 2
 R867: alpha[2,8] + beta[9,2] + gamma[8,9] - a[8,9,2,8] <= 2
 R868: alpha[2,8] + beta[9,2] + gamma[9,9] - a[9,9,2,8] <= 2
 R869: alpha[2,9] + beta[0,2] + gamma[0,0] - a[0,0,2,9] <= 2
 R870: alpha[2,9] + beta[1,2] + gamma[0,1] - a[0,1,2,9] <= 2
 R871: alpha[2,9] + beta[1,2] + gamma[1,1] - a[1,1,2,9] <= 2
 R872: alpha[2,9] + beta[2,2] + gamma[0,2] - a[0,2,2,9] <= 2
 R873: alpha[2,9] + beta[2,2] + gamma[1,2] - a[1,2,2,9] <= 2
 R874: alpha[2,9] + beta[2,2] + gamma[2,2] - a[2,2,2,9] <= 2
 R875: alpha[2,9] + beta[3,2] + gamma[0,3] - a[0,3,2,9] <= 2
 R876: alpha[2,9] + beta[3,2] + gamma[1,3] - a[1,3,2,9] <= 2
 R877: alpha[2,9] + beta[3,2] + gamma[2,3] - a[2,3,2,9] <= 2
 R878: alpha[2,9] + beta[3,2] + gamma[3,3] - a[3,3,2,9] <= 2
 R879: alpha[2,9] + beta[4,2] + gamma[0,4] - a[0,4,2,9] <= 2
 R880: alpha[2,9] + beta[4,2] + gamma[1,4] - a[1,4,2,9] <= 2
 R881: alpha[2,9] + beta[4,2] + gamma[2,4] - a[2,4,2,9] <= 2
 R882: alpha[2,9] + beta[4,2] + gamma[3,4] - a[3,4,2,9] <= 2
 R883: alpha[2,9] + beta[4,2] + gamma[4,4] - a[4,4,2,9] <= 2
 R884: alpha[2,9] + beta[5,2] + gamma[0,5] - a[0,5,2,9] <= 2
 R885: alpha[2,9] + beta[5,2] + gamma[1,5] - a[1,5,2,9] <= 2
 R886: alpha[2,9] + beta[5,2] + gamma[2,5] - a[2,5,2,9] <= 2
 R887: alpha[2,9] + beta[5,2] + gamma[3,5] - a[3,5,2,9] <= 2
 R888: alpha[2,9] + beta[5,2] + gamma[4,5] - a[4,5,2,9] <= 2
 R889: alpha[2,9] + beta[5,2] + gamma[5,5] - a[5,5,2,9] <= 2
 R890: alpha[2,9] + beta[6,2] + gamma[0,6] - a[0,6,2,9] <= 2
 R891: alpha[2,9] + beta[6,2] + gamma[1,6] - a[1,6,2,9] <= 2
 R892: alpha[2,9] + beta[6,2] + gamma[2,6] - a[2,6,2,9] <= 2
 R893: alpha[2,9] + beta[6,2] + gamma[3,6] - a[3,6,2,9] <= 2
 R894: alpha[2,9] + beta[6,2] + gamma[4,6] - a[4,6,2,9] <= 2
 R895: alpha[2,9] + beta[6,2] + gamma[5,6] - a[5,6,2,9] <= 2
 R896: alpha[2,9] + beta[6,2] + gamma[6,6] - a[6,6,2,9] <= 2
 R897: alpha[2,9] + beta[7,2] + gamma[0,7] - a[0,7,2,9] <= 2
 R898: alpha[2,9] + beta[7,2] + gamma[1,7] - a[1,7,2,9] <= 2
 R899: alpha[2,9] + beta[7,2] + gamma[2,7] - a[2,7,2,9] <= 2
 R900: alpha[2,9] + beta[7,2] + gamma[3,7] - a[3,7,2,9] <= 2
 R901: alpha[2,9] + beta[7,2] + gamma[4,7] - a[4,7,2,9] <= 2
 R902: alpha[2,9] + beta[7,2] + gamma[5,7] - a[5,7,2,9] <= 2
 R903: alpha[2,9] + beta[7,2] + gamma[6,7] - a[6,7,2,9] <= 2
 R904: alpha[2,9] + beta[7,2] + gamma[7,7] - a[7,7,2,9] <= 2
 R905: alpha[2,9] + beta[8,2] + gamma[0,8] - a[0,8,2,9] <= 2
 R906: alpha[2,9] + beta[8,2] + gamma[1,8] - a[1,8,2,9] <= 2
 R907: alpha[2,9] + beta[8,2] + gamma[2,8] - a[2,8,2,9] <= 2
 R908: alpha[2,9] + beta[8,2] + gamma[3,8] - a[3,8,2,9] <= 2
 R909: alpha[2,9] + beta[8,2] + gamma[4,8] - a[4,8,2,9] <= 2
 R910: alpha[2,9] + beta[8,2] + gamma[5,8] - a[5,8,2,9] <= 2
 R911: alpha[2,9] + beta[8,2] + gamma[6,8] - a[6,8,2,9] <= 2
 R912: alpha[2,9] + beta[8,2] + gamma[7,8] - a[7,8,2,9] <= 2
 R913: alpha[2,9] + beta[8,2] + gamma[8,8] - a[8,8,2,9] <= 2
 R914: alpha[2,9] + beta[9,2] + gamma[0,9] - a[0,9,2,9] <= 2
 R915: alpha[2,9] + beta[9,2] + gamma[1,9] - a[1,9,2,9] <= 2
 R916: alpha[2,9] + beta[9,2] + gamma[2,9] - a[2,9,2,9] <= 2
 R917: alpha[2,9] + beta[9,2] + gamma[3,9] - a[3,9,2,9] <= 2
 R918: alpha[2,9] + beta[9,2] + gamma[4,9] - a[4,9,2,9] <= 2
 R919: alpha[2,9] + beta[9,2] + gamma[5,9] - a[5,9,2,9] <= 2
 R920: alpha[2,9] + beta[9,2] + gamma[6,9] - a[6,9,2,9] <= 2
 R921: alpha[2,9] + beta[9,2] + gamma[7,9] - a[7,9,2,9] <= 2
 R922: alpha[2,9] + beta[9,2] + gamma[8,9] - a[8,9,2,9] <= 2
 R923: alpha[2,9] + beta[9,2] + gamma[9,9] - a[9,9,2,9] <= 2
 R924: alpha[3,3] + beta[0,3] + gamma[0,0] - a[0,0,3,3] <= 2
 R925: alpha[3,3] + beta[1,3] + gamma[0,1] - a[0,1,3,3] <= 2
 R926: alpha[3,3] + beta[1,3] + gamma[1,1] - a[1,1,3,3] <= 2
 R927: alpha[3,3] + beta[2,3] + gamma[0,2] - a[0,2,3,3] <= 2
 R928: alpha[3,3] + beta[2,3] + gamma[1,2] - a[1,2,3,3] <= 2
 R929: alpha[3,3] + beta[2,3] + gamma[2,2] - a[2,2,3,3] <= 2
 R930: alpha[3,3] + beta[3,3] + gamma[0,3] - a[0,3,3,3] <= 2
 R931: alpha[3,3] + beta[3,3] + gamma[1,3] - a[1,3,3,3] <= 2
 R932: alpha[3,3] + beta[3,3] + gamma[2,3] - a[2,3,3,3] <= 2
 R933: alpha[3,3] + beta[3,3] + gamma[3,3] - a[3,3,3,3] <= 2
 R934: alpha[3,3] + beta[4,3] + gamma[0,4] - a[0,4,3,3] <= 2
 R935: alpha[3,3] + beta[4,3] + gamma[1,4] - a[1,4,3,3] <= 2
 R936: alpha[3,3] + beta[4,3] + gamma[2,4] - a[2,4,3,3] <= 2
 R937: alpha[3,3] + beta[4,3] + gamma[3,4] - a[3,4,3,3] <= 2
 R938: alpha[3,3] + beta[4,3] + gamma[4,4] - a[4,4,3,3] <= 2
 R939: alpha[3,3] + beta[5,3] + gamma[0,5] - a[0,5,3,3] <= 2
 R940: alpha[3,3] + beta[5,3] + gamma[1,5] - a[1,5,3,3] <= 2
 R941: alpha[3,3] + beta[5,3] + gamma[2,5] - a[2,5,3,3] <= 2
 R942: alpha[3,3] + beta[5,3] + gamma[3,5] - a[3,5,3,3] <= 2
 R943: alpha[3,3] + beta[5,3] + gamma[4,5] - a[4,5,3,3] <= 2
 R944: alpha[3,3] + beta[5,3] + gamma[5,5] - a[5,5,3,3] <= 2
 R945: alpha[3,3] + beta[6,3] + gamma[0,6] - a[0,6,3,3] <= 2
 R946: alpha[3,3] + beta[6,3] + gamma[1,6] - a[1,6,3,3] <= 2
 R947: alpha[3,3] + beta[6,3] + gamma[2,6] - a[2,6,3,3] <= 2
 R948: alpha[3,3] + beta[6,3] + gamma[3,6] - a[3,6,3,3] <= 2
 R949: alpha[3,3] + beta[6,3] + gamma[4,6] - a[4,6,3,3] <= 2
 R950: alpha[3,3] + beta[6,3] + gamma[5,6] - a[5,6,3,3] <= 2
 R951: alpha[3,3] + beta[6,3] + gamma[6,6] - a[6,6,3,3] <= 2
 R952: alpha[3,3] + beta[7,3] + gamma[0,7] - a[0,7,3,3] <= 2
 R953: alpha[3,3] + beta[7,3] + gamma[1,7] - a[1,7,3,3] <= 2
 R954: alpha[3,3] + beta[7,3] + gamma[2,7] - a[2,7,3,3] <= 2
 R955: alpha[3,3] + beta[7,3] + gamma[3,7] - a[3,7,3,3] <= 2
 R956: alpha[3,3] + beta[7,3] + gamma[4,7] - a[4,7,3,3] <= 2
 R957: alpha[3,3] + beta[7,3] + gamma[5,7] - a[5,7,3,3] <= 2
 R958: alpha[3,3] + beta[7,3] + gamma[6,7] - a[6,7,3,3] <= 2
 R959: alpha[3,3] + beta[7,3] + gamma[7,7] - a[7,7,3,3] <= 2
 R960: alpha[3,3] + beta[8,3] + gamma[0,8] - a[0,8,3,3] <= 2
 R961: alpha[3,3] + beta[8,3] + gamma[1,8] - a[1,8,3,3] <= 2
 R962: alpha[3,3] + beta[8,3] + gamma[2,8] - a[2,8,3,3] <= 2
 R963: alpha[3,3] + beta[8,3] + gamma[3,8] - a[3,8,3,3] <= 2
 R964: alpha[3,3] + beta[8,3] + gamma[4,8] - a[4,8,3,3] <= 2
 R965: alpha[3,3] + beta[8,3] + gamma[5,8] - a[5,8,3,3] <= 2
 R966: alpha[3,3] + beta[8,3] + gamma[6,8] - a[6,8,3,3] <= 2
 R967: alpha[3,3] + beta[8,3] + gamma[7,8] - a[7,8,3,3] <= 2
 R968: alpha[3,3] + beta[8,3] + gamma[8,8] - a[8,8,3,3] <= 2
 R969: alpha[3,3] + beta[9,3] + gamma[0,9] - a[0,9,3,3] <= 2
 R970: alpha[3,3] + beta[9,3] + gamma[1,9] - a[1,9,3,3] <= 2
 R971: alpha[3,3] + beta[9,3] + gamma[2,9] - a[2,9,3,3] <= 2
 R972: alpha[3,3] + beta[9,3] + gamma[3,9] - a[3,9,3,3] <= 2
 R973: alpha[3,3] + beta[9,3] + gamma[4,9] - a[4,9,3,3] <= 2
 R974: alpha[3,3] + beta[9,3] + gamma[5,9] - a[5,9,3,3] <= 2
 R975: alpha[3,3] + beta[9,3] + gamma[6,9] - a[6,9,3,3] <= 2
 R976: alpha[3,3] + beta[9,3] + gamma[7,9] - a[7,9,3,3] <= 2
 R977: alpha[3,3] + beta[9,3] + gamma[8,9] - a[8,9,3,3] <= 2
 R978: alpha[3,3] + beta[9,3] + gamma[9,9] - a[9,9,3,3] <= 2
 R979: alpha[3,7] + beta[0,3] + gamma[0,0] - a[0,0,3,7] <= 2
 R980: alpha[3,7] + beta[1,3] + gamma[0,1] - a[0,1,3,7] <= 2
 R981: alpha[3,7] + beta[1,3] + gamma[1,1] - a[1,1,3,7] <= 2
 R982: alpha[3,7] + beta[2,3] + gamma[0,2] - a[0,2,3,7] <= 2
 R983: alpha[3,7] + beta[2,3] + gamma[1,2] - a[1,2,3,7] <= 2
 R984: alpha[3,7] + beta[2,3] + gamma[2,2] - a[2,2,3,7] <= 2
 R985: alpha[3,7] + beta[3,3] + gamma[0,3] - a[0,3,3,7] <= 2
 R986: alpha[3,7] + beta[3,3] + gamma[1,3] - a[1,3,3,7] <= 2
 R987: alpha[3,7] + beta[3,3] + gamma[2,3] - a[2,3,3,7] <= 2
 R988: alpha[3,7] + beta[3,3] + gamma[3,3] - a[3,3,3,7] <= 2
 R989: alpha[3,7] + beta[4,3] + gamma[0,4] - a[0,4,3,7] <= 2
 R990: alpha[3,7] + beta[4,3] + gamma[1,4] - a[1,4,3,7] <= 2
 R991: alpha[3,7] + beta[4,3] + gamma[2,4] - a[2,4,3,7] <= 2
 R992: alpha[3,7] + beta[4,3] + gamma[3,4] - a[3,4,3,7] <= 2
 R993: alpha[3,7] + beta[4,3] + gamma[4,4] - a[4,4,3,7] <= 2
 R994: alpha[3,7] + beta[5,3] + gamma[0,5] - a[0,5,3,7] <= 2
 R995: alpha[3,7] + beta[5,3] + gamma[1,5] - a[1,5,3,7] <= 2
 R996: alpha[3,7] + beta[5,3] + gamma[2,5] - a[2,5,3,7] <= 2
 R997: alpha[3,7] + beta[5,3] + gamma[3,5] - a[3,5,3,7] <= 2
 R998: alpha[3,7] + beta[5,3] + gamma[4,5] - a[4,5,3,7] <= 2
 R999: alpha[3,7] + beta[5,3] + gamma[5,5] - a[5,5,3,7] <= 2
 R1000: alpha[3,7] + beta[6,3] + gamma[0,6] - a[0,6,3,7] <= 2
 R1001: alpha[3,7] + beta[6,3] + gamma[1,6] - a[1,6,3,7] <= 2
 R1002: alpha[3,7] + beta[6,3] + gamma[2,6] - a[2,6,3,7] <= 2
 R1003: alpha[3,7] + beta[6,3] + gamma[3,6] - a[3,6,3,7] <= 2
 R1004: alpha[3,7] + beta[6,3] + gamma[4,6] - a[4,6,3,7] <= 2
 R1005: alpha[3,7] + beta[6,3] + gamma[5,6] - a[5,6,3,7] <= 2
 R1006: alpha[3,7] + beta[6,3] + gamma[6,6] - a[6,6,3,7] <= 2
 R1007: alpha[3,7] + beta[7,3] + gamma[0,7] - a[0,7,3,7] <= 2
 R1008: alpha[3,7] + beta[7,3] + gamma[1,7] - a[1,7,3,7] <= 2
 R1009: alpha[3,7] + beta[7,3] + gamma[2,7] - a[2,7,3,7] <= 2
 R1010: alpha[3,7] + beta[7,3] + gamma[3,7] - a[3,7,3,7] <= 2
 R1011: alpha[3,7] + beta[7,3] + gamma[4,7] - a[4,7,3,7] <= 2
 R1012: alpha[3,7] + beta[7,3] + gamma[5,7] - a[5,7,3,7] <= 2
 R1013: alpha[3,7] + beta[7,3] + gamma[6,7] - a[6,7,3,7] <= 2
 R1014: alpha[3,7] + beta[7,3] + gamma[7,7] - a[7,7,3,7] <= 2
 R1015: alpha[3,7] + beta[8,3] + gamma[0,8] - a[0,8,3,7] <= 2
 R1016: alpha[3,7] + beta[8,3] + gamma[1,8] - a[1,8,3,7] <= 2
 R1017: alpha[3,7] + beta[8,3] + gamma[2,8] - a[2,8,3,7] <= 2
 R1018: alpha[3,7] + beta[8,3] + gamma[3,8] - a[3,8,3,7] <= 2
 R1019: alpha[3,7] + beta[8,3] + gamma[4,8] - a[4,8,3,7] <= 2
 R1020: alpha[3,7] + beta[8,3] + gamma[5,8] - a[5,8,3,7] <= 2
 R1021: alpha[3,7] + beta[8,3] + gamma[6,8] - a[6,8,3,7] <= 2
 R1022: alpha[3,7] + beta[8,3] + gamma[7,8] - a[7,8,3,7] <= 2
 R1023: alpha[3,7] + beta[8,3] + gamma[8,8] - a[8,8,3,7] <= 2
 R1024: alpha[3,7] + beta[9,3] + gamma[0,9] - a[0,9,3,7] <= 2
 R1025: alpha[3,7] + beta[9,3] + gamma[1,9] - a[1,9,3,7] <= 2
 R1026: alpha[3,7] + beta[9,3] + gamma[2,9] - a[2,9,3,7] <= 2
 R1027: alpha[3,7] + beta[9,3] + gamma[3,9] - a[3,9,3,7] <= 2
 R1028: alpha[3,7] + beta[9,3] + gamma[4,9] - a[4,9,3,7] <= 2
 R1029: alpha[3,7] + beta[9,3] + gamma[5,9] - a[5,9,3,7] <= 2
 R1030: alpha[3,7] + beta[9,3] + gamma[6,9] - a[6,9,3,7] <= 2
 R1031: alpha[3,7] + beta[9,3] + gamma[7,9] - a[7,9,3,7] <= 2
 R1032: alpha[3,7] + beta[9,3] + gamma[8,9] - a[8,9,3,7] <= 2
 R1033: alpha[3,7] + beta[9,3] + gamma[9,9] - a[9,9,3,7] <= 2
 R1034: alpha[4,4] + beta[0,4] + gamma[0,0] - a[0,0,4,4] <= 2
 R1035: alpha[4,4] + beta[1,4] + gamma[0,1] - a[0,1,4,4] <= 2
 R1036: alpha[4,4] + beta[1,4] + gamma[1,1] - a[1,1,4,4] <= 2
 R1037: alpha[4,4] + beta[2,4] + gamma[0,2] - a[0,2,4,4] <= 2
 R1038: alpha[4,4] + beta[2,4] + gamma[1,2] - a[1,2,4,4] <= 2
 R1039: alpha[4,4] + beta[2,4] + gamma[2,2] - a[2,2,4,4] <= 2
 R1040: alpha[4,4] + beta[3,4] + gamma[0,3] - a[0,3,4,4] <= 2
 R1041: alpha[4,4] + beta[3,4] + gamma[1,3] - a[1,3,4,4] <= 2
 R1042: alpha[4,4] + beta[3,4] + gamma[2,3] - a[2,3,4,4] <= 2
 R1043: alpha[4,4] + beta[3,4] + gamma[3,3] - a[3,3,4,4] <= 2
 R1044: alpha[4,4] + beta[4,4] + gamma[0,4] - a[0,4,4,4] <= 2
 R1045: alpha[4,4] + beta[4,4] + gamma[1,4] - a[1,4,4,4] <= 2
 R1046: alpha[4,4] + beta[4,4] + gamma[2,4] - a[2,4,4,4] <= 2
 R1047: alpha[4,4] + beta[4,4] + gamma[3,4] - a[3,4,4,4] <= 2
 R1048: alpha[4,4] + beta[4,4] + gamma[4,4] - a[4,4,4,4] <= 2
 R1049: alpha[4,4] + beta[5,4] + gamma[0,5] - a[0,5,4,4] <= 2
 R1050: alpha[4,4] + beta[5,4] + gamma[1,5] - a[1,5,4,4] <= 2
 R1051: alpha[4,4] + beta[5,4] + gamma[2,5] - a[2,5,4,4] <= 2
 R1052: alpha[4,4] + beta[5,4] + gamma[3,5] - a[3,5,4,4] <= 2
 R1053: alpha[4,4] + beta[5,4] + gamma[4,5] - a[4,5,4,4] <= 2
 R1054: alpha[4,4] + beta[5,4] + gamma[5,5] - a[5,5,4,4] <= 2
 R1055: alpha[4,4] + beta[6,4] + gamma[0,6] - a[0,6,4,4] <= 2
 R1056: alpha[4,4] + beta[6,4] + gamma[1,6] - a[1,6,4,4] <= 2
 R1057: alpha[4,4] + beta[6,4] + gamma[2,6] - a[2,6,4,4] <= 2
 R1058: alpha[4,4] + beta[6,4] + gamma[3,6] - a[3,6,4,4] <= 2
 R1059: alpha[4,4] + beta[6,4] + gamma[4,6] - a[4,6,4,4] <= 2
 R1060: alpha[4,4] + beta[6,4] + gamma[5,6] - a[5,6,4,4] <= 2
 R1061: alpha[4,4] + beta[6,4] + gamma[6,6] - a[6,6,4,4] <= 2
 R1062: alpha[4,4] + beta[7,4] + gamma[0,7] - a[0,7,4,4] <= 2
 R1063: alpha[4,4] + beta[7,4] + gamma[1,7] - a[1,7,4,4] <= 2
 R1064: alpha[4,4] + beta[7,4] + gamma[2,7] - a[2,7,4,4] <= 2
 R1065: alpha[4,4] + beta[7,4] + gamma[3,7] - a[3,7,4,4] <= 2
 R1066: alpha[4,4] + beta[7,4] + gamma[4,7] - a[4,7,4,4] <= 2
 R1067: alpha[4,4] + beta[7,4] + gamma[5,7] - a[5,7,4,4] <= 2
 R1068: alpha[4,4] + beta[7,4] + gamma[6,7] - a[6,7,4,4] <= 2
 R1069: alpha[4,4] + beta[7,4] + gamma[7,7] - a[7,7,4,4] <= 2
 R1070: alpha[4,4] + beta[8,4] + gamma[0,8] - a[0,8,4,4] <= 2
 R1071: alpha[4,4] + beta[8,4] + gamma[1,8] - a[1,8,4,4] <= 2
 R1072: alpha[4,4] + beta[8,4] + gamma[2,8] - a[2,8,4,4] <= 2
 R1073: alpha[4,4] + beta[8,4] + gamma[3,8] - a[3,8,4,4] <= 2
 R1074: alpha[4,4] + beta[8,4] + gamma[4,8] - a[4,8,4,4] <= 2
 R1075: alpha[4,4] + beta[8,4] + gamma[5,8] - a[5,8,4,4] <= 2
 R1076: alpha[4,4] + beta[8,4] + gamma[6,8] - a[6,8,4,4] <= 2
 R1077: alpha[4,4] + beta[8,4] + gamma[7,8] - a[7,8,4,4] <= 2
 R1078: alpha[4,4] + beta[8,4] + gamma[8,8] - a[8,8,4,4] <= 2
 R1079: alpha[4,4] + beta[9,4] + gamma[0,9] - a[0,9,4,4] <= 2
 R1080: alpha[4,4] + beta[9,4] + gamma[1,9] - a[1,9,4,4] <= 2
 R1081: alpha[4,4] + beta[9,4] + gamma[2,9] - a[2,9,4,4] <= 2
 R1082: alpha[4,4] + beta[9,4] + gamma[3,9] - a[3,9,4,4] <= 2
 R1083: alpha[4,4] + beta[9,4] + gamma[4,9] - a[4,9,4,4] <= 2
 R1084: alpha[4,4] + beta[9,4] + gamma[5,9] - a[5,9,4,4] <= 2
 R1085: alpha[4,4] + beta[9,4] + gamma[6,9] - a[6,9,4,4] <= 2
 R1086: alpha[4,4] + beta[9,4] + gamma[7,9] - a[7,9,4,4] <= 2
 R1087: alpha[4,4] + beta[9,4] + gamma[8,9] - a[8,9,4,4] <= 2
 R1088: alpha[4,4] + beta[9,4] + gamma[9,9] - a[9,9,4,4] <= 2
 R1089: alpha[5,5] + beta[0,5] + gamma[0,0] - a[0,0,5,5] <= 2
 R1090: alpha[5,5] + beta[1,5] + gamma[0,1] - a[0,1,5,5] <= 2
 R1091: alpha[5,5] + beta[1,5] + gamma[1,1] - a[1,1,5,5] <= 2
 R1092: alpha[5,5] + beta[2,5] + gamma[0,2] - a[0,2,5,5] <= 2
 R1093: alpha[5,5] + beta[2,5] + gamma[1,2] - a[1,2,5,5] <= 2
 R1094: alpha[5,5] + beta[2,5] + gamma[2,2] - a[2,2,5,5] <= 2
 R1095: alpha[5,5] + beta[3,5] + gamma[0,3] - a[0,3,5,5] <= 2
 R1096: alpha[5,5] + beta[3,5] + gamma[1,3] - a[1,3,5,5] <= 2
 R1097: alpha[5,5] + beta[3,5] + gamma[2,3] - a[2,3,5,5] <= 2
 R1098: alpha[5,5] + beta[3,5] + gamma[3,3] - a[3,3,5,5] <= 2
 R1099: alpha[5,5] + beta[4,5] + gamma[0,4] - a[0,4,5,5] <= 2
 R1100: alpha[5,5] + beta[4,5] + gamma[1,4] - a[1,4,5,5] <= 2
 R1101: alpha[5,5] + beta[4,5] + gamma[2,4] - a[2,4,5,5] <= 2
 R1102: alpha[5,5] + beta[4,5] + gamma[3,4] - a[3,4,5,5] <= 2
 R1103: alpha[5,5] + beta[4,5] + gamma[4,4] - a[4,4,5,5] <= 2
 R1104: alpha[5,5] + beta[5,5] + gamma[0,5] - a[0,5,5,5] <= 2
 R1105: alpha[5,5] + beta[5,5] + gamma[1,5] - a[1,5,5,5] <= 2
 R1106: alpha[5,5] + beta[5,5] + gamma[2,5] - a[2,5,5,5] <= 2
 R1107: alpha[5,5] + beta[5,5] + gamma[3,5] - a[3,5,5,5] <= 2
 R1108: alpha[5,5] + beta[5,5] + gamma[4,5] - a[4,5,5,5] <= 2
 R1109: alpha[5,5] + beta[5,5] + gamma[5,5] - a[5,5,5,5] <= 2
 R1110: alpha[5,5] + beta[6,5] + gamma[0,6] - a[0,6,5,5] <= 2
 R1111: alpha[5,5] + beta[6,5] + gamma[1,6] - a[1,6,5,5] <= 2
 R1112: alpha[5,5] + beta[6,5] + gamma[2,6] - a[2,6,5,5] <= 2
 R1113: alpha[5,5] + beta[6,5] + gamma[3,6] - a[3,6,5,5] <= 2
 R1114: alpha[5,5] + beta[6,5] + gamma[4,6] - a[4,6,5,5] <= 2
 R1115: alpha[5,5] + beta[6,5] + gamma[5,6] - a[5,6,5,5] <= 2
 R1116: alpha[5,5] + beta[6,5] + gamma[6,6] - a[6,6,5,5] <= 2
 R1117: alpha[5,5] + beta[7,5] + gamma[0,7] - a[0,7,5,5] <= 2
 R1118: alpha[5,5] + beta[7,5] + gamma[1,7] - a[1,7,5,5] <= 2
 R1119: alpha[5,5] + beta[7,5] + gamma[2,7] - a[2,7,5,5] <= 2
 R1120: alpha[5,5] + beta[7,5] + gamma[3,7] - a[3,7,5,5] <= 2
 R1121: alpha[5,5] + beta[7,5] + gamma[4,7] - a[4,7,5,5] <= 2
 R1122: alpha[5,5] + beta[7,5] + gamma[5,7] - a[5,7,5,5] <= 2
 R1123: alpha[5,5] + beta[7,5] + gamma[6,7] - a[6,7,5,5] <= 2
 R1124: alpha[5,5] + beta[7,5] + gamma[7,7] - a[7,7,5,5] <= 2
 R1125: alpha[5,5] + beta[8,5] + gamma[0,8] - a[0,8,5,5] <= 2
 R1126: alpha[5,5] + beta[8,5] + gamma[1,8] - a[1,8,5,5] <= 2
 R1127: alpha[5,5] + beta[8,5] + gamma[2,8] - a[2,8,5,5] <= 2
 R1128: alpha[5,5] + beta[8,5] + gamma[3,8] - a[3,8,5,5] <= 2
 R1129: alpha[5,5] + beta[8,5] + gamma[4,8] - a[4,8,5,5] <= 2
 R1130: alpha[5,5] + beta[8,5] + gamma[5,8] - a[5,8,5,5] <= 2
 R1131: alpha[5,5] + beta[8,5] + gamma[6,8] - a[6,8,5,5] <= 2
 R1132: alpha[5,5] + beta[8,5] + gamma[7,8] - a[7,8,5,5] <= 2
 R1133: alpha[5,5] + beta[8,5] + gamma[8,8] - a[8,8,5,5] <= 2
 R1134: alpha[5,5] + beta[9,5] + gamma[0,9] - a[0,9,5,5] <= 2
 R1135: alpha[5,5] + beta[9,5] + gamma[1,9] - a[1,9,5,5] <= 2
 R1136: alpha[5,5] + beta[9,5] + gamma[2,9] - a[2,9,5,5] <= 2
 R1137: alpha[5,5] + beta[9,5] + gamma[3,9] - a[3,9,5,5] <= 2
 R1138: alpha[5,5] + beta[9,5] + gamma[4,9] - a[4,9,5,5] <= 2
 R1139: alpha[5,5] + beta[9,5] + gamma[5,9] - a[5,9,5,5] <= 2
 R1140: alpha[5,5] + beta[9,5] + gamma[6,9] - a[6,9,5,5] <= 2
 R1141: alpha[5,5] + beta[9,5] + gamma[7,9] - a[7,9,5,5] <= 2
 R1142: alpha[5,5] + beta[9,5] + gamma[8,9] - a[8,9,5,5] <= 2
 R1143: alpha[5,5] + beta[9,5] + gamma[9,9] - a[9,9,5,5] <= 2
 R1144: alpha[6,6] + beta[0,6] + gamma[0,0] - a[0,0,6,6] <= 2
 R1145: alpha[6,6] + beta[1,6] + gamma[0,1] - a[0,1,6,6] <= 2
 R1146: alpha[6,6] + beta[1,6] + gamma[1,1] - a[1,1,6,6] <= 2
 R1147: alpha[6,6] + beta[2,6] + gamma[0,2] - a[0,2,6,6] <= 2
 R1148: alpha[6,6] + beta[2,6] + gamma[1,2] - a[1,2,6,6] <= 2
 R1149: alpha[6,6] + beta[2,6] + gamma[2,2] - a[2,2,6,6] <= 2
 R1150: alpha[6,6] + beta[3,6] + gamma[0,3] - a[0,3,6,6] <= 2
 R1151: alpha[6,6] + beta[3,6] + gamma[1,3] - a[1,3,6,6] <= 2
 R1152: alpha[6,6] + beta[3,6] + gamma[2,3] - a[2,3,6,6] <= 2
 R1153: alpha[6,6] + beta[3,6] + gamma[3,3] - a[3,3,6,6] <= 2
 R1154: alpha[6,6] + beta[4,6] + gamma[0,4] - a[0,4,6,6] <= 2
 R1155: alpha[6,6] + beta[4,6] + gamma[1,4] - a[1,4,6,6] <= 2
 R1156: alpha[6,6] + beta[4,6] + gamma[2,4] - a[2,4,6,6] <= 2
 R1157: alpha[6,6] + beta[4,6] + gamma[3,4] - a[3,4,6,6] <= 2
 R1158: alpha[6,6] + beta[4,6] + gamma[4,4] - a[4,4,6,6] <= 2
 R1159: alpha[6,6] + beta[5,6] + gamma[0,5] - a[0,5,6,6] <= 2
 R1160: alpha[6,6] + beta[5,6] + gamma[1,5] - a[1,5,6,6] <= 2
 R1161: alpha[6,6] + beta[5,6] + gamma[2,5] - a[2,5,6,6] <= 2
 R1162: alpha[6,6] + beta[5,6] + gamma[3,5] - a[3,5,6,6] <= 2
 R1163: alpha[6,6] + beta[5,6] + gamma[4,5] - a[4,5,6,6] <= 2
 R1164: alpha[6,6] + beta[5,6] + gamma[5,5] - a[5,5,6,6] <= 2
 R1165: alpha[6,6] + beta[6,6] + gamma[0,6] - a[0,6,6,6] <= 2
 R1166: alpha[6,6] + beta[6,6] + gamma[1,6] - a[1,6,6,6] <= 2
 R1167: alpha[6,6] + beta[6,6] + gamma[2,6] - a[2,6,6,6] <= 2
 R1168: alpha[6,6] + beta[6,6] + gamma[3,6] - a[3,6,6,6] <= 2
 R1169: alpha[6,6] + beta[6,6] + gamma[4,6] - a[4,6,6,6] <= 2
 R1170: alpha[6,6] + beta[6,6] + gamma[5,6] - a[5,6,6,6] <= 2
 R1171: alpha[6,6] + beta[6,6] + gamma[6,6] - a[6,6,6,6] <= 2
 R1172: alpha[6,6] + beta[7,6] + gamma[0,7] - a[0,7,6,6] <= 2
 R1173: alpha[6,6] + beta[7,6] + gamma[1,7] - a[1,7,6,6] <= 2
 R1174: alpha[6,6] + beta[7,6] + gamma[2,7] - a[2,7,6,6] <= 2
 R1175: alpha[6,6] + beta[7,6] + gamma[3,7] - a[3,7,6,6] <= 2
 R1176: alpha[6,6] + beta[7,6] + gamma[4,7] - a[4,7,6,6] <= 2
 R1177: alpha[6,6] + beta[7,6] + gamma[5,7] - a[5,7,6,6] <= 2
 R1178: alpha[6,6] + beta[7,6] + gamma[6,7] - a[6,7,6,6] <= 2
 R1179: alpha[6,6] + beta[7,6] + gamma[7,7] - a[7,7,6,6] <= 2
 R1180: alpha[6,6] + beta[8,6] + gamma[0,8] - a[0,8,6,6] <= 2
 R1181: alpha[6,6] + beta[8,6] + gamma[1,8] - a[1,8,6,6] <= 2
 R1182: alpha[6,6] + beta[8,6] + gamma[2,8] - a[2,8,6,6] <= 2
 R1183: alpha[6,6] + beta[8,6] + gamma[3,8] - a[3,8,6,6] <= 2
 R1184: alpha[6,6] + beta[8,6] + gamma[4,8] - a[4,8,6,6] <= 2
 R1185: alpha[6,6] + beta[8,6] + gamma[5,8] - a[5,8,6,6] <= 2
 R1186: alpha[6,6] + beta[8,6] + gamma[6,8] - a[6,8,6,6] <= 2
 R1187: alpha[6,6] + beta[8,6] + gamma[7,8] - a[7,8,6,6] <= 2
 R1188: alpha[6,6] + beta[8,6] + gamma[8,8] - a[8,8,6,6] <= 2
 R1189: alpha[6,6] + beta[9,6] + gamma[0,9] - a[0,9,6,6] <= 2
 R1190: alpha[6,6] + beta[9,6] + gamma[1,9] - a[1,9,6,6] <= 2
 R1191: alpha[6,6] + beta[9,6] + gamma[2,9] - a[2,9,6,6] <= 2
 R1192: alpha[6,6] + beta[9,6] + gamma[3,9] - a[3,9,6,6] <= 2
 R1193: alpha[6,6] + beta[9,6] + gamma[4,9] - a[4,9,6,6] <= 2
 R1194: alpha[6,6] + beta[9,6] + gamma[5,9] - a[5,9,6,6] <= 2
 R1195: alpha[6,6] + beta[9,6] + gamma[6,9] - a[6,9,6,6] <= 2
 R1196: alpha[6,6] + beta[9,6] + gamma[7,9] - a[7,9,6,6] <= 2
 R1197: alpha[6,6] + beta[9,6] + gamma[8,9] - a[8,9,6,6] <= 2
 R1198: alpha[6,6] + beta[9,6] + gamma[9,9] - a[9,9,6,6] <= 2
 R1199: alpha[6,8] + beta[0,6] + gamma[0,0] - a[0,0,6,8] <= 2
 R1200: alpha[6,8] + beta[1,6] + gamma[0,1] - a[0,1,6,8] <= 2
 R1201: alpha[6,8] + beta[1,6] + gamma[1,1] - a[1,1,6,8] <= 2
 R1202: alpha[6,8] + beta[2,6] + gamma[0,2] - a[0,2,6,8] <= 2
 R1203: alpha[6,8] + beta[2,6] + gamma[1,2] - a[1,2,6,8] <= 2
 R1204: alpha[6,8] + beta[2,6] + gamma[2,2] - a[2,2,6,8] <= 2
 R1205: alpha[6,8] + beta[3,6] + gamma[0,3] - a[0,3,6,8] <= 2
 R1206: alpha[6,8] + beta[3,6] + gamma[1,3] - a[1,3,6,8] <= 2
 R1207: alpha[6,8] + beta[3,6] + gamma[2,3] - a[2,3,6,8] <= 2
 R1208: alpha[6,8] + beta[3,6] + gamma[3,3] - a[3,3,6,8] <= 2
 R1209: alpha[6,8] + beta[4,6] + gamma[0,4] - a[0,4,6,8] <= 2
 R1210: alpha[6,8] + beta[4,6] + gamma[1,4] - a[1,4,6,8] <= 2
 R1211: alpha[6,8] + beta[4,6] + gamma[2,4] - a[2,4,6,8] <= 2
 R1212: alpha[6,8] + beta[4,6] + gamma[3,4] - a[3,4,6,8] <= 2
 R1213: alpha[6,8] + beta[4,6] + gamma[4,4] - a[4,4,6,8] <= 2
 R1214: alpha[6,8] + beta[5,6] + gamma[0,5] - a[0,5,6,8] <= 2
 R1215: alpha[6,8] + beta[5,6] + gamma[1,5] - a[1,5,6,8] <= 2
 R1216: alpha[6,8] + beta[5,6] + gamma[2,5] - a[2,5,6,8] <= 2
 R1217: alpha[6,8] + beta[5,6] + gamma[3,5] - a[3,5,6,8] <= 2
 R1218: alpha[6,8] + beta[5,6] + gamma[4,5] - a[4,5,6,8] <= 2
 R1219: alpha[6,8] + beta[5,6] + gamma[5,5] - a[5,5,6,8] <= 2
 R1220: alpha[6,8] + beta[6,6] + gamma[0,6] - a[0,6,6,8] <= 2
 R1221: alpha[6,8] + beta[6,6] + gamma[1,6] - a[1,6,6,8] <= 2
 R1222: alpha[6,8] + beta[6,6] + gamma[2,6] - a[2,6,6,8] <= 2
 R1223: alpha[6,8] + beta[6,6] + gamma[3,6] - a[3,6,6,8] <= 2
 R1224: alpha[6,8] + beta[6,6] + gamma[4,6] - a[4,6,6,8] <= 2
 R1225: alpha[6,8] + beta[6,6] + gamma[5,6] - a[5,6,6,8] <= 2
 R1226: alpha[6,8] + beta[6,6] + gamma[6,6] - a[6,6,6,8] <= 2
 R1227: alpha[6,8] + beta[7,6] + gamma[0,7] - a[0,7,6,8] <= 2
 R1228: alpha[6,8] + beta[7,6] + gamma[1,7] - a[1,7,6,8] <= 2
 R1229: alpha[6,8] + beta[7,6] + gamma[2,7] - a[2,7,6,8] <= 2
 R1230: alpha[6,8] + beta[7,6] + gamma[3,7] - a[3,7,6,8] <= 2
 R1231: alpha[6,8] + beta[7,6] + gamma[4,7] - a[4,7,6,8] <= 2
 R1232: alpha[6,8] + beta[7,6] + gamma[5,7] - a[5,7,6,8] <= 2
 R1233: alpha[6,8] + beta[7,6] + gamma[6,7] - a[6,7,6,8] <= 2
 R1234: alpha[6,8] + beta[7,6] + gamma[7,7] - a[7,7,6,8] <= 2
 R1235: alpha[6,8] + beta[8,6] + gamma[0,8] - a[0,8,6,8] <= 2
 R1236: alpha[6,8] + beta[8,6] + gamma[1,8] - a[1,8,6,8] <= 2
 R1237: alpha[6,8] + beta[8,6] + gamma[2,8] - a[2,8,6,8] <= 2
 R1238: alpha[6,8] + beta[8,6] + gamma[3,8] - a[3,8,6,8] <= 2
 R1239: alpha[6,8] + beta[8,6] + gamma[4,8] - a[4,8,6,8] <= 2
 R1240: alpha[6,8] + beta[8,6] + gamma[5,8] - a[5,8,6,8] <= 2
 R1241: alpha[6,8] + beta[8,6] + gamma[6,8] - a[6,8,6,8] <= 2
 R1242: alpha[6,8] + beta[8,6] + gamma[7,8] - a[7,8,6,8] <= 2
 R1243: alpha[6,8] + beta[8,6] + gamma[8,8] - a[8,8,6,8] <= 2
 R1244: alpha[6,8] + beta[9,6] + gamma[0,9] - a[0,9,6,8] <= 2
 R1245: alpha[6,8] + beta[9,6] + gamma[1,9] - a[1,9,6,8] <= 2
 R1246: alpha[6,8] + beta[9,6] + gamma[2,9] - a[2,9,6,8] <= 2
 R1247: alpha[6,8] + beta[9,6] + gamma[3,9] - a[3,9,6,8] <= 2
 R1248: alpha[6,8] + beta[9,6] + gamma[4,9] - a[4,9,6,8] <= 2
 R1249: alpha[6,8] + beta[9,6] + gamma[5,9] - a[5,9,6,8] <= 2
 R1250: alpha[6,8] + beta[9,6] + gamma[6,9] - a[6,9,6,8] <= 2
 R1251: alpha[6,8] + beta[9,6] + gamma[7,9] - a[7,9,6,8] <= 2
 R1252: alpha[6,8] + beta[9,6] + gamma[8,9] - a[8,9,6,8] <= 2
 R1253: alpha[6,8] + beta[9,6] + gamma[9,9] - a[9,9,6,8] <= 2
 R1254: alpha[6,9] + beta[0,6] + gamma[0,0] - a[0,0,6,9] <= 2
 R1255: alpha[6,9] + beta[1,6] + gamma[0,1] - a[0,1,6,9] <= 2
 R1256: alpha[6,9] + beta[1,6] + gamma[1,1] - a[1,1,6,9] <= 2
 R1257: alpha[6,9] + beta[2,6] + gamma[0,2] - a[0,2,6,9] <= 2
 R1258: alpha[6,9] + beta[2,6] + gamma[1,2] - a[1,2,6,9] <= 2
 R1259: alpha[6,9] + beta[2,6] + gamma[2,2] - a[2,2,6,9] <= 2
 R1260: alpha[6,9] + beta[3,6] + gamma[0,3] - a[0,3,6,9] <= 2
 R1261: alpha[6,9] + beta[3,6] + gamma[1,3] - a[1,3,6,9] <= 2
 R1262: alpha[6,9] + beta[3,6] + gamma[2,3] - a[2,3,6,9] <= 2
 R1263: alpha[6,9] + beta[3,6] + gamma[3,3] - a[3,3,6,9] <= 2
 R1264: alpha[6,9] + beta[4,6] + gamma[0,4] - a[0,4,6,9] <= 2
 R1265: alpha[6,9] + beta[4,6] + gamma[1,4] - a[1,4,6,9] <= 2
 R1266: alpha[6,9] + beta[4,6] + gamma[2,4] - a[2,4,6,9] <= 2
 R1267: alpha[6,9] + beta[4,6] + gamma[3,4] - a[3,4,6,9] <= 2
 R1268: alpha[6,9] + beta[4,6] + gamma[4,4] - a[4,4,6,9] <= 2
 R1269: alpha[6,9] + beta[5,6] + gamma[0,5] - a[0,5,6,9] <= 2
 R1270: alpha[6,9] + beta[5,6] + gamma[1,5] - a[1,5,6,9] <= 2
 R1271: alpha[6,9] + beta[5,6] + gamma[2,5] - a[2,5,6,9] <= 2
 R1272: alpha[6,9] + beta[5,6] + gamma[3,5] - a[3,5,6,9] <= 2
 R1273: alpha[6,9] + beta[5,6] + gamma[4,5] - a[4,5,6,9] <= 2
 R1274: alpha[6,9] + beta[5,6] + gamma[5,5] - a[5,5,6,9] <= 2
 R1275: alpha[6,9] + beta[6,6] + gamma[0,6] - a[0,6,6,9] <= 2
 R1276: alpha[6,9] + beta[6,6] + gamma[1,6] - a[1,6,6,9] <= 2
 R1277: alpha[6,9] + beta[6,6] + gamma[2,6] - a[2,6,6,9] <= 2
 R1278: alpha[6,9] + beta[6,6] + gamma[3,6] - a[3,6,6,9] <= 2
 R1279: alpha[6,9] + beta[6,6] + gamma[4,6] - a[4,6,6,9] <= 2
 R1280: alpha[6,9] + beta[6,6] + gamma[5,6] - a[5,6,6,9] <= 2
 R1281: alpha[6,9] + beta[6,6] + gamma[6,6] - a[6,6,6,9] <= 2
 R1282: alpha[6,9] + beta[7,6] + gamma[0,7] - a[0,7,6,9] <= 2
 R1283: alpha[6,9] + beta[7,6] + gamma[1,7] - a[1,7,6,9] <= 2
 R1284: alpha[6,9] + beta[7,6] + gamma[2,7] - a[2,7,6,9] <= 2
 R1285: alpha[6,9] + beta[7,6] + gamma[3,7] - a[3,7,6,9] <= 2
 R1286: alpha[6,9] + beta[7,6] + gamma[4,7] - a[4,7,6,9] <= 2
 R1287: alpha[6,9] + beta[7,6] + gamma[5,7] - a[5,7,6,9] <= 2
 R1288: alpha[6,9] + beta[7,6] + gamma[6,7] - a[6,7,6,9] <= 2
 R1289: alpha[6,9] + beta[7,6] + gamma[7,7] - a[7,7,6,9] <= 2
 R1290: alpha[6,9] + beta[8,6] + gamma[0,8] - a[0,8,6,9] <= 2
 R1291: alpha[6,9] + beta[8,6] + gamma[1,8] - a[1,8,6,9] <= 2
 R1292: alpha[6,9] + beta[8,6] + gamma[2,8] - a[2,8,6,9] <= 2
 R1293: alpha[6,9] + beta[8,6] + gamma[3,8] - a[3,8,6,9] <= 2
 R1294: alpha[6,9] + beta[8,6] + gamma[4,8] - a[4,8,6,9] <= 2
 R1295: alpha[6,9] + beta[8,6] + gamma[5,8] - a[5,8,6,9] <= 2
 R1296: alpha[6,9] + beta[8,6] + gamma[6,8] - a[6,8,6,9] <= 2
 R1297: alpha[6,9] + beta[8,6] + gamma[7,8] - a[7,8,6,9] <= 2
 R1298: alpha[6,9] + beta[8,6] + gamma[8,8] - a[8,8,6,9] <= 2
 R1299: alpha[6,9] + beta[9,6] + gamma[0,9] - a[0,9,6,9] <= 2
 R1300: alpha[6,9] + beta[9,6] + gamma[1,9] - a[1,9,6,9] <= 2
 R1301: alpha[6,9] + beta[9,6] + gamma[2,9] - a[2,9,6,9] <= 2
 R1302: alpha[6,9] + beta[9,6] + gamma[3,9] - a[3,9,6,9] <= 2
 R1303: alpha[6,9] + beta[9,6] + gamma[4,9] - a[4,9,6,9] <= 2
 R1304: alpha[6,9] + beta[9,6] + gamma[5,9] - a[5,9,6,9] <= 2
 R1305: alpha[6,9] + beta[9,6] + gamma[6,9] - a[6,9,6,9] <= 2
 R1306: alpha[6,9] + beta[9,6] + gamma[7,9] - a[7,9,6,9] <= 2
 R1307: alpha[6,9] + beta[9,6] + gamma[8,9] - a[8,9,6,9] <= 2
 R1308: alpha[6,9] + beta[9,6] + gamma[9,9] - a[9,9,6,9] <= 2
 R1309: alpha[7,7] + beta[0,7] + gamma[0,0] - a[0,0,7,7] <= 2
 R1310: alpha[7,7] + beta[1,7] + gamma[0,1] - a[0,1,7,7] <= 2
 R1311: alpha[7,7] + beta[1,7] + gamma[1,1] - a[1,1,7,7] <= 2
 R1312: alpha[7,7] + beta[2,7] + gamma[0,2] - a[0,2,7,7] <= 2
 R1313: alpha[7,7] + beta[2,7] + gamma[1,2] - a[1,2,7,7] <= 2
 R1314: alpha[7,7] + beta[2,7] + gamma[2,2] - a[2,2,7,7] <= 2
 R1315: alpha[7,7] + beta[3,7] + gamma[0,3] - a[0,3,7,7] <= 2
 R1316: alpha[7,7] + beta[3,7] + gamma[1,3] - a[1,3,7,7] <= 2
 R1317: alpha[7,7] + beta[3,7] + gamma[2,3] - a[2,3,7,7] <= 2
 R1318: alpha[7,7] + beta[3,7] + gamma[3,3] - a[3,3,7,7] <= 2
 R1319: alpha[7,7] + beta[4,7] + gamma[0,4] - a[0,4,7,7] <= 2
 R1320: alpha[7,7] + beta[4,7] + gamma[1,4] - a[1,4,7,7] <= 2
 R1321: alpha[7,7] + beta[4,7] + gamma[2,4] - a[2,4,7,7] <= 2
 R1322: alpha[7,7] + beta[4,7] + gamma[3,4] - a[3,4,7,7] <= 2
 R1323: alpha[7,7] + beta[4,7] + gamma[4,4] - a[4,4,7,7] <= 2
 R1324: alpha[7,7] + beta[5,7] + gamma[0,5] - a[0,5,7,7] <= 2
 R1325: alpha[7,7] + beta[5,7] + gamma[1,5] - a[1,5,7,7] <= 2
 R1326: alpha[7,7] + beta[5,7] + gamma[2,5] - a[2,5,7,7] <= 2
 R1327: alpha[7,7] + beta[5,7] + gamma[3,5] - a[3,5,7,7] <= 2
 R1328: alpha[7,7] + beta[5,7] + gamma[4,5] - a[4,5,7,7] <= 2
 R1329: alpha[7,7] + beta[5,7] + gamma[5,5] - a[5,5,7,7] <= 2
 R1330: alpha[7,7] + beta[6,7] + gamma[0,6] - a[0,6,7,7] <= 2
 R1331: alpha[7,7] + beta[6,7] + gamma[1,6] - a[1,6,7,7] <= 2
 R1332: alpha[7,7] + beta[6,7] + gamma[2,6] - a[2,6,7,7] <= 2
 R1333: alpha[7,7] + beta[6,7] + gamma[3,6] - a[3,6,7,7] <= 2
 R1334: alpha[7,7] + beta[6,7] + gamma[4,6] - a[4,6,7,7] <= 2
 R1335: alpha[7,7] + beta[6,7] + gamma[5,6] - a[5,6,7,7] <= 2
 R1336: alpha[7,7] + beta[6,7] + gamma[6,6] - a[6,6,7,7] <= 2
 R1337: alpha[7,7] + beta[7,7] + gamma[0,7] - a[0,7,7,7] <= 2
 R1338: alpha[7,7] + beta[7,7] + gamma[1,7] - a[1,7,7,7] <= 2
 R1339: alpha[7,7] + beta[7,7] + gamma[2,7] - a[2,7,7,7] <= 2
 R1340: alpha[7,7] + beta[7,7] + gamma[3,7] - a[3,7,7,7] <= 2
 R1341: alpha[7,7] + beta[7,7] + gamma[4,7] - a[4,7,7,7] <= 2
 R1342: alpha[7,7] + beta[7,7] + gamma[5,7] - a[5,7,7,7] <= 2
 R1343: alpha[7,7] + beta[7,7] + gamma[6,7] - a[6,7,7,7] <= 2
 R1344: alpha[7,7] + beta[7,7] + gamma[7,7] - a[7,7,7,7] <= 2
 R1345: alpha[7,7] + beta[8,7] + gamma[0,8] - a[0,8,7,7] <= 2
 R1346: alpha[7,7] + beta[8,7] + gamma[1,8] - a[1,8,7,7] <= 2
 R1347: alpha[7,7] + beta[8,7] + gamma[2,8] - a[2,8,7,7] <= 2
 R1348: alpha[7,7] + beta[8,7] + gamma[3,8] - a[3,8,7,7] <= 2
 R1349: alpha[7,7] + beta[8,7] + gamma[4,8] - a[4,8,7,7] <= 2
 R1350: alpha[7,7] + beta[8,7] + gamma[5,8] - a[5,8,7,7] <= 2
 R1351: alpha[7,7] + beta[8,7] + gamma[6,8] - a[6,8,7,7] <= 2
 R1352: alpha[7,7] + beta[8,7] + gamma[7,8] - a[7,8,7,7] <= 2
 R1353: alpha[7,7] + beta[8,7] + gamma[8,8] - a[8,8,7,7] <= 2
 R1354: alpha[7,7] + beta[9,7] + gamma[0,9] - a[0,9,7,7] <= 2
 R1355: alpha[7,7] + beta[9,7] + gamma[1,9] - a[1,9,7,7] <= 2
 R1356: alpha[7,7] + beta[9,7] + gamma[2,9] - a[2,9,7,7] <= 2
 R1357: alpha[7,7] + beta[9,7] + gamma[3,9] - a[3,9,7,7] <= 2
 R1358: alpha[7,7] + beta[9,7] + gamma[4,9] - a[4,9,7,7] <= 2
 R1359: alpha[7,7] + beta[9,7] + gamma[5,9] - a[5,9,7,7] <= 2
 R1360: alpha[7,7] + beta[9,7] + gamma[6,9] - a[6,9,7,7] <= 2
 R1361: alpha[7,7] + beta[9,7] + gamma[7,9] - a[7,9,7,7] <= 2
 R1362: alpha[7,7] + beta[9,7] + gamma[8,9] - a[8,9,7,7] <= 2
 R1363: alpha[7,7] + beta[9,7] + gamma[9,9] - a[9,9,7,7] <= 2
 R1364: alpha[8,8] + beta[0,8] + gamma[0,0] - a[0,0,8,8] <= 2
 R1365: alpha[8,8] + beta[1,8] + gamma[0,1] - a[0,1,8,8] <= 2
 R1366: alpha[8,8] + beta[1,8] + gamma[1,1] - a[1,1,8,8] <= 2
 R1367: alpha[8,8] + beta[2,8] + gamma[0,2] - a[0,2,8,8] <= 2
 R1368: alpha[8,8] + beta[2,8] + gamma[1,2] - a[1,2,8,8] <= 2
 R1369: alpha[8,8] + beta[2,8] + gamma[2,2] - a[2,2,8,8] <= 2
 R1370: alpha[8,8] + beta[3,8] + gamma[0,3] - a[0,3,8,8] <= 2
 R1371: alpha[8,8] + beta[3,8] + gamma[1,3] - a[1,3,8,8] <= 2
 R1372: alpha[8,8] + beta[3,8] + gamma[2,3] - a[2,3,8,8] <= 2
 R1373: alpha[8,8] + beta[3,8] + gamma[3,3] - a[3,3,8,8] <= 2
 R1374: alpha[8,8] + beta[4,8] + gamma[0,4] - a[0,4,8,8] <= 2
 R1375: alpha[8,8] + beta[4,8] + gamma[1,4] - a[1,4,8,8] <= 2
 R1376: alpha[8,8] + beta[4,8] + gamma[2,4] - a[2,4,8,8] <= 2
 R1377: alpha[8,8] + beta[4,8] + gamma[3,4] - a[3,4,8,8] <= 2
 R1378: alpha[8,8] + beta[4,8] + gamma[4,4] - a[4,4,8,8] <= 2
 R1379: alpha[8,8] + beta[5,8] + gamma[0,5] - a[0,5,8,8] <= 2
 R1380: alpha[8,8] + beta[5,8] + gamma[1,5] - a[1,5,8,8] <= 2
 R1381: alpha[8,8] + beta[5,8] + gamma[2,5] - a[2,5,8,8] <= 2
 R1382: alpha[8,8] + beta[5,8] + gamma[3,5] - a[3,5,8,8] <= 2
 R1383: alpha[8,8] + beta[5,8] + gamma[4,5] - a[4,5,8,8] <= 2
 R1384: alpha[8,8] + beta[5,8] + gamma[5,5] - a[5,5,8,8] <= 2
 R1385: alpha[8,8] + beta[6,8] + gamma[0,6] - a[0,6,8,8] <= 2
 R1386: alpha[8,8] + beta[6,8] + gamma[1,6] - a[1,6,8,8] <= 2
 R1387: alpha[8,8] + beta[6,8] + gamma[2,6] - a[2,6,8,8] <= 2
 R1388: alpha[8,8] + beta[6,8] + gamma[3,6] - a[3,6,8,8] <= 2
 R1389: alpha[8,8] + beta[6,8] + gamma[4,6] - a[4,6,8,8] <= 2
 R1390: alpha[8,8] + beta[6,8] + gamma[5,6] - a[5,6,8,8] <= 2
 R1391: alpha[8,8] + beta[6,8] + gamma[6,6] - a[6,6,8,8] <= 2
 R1392: alpha[8,8] + beta[7,8] + gamma[0,7] - a[0,7,8,8] <= 2
 R1393: alpha[8,8] + beta[7,8] + gamma[1,7] - a[1,7,8,8] <= 2
 R1394: alpha[8,8] + beta[7,8] + gamma[2,7] - a[2,7,8,8] <= 2
 R1395: alpha[8,8] + beta[7,8] + gamma[3,7] - a[3,7,8,8] <= 2
 R1396: alpha[8,8] + beta[7,8] + gamma[4,7] - a[4,7,8,8] <= 2
 R1397: alpha[8,8] + beta[7,8] + gamma[5,7] - a[5,7,8,8] <= 2
 R1398: alpha[8,8] + beta[7,8] + gamma[6,7] - a[6,7,8,8] <= 2
 R1399: alpha[8,8] + beta[7,8] + gamma[7,7] - a[7,7,8,8] <= 2
 R1400: alpha[8,8] + beta[8,8] + gamma[0,8] - a[0,8,8,8] <= 2
 R1401: alpha[8,8] + beta[8,8] + gamma[1,8] - a[1,8,8,8] <= 2
 R1402: alpha[8,8] + beta[8,8] + gamma[2,8] - a[2,8,8,8] <= 2
 R1403: alpha[8,8] + beta[8,8] + gamma[3,8] - a[3,8,8,8] <= 2
 R1404: alpha[8,8] + beta[8,8] + gamma[4,8] - a[4,8,8,8] <= 2
 R1405: alpha[8,8] + beta[8,8] + gamma[5,8] - a[5,8,8,8] <= 2
 R1406: alpha[8,8] + beta[8,8] + gamma[6,8] - a[6,8,8,8] <= 2
 R1407: alpha[8,8] + beta[8,8] + gamma[7,8] - a[7,8,8,8] <= 2
 R1408: alpha[8,8] + beta[8,8] + gamma[8,8] - a[8,8,8,8] <= 2
 R1409: alpha[8,8] + beta[9,8] + gamma[0,9] - a[0,9,8,8] <= 2
 R1410: alpha[8,8] + beta[9,8] + gamma[1,9] - a[1,9,8,8] <= 2
 R1411: alpha[8,8] + beta[9,8] + gamma[2,9] - a[2,9,8,8] <= 2
 R1412: alpha[8,8] + beta[9,8] + gamma[3,9] - a[3,9,8,8] <= 2
 R1413: alpha[8,8] + beta[9,8] + gamma[4,9] - a[4,9,8,8] <= 2
 R1414: alpha[8,8] + beta[9,8] + gamma[5,9] - a[5,9,8,8] <= 2
 R1415: alpha[8,8] + beta[9,8] + gamma[6,9] - a[6,9,8,8] <= 2
 R1416: alpha[8,8] + beta[9,8] + gamma[7,9] - a[7,9,8,8] <= 2
 R1417: alpha[8,8] + beta[9,8] + gamma[8,9] - a[8,9,8,8] <= 2
 R1418: alpha[8,8] + beta[9,8] + gamma[9,9] - a[9,9,8,8] <= 2
 R1419: alpha[8,9] + beta[0,8] + gamma[0,0] - a[0,0,8,9] <= 2
 R1420: alpha[8,9] + beta[1,8] + gamma[0,1] - a[0,1,8,9] <= 2
 R1421: alpha[8,9] + beta[1,8] + gamma[1,1] - a[1,1,8,9] <= 2
 R1422: alpha[8,9] + beta[2,8] + gamma[0,2] - a[0,2,8,9] <= 2
 R1423: alpha[8,9] + beta[2,8] + gamma[1,2] - a[1,2,8,9] <= 2
 R1424: alpha[8,9] + beta[2,8] + gamma[2,2] - a[2,2,8,9] <= 2
 R1425: alpha[8,9] + beta[3,8] + gamma[0,3] - a[0,3,8,9] <= 2
 R1426: alpha[8,9] + beta[3,8] + gamma[1,3] - a[1,3,8,9] <= 2
 R1427: alpha[8,9] + beta[3,8] + gamma[2,3] - a[2,3,8,9] <= 2
 R1428: alpha[8,9] + beta[3,8] + gamma[3,3] - a[3,3,8,9] <= 2
 R1429: alpha[8,9] + beta[4,8] + gamma[0,4] - a[0,4,8,9] <= 2
 R1430: alpha[8,9] + beta[4,8] + gamma[1,4] - a[1,4,8,9] <= 2
 R1431: alpha[8,9] + beta[4,8] + gamma[2,4] - a[2,4,8,9] <= 2
 R1432: alpha[8,9] + beta[4,8] + gamma[3,4] - a[3,4,8,9] <= 2
 R1433: alpha[8,9] + beta[4,8] + gamma[4,4] - a[4,4,8,9] <= 2
 R1434: alpha[8,9] + beta[5,8] + gamma[0,5] - a[0,5,8,9] <= 2
 R1435: alpha[8,9] + beta[5,8] + gamma[1,5] - a[1,5,8,9] <= 2
 R1436: alpha[8,9] + beta[5,8] + gamma[2,5] - a[2,5,8,9] <= 2
 R1437: alpha[8,9] + beta[5,8] + gamma[3,5] - a[3,5,8,9] <= 2
 R1438: alpha[8,9] + beta[5,8] + gamma[4,5] - a[4,5,8,9] <= 2
 R1439: alpha[8,9] + beta[5,8] + gamma[5,5] - a[5,5,8,9] <= 2
 R1440: alpha[8,9] + beta[6,8] + gamma[0,6] - a[0,6,8,9] <= 2
 R1441: alpha[8,9] + beta[6,8] + gamma[1,6] - a[1,6,8,9] <= 2
 R1442: alpha[8,9] + beta[6,8] + gamma[2,6] - a[2,6,8,9] <= 2
 R1443: alpha[8,9] + beta[6,8] + gamma[3,6] - a[3,6,8,9] <= 2
 R1444: alpha[8,9] + beta[6,8] + gamma[4,6] - a[4,6,8,9] <= 2
 R1445: alpha[8,9] + beta[6,8] + gamma[5,6] - a[5,6,8,9] <= 2
 R1446: alpha[8,9] + beta[6,8] + gamma[6,6] - a[6,6,8,9] <= 2
 R1447: alpha[8,9] + beta[7,8] + gamma[0,7] - a[0,7,8,9] <= 2
 R1448: alpha[8,9] + beta[7,8] + gamma[1,7] - a[1,7,8,9] <= 2
 R1449: alpha[8,9] + beta[7,8] + gamma[2,7] - a[2,7,8,9] <= 2
 R1450: alpha[8,9] + beta[7,8] + gamma[3,7] - a[3,7,8,9] <= 2
 R1451: alpha[8,9] + beta[7,8] + gamma[4,7] - a[4,7,8,9] <= 2
 R1452: alpha[8,9] + beta[7,8] + gamma[5,7] - a[5,7,8,9] <= 2
 R1453: alpha[8,9] + beta[7,8] + gamma[6,7] - a[6,7,8,9] <= 2
 R1454: alpha[8,9] + beta[7,8] + gamma[7,7] - a[7,7,8,9] <= 2
 R1455: alpha[8,9] + beta[8,8] + gamma[0,8] - a[0,8,8,9] <= 2
 R1456: alpha[8,9] + beta[8,8] + gamma[1,8] - a[1,8,8,9] <= 2
 R1457: alpha[8,9] + beta[8,8] + gamma[2,8] - a[2,8,8,9] <= 2
 R1458: alpha[8,9] + beta[8,8] + gamma[3,8] - a[3,8,8,9] <= 2
 R1459: alpha[8,9] + beta[8,8] + gamma[4,8] - a[4,8,8,9] <= 2
 R1460: alpha[8,9] + beta[8,8] + gamma[5,8] - a[5,8,8,9] <= 2
 R1461: alpha[8,9] + beta[8,8] + gamma[6,8] - a[6,8,8,9] <= 2
 R1462: alpha[8,9] + beta[8,8] + gamma[7,8] - a[7,8,8,9] <= 2
 R1463: alpha[8,9] + beta[8,8] + gamma[8,8] - a[8,8,8,9] <= 2
 R1464: alpha[8,9] + beta[9,8] + gamma[0,9] - a[0,9,8,9] <= 2
 R1465: alpha[8,9] + beta[9,8] + gamma[1,9] - a[1,9,8,9] <= 2
 R1466: alpha[8,9] + beta[9,8] + gamma[2,9] - a[2,9,8,9] <= 2
 R1467: alpha[8,9] + beta[9,8] + gamma[3,9] - a[3,9,8,9] <= 2
 R1468: alpha[8,9] + beta[9,8] + gamma[4,9] - a[4,9,8,9] <= 2
 R1469: alpha[8,9] + beta[9,8] + gamma[5,9] - a[5,9,8,9] <= 2
 R1470: alpha[8,9] + beta[9,8] + gamma[6,9] - a[6,9,8,9] <= 2
 R1471: alpha[8,9] + beta[9,8] + gamma[7,9] - a[7,9,8,9] <= 2
 R1472: alpha[8,9] + beta[9,8] + gamma[8,9] - a[8,9,8,9] <= 2
 R1473: alpha[8,9] + beta[9,8] + gamma[9,9] - a[9,9,8,9] <= 2
 R1474: alpha[9,9] + beta[0,9] + gamma[0,0] - a[0,0,9,9] <= 2
 R1475: alpha[9,9] + beta[1,9] + gamma[0,1] - a[0,1,9,9] <= 2
 R1476: alpha[9,9] + beta[1,9] + gamma[1,1] - a[1,1,9,9] <= 2
 R1477: alpha[9,9] + beta[2,9] + gamma[0,2] - a[0,2,9,9] <= 2
 R1478: alpha[9,9] + beta[2,9] + gamma[1,2] - a[1,2,9,9] <= 2
 R1479: alpha[9,9] + beta[2,9] + gamma[2,2] - a[2,2,9,9] <= 2
 R1480: alpha[9,9] + beta[3,9] + gamma[0,3] - a[0,3,9,9] <= 2
 R1481: alpha[9,9] + beta[3,9] + gamma[1,3] - a[1,3,9,9] <= 2
 R1482: alpha[9,9] + beta[3,9] + gamma[2,3] - a[2,3,9,9] <= 2
 R1483: alpha[9,9] + beta[3,9] + gamma[3,3] - a[3,3,9,9] <= 2
 R1484: alpha[9,9] + beta[4,9] + gamma[0,4] - a[0,4,9,9] <= 2
 R1485: alpha[9,9] + beta[4,9] + gamma[1,4] - a[1,4,9,9] <= 2
 R1486: alpha[9,9] + beta[4,9] + gamma[2,4] - a[2,4,9,9] <= 2
 R1487: alpha[9,9] + beta[4,9] + gamma[3,4] - a[3,4,9,9] <= 2
 R1488: alpha[9,9] + beta[4,9] + gamma[4,4] - a[4,4,9,9] <= 2
 R1489: alpha[9,9] + beta[5,9] + gamma[0,5] - a[0,5,9,9] <= 2
 R1490: alpha[9,9] + beta[5,9] + gamma[1,5] - a[1,5,9,9] <= 2
 R1491: alpha[9,9] + beta[5,9] + gamma[2,5] - a[2,5,9,9] <= 2
 R1492: alpha[9,9] + beta[5,9] + gamma[3,5] - a[3,5,9,9] <= 2
 R1493: alpha[9,9] + beta[5,9] + gamma[4,5] - a[4,5,9,9] <= 2
 R1494: alpha[9,9] + beta[5,9] + gamma[5,5] - a[5,5,9,9] <= 2
 R1495: alpha[9,9] + beta[6,9] + gamma[0,6] - a[0,6,9,9] <= 2
 R1496: alpha[9,9] + beta[6,9] + gamma[1,6] - a[1,6,9,9] <= 2
 R1497: alpha[9,9] + beta[6,9] + gamma[2,6] - a[2,6,9,9] <= 2
 R1498: alpha[9,9] + beta[6,9] + gamma[3,6] - a[3,6,9,9] <= 2
 R1499: alpha[9,9] + beta[6,9] + gamma[4,6] - a[4,6,9,9] <= 2
 R1500: alpha[9,9] + beta[6,9] + gamma[5,6] - a[5,6,9,9] <= 2
 R1501: alpha[9,9] + beta[6,9] + gamma[6,6] - a[6,6,9,9] <= 2
 R1502: alpha[9,9] + beta[7,9] + gamma[0,7] - a[0,7,9,9] <= 2
 R1503: alpha[9,9] + beta[7,9] + gamma[1,7] - a[1,7,9,9] <= 2
 R1504: alpha[9,9] + beta[7,9] + gamma[2,7] - a[2,7,9,9] <= 2
 R1505: alpha[9,9] + beta[7,9] + gamma[3,7] - a[3,7,9,9] <= 2
 R1506: alpha[9,9] + beta[7,9] + gamma[4,7] - a[4,7,9,9] <= 2
 R1507: alpha[9,9] + beta[7,9] + gamma[5,7] - a[5,7,9,9] <= 2
 R1508: alpha[9,9] + beta[7,9] + gamma[6,7] - a[6,7,9,9] <= 2
 R1509: alpha[9,9] + beta[7,9] + gamma[7,7] - a[7,7,9,9] <= 2
 R1510: alpha[9,9] + beta[8,9] + gamma[0,8] - a[0,8,9,9] <= 2
 R1511: alpha[9,9] + beta[8,9] + gamma[1,8] - a[1,8,9,9] <= 2
 R1512: alpha[9,9] + beta[8,9] + gamma[2,8] - a[2,8,9,9] <= 2
 R1513: alpha[9,9] + beta[8,9] + gamma[3,8] - a[3,8,9,9] <= 2
 R1514: alpha[9,9] + beta[8,9] + gamma[4,8] - a[4,8,9,9] <= 2
 R1515: alpha[9,9] + beta[8,9] + gamma[5,8] - a[5,8,9,9] <= 2
 R1516: alpha[9,9] + beta[8,9] + gamma[6,8] - a[6,8,9,9] <= 2
 R1517: alpha[9,9] + beta[8,9] + gamma[7,8] - a[7,8,9,9] <= 2
 R1518: alpha[9,9] + beta[8,9] + gamma[8,8] - a[8,8,9,9] <= 2
 R1519: alpha[9,9] + beta[9,9] + gamma[0,9] - a[0,9,9,9] <= 2
 R1520: alpha[9,9] + beta[9,9] + gamma[1,9] - a[1,9,9,9] <= 2
 R1521: alpha[9,9] + beta[9,9] + gamma[2,9] - a[2,9,9,9] <= 2
 R1522: alpha[9,9] + beta[9,9] + gamma[3,9] - a[3,9,9,9] <= 2
 R1523: alpha[9,9] + beta[9,9] + gamma[4,9] - a[4,9,9,9] <= 2
 R1524: alpha[9,9] + beta[9,9] + gamma[5,9] - a[5,9,9,9] <= 2
 R1525: alpha[9,9] + beta[9,9] + gamma[6,9] - a[6,9,9,9] <= 2
 R1526: alpha[9,9] + beta[9,9] + gamma[7,9] - a[7,9,9,9] <= 2
 R1527: alpha[9,9] + beta[9,9] + gamma[8,9] - a[8,9,9,9] <= 2
 R1528: alpha[9,9] + beta[9,9] + gamma[9,9] - a[9,9,9,9] <= 2
 R1529: alpha[0,0] + beta[0,0] + gamma[0,0] - 3 a[0,0,0,0] >= 0
 R1530: alpha[0,0] + beta[1,0] + gamma[0,1] - 3 a[0,1,0,0] >= 0
 R1531: alpha[0,0] + beta[1,0] + gamma[1,1] - 3 a[1,1,0,0] >= 0
 R1532: alpha[0,0] + beta[2,0] + gamma[0,2] - 3 a[0,2,0,0] >= 0
 R1533: alpha[0,0] + beta[2,0] + gamma[1,2] - 3 a[1,2,0,0] >= 0
 R1534: alpha[0,0] + beta[2,0] + gamma[2,2] - 3 a[2,2,0,0] >= 0
 R1535: alpha[0,0] + beta[3,0] + gamma[0,3] - 3 a[0,3,0,0] >= 0
 R1536: alpha[0,0] + beta[3,0] + gamma[1,3] - 3 a[1,3,0,0] >= 0
 R1537: alpha[0,0] + beta[3,0] + gamma[2,3] - 3 a[2,3,0,0] >= 0
 R1538: alpha[0,0] + beta[3,0] + gamma[3,3] - 3 a[3,3,0,0] >= 0
 R1539: alpha[0,0] + beta[4,0] + gamma[0,4] - 3 a[0,4,0,0] >= 0
 R1540: alpha[0,0] + beta[4,0] + gamma[1,4] - 3 a[1,4,0,0] >= 0
 R1541: alpha[0,0] + beta[4,0] + gamma[2,4] - 3 a[2,4,0,0] >= 0
 R1542: alpha[0,0] + beta[4,0] + gamma[3,4] - 3 a[3,4,0,0] >= 0
 R1543: alpha[0,0] + beta[4,0] + gamma[4,4] - 3 a[4,4,0,0] >= 0
 R1544: alpha[0,0] + beta[5,0] + gamma[0,5] - 3 a[0,5,0,0] >= 0
 R1545: alpha[0,0] + beta[5,0] + gamma[1,5] - 3 a[1,5,0,0] >= 0
 R1546: alpha[0,0] + beta[5,0] + gamma[2,5] - 3 a[2,5,0,0] >= 0
 R1547: alpha[0,0] + beta[5,0] + gamma[3,5] - 3 a[3,5,0,0] >= 0
 R1548: alpha[0,0] + beta[5,0] + gamma[4,5] - 3 a[4,5,0,0] >= 0
 R1549: alpha[0,0] + beta[5,0] + gamma[5,5] - 3 a[5,5,0,0] >= 0
 R1550: alpha[0,0] + beta[6,0] + gamma[0,6] - 3 a[0,6,0,0] >= 0
 R1551: alpha[0,0] + beta[6,0] + gamma[1,6] - 3 a[1,6,0,0] >= 0
 R1552: alpha[0,0] + beta[6,0] + gamma[2,6] - 3 a[2,6,0,0] >= 0
 R1553: alpha[0,0] + beta[6,0] + gamma[3,6] - 3 a[3,6,0,0] >= 0
 R1554: alpha[0,0] + beta[6,0] + gamma[4,6] - 3 a[4,6,0,0] >= 0
 R1555: alpha[0,0] + beta[6,0] + gamma[5,6] - 3 a[5,6,0,0] >= 0
 R1556: alpha[0,0] + beta[6,0] + gamma[6,6] - 3 a[6,6,0,0] >= 0
 R1557: alpha[0,0] + beta[7,0] + gamma[0,7] - 3 a[0,7,0,0] >= 0
 R1558: alpha[0,0] + beta[7,0] + gamma[1,7] - 3 a[1,7,0,0] >= 0
 R1559: alpha[0,0] + beta[7,0] + gamma[2,7] - 3 a[2,7,0,0] >= 0
 R1560: alpha[0,0] + beta[7,0] + gamma[3,7] - 3 a[3,7,0,0] >= 0
 R1561: alpha[0,0] + beta[7,0] + gamma[4,7] - 3 a[4,7,0,0] >= 0
 R1562: alpha[0,0] + beta[7,0] + gamma[5,7] - 3 a[5,7,0,0] >= 0
 R1563: alpha[0,0] + beta[7,0] + gamma[6,7] - 3 a[6,7,0,0] >= 0
 R1564: alpha[0,0] + beta[7,0] + gamma[7,7] - 3 a[7,7,0,0] >= 0
 R1565: alpha[0,0] + beta[8,0] + gamma[0,8] - 3 a[0,8,0,0] >= 0
 R1566: alpha[0,0] + beta[8,0] + gamma[1,8] - 3 a[1,8,0,0] >= 0
 R1567: alpha[0,0] + beta[8,0] + gamma[2,8] - 3 a[2,8,0,0] >= 0
 R1568: alpha[0,0] + beta[8,0] + gamma[3,8] - 3 a[3,8,0,0] >= 0
 R1569: alpha[0,0] + beta[8,0] + gamma[4,8] - 3 a[4,8,0,0] >= 0
 R1570: alpha[0,0] + beta[8,0] + gamma[5,8] - 3 a[5,8,0,0] >= 0
 R1571: alpha[0,0] + beta[8,0] + gamma[6,8] - 3 a[6,8,0,0] >= 0
 R1572: alpha[0,0] + beta[8,0] + gamma[7,8] - 3 a[7,8,0,0] >= 0
 R1573: alpha[0,0] + beta[8,0] + gamma[8,8] - 3 a[8,8,0,0] >= 0
 R1574: alpha[0,0] + beta[9,0] + gamma[0,9] - 3 a[0,9,0,0] >= 0
 R1575: alpha[0,0] + beta[9,0] + gamma[1,9] - 3 a[1,9,0,0] >= 0
 R1576: alpha[0,0] + beta[9,0] + gamma[2,9] - 3 a[2,9,0,0] >= 0
 R1577: alpha[0,0] + beta[9,0] + gamma[3,9] - 3 a[3,9,0,0] >= 0
 R1578: alpha[0,0] + beta[9,0] + gamma[4,9] - 3 a[4,9,0,0] >= 0
 R1579: alpha[0,0] + beta[9,0] + gamma[5,9] - 3 a[5,9,0,0] >= 0
 R1580: alpha[0,0] + beta[9,0] + gamma[6,9] - 3 a[6,9,0,0] >= 0
 R1581: alpha[0,0] + beta[9,0] + gamma[7,9] - 3 a[7,9,0,0] >= 0
 R1582: alpha[0,0] + beta[9,0] + gamma[8,9] - 3 a[8,9,0,0] >= 0
 R1583: alpha[0,0] + beta[9,0] + gamma[9,9] - 3 a[9,9,0,0] >= 0
 R1584: alpha[0,2] + beta[0,0] + gamma[0,0] - 3 a[0,0,0,2] >= 0
 R1585: alpha[0,2] + beta[1,0] + gamma[0,1] - 3 a[0,1,0,2] >= 0
 R1586: alpha[0,2] + beta[1,0] + gamma[1,1] - 3 a[1,1,0,2] >= 0
 R1587: alpha[0,2] + beta[2,0] + gamma[0,2] - 3 a[0,2,0,2] >= 0
 R1588: alpha[0,2] + beta[2,0] + gamma[1,2] - 3 a[1,2,0,2] >= 0
 R1589: alpha[0,2] + beta[2,0] + gamma[2,2] - 3 a[2,2,0,2] >= 0
 R1590: alpha[0,2] + beta[3,0] + gamma[0,3] - 3 a[0,3,0,2] >= 0
 R1591: alpha[0,2] + beta[3,0] + gamma[1,3] - 3 a[1,3,0,2] >= 0
 R1592: alpha[0,2] + beta[3,0] + gamma[2,3] - 3 a[2,3,0,2] >= 0
 R1593: alpha[0,2] + beta[3,0] + gamma[3,3] - 3 a[3,3,0,2] >= 0
 R1594: alpha[0,2] + beta[4,0] + gamma[0,4] - 3 a[0,4,0,2] >= 0
 R1595: alpha[0,2] + beta[4,0] + gamma[1,4] - 3 a[1,4,0,2] >= 0
 R1596: alpha[0,2] + beta[4,0] + gamma[2,4] - 3 a[2,4,0,2] >= 0
 R1597: alpha[0,2] + beta[4,0] + gamma[3,4] - 3 a[3,4,0,2] >= 0
 R1598: alpha[0,2] + beta[4,0] + gamma[4,4] - 3 a[4,4,0,2] >= 0
 R1599: alpha[0,2] + beta[5,0] + gamma[0,5] - 3 a[0,5,0,2] >= 0
 R1600: alpha[0,2] + beta[5,0] + gamma[1,5] - 3 a[1,5,0,2] >= 0
 R1601: alpha[0,2] + beta[5,0] + gamma[2,5] - 3 a[2,5,0,2] >= 0
 R1602: alpha[0,2] + beta[5,0] + gamma[3,5] - 3 a[3,5,0,2] >= 0
 R1603: alpha[0,2] + beta[5,0] + gamma[4,5] - 3 a[4,5,0,2] >= 0
 R1604: alpha[0,2] + beta[5,0] + gamma[5,5] - 3 a[5,5,0,2] >= 0
 R1605: alpha[0,2] + beta[6,0] + gamma[0,6] - 3 a[0,6,0,2] >= 0
 R1606: alpha[0,2] + beta[6,0] + gamma[1,6] - 3 a[1,6,0,2] >= 0
 R1607: alpha[0,2] + beta[6,0] + gamma[2,6] - 3 a[2,6,0,2] >= 0
 R1608: alpha[0,2] + beta[6,0] + gamma[3,6] - 3 a[3,6,0,2] >= 0
 R1609: alpha[0,2] + beta[6,0] + gamma[4,6] - 3 a[4,6,0,2] >= 0
 R1610: alpha[0,2] + beta[6,0] + gamma[5,6] - 3 a[5,6,0,2] >= 0
 R1611: alpha[0,2] + beta[6,0] + gamma[6,6] - 3 a[6,6,0,2] >= 0
 R1612: alpha[0,2] + beta[7,0] + gamma[0,7] - 3 a[0,7,0,2] >= 0
 R1613: alpha[0,2] + beta[7,0] + gamma[1,7] - 3 a[1,7,0,2] >= 0
 R1614: alpha[0,2] + beta[7,0] + gamma[2,7] - 3 a[2,7,0,2] >= 0
 R1615: alpha[0,2] + beta[7,0] + gamma[3,7] - 3 a[3,7,0,2] >= 0
 R1616: alpha[0,2] + beta[7,0] + gamma[4,7] - 3 a[4,7,0,2] >= 0
 R1617: alpha[0,2] + beta[7,0] + gamma[5,7] - 3 a[5,7,0,2] >= 0
 R1618: alpha[0,2] + beta[7,0] + gamma[6,7] - 3 a[6,7,0,2] >= 0
 R1619: alpha[0,2] + beta[7,0] + gamma[7,7] - 3 a[7,7,0,2] >= 0
 R1620: alpha[0,2] + beta[8,0] + gamma[0,8] - 3 a[0,8,0,2] >= 0
 R1621: alpha[0,2] + beta[8,0] + gamma[1,8] - 3 a[1,8,0,2] >= 0
 R1622: alpha[0,2] + beta[8,0] + gamma[2,8] - 3 a[2,8,0,2] >= 0
 R1623: alpha[0,2] + beta[8,0] + gamma[3,8] - 3 a[3,8,0,2] >= 0
 R1624: alpha[0,2] + beta[8,0] + gamma[4,8] - 3 a[4,8,0,2] >= 0
 R1625: alpha[0,2] + beta[8,0] + gamma[5,8] - 3 a[5,8,0,2] >= 0
 R1626: alpha[0,2] + beta[8,0] + gamma[6,8] - 3 a[6,8,0,2] >= 0
 R1627: alpha[0,2] + beta[8,0] + gamma[7,8] - 3 a[7,8,0,2] >= 0
 R1628: alpha[0,2] + beta[8,0] + gamma[8,8] - 3 a[8,8,0,2] >= 0
 R1629: alpha[0,2] + beta[9,0] + gamma[0,9] - 3 a[0,9,0,2] >= 0
 R1630: alpha[0,2] + beta[9,0] + gamma[1,9] - 3 a[1,9,0,2] >= 0
 R1631: alpha[0,2] + beta[9,0] + gamma[2,9] - 3 a[2,9,0,2] >= 0
 R1632: alpha[0,2] + beta[9,0] + gamma[3,9] - 3 a[3,9,0,2] >= 0
 R1633: alpha[0,2] + beta[9,0] + gamma[4,9] - 3 a[4,9,0,2] >= 0
 R1634: alpha[0,2] + beta[9,0] + gamma[5,9] - 3 a[5,9,0,2] >= 0
 R1635: alpha[0,2] + beta[9,0] + gamma[6,9] - 3 a[6,9,0,2] >= 0
 R1636: alpha[0,2] + beta[9,0] + gamma[7,9] - 3 a[7,9,0,2] >= 0
 R1637: alpha[0,2] + beta[9,0] + gamma[8,9] - 3 a[8,9,0,2] >= 0
 R1638: alpha[0,2] + beta[9,0] + gamma[9,9] - 3 a[9,9,0,2] >= 0
 R1639: alpha[0,6] + beta[0,0] + gamma[0,0] - 3 a[0,0,0,6] >= 0
 R1640: alpha[0,6] + beta[1,0] + gamma[0,1] - 3 a[0,1,0,6] >= 0
 R1641: alpha[0,6] + beta[1,0] + gamma[1,1] - 3 a[1,1,0,6] >= 0
 R1642: alpha[0,6] + beta[2,0] + gamma[0,2] - 3 a[0,2,0,6] >= 0
 R1643: alpha[0,6] + beta[2,0] + gamma[1,2] - 3 a[1,2,0,6] >= 0
 R1644: alpha[0,6] + beta[2,0] + gamma[2,2] - 3 a[2,2,0,6] >= 0
 R1645: alpha[0,6] + beta[3,0] + gamma[0,3] - 3 a[0,3,0,6] >= 0
 R1646: alpha[0,6] + beta[3,0] + gamma[1,3] - 3 a[1,3,0,6] >= 0
 R1647: alpha[0,6] + beta[3,0] + gamma[2,3] - 3 a[2,3,0,6] >= 0
 R1648: alpha[0,6] + beta[3,0] + gamma[3,3] - 3 a[3,3,0,6] >= 0
 R1649: alpha[0,6] + beta[4,0] + gamma[0,4] - 3 a[0,4,0,6] >= 0
 R1650: alpha[0,6] + beta[4,0] + gamma[1,4] - 3 a[1,4,0,6] >= 0
 R1651: alpha[0,6] + beta[4,0] + gamma[2,4] - 3 a[2,4,0,6] >= 0
 R1652: alpha[0,6] + beta[4,0] + gamma[3,4] - 3 a[3,4,0,6] >= 0
 R1653: alpha[0,6] + beta[4,0] + gamma[4,4] - 3 a[4,4,0,6] >= 0
 R1654: alpha[0,6] + beta[5,0] + gamma[0,5] - 3 a[0,5,0,6] >= 0
 R1655: alpha[0,6] + beta[5,0] + gamma[1,5] - 3 a[1,5,0,6] >= 0
 R1656: alpha[0,6] + beta[5,0] + gamma[2,5] - 3 a[2,5,0,6] >= 0
 R1657: alpha[0,6] + beta[5,0] + gamma[3,5] - 3 a[3,5,0,6] >= 0
 R1658: alpha[0,6] + beta[5,0] + gamma[4,5] - 3 a[4,5,0,6] >= 0
 R1659: alpha[0,6] + beta[5,0] + gamma[5,5] - 3 a[5,5,0,6] >= 0
 R1660: alpha[0,6] + beta[6,0] + gamma[0,6] - 3 a[0,6,0,6] >= 0
 R1661: alpha[0,6] + beta[6,0] + gamma[1,6] - 3 a[1,6,0,6] >= 0
 R1662: alpha[0,6] + beta[6,0] + gamma[2,6] - 3 a[2,6,0,6] >= 0
 R1663: alpha[0,6] + beta[6,0] + gamma[3,6] - 3 a[3,6,0,6] >= 0
 R1664: alpha[0,6] + beta[6,0] + gamma[4,6] - 3 a[4,6,0,6] >= 0
 R1665: alpha[0,6] + beta[6,0] + gamma[5,6] - 3 a[5,6,0,6] >= 0
 R1666: alpha[0,6] + beta[6,0] + gamma[6,6] - 3 a[6,6,0,6] >= 0
 R1667: alpha[0,6] + beta[7,0] + gamma[0,7] - 3 a[0,7,0,6] >= 0
 R1668: alpha[0,6] + beta[7,0] + gamma[1,7] - 3 a[1,7,0,6] >= 0
 R1669: alpha[0,6] + beta[7,0] + gamma[2,7] - 3 a[2,7,0,6] >= 0
 R1670: alpha[0,6] + beta[7,0] + gamma[3,7] - 3 a[3,7,0,6] >= 0
 R1671: alpha[0,6] + beta[7,0] + gamma[4,7] - 3 a[4,7,0,6] >= 0
 R1672: alpha[0,6] + beta[7,0] + gamma[5,7] - 3 a[5,7,0,6] >= 0
 R1673: alpha[0,6] + beta[7,0] + gamma[6,7] - 3 a[6,7,0,6] >= 0
 R1674: alpha[0,6] + beta[7,0] + gamma[7,7] - 3 a[7,7,0,6] >= 0
 R1675: alpha[0,6] + beta[8,0] + gamma[0,8] - 3 a[0,8,0,6] >= 0
 R1676: alpha[0,6] + beta[8,0] + gamma[1,8] - 3 a[1,8,0,6] >= 0
 R1677: alpha[0,6] + beta[8,0] + gamma[2,8] - 3 a[2,8,0,6] >= 0
 R1678: alpha[0,6] + beta[8,0] + gamma[3,8] - 3 a[3,8,0,6] >= 0
 R1679: alpha[0,6] + beta[8,0] + gamma[4,8] - 3 a[4,8,0,6] >= 0
 R1680: alpha[0,6] + beta[8,0] + gamma[5,8] - 3 a[5,8,0,6] >= 0
 R1681: alpha[0,6] + beta[8,0] + gamma[6,8] - 3 a[6,8,0,6] >= 0
 R1682: alpha[0,6] + beta[8,0] + gamma[7,8] - 3 a[7,8,0,6] >= 0
 R1683: alpha[0,6] + beta[8,0] + gamma[8,8] - 3 a[8,8,0,6] >= 0
 R1684: alpha[0,6] + beta[9,0] + gamma[0,9] - 3 a[0,9,0,6] >= 0
 R1685: alpha[0,6] + beta[9,0] + gamma[1,9] - 3 a[1,9,0,6] >= 0
 R1686: alpha[0,6] + beta[9,0] + gamma[2,9] - 3 a[2,9,0,6] >= 0
 R1687: alpha[0,6] + beta[9,0] + gamma[3,9] - 3 a[3,9,0,6] >= 0
 R1688: alpha[0,6] + beta[9,0] + gamma[4,9] - 3 a[4,9,0,6] >= 0
 R1689: alpha[0,6] + beta[9,0] + gamma[5,9] - 3 a[5,9,0,6] >= 0
 R1690: alpha[0,6] + beta[9,0] + gamma[6,9] - 3 a[6,9,0,6] >= 0
 R1691: alpha[0,6] + beta[9,0] + gamma[7,9] - 3 a[7,9,0,6] >= 0
 R1692: alpha[0,6] + beta[9,0] + gamma[8,9] - 3 a[8,9,0,6] >= 0
 R1693: alpha[0,6] + beta[9,0] + gamma[9,9] - 3 a[9,9,0,6] >= 0
 R1694: alpha[0,8] + beta[0,0] + gamma[0,0] - 3 a[0,0,0,8] >= 0
 R1695: alpha[0,8] + beta[1,0] + gamma[0,1] - 3 a[0,1,0,8] >= 0
 R1696: alpha[0,8] + beta[1,0] + gamma[1,1] - 3 a[1,1,0,8] >= 0
 R1697: alpha[0,8] + beta[2,0] + gamma[0,2] - 3 a[0,2,0,8] >= 0
 R1698: alpha[0,8] + beta[2,0] + gamma[1,2] - 3 a[1,2,0,8] >= 0
 R1699: alpha[0,8] + beta[2,0] + gamma[2,2] - 3 a[2,2,0,8] >= 0
 R1700: alpha[0,8] + beta[3,0] + gamma[0,3] - 3 a[0,3,0,8] >= 0
 R1701: alpha[0,8] + beta[3,0] + gamma[1,3] - 3 a[1,3,0,8] >= 0
 R1702: alpha[0,8] + beta[3,0] + gamma[2,3] - 3 a[2,3,0,8] >= 0
 R1703: alpha[0,8] + beta[3,0] + gamma[3,3] - 3 a[3,3,0,8] >= 0
 R1704: alpha[0,8] + beta[4,0] + gamma[0,4] - 3 a[0,4,0,8] >= 0
 R1705: alpha[0,8] + beta[4,0] + gamma[1,4] - 3 a[1,4,0,8] >= 0
 R1706: alpha[0,8] + beta[4,0] + gamma[2,4] - 3 a[2,4,0,8] >= 0
 R1707: alpha[0,8] + beta[4,0] + gamma[3,4] - 3 a[3,4,0,8] >= 0
 R1708: alpha[0,8] + beta[4,0] + gamma[4,4] - 3 a[4,4,0,8] >= 0
 R1709: alpha[0,8] + beta[5,0] + gamma[0,5] - 3 a[0,5,0,8] >= 0
 R1710: alpha[0,8] + beta[5,0] + gamma[1,5] - 3 a[1,5,0,8] >= 0
 R1711: alpha[0,8] + beta[5,0] + gamma[2,5] - 3 a[2,5,0,8] >= 0
 R1712: alpha[0,8] + beta[5,0] + gamma[3,5] - 3 a[3,5,0,8] >= 0
 R1713: alpha[0,8] + beta[5,0] + gamma[4,5] - 3 a[4,5,0,8] >= 0
 R1714: alpha[0,8] + beta[5,0] + gamma[5,5] - 3 a[5,5,0,8] >= 0
 R1715: alpha[0,8] + beta[6,0] + gamma[0,6] - 3 a[0,6,0,8] >= 0
 R1716: alpha[0,8] + beta[6,0] + gamma[1,6] - 3 a[1,6,0,8] >= 0
 R1717: alpha[0,8] + beta[6,0] + gamma[2,6] - 3 a[2,6,0,8] >= 0
 R1718: alpha[0,8] + beta[6,0] + gamma[3,6] - 3 a[3,6,0,8] >= 0
 R1719: alpha[0,8] + beta[6,0] + gamma[4,6] - 3 a[4,6,0,8] >= 0
 R1720: alpha[0,8] + beta[6,0] + gamma[5,6] - 3 a[5,6,0,8] >= 0
 R1721: alpha[0,8] + beta[6,0] + gamma[6,6] - 3 a[6,6,0,8] >= 0
 R1722: alpha[0,8] + beta[7,0] + gamma[0,7] - 3 a[0,7,0,8] >= 0
 R1723: alpha[0,8] + beta[7,0] + gamma[1,7] - 3 a[1,7,0,8] >= 0
 R1724: alpha[0,8] + beta[7,0] + gamma[2,7] - 3 a[2,7,0,8] >= 0
 R1725: alpha[0,8] + beta[7,0] + gamma[3,7] - 3 a[3,7,0,8] >= 0
 R1726: alpha[0,8] + beta[7,0] + gamma[4,7] - 3 a[4,7,0,8] >= 0
 R1727: alpha[0,8] + beta[7,0] + gamma[5,7] - 3 a[5,7,0,8] >= 0
 R1728: alpha[0,8] + beta[7,0] + gamma[6,7] - 3 a[6,7,0,8] >= 0
 R1729: alpha[0,8] + beta[7,0] + gamma[7,7] - 3 a[7,7,0,8] >= 0
 R1730: alpha[0,8] + beta[8,0] + gamma[0,8] - 3 a[0,8,0,8] >= 0
 R1731: alpha[0,8] + beta[8,0] + gamma[1,8] - 3 a[1,8,0,8] >= 0
 R1732: alpha[0,8] + beta[8,0] + gamma[2,8] - 3 a[2,8,0,8] >= 0
 R1733: alpha[0,8] + beta[8,0] + gamma[3,8] - 3 a[3,8,0,8] >= 0
 R1734: alpha[0,8] + beta[8,0] + gamma[4,8] - 3 a[4,8,0,8] >= 0
 R1735: alpha[0,8] + beta[8,0] + gamma[5,8] - 3 a[5,8,0,8] >= 0
 R1736: alpha[0,8] + beta[8,0] + gamma[6,8] - 3 a[6,8,0,8] >= 0
 R1737: alpha[0,8] + beta[8,0] + gamma[7,8] - 3 a[7,8,0,8] >= 0
 R1738: alpha[0,8] + beta[8,0] + gamma[8,8] - 3 a[8,8,0,8] >= 0
 R1739: alpha[0,8] + beta[9,0] + gamma[0,9] - 3 a[0,9,0,8] >= 0
 R1740: alpha[0,8] + beta[9,0] + gamma[1,9] - 3 a[1,9,0,8] >= 0
 R1741: alpha[0,8] + beta[9,0] + gamma[2,9] - 3 a[2,9,0,8] >= 0
 R1742: alpha[0,8] + beta[9,0] + gamma[3,9] - 3 a[3,9,0,8] >= 0
 R1743: alpha[0,8] + beta[9,0] + gamma[4,9] - 3 a[4,9,0,8] >= 0
 R1744: alpha[0,8] + beta[9,0] + gamma[5,9] - 3 a[5,9,0,8] >= 0
 R1745: alpha[0,8] + beta[9,0] + gamma[6,9] - 3 a[6,9,0,8] >= 0
 R1746: alpha[0,8] + beta[9,0] + gamma[7,9] - 3 a[7,9,0,8] >= 0
 R1747: alpha[0,8] + beta[9,0] + gamma[8,9] - 3 a[8,9,0,8] >= 0
 R1748: alpha[0,8] + beta[9,0] + gamma[9,9] - 3 a[9,9,0,8] >= 0
 R1749: alpha[0,9] + beta[0,0] + gamma[0,0] - 3 a[0,0,0,9] >= 0
 R1750: alpha[0,9] + beta[1,0] + gamma[0,1] - 3 a[0,1,0,9] >= 0
 R1751: alpha[0,9] + beta[1,0] + gamma[1,1] - 3 a[1,1,0,9] >= 0
 R1752: alpha[0,9] + beta[2,0] + gamma[0,2] - 3 a[0,2,0,9] >= 0
 R1753: alpha[0,9] + beta[2,0] + gamma[1,2] - 3 a[1,2,0,9] >= 0
 R1754: alpha[0,9] + beta[2,0] + gamma[2,2] - 3 a[2,2,0,9] >= 0
 R1755: alpha[0,9] + beta[3,0] + gamma[0,3] - 3 a[0,3,0,9] >= 0
 R1756: alpha[0,9] + beta[3,0] + gamma[1,3] - 3 a[1,3,0,9] >= 0
 R1757: alpha[0,9] + beta[3,0] + gamma[2,3] - 3 a[2,3,0,9] >= 0
 R1758: alpha[0,9] + beta[3,0] + gamma[3,3] - 3 a[3,3,0,9] >= 0
 R1759: alpha[0,9] + beta[4,0] + gamma[0,4] - 3 a[0,4,0,9] >= 0
 R1760: alpha[0,9] + beta[4,0] + gamma[1,4] - 3 a[1,4,0,9] >= 0
 R1761: alpha[0,9] + beta[4,0] + gamma[2,4] - 3 a[2,4,0,9] >= 0
 R1762: alpha[0,9] + beta[4,0] + gamma[3,4] - 3 a[3,4,0,9] >= 0
 R1763: alpha[0,9] + beta[4,0] + gamma[4,4] - 3 a[4,4,0,9] >= 0
 R1764: alpha[0,9] + beta[5,0] + gamma[0,5] - 3 a[0,5,0,9] >= 0
 R1765: alpha[0,9] + beta[5,0] + gamma[1,5] - 3 a[1,5,0,9] >= 0
 R1766: alpha[0,9] + beta[5,0] + gamma[2,5] - 3 a[2,5,0,9] >= 0
 R1767: alpha[0,9] + beta[5,0] + gamma[3,5] - 3 a[3,5,0,9] >= 0
 R1768: alpha[0,9] + beta[5,0] + gamma[4,5] - 3 a[4,5,0,9] >= 0
 R1769: alpha[0,9] + beta[5,0] + gamma[5,5] - 3 a[5,5,0,9] >= 0
 R1770: alpha[0,9] + beta[6,0] + gamma[0,6] - 3 a[0,6,0,9] >= 0
 R1771: alpha[0,9] + beta[6,0] + gamma[1,6] - 3 a[1,6,0,9] >= 0
 R1772: alpha[0,9] + beta[6,0] + gamma[2,6] - 3 a[2,6,0,9] >= 0
 R1773: alpha[0,9] + beta[6,0] + gamma[3,6] - 3 a[3,6,0,9] >= 0
 R1774: alpha[0,9] + beta[6,0] + gamma[4,6] - 3 a[4,6,0,9] >= 0
 R1775: alpha[0,9] + beta[6,0] + gamma[5,6] - 3 a[5,6,0,9] >= 0
 R1776: alpha[0,9] + beta[6,0] + gamma[6,6] - 3 a[6,6,0,9] >= 0
 R1777: alpha[0,9] + beta[7,0] + gamma[0,7] - 3 a[0,7,0,9] >= 0
 R1778: alpha[0,9] + beta[7,0] + gamma[1,7] - 3 a[1,7,0,9] >= 0
 R1779: alpha[0,9] + beta[7,0] + gamma[2,7] - 3 a[2,7,0,9] >= 0
 R1780: alpha[0,9] + beta[7,0] + gamma[3,7] - 3 a[3,7,0,9] >= 0
 R1781: alpha[0,9] + beta[7,0] + gamma[4,7] - 3 a[4,7,0,9] >= 0
 R1782: alpha[0,9] + beta[7,0] + gamma[5,7] - 3 a[5,7,0,9] >= 0
 R1783: alpha[0,9] + beta[7,0] + gamma[6,7] - 3 a[6,7,0,9] >= 0
 R1784: alpha[0,9] + beta[7,0] + gamma[7,7] - 3 a[7,7,0,9] >= 0
 R1785: alpha[0,9] + beta[8,0] + gamma[0,8] - 3 a[0,8,0,9] >= 0
 R1786: alpha[0,9] + beta[8,0] + gamma[1,8] - 3 a[1,8,0,9] >= 0
 R1787: alpha[0,9] + beta[8,0] + gamma[2,8] - 3 a[2,8,0,9] >= 0
 R1788: alpha[0,9] + beta[8,0] + gamma[3,8] - 3 a[3,8,0,9] >= 0
 R1789: alpha[0,9] + beta[8,0] + gamma[4,8] - 3 a[4,8,0,9] >= 0
 R1790: alpha[0,9] + beta[8,0] + gamma[5,8] - 3 a[5,8,0,9] >= 0
 R1791: alpha[0,9] + beta[8,0] + gamma[6,8] - 3 a[6,8,0,9] >= 0
 R1792: alpha[0,9] + beta[8,0] + gamma[7,8] - 3 a[7,8,0,9] >= 0
 R1793: alpha[0,9] + beta[8,0] + gamma[8,8] - 3 a[8,8,0,9] >= 0
 R1794: alpha[0,9] + beta[9,0] + gamma[0,9] - 3 a[0,9,0,9] >= 0
 R1795: alpha[0,9] + beta[9,0] + gamma[1,9] - 3 a[1,9,0,9] >= 0
 R1796: alpha[0,9] + beta[9,0] + gamma[2,9] - 3 a[2,9,0,9] >= 0
 R1797: alpha[0,9] + beta[9,0] + gamma[3,9] - 3 a[3,9,0,9] >= 0
 R1798: alpha[0,9] + beta[9,0] + gamma[4,9] - 3 a[4,9,0,9] >= 0
 R1799: alpha[0,9] + beta[9,0] + gamma[5,9] - 3 a[5,9,0,9] >= 0
 R1800: alpha[0,9] + beta[9,0] + gamma[6,9] - 3 a[6,9,0,9] >= 0
 R1801: alpha[0,9] + beta[9,0] + gamma[7,9] - 3 a[7,9,0,9] >= 0
 R1802: alpha[0,9] + beta[9,0] + gamma[8,9] - 3 a[8,9,0,9] >= 0
 R1803: alpha[0,9] + beta[9,0] + gamma[9,9] - 3 a[9,9,0,9] >= 0
 R1804: alpha[1,1] + beta[0,1] + gamma[0,0] - 3 a[0,0,1,1] >= 0
 R1805: alpha[1,1] + beta[1,1] + gamma[0,1] - 3 a[0,1,1,1] >= 0
 R1806: alpha[1,1] + beta[1,1] + gamma[1,1] - 3 a[1,1,1,1] >= 0
 R1807: alpha[1,1] + beta[2,1] + gamma[0,2] - 3 a[0,2,1,1] >= 0
 R1808: alpha[1,1] + beta[2,1] + gamma[1,2] - 3 a[1,2,1,1] >= 0
 R1809: alpha[1,1] + beta[2,1] + gamma[2,2] - 3 a[2,2,1,1] >= 0
 R1810: alpha[1,1] + beta[3,1] + gamma[0,3] - 3 a[0,3,1,1] >= 0
 R1811: alpha[1,1] + beta[3,1] + gamma[1,3] - 3 a[1,3,1,1] >= 0
 R1812: alpha[1,1] + beta[3,1] + gamma[2,3] - 3 a[2,3,1,1] >= 0
 R1813: alpha[1,1] + beta[3,1] + gamma[3,3] - 3 a[3,3,1,1] >= 0
 R1814: alpha[1,1] + beta[4,1] + gamma[0,4] - 3 a[0,4,1,1] >= 0
 R1815: alpha[1,1] + beta[4,1] + gamma[1,4] - 3 a[1,4,1,1] >= 0
 R1816: alpha[1,1] + beta[4,1] + gamma[2,4] - 3 a[2,4,1,1] >= 0
 R1817: alpha[1,1] + beta[4,1] + gamma[3,4] - 3 a[3,4,1,1] >= 0
 R1818: alpha[1,1] + beta[4,1] + gamma[4,4] - 3 a[4,4,1,1] >= 0
 R1819: alpha[1,1] + beta[5,1] + gamma[0,5] - 3 a[0,5,1,1] >= 0
 R1820: alpha[1,1] + beta[5,1] + gamma[1,5] - 3 a[1,5,1,1] >= 0
 R1821: alpha[1,1] + beta[5,1] + gamma[2,5] - 3 a[2,5,1,1] >= 0
 R1822: alpha[1,1] + beta[5,1] + gamma[3,5] - 3 a[3,5,1,1] >= 0
 R1823: alpha[1,1] + beta[5,1] + gamma[4,5] - 3 a[4,5,1,1] >= 0
 R1824: alpha[1,1] + beta[5,1] + gamma[5,5] - 3 a[5,5,1,1] >= 0
 R1825: alpha[1,1] + beta[6,1] + gamma[0,6] - 3 a[0,6,1,1] >= 0
 R1826: alpha[1,1] + beta[6,1] + gamma[1,6] - 3 a[1,6,1,1] >= 0
 R1827: alpha[1,1] + beta[6,1] + gamma[2,6] - 3 a[2,6,1,1] >= 0
 R1828: alpha[1,1] + beta[6,1] + gamma[3,6] - 3 a[3,6,1,1] >= 0
 R1829: alpha[1,1] + beta[6,1] + gamma[4,6] - 3 a[4,6,1,1] >= 0
 R1830: alpha[1,1] + beta[6,1] + gamma[5,6] - 3 a[5,6,1,1] >= 0
 R1831: alpha[1,1] + beta[6,1] + gamma[6,6] - 3 a[6,6,1,1] >= 0
 R1832: alpha[1,1] + beta[7,1] + gamma[0,7] - 3 a[0,7,1,1] >= 0
 R1833: alpha[1,1] + beta[7,1] + gamma[1,7] - 3 a[1,7,1,1] >= 0
 R1834: alpha[1,1] + beta[7,1] + gamma[2,7] - 3 a[2,7,1,1] >= 0
 R1835: alpha[1,1] + beta[7,1] + gamma[3,7] - 3 a[3,7,1,1] >= 0
 R1836: alpha[1,1] + beta[7,1] + gamma[4,7] - 3 a[4,7,1,1] >= 0
 R1837: alpha[1,1] + beta[7,1] + gamma[5,7] - 3 a[5,7,1,1] >= 0
 R1838: alpha[1,1] + beta[7,1] + gamma[6,7] - 3 a[6,7,1,1] >= 0
 R1839: alpha[1,1] + beta[7,1] + gamma[7,7] - 3 a[7,7,1,1] >= 0
 R1840: alpha[1,1] + beta[8,1] + gamma[0,8] - 3 a[0,8,1,1] >= 0
 R1841: alpha[1,1] + beta[8,1] + gamma[1,8] - 3 a[1,8,1,1] >= 0
 R1842: alpha[1,1] + beta[8,1] + gamma[2,8] - 3 a[2,8,1,1] >= 0
 R1843: alpha[1,1] + beta[8,1] + gamma[3,8] - 3 a[3,8,1,1] >= 0
 R1844: alpha[1,1] + beta[8,1] + gamma[4,8] - 3 a[4,8,1,1] >= 0
 R1845: alpha[1,1] + beta[8,1] + gamma[5,8] - 3 a[5,8,1,1] >= 0
 R1846: alpha[1,1] + beta[8,1] + gamma[6,8] - 3 a[6,8,1,1] >= 0
 R1847: alpha[1,1] + beta[8,1] + gamma[7,8] - 3 a[7,8,1,1] >= 0
 R1848: alpha[1,1] + beta[8,1] + gamma[8,8] - 3 a[8,8,1,1] >= 0
 R1849: alpha[1,1] + beta[9,1] + gamma[0,9] - 3 a[0,9,1,1] >= 0
 R1850: alpha[1,1] + beta[9,1] + gamma[1,9] - 3 a[1,9,1,1] >= 0
 R1851: alpha[1,1] + beta[9,1] + gamma[2,9] - 3 a[2,9,1,1] >= 0
 R1852: alpha[1,1] + beta[9,1] + gamma[3,9] - 3 a[3,9,1,1] >= 0
 R1853: alpha[1,1] + beta[9,1] + gamma[4,9] - 3 a[4,9,1,1] >= 0
 R1854: alpha[1,1] + beta[9,1] + gamma[5,9] - 3 a[5,9,1,1] >= 0
 R1855: alpha[1,1] + beta[9,1] + gamma[6,9] - 3 a[6,9,1,1] >= 0
 R1856: alpha[1,1] + beta[9,1] + gamma[7,9] - 3 a[7,9,1,1] >= 0
 R1857: alpha[1,1] + beta[9,1] + gamma[8,9] - 3 a[8,9,1,1] >= 0
 R1858: alpha[1,1] + beta[9,1] + gamma[9,9] - 3 a[9,9,1,1] >= 0
 R1859: alpha[1,3] + beta[0,1] + gamma[0,0] - 3 a[0,0,1,3] >= 0
 R1860: alpha[1,3] + beta[1,1] + gamma[0,1] - 3 a[0,1,1,3] >= 0
 R1861: alpha[1,3] + beta[1,1] + gamma[1,1] - 3 a[1,1,1,3] >= 0
 R1862: alpha[1,3] + beta[2,1] + gamma[0,2] - 3 a[0,2,1,3] >= 0
 R1863: alpha[1,3] + beta[2,1] + gamma[1,2] - 3 a[1,2,1,3] >= 0
 R1864: alpha[1,3] + beta[2,1] + gamma[2,2] - 3 a[2,2,1,3] >= 0
 R1865: alpha[1,3] + beta[3,1] + gamma[0,3] - 3 a[0,3,1,3] >= 0
 R1866: alpha[1,3] + beta[3,1] + gamma[1,3] - 3 a[1,3,1,3] >= 0
 R1867: alpha[1,3] + beta[3,1] + gamma[2,3] - 3 a[2,3,1,3] >= 0
 R1868: alpha[1,3] + beta[3,1] + gamma[3,3] - 3 a[3,3,1,3] >= 0
 R1869: alpha[1,3] + beta[4,1] + gamma[0,4] - 3 a[0,4,1,3] >= 0
 R1870: alpha[1,3] + beta[4,1] + gamma[1,4] - 3 a[1,4,1,3] >= 0
 R1871: alpha[1,3] + beta[4,1] + gamma[2,4] - 3 a[2,4,1,3] >= 0
 R1872: alpha[1,3] + beta[4,1] + gamma[3,4] - 3 a[3,4,1,3] >= 0
 R1873: alpha[1,3] + beta[4,1] + gamma[4,4] - 3 a[4,4,1,3] >= 0
 R1874: alpha[1,3] + beta[5,1] + gamma[0,5] - 3 a[0,5,1,3] >= 0
 R1875: alpha[1,3] + beta[5,1] + gamma[1,5] - 3 a[1,5,1,3] >= 0
 R1876: alpha[1,3] + beta[5,1] + gamma[2,5] - 3 a[2,5,1,3] >= 0
 R1877: alpha[1,3] + beta[5,1] + gamma[3,5] - 3 a[3,5,1,3] >= 0
 R1878: alpha[1,3] + beta[5,1] + gamma[4,5] - 3 a[4,5,1,3] >= 0
 R1879: alpha[1,3] + beta[5,1] + gamma[5,5] - 3 a[5,5,1,3] >= 0
 R1880: alpha[1,3] + beta[6,1] + gamma[0,6] - 3 a[0,6,1,3] >= 0
 R1881: alpha[1,3] + beta[6,1] + gamma[1,6] - 3 a[1,6,1,3] >= 0
 R1882: alpha[1,3] + beta[6,1] + gamma[2,6] - 3 a[2,6,1,3] >= 0
 R1883: alpha[1,3] + beta[6,1] + gamma[3,6] - 3 a[3,6,1,3] >= 0
 R1884: alpha[1,3] + beta[6,1] + gamma[4,6] - 3 a[4,6,1,3] >= 0
 R1885: alpha[1,3] + beta[6,1] + gamma[5,6] - 3 a[5,6,1,3] >= 0
 R1886: alpha[1,3] + beta[6,1] + gamma[6,6] - 3 a[6,6,1,3] >= 0
 R1887: alpha[1,3] + beta[7,1] + gamma[0,7] - 3 a[0,7,1,3] >= 0
 R1888: alpha[1,3] + beta[7,1] + gamma[1,7] - 3 a[1,7,1,3] >= 0
 R1889: alpha[1,3] + beta[7,1] + gamma[2,7] - 3 a[2,7,1,3] >= 0
 R1890: alpha[1,3] + beta[7,1] + gamma[3,7] - 3 a[3,7,1,3] >= 0
 R1891: alpha[1,3] + beta[7,1] + gamma[4,7] - 3 a[4,7,1,3] >= 0
 R1892: alpha[1,3] + beta[7,1] + gamma[5,7] - 3 a[5,7,1,3] >= 0
 R1893: alpha[1,3] + beta[7,1] + gamma[6,7] - 3 a[6,7,1,3] >= 0
 R1894: alpha[1,3] + beta[7,1] + gamma[7,7] - 3 a[7,7,1,3] >= 0
 R1895: alpha[1,3] + beta[8,1] + gamma[0,8] - 3 a[0,8,1,3] >= 0
 R1896: alpha[1,3] + beta[8,1] + gamma[1,8] - 3 a[1,8,1,3] >= 0
 R1897: alpha[1,3] + beta[8,1] + gamma[2,8] - 3 a[2,8,1,3] >= 0
 R1898: alpha[1,3] + beta[8,1] + gamma[3,8] - 3 a[3,8,1,3] >= 0
 R1899: alpha[1,3] + beta[8,1] + gamma[4,8] - 3 a[4,8,1,3] >= 0
 R1900: alpha[1,3] + beta[8,1] + gamma[5,8] - 3 a[5,8,1,3] >= 0
 R1901: alpha[1,3] + beta[8,1] + gamma[6,8] - 3 a[6,8,1,3] >= 0
 R1902: alpha[1,3] + beta[8,1] + gamma[7,8] - 3 a[7,8,1,3] >= 0
 R1903: alpha[1,3] + beta[8,1] + gamma[8,8] - 3 a[8,8,1,3] >= 0
 R1904: alpha[1,3] + beta[9,1] + gamma[0,9] - 3 a[0,9,1,3] >= 0
 R1905: alpha[1,3] + beta[9,1] + gamma[1,9] - 3 a[1,9,1,3] >= 0
 R1906: alpha[1,3] + beta[9,1] + gamma[2,9] - 3 a[2,9,1,3] >= 0
 R1907: alpha[1,3] + beta[9,1] + gamma[3,9] - 3 a[3,9,1,3] >= 0
 R1908: alpha[1,3] + beta[9,1] + gamma[4,9] - 3 a[4,9,1,3] >= 0
 R1909: alpha[1,3] + beta[9,1] + gamma[5,9] - 3 a[5,9,1,3] >= 0
 R1910: alpha[1,3] + beta[9,1] + gamma[6,9] - 3 a[6,9,1,3] >= 0
 R1911: alpha[1,3] + beta[9,1] + gamma[7,9] - 3 a[7,9,1,3] >= 0
 R1912: alpha[1,3] + beta[9,1] + gamma[8,9] - 3 a[8,9,1,3] >= 0
 R1913: alpha[1,3] + beta[9,1] + gamma[9,9] - 3 a[9,9,1,3] >= 0
 R1914: alpha[1,7] + beta[0,1] + gamma[0,0] - 3 a[0,0,1,7] >= 0
 R1915: alpha[1,7] + beta[1,1] + gamma[0,1] - 3 a[0,1,1,7] >= 0
 R1916: alpha[1,7] + beta[1,1] + gamma[1,1] - 3 a[1,1,1,7] >= 0
 R1917: alpha[1,7] + beta[2,1] + gamma[0,2] - 3 a[0,2,1,7] >= 0
 R1918: alpha[1,7] + beta[2,1] + gamma[1,2] - 3 a[1,2,1,7] >= 0
 R1919: alpha[1,7] + beta[2,1] + gamma[2,2] - 3 a[2,2,1,7] >= 0
 R1920: alpha[1,7] + beta[3,1] + gamma[0,3] - 3 a[0,3,1,7] >= 0
 R1921: alpha[1,7] + beta[3,1] + gamma[1,3] - 3 a[1,3,1,7] >= 0
 R1922: alpha[1,7] + beta[3,1] + gamma[2,3] - 3 a[2,3,1,7] >= 0
 R1923: alpha[1,7] + beta[3,1] + gamma[3,3] - 3 a[3,3,1,7] >= 0
 R1924: alpha[1,7] + beta[4,1] + gamma[0,4] - 3 a[0,4,1,7] >= 0
 R1925: alpha[1,7] + beta[4,1] + gamma[1,4] - 3 a[1,4,1,7] >= 0
 R1926: alpha[1,7] + beta[4,1] + gamma[2,4] - 3 a[2,4,1,7] >= 0
 R1927: alpha[1,7] + beta[4,1] + gamma[3,4] - 3 a[3,4,1,7] >= 0
 R1928: alpha[1,7] + beta[4,1] + gamma[4,4] - 3 a[4,4,1,7] >= 0
 R1929: alpha[1,7] + beta[5,1] + gamma[0,5] - 3 a[0,5,1,7] >= 0
 R1930: alpha[1,7] + beta[5,1] + gamma[1,5] - 3 a[1,5,1,7] >= 0
 R1931: alpha[1,7] + beta[5,1] + gamma[2,5] - 3 a[2,5,1,7] >= 0
 R1932: alpha[1,7] + beta[5,1] + gamma[3,5] - 3 a[3,5,1,7] >= 0
 R1933: alpha[1,7] + beta[5,1] + gamma[4,5] - 3 a[4,5,1,7] >= 0
 R1934: alpha[1,7] + beta[5,1] + gamma[5,5] - 3 a[5,5,1,7] >= 0
 R1935: alpha[1,7] + beta[6,1] + gamma[0,6] - 3 a[0,6,1,7] >= 0
 R1936: alpha[1,7] + beta[6,1] + gamma[1,6] - 3 a[1,6,1,7] >= 0
 R1937: alpha[1,7] + beta[6,1] + gamma[2,6] - 3 a[2,6,1,7] >= 0
 R1938: alpha[1,7] + beta[6,1] + gamma[3,6] - 3 a[3,6,1,7] >= 0
 R1939: alpha[1,7] + beta[6,1] + gamma[4,6] - 3 a[4,6,1,7] >= 0
 R1940: alpha[1,7] + beta[6,1] + gamma[5,6] - 3 a[5,6,1,7] >= 0
 R1941: alpha[1,7] + beta[6,1] + gamma[6,6] - 3 a[6,6,1,7] >= 0
 R1942: alpha[1,7] + beta[7,1] + gamma[0,7] - 3 a[0,7,1,7] >= 0
 R1943: alpha[1,7] + beta[7,1] + gamma[1,7] - 3 a[1,7,1,7] >= 0
 R1944: alpha[1,7] + beta[7,1] + gamma[2,7] - 3 a[2,7,1,7] >= 0
 R1945: alpha[1,7] + beta[7,1] + gamma[3,7] - 3 a[3,7,1,7] >= 0
 R1946: alpha[1,7] + beta[7,1] + gamma[4,7] - 3 a[4,7,1,7] >= 0
 R1947: alpha[1,7] + beta[7,1] + gamma[5,7] - 3 a[5,7,1,7] >= 0
 R1948: alpha[1,7] + beta[7,1] + gamma[6,7] - 3 a[6,7,1,7] >= 0
 R1949: alpha[1,7] + beta[7,1] + gamma[7,7] - 3 a[7,7,1,7] >= 0
 R1950: alpha[1,7] + beta[8,1] + gamma[0,8] - 3 a[0,8,1,7] >= 0
 R1951: alpha[1,7] + beta[8,1] + gamma[1,8] - 3 a[1,8,1,7] >= 0
 R1952: alpha[1,7] + beta[8,1] + gamma[2,8] - 3 a[2,8,1,7] >= 0
 R1953: alpha[1,7] + beta[8,1] + gamma[3,8] - 3 a[3,8,1,7] >= 0
 R1954: alpha[1,7] + beta[8,1] + gamma[4,8] - 3 a[4,8,1,7] >= 0
 R1955: alpha[1,7] + beta[8,1] + gamma[5,8] - 3 a[5,8,1,7] >= 0
 R1956: alpha[1,7] + beta[8,1] + gamma[6,8] - 3 a[6,8,1,7] >= 0
 R1957: alpha[1,7] + beta[8,1] + gamma[7,8] - 3 a[7,8,1,7] >= 0
 R1958: alpha[1,7] + beta[8,1] + gamma[8,8] - 3 a[8,8,1,7] >= 0
 R1959: alpha[1,7] + beta[9,1] + gamma[0,9] - 3 a[0,9,1,7] >= 0
 R1960: alpha[1,7] + beta[9,1] + gamma[1,9] - 3 a[1,9,1,7] >= 0
 R1961: alpha[1,7] + beta[9,1] + gamma[2,9] - 3 a[2,9,1,7] >= 0
 R1962: alpha[1,7] + beta[9,1] + gamma[3,9] - 3 a[3,9,1,7] >= 0
 R1963: alpha[1,7] + beta[9,1] + gamma[4,9] - 3 a[4,9,1,7] >= 0
 R1964: alpha[1,7] + beta[9,1] + gamma[5,9] - 3 a[5,9,1,7] >= 0
 R1965: alpha[1,7] + beta[9,1] + gamma[6,9] - 3 a[6,9,1,7] >= 0
 R1966: alpha[1,7] + beta[9,1] + gamma[7,9] - 3 a[7,9,1,7] >= 0
 R1967: alpha[1,7] + beta[9,1] + gamma[8,9] - 3 a[8,9,1,7] >= 0
 R1968: alpha[1,7] + beta[9,1] + gamma[9,9] - 3 a[9,9,1,7] >= 0
 R1969: alpha[2,2] + beta[0,2] + gamma[0,0] - 3 a[0,0,2,2] >= 0
 R1970: alpha[2,2] + beta[1,2] + gamma[0,1] - 3 a[0,1,2,2] >= 0
 R1971: alpha[2,2] + beta[1,2] + gamma[1,1] - 3 a[1,1,2,2] >= 0
 R1972: alpha[2,2] + beta[2,2] + gamma[0,2] - 3 a[0,2,2,2] >= 0
 R1973: alpha[2,2] + beta[2,2] + gamma[1,2] - 3 a[1,2,2,2] >= 0
 R1974: alpha[2,2] + beta[2,2] + gamma[2,2] - 3 a[2,2,2,2] >= 0
 R1975: alpha[2,2] + beta[3,2] + gamma[0,3] - 3 a[0,3,2,2] >= 0
 R1976: alpha[2,2] + beta[3,2] + gamma[1,3] - 3 a[1,3,2,2] >= 0
 R1977: alpha[2,2] + beta[3,2] + gamma[2,3] - 3 a[2,3,2,2] >= 0
 R1978: alpha[2,2] + beta[3,2] + gamma[3,3] - 3 a[3,3,2,2] >= 0
 R1979: alpha[2,2] + beta[4,2] + gamma[0,4] - 3 a[0,4,2,2] >= 0
 R1980: alpha[2,2] + beta[4,2] + gamma[1,4] - 3 a[1,4,2,2] >= 0
 R1981: alpha[2,2] + beta[4,2] + gamma[2,4] - 3 a[2,4,2,2] >= 0
 R1982: alpha[2,2] + beta[4,2] + gamma[3,4] - 3 a[3,4,2,2] >= 0
 R1983: alpha[2,2] + beta[4,2] + gamma[4,4] - 3 a[4,4,2,2] >= 0
 R1984: alpha[2,2] + beta[5,2] + gamma[0,5] - 3 a[0,5,2,2] >= 0
 R1985: alpha[2,2] + beta[5,2] + gamma[1,5] - 3 a[1,5,2,2] >= 0
 R1986: alpha[2,2] + beta[5,2] + gamma[2,5] - 3 a[2,5,2,2] >= 0
 R1987: alpha[2,2] + beta[5,2] + gamma[3,5] - 3 a[3,5,2,2] >= 0
 R1988: alpha[2,2] + beta[5,2] + gamma[4,5] - 3 a[4,5,2,2] >= 0
 R1989: alpha[2,2] + beta[5,2] + gamma[5,5] - 3 a[5,5,2,2] >= 0
 R1990: alpha[2,2] + beta[6,2] + gamma[0,6] - 3 a[0,6,2,2] >= 0
 R1991: alpha[2,2] + beta[6,2] + gamma[1,6] - 3 a[1,6,2,2] >= 0
 R1992: alpha[2,2] + beta[6,2] + gamma[2,6] - 3 a[2,6,2,2] >= 0
 R1993: alpha[2,2] + beta[6,2] + gamma[3,6] - 3 a[3,6,2,2] >= 0
 R1994: alpha[2,2] + beta[6,2] + gamma[4,6] - 3 a[4,6,2,2] >= 0
 R1995: alpha[2,2] + beta[6,2] + gamma[5,6] - 3 a[5,6,2,2] >= 0
 R1996: alpha[2,2] + beta[6,2] + gamma[6,6] - 3 a[6,6,2,2] >= 0
 R1997: alpha[2,2] + beta[7,2] + gamma[0,7] - 3 a[0,7,2,2] >= 0
 R1998: alpha[2,2] + beta[7,2] + gamma[1,7] - 3 a[1,7,2,2] >= 0
 R1999: alpha[2,2] + beta[7,2] + gamma[2,7] - 3 a[2,7,2,2] >= 0
 R2000: alpha[2,2] + beta[7,2] + gamma[3,7] - 3 a[3,7,2,2] >= 0
 R2001: alpha[2,2] + beta[7,2] + gamma[4,7] - 3 a[4,7,2,2] >= 0
 R2002: alpha[2,2] + beta[7,2] + gamma[5,7] - 3 a[5,7,2,2] >= 0
 R2003: alpha[2,2] + beta[7,2] + gamma[6,7] - 3 a[6,7,2,2] >= 0
 R2004: alpha[2,2] + beta[7,2] + gamma[7,7] - 3 a[7,7,2,2] >= 0
 R2005: alpha[2,2] + beta[8,2] + gamma[0,8] - 3 a[0,8,2,2] >= 0
 R2006: alpha[2,2] + beta[8,2] + gamma[1,8] - 3 a[1,8,2,2] >= 0
 R2007: alpha[2,2] + beta[8,2] + gamma[2,8] - 3 a[2,8,2,2] >= 0
 R2008: alpha[2,2] + beta[8,2] + gamma[3,8] - 3 a[3,8,2,2] >= 0
 R2009: alpha[2,2] + beta[8,2] + gamma[4,8] - 3 a[4,8,2,2] >= 0
 R2010: alpha[2,2] + beta[8,2] + gamma[5,8] - 3 a[5,8,2,2] >= 0
 R2011: alpha[2,2] + beta[8,2] + gamma[6,8] - 3 a[6,8,2,2] >= 0
 R2012: alpha[2,2] + beta[8,2] + gamma[7,8] - 3 a[7,8,2,2] >= 0
 R2013: alpha[2,2] + beta[8,2] + gamma[8,8] - 3 a[8,8,2,2] >= 0
 R2014: alpha[2,2] + beta[9,2] + gamma[0,9] - 3 a[0,9,2,2] >= 0
 R2015: alpha[2,2] + beta[9,2] + gamma[1,9] - 3 a[1,9,2,2] >= 0
 R2016: alpha[2,2] + beta[9,2] + gamma[2,9] - 3 a[2,9,2,2] >= 0
 R2017: alpha[2,2] + beta[9,2] + gamma[3,9] - 3 a[3,9,2,2] >= 0
 R2018: alpha[2,2] + beta[9,2] + gamma[4,9] - 3 a[4,9,2,2] >= 0
 R2019: alpha[2,2] + beta[9,2] + gamma[5,9] - 3 a[5,9,2,2] >= 0
 R2020: alpha[2,2] + beta[9,2] + gamma[6,9] - 3 a[6,9,2,2] >= 0
 R2021: alpha[2,2] + beta[9,2] + gamma[7,9] - 3 a[7,9,2,2] >= 0
 R2022: alpha[2,2] + beta[9,2] + gamma[8,9] - 3 a[8,9,2,2] >= 0
 R2023: alpha[2,2] + beta[9,2] + gamma[9,9] - 3 a[9,9,2,2] >= 0
 R2024: alpha[2,6] + beta[0,2] + gamma[0,0] - 3 a[0,0,2,6] >= 0
 R2025: alpha[2,6] + beta[1,2] + gamma[0,1] - 3 a[0,1,2,6] >= 0
 R2026: alpha[2,6] + beta[1,2] + gamma[1,1] - 3 a[1,1,2,6] >= 0
 R2027: alpha[2,6] + beta[2,2] + gamma[0,2] - 3 a[0,2,2,6] >= 0
 R2028: alpha[2,6] + beta[2,2] + gamma[1,2] - 3 a[1,2,2,6] >= 0
 R2029: alpha[2,6] + beta[2,2] + gamma[2,2] - 3 a[2,2,2,6] >= 0
 R2030: alpha[2,6] + beta[3,2] + gamma[0,3] - 3 a[0,3,2,6] >= 0
 R2031: alpha[2,6] + beta[3,2] + gamma[1,3] - 3 a[1,3,2,6] >= 0
 R2032: alpha[2,6] + beta[3,2] + gamma[2,3] - 3 a[2,3,2,6] >= 0
 R2033: alpha[2,6] + beta[3,2] + gamma[3,3] - 3 a[3,3,2,6] >= 0
 R2034: alpha[2,6] + beta[4,2] + gamma[0,4] - 3 a[0,4,2,6] >= 0
 R2035: alpha[2,6] + beta[4,2] + gamma[1,4] - 3 a[1,4,2,6] >= 0
 R2036: alpha[2,6] + beta[4,2] + gamma[2,4] - 3 a[2,4,2,6] >= 0
 R2037: alpha[2,6] + beta[4,2] + gamma[3,4] - 3 a[3,4,2,6] >= 0
 R2038: alpha[2,6] + beta[4,2] + gamma[4,4] - 3 a[4,4,2,6] >= 0
 R2039: alpha[2,6] + beta[5,2] + gamma[0,5] - 3 a[0,5,2,6] >= 0
 R2040: alpha[2,6] + beta[5,2] + gamma[1,5] - 3 a[1,5,2,6] >= 0
 R2041: alpha[2,6] + beta[5,2] + gamma[2,5] - 3 a[2,5,2,6] >= 0
 R2042: alpha[2,6] + beta[5,2] + gamma[3,5] - 3 a[3,5,2,6] >= 0
 R2043: alpha[2,6] + beta[5,2] + gamma[4,5] - 3 a[4,5,2,6] >= 0
 R2044: alpha[2,6] + beta[5,2] + gamma[5,5] - 3 a[5,5,2,6] >= 0
 R2045: alpha[2,6] + beta[6,2] + gamma[0,6] - 3 a[0,6,2,6] >= 0
 R2046: alpha[2,6] + beta[6,2] + gamma[1,6] - 3 a[1,6,2,6] >= 0
 R2047: alpha[2,6] + beta[6,2] + gamma[2,6] - 3 a[2,6,2,6] >= 0
 R2048: alpha[2,6] + beta[6,2] + gamma[3,6] - 3 a[3,6,2,6] >= 0
 R2049: alpha[2,6] + beta[6,2] + gamma[4,6] - 3 a[4,6,2,6] >= 0
 R2050: alpha[2,6] + beta[6,2] + gamma[5,6] - 3 a[5,6,2,6] >= 0
 R2051: alpha[2,6] + beta[6,2] + gamma[6,6] - 3 a[6,6,2,6] >= 0
 R2052: alpha[2,6] + beta[7,2] + gamma[0,7] - 3 a[0,7,2,6] >= 0
 R2053: alpha[2,6] + beta[7,2] + gamma[1,7] - 3 a[1,7,2,6] >= 0
 R2054: alpha[2,6] + beta[7,2] + gamma[2,7] - 3 a[2,7,2,6] >= 0
 R2055: alpha[2,6] + beta[7,2] + gamma[3,7] - 3 a[3,7,2,6] >= 0
 R2056: alpha[2,6] + beta[7,2] + gamma[4,7] - 3 a[4,7,2,6] >= 0
 R2057: alpha[2,6] + beta[7,2] + gamma[5,7] - 3 a[5,7,2,6] >= 0
 R2058: alpha[2,6] + beta[7,2] + gamma[6,7] - 3 a[6,7,2,6] >= 0
 R2059: alpha[2,6] + beta[7,2] + gamma[7,7] - 3 a[7,7,2,6] >= 0
 R2060: alpha[2,6] + beta[8,2] + gamma[0,8] - 3 a[0,8,2,6] >= 0
 R2061: alpha[2,6] + beta[8,2] + gamma[1,8] - 3 a[1,8,2,6] >= 0
 R2062: alpha[2,6] + beta[8,2] + gamma[2,8] - 3 a[2,8,2,6] >= 0
 R2063: alpha[2,6] + beta[8,2] + gamma[3,8] - 3 a[3,8,2,6] >= 0
 R2064: alpha[2,6] + beta[8,2] + gamma[4,8] - 3 a[4,8,2,6] >= 0
 R2065: alpha[2,6] + beta[8,2] + gamma[5,8] - 3 a[5,8,2,6] >= 0
 R2066: alpha[2,6] + beta[8,2] + gamma[6,8] - 3 a[6,8,2,6] >= 0
 R2067: alpha[2,6] + beta[8,2] + gamma[7,8] - 3 a[7,8,2,6] >= 0
 R2068: alpha[2,6] + beta[8,2] + gamma[8,8] - 3 a[8,8,2,6] >= 0
 R2069: alpha[2,6] + beta[9,2] + gamma[0,9] - 3 a[0,9,2,6] >= 0
 R2070: alpha[2,6] + beta[9,2] + gamma[1,9] - 3 a[1,9,2,6] >= 0
 R2071: alpha[2,6] + beta[9,2] + gamma[2,9] - 3 a[2,9,2,6] >= 0
 R2072: alpha[2,6] + beta[9,2] + gamma[3,9] - 3 a[3,9,2,6] >= 0
 R2073: alpha[2,6] + beta[9,2] + gamma[4,9] - 3 a[4,9,2,6] >= 0
 R2074: alpha[2,6] + beta[9,2] + gamma[5,9] - 3 a[5,9,2,6] >= 0
 R2075: alpha[2,6] + beta[9,2] + gamma[6,9] - 3 a[6,9,2,6] >= 0
 R2076: alpha[2,6] + beta[9,2] + gamma[7,9] - 3 a[7,9,2,6] >= 0
 R2077: alpha[2,6] + beta[9,2] + gamma[8,9] - 3 a[8,9,2,6] >= 0
 R2078: alpha[2,6] + beta[9,2] + gamma[9,9] - 3 a[9,9,2,6] >= 0
 R2079: alpha[2,8] + beta[0,2] + gamma[0,0] - 3 a[0,0,2,8] >= 0
 R2080: alpha[2,8] + beta[1,2] + gamma[0,1] - 3 a[0,1,2,8] >= 0
 R2081: alpha[2,8] + beta[1,2] + gamma[1,1] - 3 a[1,1,2,8] >= 0
 R2082: alpha[2,8] + beta[2,2] + gamma[0,2] - 3 a[0,2,2,8] >= 0
 R2083: alpha[2,8] + beta[2,2] + gamma[1,2] - 3 a[1,2,2,8] >= 0
 R2084: alpha[2,8] + beta[2,2] + gamma[2,2] - 3 a[2,2,2,8] >= 0
 R2085: alpha[2,8] + beta[3,2] + gamma[0,3] - 3 a[0,3,2,8] >= 0
 R2086: alpha[2,8] + beta[3,2] + gamma[1,3] - 3 a[1,3,2,8] >= 0
 R2087: alpha[2,8] + beta[3,2] + gamma[2,3] - 3 a[2,3,2,8] >= 0
 R2088: alpha[2,8] + beta[3,2] + gamma[3,3] - 3 a[3,3,2,8] >= 0
 R2089: alpha[2,8] + beta[4,2] + gamma[0,4] - 3 a[0,4,2,8] >= 0
 R2090: alpha[2,8] + beta[4,2] + gamma[1,4] - 3 a[1,4,2,8] >= 0
 R2091: alpha[2,8] + beta[4,2] + gamma[2,4] - 3 a[2,4,2,8] >= 0
 R2092: alpha[2,8] + beta[4,2] + gamma[3,4] - 3 a[3,4,2,8] >= 0
 R2093: alpha[2,8] + beta[4,2] + gamma[4,4] - 3 a[4,4,2,8] >= 0
 R2094: alpha[2,8] + beta[5,2] + gamma[0,5] - 3 a[0,5,2,8] >= 0
 R2095: alpha[2,8] + beta[5,2] + gamma[1,5] - 3 a[1,5,2,8] >= 0
 R2096: alpha[2,8] + beta[5,2] + gamma[2,5] - 3 a[2,5,2,8] >= 0
 R2097: alpha[2,8] + beta[5,2] + gamma[3,5] - 3 a[3,5,2,8] >= 0
 R2098: alpha[2,8] + beta[5,2] + gamma[4,5] - 3 a[4,5,2,8] >= 0
 R2099: alpha[2,8] + beta[5,2] + gamma[5,5] - 3 a[5,5,2,8] >= 0
 R2100: alpha[2,8] + beta[6,2] + gamma[0,6] - 3 a[0,6,2,8] >= 0
 R2101: alpha[2,8] + beta[6,2] + gamma[1,6] - 3 a[1,6,2,8] >= 0
 R2102: alpha[2,8] + beta[6,2] + gamma[2,6] - 3 a[2,6,2,8] >= 0
 R2103: alpha[2,8] + beta[6,2] + gamma[3,6] - 3 a[3,6,2,8] >= 0
 R2104: alpha[2,8] + beta[6,2] + gamma[4,6] - 3 a[4,6,2,8] >= 0
 R2105: alpha[2,8] + beta[6,2] + gamma[5,6] - 3 a[5,6,2,8] >= 0
 R2106: alpha[2,8] + beta[6,2] + gamma[6,6] - 3 a[6,6,2,8] >= 0
 R2107: alpha[2,8] + beta[7,2] + gamma[0,7] - 3 a[0,7,2,8] >= 0
 R2108: alpha[2,8] + beta[7,2] + gamma[1,7] - 3 a[1,7,2,8] >= 0
 R2109: alpha[2,8] + beta[7,2] + gamma[2,7] - 3 a[2,7,2,8] >= 0
 R2110: alpha[2,8] + beta[7,2] + gamma[3,7] - 3 a[3,7,2,8] >= 0
 R2111: alpha[2,8] + beta[7,2] + gamma[4,7] - 3 a[4,7,2,8] >= 0
 R2112: alpha[2,8] + beta[7,2] + gamma[5,7] - 3 a[5,7,2,8] >= 0
 R2113: alpha[2,8] + beta[7,2] + gamma[6,7] - 3 a[6,7,2,8] >= 0
 R2114: alpha[2,8] + beta[7,2] + gamma[7,7] - 3 a[7,7,2,8] >= 0
 R2115: alpha[2,8] + beta[8,2] + gamma[0,8] - 3 a[0,8,2,8] >= 0
 R2116: alpha[2,8] + beta[8,2] + gamma[1,8] - 3 a[1,8,2,8] >= 0
 R2117: alpha[2,8] + beta[8,2] + gamma[2,8] - 3 a[2,8,2,8] >= 0
 R2118: alpha[2,8] + beta[8,2] + gamma[3,8] - 3 a[3,8,2,8] >= 0
 R2119: alpha[2,8] + beta[8,2] + gamma[4,8] - 3 a[4,8,2,8] >= 0
 R2120: alpha[2,8] + beta[8,2] + gamma[5,8] - 3 a[5,8,2,8] >= 0
 R2121: alpha[2,8] + beta[8,2] + gamma[6,8] - 3 a[6,8,2,8] >= 0
 R2122: alpha[2,8] + beta[8,2] + gamma[7,8] - 3 a[7,8,2,8] >= 0
 R2123: alpha[2,8] + beta[8,2] + gamma[8,8] - 3 a[8,8,2,8] >= 0
 R2124: alpha[2,8] + beta[9,2] + gamma[0,9] - 3 a[0,9,2,8] >= 0
 R2125: alpha[2,8] + beta[9,2] + gamma[1,9] - 3 a[1,9,2,8] >= 0
 R2126: alpha[2,8] + beta[9,2] + gamma[2,9] - 3 a[2,9,2,8] >= 0
 R2127: alpha[2,8] + beta[9,2] + gamma[3,9] - 3 a[3,9,2,8] >= 0
 R2128: alpha[2,8] + beta[9,2] + gamma[4,9] - 3 a[4,9,2,8] >= 0
 R2129: alpha[2,8] + beta[9,2] + gamma[5,9] - 3 a[5,9,2,8] >= 0
 R2130: alpha[2,8] + beta[9,2] + gamma[6,9] - 3 a[6,9,2,8] >= 0
 R2131: alpha[2,8] + beta[9,2] + gamma[7,9] - 3 a[7,9,2,8] >= 0
 R2132: alpha[2,8] + beta[9,2] + gamma[8,9] - 3 a[8,9,2,8] >= 0
 R2133: alpha[2,8] + beta[9,2] + gamma[9,9] - 3 a[9,9,2,8] >= 0
 R2134: alpha[2,9] + beta[0,2] + gamma[0,0] - 3 a[0,0,2,9] >= 0
 R2135: alpha[2,9] + beta[1,2] + gamma[0,1] - 3 a[0,1,2,9] >= 0
 R2136: alpha[2,9] + beta[1,2] + gamma[1,1] - 3 a[1,1,2,9] >= 0
 R2137: alpha[2,9] + beta[2,2] + gamma[0,2] - 3 a[0,2,2,9] >= 0
 R2138: alpha[2,9] + beta[2,2] + gamma[1,2] - 3 a[1,2,2,9] >= 0
 R2139: alpha[2,9] + beta[2,2] + gamma[2,2] - 3 a[2,2,2,9] >= 0
 R2140: alpha[2,9] + beta[3,2] + gamma[0,3] - 3 a[0,3,2,9] >= 0
 R2141: alpha[2,9] + beta[3,2] + gamma[1,3] - 3 a[1,3,2,9] >= 0
 R2142: alpha[2,9] + beta[3,2] + gamma[2,3] - 3 a[2,3,2,9] >= 0
 R2143: alpha[2,9] + beta[3,2] + gamma[3,3] - 3 a[3,3,2,9] >= 0
 R2144: alpha[2,9] + beta[4,2] + gamma[0,4] - 3 a[0,4,2,9] >= 0
 R2145: alpha[2,9] + beta[4,2] + gamma[1,4] - 3 a[1,4,2,9] >= 0
 R2146: alpha[2,9] + beta[4,2] + gamma[2,4] - 3 a[2,4,2,9] >= 0
 R2147: alpha[2,9] + beta[4,2] + gamma[3,4] - 3 a[3,4,2,9] >= 0
 R2148: alpha[2,9] + beta[4,2] + gamma[4,4] - 3 a[4,4,2,9] >= 0
 R2149: alpha[2,9] + beta[5,2] + gamma[0,5] - 3 a[0,5,2,9] >= 0
 R2150: alpha[2,9] + beta[5,2] + gamma[1,5] - 3 a[1,5,2,9] >= 0
 R2151: alpha[2,9] + beta[5,2] + gamma[2,5] - 3 a[2,5,2,9] >= 0
 R2152: alpha[2,9] + beta[5,2] + gamma[3,5] - 3 a[3,5,2,9] >= 0
 R2153: alpha[2,9] + beta[5,2] + gamma[4,5] - 3 a[4,5,2,9] >= 0
 R2154: alpha[2,9] + beta[5,2] + gamma[5,5] - 3 a[5,5,2,9] >= 0
 R2155: alpha[2,9] + beta[6,2] + gamma[0,6] - 3 a[0,6,2,9] >= 0
 R2156: alpha[2,9] + beta[6,2] + gamma[1,6] - 3 a[1,6,2,9] >= 0
 R2157: alpha[2,9] + beta[6,2] + gamma[2,6] - 3 a[2,6,2,9] >= 0
 R2158: alpha[2,9] + beta[6,2] + gamma[3,6] - 3 a[3,6,2,9] >= 0
 R2159: alpha[2,9] + beta[6,2] + gamma[4,6] - 3 a[4,6,2,9] >= 0
 R2160: alpha[2,9] + beta[6,2] + gamma[5,6] - 3 a[5,6,2,9] >= 0
 R2161: alpha[2,9] + beta[6,2] + gamma[6,6] - 3 a[6,6,2,9] >= 0
 R2162: alpha[2,9] + beta[7,2] + gamma[0,7] - 3 a[0,7,2,9] >= 0
 R2163: alpha[2,9] + beta[7,2] + gamma[1,7] - 3 a[1,7,2,9] >= 0
 R2164: alpha[2,9] + beta[7,2] + gamma[2,7] - 3 a[2,7,2,9] >= 0
 R2165: alpha[2,9] + beta[7,2] + gamma[3,7] - 3 a[3,7,2,9] >= 0
 R2166: alpha[2,9] + beta[7,2] + gamma[4,7] - 3 a[4,7,2,9] >= 0
 R2167: alpha[2,9] + beta[7,2] + gamma[5,7] - 3 a[5,7,2,9] >= 0
 R2168: alpha[2,9] + beta[7,2] + gamma[6,7] - 3 a[6,7,2,9] >= 0
 R2169: alpha[2,9] + beta[7,2] + gamma[7,7] - 3 a[7,7,2,9] >= 0
 R2170: alpha[2,9] + beta[8,2] + gamma[0,8] - 3 a[0,8,2,9] >= 0
 R2171: alpha[2,9] + beta[8,2] + gamma[1,8] - 3 a[1,8,2,9] >= 0
 R2172: alpha[2,9] + beta[8,2] + gamma[2,8] - 3 a[2,8,2,9] >= 0
 R2173: alpha[2,9] + beta[8,2] + gamma[3,8] - 3 a[3,8,2,9] >= 0
 R2174: alpha[2,9] + beta[8,2] + gamma[4,8] - 3 a[4,8,2,9] >= 0
 R2175: alpha[2,9] + beta[8,2] + gamma[5,8] - 3 a[5,8,2,9] >= 0
 R2176: alpha[2,9] + beta[8,2] + gamma[6,8] - 3 a[6,8,2,9] >= 0
 R2177: alpha[2,9] + beta[8,2] + gamma[7,8] - 3 a[7,8,2,9] >= 0
 R2178: alpha[2,9] + beta[8,2] + gamma[8,8] - 3 a[8,8,2,9] >= 0
 R2179: alpha[2,9] + beta[9,2] + gamma[0,9] - 3 a[0,9,2,9] >= 0
 R2180: alpha[2,9] + beta[9,2] + gamma[1,9] - 3 a[1,9,2,9] >= 0
 R2181: alpha[2,9] + beta[9,2] + gamma[2,9] - 3 a[2,9,2,9] >= 0
 R2182: alpha[2,9] + beta[9,2] + gamma[3,9] - 3 a[3,9,2,9] >= 0
 R2183: alpha[2,9] + beta[9,2] + gamma[4,9] - 3 a[4,9,2,9] >= 0
 R2184: alpha[2,9] + beta[9,2] + gamma[5,9] - 3 a[5,9,2,9] >= 0
 R2185: alpha[2,9] + beta[9,2] + gamma[6,9] - 3 a[6,9,2,9] >= 0
 R2186: alpha[2,9] + beta[9,2] + gamma[7,9] - 3 a[7,9,2,9] >= 0
 R2187: alpha[2,9] + beta[9,2] + gamma[8,9] - 3 a[8,9,2,9] >= 0
 R2188: alpha[2,9] + beta[9,2] + gamma[9,9] - 3 a[9,9,2,9] >= 0
 R2189: alpha[3,3] + beta[0,3] + gamma[0,0] - 3 a[0,0,3,3] >= 0
 R2190: alpha[3,3] + beta[1,3] + gamma[0,1] - 3 a[0,1,3,3] >= 0
 R2191: alpha[3,3] + beta[1,3] + gamma[1,1] - 3 a[1,1,3,3] >= 0
 R2192: alpha[3,3] + beta[2,3] + gamma[0,2] - 3 a[0,2,3,3] >= 0
 R2193: alpha[3,3] + beta[2,3] + gamma[1,2] - 3 a[1,2,3,3] >= 0
 R2194: alpha[3,3] + beta[2,3] + gamma[2,2] - 3 a[2,2,3,3] >= 0
 R2195: alpha[3,3] + beta[3,3] + gamma[0,3] - 3 a[0,3,3,3] >= 0
 R2196: alpha[3,3] + beta[3,3] + gamma[1,3] - 3 a[1,3,3,3] >= 0
 R2197: alpha[3,3] + beta[3,3] + gamma[2,3] - 3 a[2,3,3,3] >= 0
 R2198: alpha[3,3] + beta[3,3] + gamma[3,3] - 3 a[3,3,3,3] >= 0
 R2199: alpha[3,3] + beta[4,3] + gamma[0,4] - 3 a[0,4,3,3] >= 0
 R2200: alpha[3,3] + beta[4,3] + gamma[1,4] - 3 a[1,4,3,3] >= 0
 R2201: alpha[3,3] + beta[4,3] + gamma[2,4] - 3 a[2,4,3,3] >= 0
 R2202: alpha[3,3] + beta[4,3] + gamma[3,4] - 3 a[3,4,3,3] >= 0
 R2203: alpha[3,3] + beta[4,3] + gamma[4,4] - 3 a[4,4,3,3] >= 0
 R2204: alpha[3,3] + beta[5,3] + gamma[0,5] - 3 a[0,5,3,3] >= 0
 R2205: alpha[3,3] + beta[5,3] + gamma[1,5] - 3 a[1,5,3,3] >= 0
 R2206: alpha[3,3] + beta[5,3] + gamma[2,5] - 3 a[2,5,3,3] >= 0
 R2207: alpha[3,3] + beta[5,3] + gamma[3,5] - 3 a[3,5,3,3] >= 0
 R2208: alpha[3,3] + beta[5,3] + gamma[4,5] - 3 a[4,5,3,3] >= 0
 R2209: alpha[3,3] + beta[5,3] + gamma[5,5] - 3 a[5,5,3,3] >= 0
 R2210: alpha[3,3] + beta[6,3] + gamma[0,6] - 3 a[0,6,3,3] >= 0
 R2211: alpha[3,3] + beta[6,3] + gamma[1,6] - 3 a[1,6,3,3] >= 0
 R2212: alpha[3,3] + beta[6,3] + gamma[2,6] - 3 a[2,6,3,3] >= 0
 R2213: alpha[3,3] + beta[6,3] + gamma[3,6] - 3 a[3,6,3,3] >= 0
 R2214: alpha[3,3] + beta[6,3] + gamma[4,6] - 3 a[4,6,3,3] >= 0
 R2215: alpha[3,3] + beta[6,3] + gamma[5,6] - 3 a[5,6,3,3] >= 0
 R2216: alpha[3,3] + beta[6,3] + gamma[6,6] - 3 a[6,6,3,3] >= 0
 R2217: alpha[3,3] + beta[7,3] + gamma[0,7] - 3 a[0,7,3,3] >= 0
 R2218: alpha[3,3] + beta[7,3] + gamma[1,7] - 3 a[1,7,3,3] >= 0
 R2219: alpha[3,3] + beta[7,3] + gamma[2,7] - 3 a[2,7,3,3] >= 0
 R2220: alpha[3,3] + beta[7,3] + gamma[3,7] - 3 a[3,7,3,3] >= 0
 R2221: alpha[3,3] + beta[7,3] + gamma[4,7] - 3 a[4,7,3,3] >= 0
 R2222: alpha[3,3] + beta[7,3] + gamma[5,7] - 3 a[5,7,3,3] >= 0
 R2223: alpha[3,3] + beta[7,3] + gamma[6,7] - 3 a[6,7,3,3] >= 0
 R2224: alpha[3,3] + beta[7,3] + gamma[7,7] - 3 a[7,7,3,3] >= 0
 R2225: alpha[3,3] + beta[8,3] + gamma[0,8] - 3 a[0,8,3,3] >= 0
 R2226: alpha[3,3] + beta[8,3] + gamma[1,8] - 3 a[1,8,3,3] >= 0
 R2227: alpha[3,3] + beta[8,3] + gamma[2,8] - 3 a[2,8,3,3] >= 0
 R2228: alpha[3,3] + beta[8,3] + gamma[3,8] - 3 a[3,8,3,3] >= 0
 R2229: alpha[3,3] + beta[8,3] + gamma[4,8] - 3 a[4,8,3,3] >= 0
 R2230: alpha[3,3] + beta[8,3] + gamma[5,8] - 3 a[5,8,3,3] >= 0
 R2231: alpha[3,3] + beta[8,3] + gamma[6,8] - 3 a[6,8,3,3] >= 0
 R2232: alpha[3,3] + beta[8,3] + gamma[7,8] - 3 a[7,8,3,3] >= 0
 R2233: alpha[3,3] + beta[8,3] + gamma[8,8] - 3 a[8,8,3,3] >= 0
 R2234: alpha[3,3] + beta[9,3] + gamma[0,9] - 3 a[0,9,3,3] >= 0
 R2235: alpha[3,3] + beta[9,3] + gamma[1,9] - 3 a[1,9,3,3] >= 0
 R2236: alpha[3,3] + beta[9,3] + gamma[2,9] - 3 a[2,9,3,3] >= 0
 R2237: alpha[3,3] + beta[9,3] + gamma[3,9] - 3 a[3,9,3,3] >= 0
 R2238: alpha[3,3] + beta[9,3] + gamma[4,9] - 3 a[4,9,3,3] >= 0
 R2239: alpha[3,3] + beta[9,3] + gamma[5,9] - 3 a[5,9,3,3] >= 0
 R2240: alpha[3,3] + beta[9,3] + gamma[6,9] - 3 a[6,9,3,3] >= 0
 R2241: alpha[3,3] + beta[9,3] + gamma[7,9] - 3 a[7,9,3,3] >= 0
 R2242: alpha[3,3] + beta[9,3] + gamma[8,9] - 3 a[8,9,3,3] >= 0
 R2243: alpha[3,3] + beta[9,3] + gamma[9,9] - 3 a[9,9,3,3] >= 0
 R2244: alpha[3,7] + beta[0,3] + gamma[0,0] - 3 a[0,0,3,7] >= 0
 R2245: alpha[3,7] + beta[1,3] + gamma[0,1] - 3 a[0,1,3,7] >= 0
 R2246: alpha[3,7] + beta[1,3] + gamma[1,1] - 3 a[1,1,3,7] >= 0
 R2247: alpha[3,7] + beta[2,3] + gamma[0,2] - 3 a[0,2,3,7] >= 0
 R2248: alpha[3,7] + beta[2,3] + gamma[1,2] - 3 a[1,2,3,7] >= 0
 R2249: alpha[3,7] + beta[2,3] + gamma[2,2] - 3 a[2,2,3,7] >= 0
 R2250: alpha[3,7] + beta[3,3] + gamma[0,3] - 3 a[0,3,3,7] >= 0
 R2251: alpha[3,7] + beta[3,3] + gamma[1,3] - 3 a[1,3,3,7] >= 0
 R2252: alpha[3,7] + beta[3,3] + gamma[2,3] - 3 a[2,3,3,7] >= 0
 R2253: alpha[3,7] + beta[3,3] + gamma[3,3] - 3 a[3,3,3,7] >= 0
 R2254: alpha[3,7] + beta[4,3] + gamma[0,4] - 3 a[0,4,3,7] >= 0
 R2255: alpha[3,7] + beta[4,3] + gamma[1,4] - 3 a[1,4,3,7] >= 0
 R2256: alpha[3,7] + beta[4,3] + gamma[2,4] - 3 a[2,4,3,7] >= 0
 R2257: alpha[3,7] + beta[4,3] + gamma[3,4] - 3 a[3,4,3,7] >= 0
 R2258: alpha[3,7] + beta[4,3] + gamma[4,4] - 3 a[4,4,3,7] >= 0
 R2259: alpha[3,7] + beta[5,3] + gamma[0,5] - 3 a[0,5,3,7] >= 0
 R2260: alpha[3,7] + beta[5,3] + gamma[1,5] - 3 a[1,5,3,7] >= 0
 R2261: alpha[3,7] + beta[5,3] + gamma[2,5] - 3 a[2,5,3,7] >= 0
 R2262: alpha[3,7] + beta[5,3] + gamma[3,5] - 3 a[3,5,3,7] >= 0
 R2263: alpha[3,7] + beta[5,3] + gamma[4,5] - 3 a[4,5,3,7] >= 0
 R2264: alpha[3,7] + beta[5,3] + gamma[5,5] - 3 a[5,5,3,7] >= 0
 R2265: alpha[3,7] + beta[6,3] + gamma[0,6] - 3 a[0,6,3,7] >= 0
 R2266: alpha[3,7] + beta[6,3] + gamma[1,6] - 3 a[1,6,3,7] >= 0
 R2267: alpha[3,7] + beta[6,3] + gamma[2,6] - 3 a[2,6,3,7] >= 0
 R2268: alpha[3,7] + beta[6,3] + gamma[3,6] - 3 a[3,6,3,7] >= 0
 R2269: alpha[3,7] + beta[6,3] + gamma[4,6] - 3 a[4,6,3,7] >= 0
 R2270: alpha[3,7] + beta[6,3] + gamma[5,6] - 3 a[5,6,3,7] >= 0
 R2271: alpha[3,7] + beta[6,3] + gamma[6,6] - 3 a[6,6,3,7] >= 0
 R2272: alpha[3,7] + beta[7,3] + gamma[0,7] - 3 a[0,7,3,7] >= 0
 R2273: alpha[3,7] + beta[7,3] + gamma[1,7] - 3 a[1,7,3,7] >= 0
 R2274: alpha[3,7] + beta[7,3] + gamma[2,7] - 3 a[2,7,3,7] >= 0
 R2275: alpha[3,7] + beta[7,3] + gamma[3,7] - 3 a[3,7,3,7] >= 0
 R2276: alpha[3,7] + beta[7,3] + gamma[4,7] - 3 a[4,7,3,7] >= 0
 R2277: alpha[3,7] + beta[7,3] + gamma[5,7] - 3 a[5,7,3,7] >= 0
 R2278: alpha[3,7] + beta[7,3] + gamma[6,7] - 3 a[6,7,3,7] >= 0
 R2279: alpha[3,7] + beta[7,3] + gamma[7,7] - 3 a[7,7,3,7] >= 0
 R2280: alpha[3,7] + beta[8,3] + gamma[0,8] - 3 a[0,8,3,7] >= 0
 R2281: alpha[3,7] + beta[8,3] + gamma[1,8] - 3 a[1,8,3,7] >= 0
 R2282: alpha[3,7] + beta[8,3] + gamma[2,8] - 3 a[2,8,3,7] >= 0
 R2283: alpha[3,7] + beta[8,3] + gamma[3,8] - 3 a[3,8,3,7] >= 0
 R2284: alpha[3,7] + beta[8,3] + gamma[4,8] - 3 a[4,8,3,7] >= 0
 R2285: alpha[3,7] + beta[8,3] + gamma[5,8] - 3 a[5,8,3,7] >= 0
 R2286: alpha[3,7] + beta[8,3] + gamma[6,8] - 3 a[6,8,3,7] >= 0
 R2287: alpha[3,7] + beta[8,3] + gamma[7,8] - 3 a[7,8,3,7] >= 0
 R2288: alpha[3,7] + beta[8,3] + gamma[8,8] - 3 a[8,8,3,7] >= 0
 R2289: alpha[3,7] + beta[9,3] + gamma[0,9] - 3 a[0,9,3,7] >= 0
 R2290: alpha[3,7] + beta[9,3] + gamma[1,9] - 3 a[1,9,3,7] >= 0
 R2291: alpha[3,7] + beta[9,3] + gamma[2,9] - 3 a[2,9,3,7] >= 0
 R2292: alpha[3,7] + beta[9,3] + gamma[3,9] - 3 a[3,9,3,7] >= 0
 R2293: alpha[3,7] + beta[9,3] + gamma[4,9] - 3 a[4,9,3,7] >= 0
 R2294: alpha[3,7] + beta[9,3] + gamma[5,9] - 3 a[5,9,3,7] >= 0
 R2295: alpha[3,7] + beta[9,3] + gamma[6,9] - 3 a[6,9,3,7] >= 0
 R2296: alpha[3,7] + beta[9,3] + gamma[7,9] - 3 a[7,9,3,7] >= 0
 R2297: alpha[3,7] + beta[9,3] + gamma[8,9] - 3 a[8,9,3,7] >= 0
 R2298: alpha[3,7] + beta[9,3] + gamma[9,9] - 3 a[9,9,3,7] >= 0
 R2299: alpha[4,4] + beta[0,4] + gamma[0,0] - 3 a[0,0,4,4] >= 0
 R2300: alpha[4,4] + beta[1,4] + gamma[0,1] - 3 a[0,1,4,4] >= 0
 R2301: alpha[4,4] + beta[1,4] + gamma[1,1] - 3 a[1,1,4,4] >= 0
 R2302: alpha[4,4] + beta[2,4] + gamma[0,2] - 3 a[0,2,4,4] >= 0
 R2303: alpha[4,4] + beta[2,4] + gamma[1,2] - 3 a[1,2,4,4] >= 0
 R2304: alpha[4,4] + beta[2,4] + gamma[2,2] - 3 a[2,2,4,4] >= 0
 R2305: alpha[4,4] + beta[3,4] + gamma[0,3] - 3 a[0,3,4,4] >= 0
 R2306: alpha[4,4] + beta[3,4] + gamma[1,3] - 3 a[1,3,4,4] >= 0
 R2307: alpha[4,4] + beta[3,4] + gamma[2,3] - 3 a[2,3,4,4] >= 0
 R2308: alpha[4,4] + beta[3,4] + gamma[3,3] - 3 a[3,3,4,4] >= 0
 R2309: alpha[4,4] + beta[4,4] + gamma[0,4] - 3 a[0,4,4,4] >= 0
 R2310: alpha[4,4] + beta[4,4] + gamma[1,4] - 3 a[1,4,4,4] >= 0
 R2311: alpha[4,4] + beta[4,4] + gamma[2,4] - 3 a[2,4,4,4] >= 0
 R2312: alpha[4,4] + beta[4,4] + gamma[3,4] - 3 a[3,4,4,4] >= 0
 R2313: alpha[4,4] + beta[4,4] + gamma[4,4] - 3 a[4,4,4,4] >= 0
 R2314: alpha[4,4] + beta[5,4] + gamma[0,5] - 3 a[0,5,4,4] >= 0
 R2315: alpha[4,4] + beta[5,4] + gamma[1,5] - 3 a[1,5,4,4] >= 0
 R2316: alpha[4,4] + beta[5,4] + gamma[2,5] - 3 a[2,5,4,4] >= 0
 R2317: alpha[4,4] + beta[5,4] + gamma[3,5] - 3 a[3,5,4,4] >= 0
 R2318: alpha[4,4] + beta[5,4] + gamma[4,5] - 3 a[4,5,4,4] >= 0
 R2319: alpha[4,4] + beta[5,4] + gamma[5,5] - 3 a[5,5,4,4] >= 0
 R2320: alpha[4,4] + beta[6,4] + gamma[0,6] - 3 a[0,6,4,4] >= 0
 R2321: alpha[4,4] + beta[6,4] + gamma[1,6] - 3 a[1,6,4,4] >= 0
 R2322: alpha[4,4] + beta[6,4] + gamma[2,6] - 3 a[2,6,4,4] >= 0
 R2323: alpha[4,4] + beta[6,4] + gamma[3,6] - 3 a[3,6,4,4] >= 0
 R2324: alpha[4,4] + beta[6,4] + gamma[4,6] - 3 a[4,6,4,4] >= 0
 R2325: alpha[4,4] + beta[6,4] + gamma[5,6] - 3 a[5,6,4,4] >= 0
 R2326: alpha[4,4] + beta[6,4] + gamma[6,6] - 3 a[6,6,4,4] >= 0
 R2327: alpha[4,4] + beta[7,4] + gamma[0,7] - 3 a[0,7,4,4] >= 0
 R2328: alpha[4,4] + beta[7,4] + gamma[1,7] - 3 a[1,7,4,4] >= 0
 R2329: alpha[4,4] + beta[7,4] + gamma[2,7] - 3 a[2,7,4,4] >= 0
 R2330: alpha[4,4] + beta[7,4] + gamma[3,7] - 3 a[3,7,4,4] >= 0
 R2331: alpha[4,4] + beta[7,4] + gamma[4,7] - 3 a[4,7,4,4] >= 0
 R2332: alpha[4,4] + beta[7,4] + gamma[5,7] - 3 a[5,7,4,4] >= 0
 R2333: alpha[4,4] + beta[7,4] + gamma[6,7] - 3 a[6,7,4,4] >= 0
 R2334: alpha[4,4] + beta[7,4] + gamma[7,7] - 3 a[7,7,4,4] >= 0
 R2335: alpha[4,4] + beta[8,4] + gamma[0,8] - 3 a[0,8,4,4] >= 0
 R2336: alpha[4,4] + beta[8,4] + gamma[1,8] - 3 a[1,8,4,4] >= 0
 R2337: alpha[4,4] + beta[8,4] + gamma[2,8] - 3 a[2,8,4,4] >= 0
 R2338: alpha[4,4] + beta[8,4] + gamma[3,8] - 3 a[3,8,4,4] >= 0
 R2339: alpha[4,4] + beta[8,4] + gamma[4,8] - 3 a[4,8,4,4] >= 0
 R2340: alpha[4,4] + beta[8,4] + gamma[5,8] - 3 a[5,8,4,4] >= 0
 R2341: alpha[4,4] + beta[8,4] + gamma[6,8] - 3 a[6,8,4,4] >= 0
 R2342: alpha[4,4] + beta[8,4] + gamma[7,8] - 3 a[7,8,4,4] >= 0
 R2343: alpha[4,4] + beta[8,4] + gamma[8,8] - 3 a[8,8,4,4] >= 0
 R2344: alpha[4,4] + beta[9,4] + gamma[0,9] - 3 a[0,9,4,4] >= 0
 R2345: alpha[4,4] + beta[9,4] + gamma[1,9] - 3 a[1,9,4,4] >= 0
 R2346: alpha[4,4] + beta[9,4] + gamma[2,9] - 3 a[2,9,4,4] >= 0
 R2347: alpha[4,4] + beta[9,4] + gamma[3,9] - 3 a[3,9,4,4] >= 0
 R2348: alpha[4,4] + beta[9,4] + gamma[4,9] - 3 a[4,9,4,4] >= 0
 R2349: alpha[4,4] + beta[9,4] + gamma[5,9] - 3 a[5,9,4,4] >= 0
 R2350: alpha[4,4] + beta[9,4] + gamma[6,9] - 3 a[6,9,4,4] >= 0
 R2351: alpha[4,4] + beta[9,4] + gamma[7,9] - 3 a[7,9,4,4] >= 0
 R2352: alpha[4,4] + beta[9,4] + gamma[8,9] - 3 a[8,9,4,4] >= 0
 R2353: alpha[4,4] + beta[9,4] + gamma[9,9] - 3 a[9,9,4,4] >= 0
 R2354: alpha[5,5] + beta[0,5] + gamma[0,0] - 3 a[0,0,5,5] >= 0
 R2355: alpha[5,5] + beta[1,5] + gamma[0,1] - 3 a[0,1,5,5] >= 0
 R2356: alpha[5,5] + beta[1,5] + gamma[1,1] - 3 a[1,1,5,5] >= 0
 R2357: alpha[5,5] + beta[2,5] + gamma[0,2] - 3 a[0,2,5,5] >= 0
 R2358: alpha[5,5] + beta[2,5] + gamma[1,2] - 3 a[1,2,5,5] >= 0
 R2359: alpha[5,5] + beta[2,5] + gamma[2,2] - 3 a[2,2,5,5] >= 0
 R2360: alpha[5,5] + beta[3,5] + gamma[0,3] - 3 a[0,3,5,5] >= 0
 R2361: alpha[5,5] + beta[3,5] + gamma[1,3] - 3 a[1,3,5,5] >= 0
 R2362: alpha[5,5] + beta[3,5] + gamma[2,3] - 3 a[2,3,5,5] >= 0
 R2363: alpha[5,5] + beta[3,5] + gamma[3,3] - 3 a[3,3,5,5] >= 0
 R2364: alpha[5,5] + beta[4,5] + gamma[0,4] - 3 a[0,4,5,5] >= 0
 R2365: alpha[5,5] + beta[4,5] + gamma[1,4] - 3 a[1,4,5,5] >= 0
 R2366: alpha[5,5] + beta[4,5] + gamma[2,4] - 3 a[2,4,5,5] >= 0
 R2367: alpha[5,5] + beta[4,5] + gamma[3,4] - 3 a[3,4,5,5] >= 0
 R2368: alpha[5,5] + beta[4,5] + gamma[4,4] - 3 a[4,4,5,5] >= 0
 R2369: alpha[5,5] + beta[5,5] + gamma[0,5] - 3 a[0,5,5,5] >= 0
 R2370: alpha[5,5] + beta[5,5] + gamma[1,5] - 3 a[1,5,5,5] >= 0
 R2371: alpha[5,5] + beta[5,5] + gamma[2,5] - 3 a[2,5,5,5] >= 0
 R2372: alpha[5,5] + beta[5,5] + gamma[3,5] - 3 a[3,5,5,5] >= 0
 R2373: alpha[5,5] + beta[5,5] + gamma[4,5] - 3 a[4,5,5,5] >= 0
 R2374: alpha[5,5] + beta[5,5] + gamma[5,5] - 3 a[5,5,5,5] >= 0
 R2375: alpha[5,5] + beta[6,5] + gamma[0,6] - 3 a[0,6,5,5] >= 0
 R2376: alpha[5,5] + beta[6,5] + gamma[1,6] - 3 a[1,6,5,5] >= 0
 R2377: alpha[5,5] + beta[6,5] + gamma[2,6] - 3 a[2,6,5,5] >= 0
 R2378: alpha[5,5] + beta[6,5] + gamma[3,6] - 3 a[3,6,5,5] >= 0
 R2379: alpha[5,5] + beta[6,5] + gamma[4,6] - 3 a[4,6,5,5] >= 0
 R2380: alpha[5,5] + beta[6,5] + gamma[5,6] - 3 a[5,6,5,5] >= 0
 R2381: alpha[5,5] + beta[6,5] + gamma[6,6] - 3 a[6,6,5,5] >= 0
 R2382: alpha[5,5] + beta[7,5] + gamma[0,7] - 3 a[0,7,5,5] >= 0
 R2383: alpha[5,5] + beta[7,5] + gamma[1,7] - 3 a[1,7,5,5] >= 0
 R2384: alpha[5,5] + beta[7,5] + gamma[2,7] - 3 a[2,7,5,5] >= 0
 R2385: alpha[5,5] + beta[7,5] + gamma[3,7] - 3 a[3,7,5,5] >= 0
 R2386: alpha[5,5] + beta[7,5] + gamma[4,7] - 3 a[4,7,5,5] >= 0
 R2387: alpha[5,5] + beta[7,5] + gamma[5,7] - 3 a[5,7,5,5] >= 0
 R2388: alpha[5,5] + beta[7,5] + gamma[6,7] - 3 a[6,7,5,5] >= 0
 R2389: alpha[5,5] + beta[7,5] + gamma[7,7] - 3 a[7,7,5,5] >= 0
 R2390: alpha[5,5] + beta[8,5] + gamma[0,8] - 3 a[0,8,5,5] >= 0
 R2391: alpha[5,5] + beta[8,5] + gamma[1,8] - 3 a[1,8,5,5] >= 0
 R2392: alpha[5,5] + beta[8,5] + gamma[2,8] - 3 a[2,8,5,5] >= 0
 R2393: alpha[5,5] + beta[8,5] + gamma[3,8] - 3 a[3,8,5,5] >= 0
 R2394: alpha[5,5] + beta[8,5] + gamma[4,8] - 3 a[4,8,5,5] >= 0
 R2395: alpha[5,5] + beta[8,5] + gamma[5,8] - 3 a[5,8,5,5] >= 0
 R2396: alpha[5,5] + beta[8,5] + gamma[6,8] - 3 a[6,8,5,5] >= 0
 R2397: alpha[5,5] + beta[8,5] + gamma[7,8] - 3 a[7,8,5,5] >= 0
 R2398: alpha[5,5] + beta[8,5] + gamma[8,8] - 3 a[8,8,5,5] >= 0
 R2399: alpha[5,5] + beta[9,5] + gamma[0,9] - 3 a[0,9,5,5] >= 0
 R2400: alpha[5,5] + beta[9,5] + gamma[1,9] - 3 a[1,9,5,5] >= 0
 R2401: alpha[5,5] + beta[9,5] + gamma[2,9] - 3 a[2,9,5,5] >= 0
 R2402: alpha[5,5] + beta[9,5] + gamma[3,9] - 3 a[3,9,5,5] >= 0
 R2403: alpha[5,5] + beta[9,5] + gamma[4,9] - 3 a[4,9,5,5] >= 0
 R2404: alpha[5,5] + beta[9,5] + gamma[5,9] - 3 a[5,9,5,5] >= 0
 R2405: alpha[5,5] + beta[9,5] + gamma[6,9] - 3 a[6,9,5,5] >= 0
 R2406: alpha[5,5] + beta[9,5] + gamma[7,9] - 3 a[7,9,5,5] >= 0
 R2407: alpha[5,5] + beta[9,5] + gamma[8,9] - 3 a[8,9,5,5] >= 0
 R2408: alpha[5,5] + beta[9,5] + gamma[9,9] - 3 a[9,9,5,5] >= 0
 R2409: alpha[6,6] + beta[0,6] + gamma[0,0] - 3 a[0,0,6,6] >= 0
 R2410: alpha[6,6] + beta[1,6] + gamma[0,1] - 3 a[0,1,6,6] >= 0
 R2411: alpha[6,6] + beta[1,6] + gamma[1,1] - 3 a[1,1,6,6] >= 0
 R2412: alpha[6,6] + beta[2,6] + gamma[0,2] - 3 a[0,2,6,6] >= 0
 R2413: alpha[6,6] + beta[2,6] + gamma[1,2] - 3 a[1,2,6,6] >= 0
 R2414: alpha[6,6] + beta[2,6] + gamma[2,2] - 3 a[2,2,6,6] >= 0
 R2415: alpha[6,6] + beta[3,6] + gamma[0,3] - 3 a[0,3,6,6] >= 0
 R2416: alpha[6,6] + beta[3,6] + gamma[1,3] - 3 a[1,3,6,6] >= 0
 R2417: alpha[6,6] + beta[3,6] + gamma[2,3] - 3 a[2,3,6,6] >= 0
 R2418: alpha[6,6] + beta[3,6] + gamma[3,3] - 3 a[3,3,6,6] >= 0
 R2419: alpha[6,6] + beta[4,6] + gamma[0,4] - 3 a[0,4,6,6] >= 0
 R2420: alpha[6,6] + beta[4,6] + gamma[1,4] - 3 a[1,4,6,6] >= 0
 R2421: alpha[6,6] + beta[4,6] + gamma[2,4] - 3 a[2,4,6,6] >= 0
 R2422: alpha[6,6] + beta[4,6] + gamma[3,4] - 3 a[3,4,6,6] >= 0
 R2423: alpha[6,6] + beta[4,6] + gamma[4,4] - 3 a[4,4,6,6] >= 0
 R2424: alpha[6,6] + beta[5,6] + gamma[0,5] - 3 a[0,5,6,6] >= 0
 R2425: alpha[6,6] + beta[5,6] + gamma[1,5] - 3 a[1,5,6,6] >= 0
 R2426: alpha[6,6] + beta[5,6] + gamma[2,5] - 3 a[2,5,6,6] >= 0
 R2427: alpha[6,6] + beta[5,6] + gamma[3,5] - 3 a[3,5,6,6] >= 0
 R2428: alpha[6,6] + beta[5,6] + gamma[4,5] - 3 a[4,5,6,6] >= 0
 R2429: alpha[6,6] + beta[5,6] + gamma[5,5] - 3 a[5,5,6,6] >= 0
 R2430: alpha[6,6] + beta[6,6] + gamma[0,6] - 3 a[0,6,6,6] >= 0
 R2431: alpha[6,6] + beta[6,6] + gamma[1,6] - 3 a[1,6,6,6] >= 0
 R2432: alpha[6,6] + beta[6,6] + gamma[2,6] - 3 a[2,6,6,6] >= 0
 R2433: alpha[6,6] + beta[6,6] + gamma[3,6] - 3 a[3,6,6,6] >= 0
 R2434: alpha[6,6] + beta[6,6] + gamma[4,6] - 3 a[4,6,6,6] >= 0
 R2435: alpha[6,6] + beta[6,6] + gamma[5,6] - 3 a[5,6,6,6] >= 0
 R2436: alpha[6,6] + beta[6,6] + gamma[6,6] - 3 a[6,6,6,6] >= 0
 R2437: alpha[6,6] + beta[7,6] + gamma[0,7] - 3 a[0,7,6,6] >= 0
 R2438: alpha[6,6] + beta[7,6] + gamma[1,7] - 3 a[1,7,6,6] >= 0
 R2439: alpha[6,6] + beta[7,6] + gamma[2,7] - 3 a[2,7,6,6] >= 0
 R2440: alpha[6,6] + beta[7,6] + gamma[3,7] - 3 a[3,7,6,6] >= 0
 R2441: alpha[6,6] + beta[7,6] + gamma[4,7] - 3 a[4,7,6,6] >= 0
 R2442: alpha[6,6] + beta[7,6] + gamma[5,7] - 3 a[5,7,6,6] >= 0
 R2443: alpha[6,6] + beta[7,6] + gamma[6,7] - 3 a[6,7,6,6] >= 0
 R2444: alpha[6,6] + beta[7,6] + gamma[7,7] - 3 a[7,7,6,6] >= 0
 R2445: alpha[6,6] + beta[8,6] + gamma[0,8] - 3 a[0,8,6,6] >= 0
 R2446: alpha[6,6] + beta[8,6] + gamma[1,8] - 3 a[1,8,6,6] >= 0
 R2447: alpha[6,6] + beta[8,6] + gamma[2,8] - 3 a[2,8,6,6] >= 0
 R2448: alpha[6,6] + beta[8,6] + gamma[3,8] - 3 a[3,8,6,6] >= 0
 R2449: alpha[6,6] + beta[8,6] + gamma[4,8] - 3 a[4,8,6,6] >= 0
 R2450: alpha[6,6] + beta[8,6] + gamma[5,8] - 3 a[5,8,6,6] >= 0
 R2451: alpha[6,6] + beta[8,6] + gamma[6,8] - 3 a[6,8,6,6] >= 0
 R2452: alpha[6,6] + beta[8,6] + gamma[7,8] - 3 a[7,8,6,6] >= 0
 R2453: alpha[6,6] + beta[8,6] + gamma[8,8] - 3 a[8,8,6,6] >= 0
 R2454: alpha[6,6] + beta[9,6] + gamma[0,9] - 3 a[0,9,6,6] >= 0
 R2455: alpha[6,6] + beta[9,6] + gamma[1,9] - 3 a[1,9,6,6] >= 0
 R2456: alpha[6,6] + beta[9,6] + gamma[2,9] - 3 a[2,9,6,6] >= 0
 R2457: alpha[6,6] + beta[9,6] + gamma[3,9] - 3 a[3,9,6,6] >= 0
 R2458: alpha[6,6] + beta[9,6] + gamma[4,9] - 3 a[4,9,6,6] >= 0
 R2459: alpha[6,6] + beta[9,6] + gamma[5,9] - 3 a[5,9,6,6] >= 0
 R2460: alpha[6,6] + beta[9,6] + gamma[6,9] - 3 a[6,9,6,6] >= 0
 R2461: alpha[6,6] + beta[9,6] + gamma[7,9] - 3 a[7,9,6,6] >= 0
 R2462: alpha[6,6] + beta[9,6] + gamma[8,9] - 3 a[8,9,6,6] >= 0
 R2463: alpha[6,6] + beta[9,6] + gamma[9,9] - 3 a[9,9,6,6] >= 0
 R2464: alpha[6,8] + beta[0,6] + gamma[0,0] - 3 a[0,0,6,8] >= 0
 R2465: alpha[6,8] + beta[1,6] + gamma[0,1] - 3 a[0,1,6,8] >= 0
 R2466: alpha[6,8] + beta[1,6] + gamma[1,1] - 3 a[1,1,6,8] >= 0
 R2467: alpha[6,8] + beta[2,6] + gamma[0,2] - 3 a[0,2,6,8] >= 0
 R2468: alpha[6,8] + beta[2,6] + gamma[1,2] - 3 a[1,2,6,8] >= 0
 R2469: alpha[6,8] + beta[2,6] + gamma[2,2] - 3 a[2,2,6,8] >= 0
 R2470: alpha[6,8] + beta[3,6] + gamma[0,3] - 3 a[0,3,6,8] >= 0
 R2471: alpha[6,8] + beta[3,6] + gamma[1,3] - 3 a[1,3,6,8] >= 0
 R2472: alpha[6,8] + beta[3,6] + gamma[2,3] - 3 a[2,3,6,8] >= 0
 R2473: alpha[6,8] + beta[3,6] + gamma[3,3] - 3 a[3,3,6,8] >= 0
 R2474: alpha[6,8] + beta[4,6] + gamma[0,4] - 3 a[0,4,6,8] >= 0
 R2475: alpha[6,8] + beta[4,6] + gamma[1,4] - 3 a[1,4,6,8] >= 0
 R2476: alpha[6,8] + beta[4,6] + gamma[2,4] - 3 a[2,4,6,8] >= 0
 R2477: alpha[6,8] + beta[4,6] + gamma[3,4] - 3 a[3,4,6,8] >= 0
 R2478: alpha[6,8] + beta[4,6] + gamma[4,4] - 3 a[4,4,6,8] >= 0
 R2479: alpha[6,8] + beta[5,6] + gamma[0,5] - 3 a[0,5,6,8] >= 0
 R2480: alpha[6,8] + beta[5,6] + gamma[1,5] - 3 a[1,5,6,8] >= 0
 R2481: alpha[6,8] + beta[5,6] + gamma[2,5] - 3 a[2,5,6,8] >= 0
 R2482: alpha[6,8] + beta[5,6] + gamma[3,5] - 3 a[3,5,6,8] >= 0
 R2483: alpha[6,8] + beta[5,6] + gamma[4,5] - 3 a[4,5,6,8] >= 0
 R2484: alpha[6,8] + beta[5,6] + gamma[5,5] - 3 a[5,5,6,8] >= 0
 R2485: alpha[6,8] + beta[6,6] + gamma[0,6] - 3 a[0,6,6,8] >= 0
 R2486: alpha[6,8] + beta[6,6] + gamma[1,6] - 3 a[1,6,6,8] >= 0
 R2487: alpha[6,8] + beta[6,6] + gamma[2,6] - 3 a[2,6,6,8] >= 0
 R2488: alpha[6,8] + beta[6,6] + gamma[3,6] - 3 a[3,6,6,8] >= 0
 R2489: alpha[6,8] + beta[6,6] + gamma[4,6] - 3 a[4,6,6,8] >= 0
 R2490: alpha[6,8] + beta[6,6] + gamma[5,6] - 3 a[5,6,6,8] >= 0
 R2491: alpha[6,8] + beta[6,6] + gamma[6,6] - 3 a[6,6,6,8] >= 0
 R2492: alpha[6,8] + beta[7,6] + gamma[0,7] - 3 a[0,7,6,8] >= 0
 R2493: alpha[6,8] + beta[7,6] + gamma[1,7] - 3 a[1,7,6,8] >= 0
 R2494: alpha[6,8] + beta[7,6] + gamma[2,7] - 3 a[2,7,6,8] >= 0
 R2495: alpha[6,8] + beta[7,6] + gamma[3,7] - 3 a[3,7,6,8] >= 0
 R2496: alpha[6,8] + beta[7,6] + gamma[4,7] - 3 a[4,7,6,8] >= 0
 R2497: alpha[6,8] + beta[7,6] + gamma[5,7] - 3 a[5,7,6,8] >= 0
 R2498: alpha[6,8] + beta[7,6] + gamma[6,7] - 3 a[6,7,6,8] >= 0
 R2499: alpha[6,8] + beta[7,6] + gamma[7,7] - 3 a[7,7,6,8] >= 0
 R2500: alpha[6,8] + beta[8,6] + gamma[0,8] - 3 a[0,8,6,8] >= 0
 R2501: alpha[6,8] + beta[8,6] + gamma[1,8] - 3 a[1,8,6,8] >= 0
 R2502: alpha[6,8] + beta[8,6] + gamma[2,8] - 3 a[2,8,6,8] >= 0
 R2503: alpha[6,8] + beta[8,6] + gamma[3,8] - 3 a[3,8,6,8] >= 0
 R2504: alpha[6,8] + beta[8,6] + gamma[4,8] - 3 a[4,8,6,8] >= 0
 R2505: alpha[6,8] + beta[8,6] + gamma[5,8] - 3 a[5,8,6,8] >= 0
 R2506: alpha[6,8] + beta[8,6] + gamma[6,8] - 3 a[6,8,6,8] >= 0
 R2507: alpha[6,8] + beta[8,6] + gamma[7,8] - 3 a[7,8,6,8] >= 0
 R2508: alpha[6,8] + beta[8,6] + gamma[8,8] - 3 a[8,8,6,8] >= 0
 R2509: alpha[6,8] + beta[9,6] + gamma[0,9] - 3 a[0,9,6,8] >= 0
 R2510: alpha[6,8] + beta[9,6] + gamma[1,9] - 3 a[1,9,6,8] >= 0
 R2511: alpha[6,8] + beta[9,6] + gamma[2,9] - 3 a[2,9,6,8] >= 0
 R2512: alpha[6,8] + beta[9,6] + gamma[3,9] - 3 a[3,9,6,8] >= 0
 R2513: alpha[6,8] + beta[9,6] + gamma[4,9] - 3 a[4,9,6,8] >= 0
 R2514: alpha[6,8] + beta[9,6] + gamma[5,9] - 3 a[5,9,6,8] >= 0
 R2515: alpha[6,8] + beta[9,6] + gamma[6,9] - 3 a[6,9,6,8] >= 0
 R2516: alpha[6,8] + beta[9,6] + gamma[7,9] - 3 a[7,9,6,8] >= 0
 R2517: alpha[6,8] + beta[9,6] + gamma[8,9] - 3 a[8,9,6,8] >= 0
 R2518: alpha[6,8] + beta[9,6] + gamma[9,9] - 3 a[9,9,6,8] >= 0
 R2519: alpha[6,9] + beta[0,6] + gamma[0,0] - 3 a[0,0,6,9] >= 0
 R2520: alpha[6,9] + beta[1,6] + gamma[0,1] - 3 a[0,1,6,9] >= 0
 R2521: alpha[6,9] + beta[1,6] + gamma[1,1] - 3 a[1,1,6,9] >= 0
 R2522: alpha[6,9] + beta[2,6] + gamma[0,2] - 3 a[0,2,6,9] >= 0
 R2523: alpha[6,9] + beta[2,6] + gamma[1,2] - 3 a[1,2,6,9] >= 0
 R2524: alpha[6,9] + beta[2,6] + gamma[2,2] - 3 a[2,2,6,9] >= 0
 R2525: alpha[6,9] + beta[3,6] + gamma[0,3] - 3 a[0,3,6,9] >= 0
 R2526: alpha[6,9] + beta[3,6] + gamma[1,3] - 3 a[1,3,6,9] >= 0
 R2527: alpha[6,9] + beta[3,6] + gamma[2,3] - 3 a[2,3,6,9] >= 0
 R2528: alpha[6,9] + beta[3,6] + gamma[3,3] - 3 a[3,3,6,9] >= 0
 R2529: alpha[6,9] + beta[4,6] + gamma[0,4] - 3 a[0,4,6,9] >= 0
 R2530: alpha[6,9] + beta[4,6] + gamma[1,4] - 3 a[1,4,6,9] >= 0
 R2531: alpha[6,9] + beta[4,6] + gamma[2,4] - 3 a[2,4,6,9] >= 0
 R2532: alpha[6,9] + beta[4,6] + gamma[3,4] - 3 a[3,4,6,9] >= 0
 R2533: alpha[6,9] + beta[4,6] + gamma[4,4] - 3 a[4,4,6,9] >= 0
 R2534: alpha[6,9] + beta[5,6] + gamma[0,5] - 3 a[0,5,6,9] >= 0
 R2535: alpha[6,9] + beta[5,6] + gamma[1,5] - 3 a[1,5,6,9] >= 0
 R2536: alpha[6,9] + beta[5,6] + gamma[2,5] - 3 a[2,5,6,9] >= 0
 R2537: alpha[6,9] + beta[5,6] + gamma[3,5] - 3 a[3,5,6,9] >= 0
 R2538: alpha[6,9] + beta[5,6] + gamma[4,5] - 3 a[4,5,6,9] >= 0
 R2539: alpha[6,9] + beta[5,6] + gamma[5,5] - 3 a[5,5,6,9] >= 0
 R2540: alpha[6,9] + beta[6,6] + gamma[0,6] - 3 a[0,6,6,9] >= 0
 R2541: alpha[6,9] + beta[6,6] + gamma[1,6] - 3 a[1,6,6,9] >= 0
 R2542: alpha[6,9] + beta[6,6] + gamma[2,6] - 3 a[2,6,6,9] >= 0
 R2543: alpha[6,9] + beta[6,6] + gamma[3,6] - 3 a[3,6,6,9] >= 0
 R2544: alpha[6,9] + beta[6,6] + gamma[4,6] - 3 a[4,6,6,9] >= 0
 R2545: alpha[6,9] + beta[6,6] + gamma[5,6] - 3 a[5,6,6,9] >= 0
 R2546: alpha[6,9] + beta[6,6] + gamma[6,6] - 3 a[6,6,6,9] >= 0
 R2547: alpha[6,9] + beta[7,6] + gamma[0,7] - 3 a[0,7,6,9] >= 0
 R2548: alpha[6,9] + beta[7,6] + gamma[1,7] - 3 a[1,7,6,9] >= 0
 R2549: alpha[6,9] + beta[7,6] + gamma[2,7] - 3 a[2,7,6,9] >= 0
 R2550: alpha[6,9] + beta[7,6] + gamma[3,7] - 3 a[3,7,6,9] >= 0
 R2551: alpha[6,9] + beta[7,6] + gamma[4,7] - 3 a[4,7,6,9] >= 0
 R2552: alpha[6,9] + beta[7,6] + gamma[5,7] - 3 a[5,7,6,9] >= 0
 R2553: alpha[6,9] + beta[7,6] + gamma[6,7] - 3 a[6,7,6,9] >= 0
 R2554: alpha[6,9] + beta[7,6] + gamma[7,7] - 3 a[7,7,6,9] >= 0
 R2555: alpha[6,9] + beta[8,6] + gamma[0,8] - 3 a[0,8,6,9] >= 0
 R2556: alpha[6,9] + beta[8,6] + gamma[1,8] - 3 a[1,8,6,9] >= 0
 R2557: alpha[6,9] + beta[8,6] + gamma[2,8] - 3 a[2,8,6,9] >= 0
 R2558: alpha[6,9] + beta[8,6] + gamma[3,8] - 3 a[3,8,6,9] >= 0
 R2559: alpha[6,9] + beta[8,6] + gamma[4,8] - 3 a[4,8,6,9] >= 0
 R2560: alpha[6,9] + beta[8,6] + gamma[5,8] - 3 a[5,8,6,9] >= 0
 R2561: alpha[6,9] + beta[8,6] + gamma[6,8] - 3 a[6,8,6,9] >= 0
 R2562: alpha[6,9] + beta[8,6] + gamma[7,8] - 3 a[7,8,6,9] >= 0
 R2563: alpha[6,9] + beta[8,6] + gamma[8,8] - 3 a[8,8,6,9] >= 0
 R2564: alpha[6,9] + beta[9,6] + gamma[0,9] - 3 a[0,9,6,9] >= 0
 R2565: alpha[6,9] + beta[9,6] + gamma[1,9] - 3 a[1,9,6,9] >= 0
 R2566: alpha[6,9] + beta[9,6] + gamma[2,9] - 3 a[2,9,6,9] >= 0
 R2567: alpha[6,9] + beta[9,6] + gamma[3,9] - 3 a[3,9,6,9] >= 0
 R2568: alpha[6,9] + beta[9,6] + gamma[4,9] - 3 a[4,9,6,9] >= 0
 R2569: alpha[6,9] + beta[9,6] + gamma[5,9] - 3 a[5,9,6,9] >= 0
 R2570: alpha[6,9] + beta[9,6] + gamma[6,9] - 3 a[6,9,6,9] >= 0
 R2571: alpha[6,9] + beta[9,6] + gamma[7,9] - 3 a[7,9,6,9] >= 0
 R2572: alpha[6,9] + beta[9,6] + gamma[8,9] - 3 a[8,9,6,9] >= 0
 R2573: alpha[6,9] + beta[9,6] + gamma[9,9] - 3 a[9,9,6,9] >= 0
 R2574: alpha[7,7] + beta[0,7] + gamma[0,0] - 3 a[0,0,7,7] >= 0
 R2575: alpha[7,7] + beta[1,7] + gamma[0,1] - 3 a[0,1,7,7] >= 0
 R2576: alpha[7,7] + beta[1,7] + gamma[1,1] - 3 a[1,1,7,7] >= 0
 R2577: alpha[7,7] + beta[2,7] + gamma[0,2] - 3 a[0,2,7,7] >= 0
 R2578: alpha[7,7] + beta[2,7] + gamma[1,2] - 3 a[1,2,7,7] >= 0
 R2579: alpha[7,7] + beta[2,7] + gamma[2,2] - 3 a[2,2,7,7] >= 0
 R2580: alpha[7,7] + beta[3,7] + gamma[0,3] - 3 a[0,3,7,7] >= 0
 R2581: alpha[7,7] + beta[3,7] + gamma[1,3] - 3 a[1,3,7,7] >= 0
 R2582: alpha[7,7] + beta[3,7] + gamma[2,3] - 3 a[2,3,7,7] >= 0
 R2583: alpha[7,7] + beta[3,7] + gamma[3,3] - 3 a[3,3,7,7] >= 0
 R2584: alpha[7,7] + beta[4,7] + gamma[0,4] - 3 a[0,4,7,7] >= 0
 R2585: alpha[7,7] + beta[4,7] + gamma[1,4] - 3 a[1,4,7,7] >= 0
 R2586: alpha[7,7] + beta[4,7] + gamma[2,4] - 3 a[2,4,7,7] >= 0
 R2587: alpha[7,7] + beta[4,7] + gamma[3,4] - 3 a[3,4,7,7] >= 0
 R2588: alpha[7,7] + beta[4,7] + gamma[4,4] - 3 a[4,4,7,7] >= 0
 R2589: alpha[7,7] + beta[5,7] + gamma[0,5] - 3 a[0,5,7,7] >= 0
 R2590: alpha[7,7] + beta[5,7] + gamma[1,5] - 3 a[1,5,7,7] >= 0
 R2591: alpha[7,7] + beta[5,7] + gamma[2,5] - 3 a[2,5,7,7] >= 0
 R2592: alpha[7,7] + beta[5,7] + gamma[3,5] - 3 a[3,5,7,7] >= 0
 R2593: alpha[7,7] + beta[5,7] + gamma[4,5] - 3 a[4,5,7,7] >= 0
 R2594: alpha[7,7] + beta[5,7] + gamma[5,5] - 3 a[5,5,7,7] >= 0
 R2595: alpha[7,7] + beta[6,7] + gamma[0,6] - 3 a[0,6,7,7] >= 0
 R2596: alpha[7,7] + beta[6,7] + gamma[1,6] - 3 a[1,6,7,7] >= 0
 R2597: alpha[7,7] + beta[6,7] + gamma[2,6] - 3 a[2,6,7,7] >= 0
 R2598: alpha[7,7] + beta[6,7] + gamma[3,6] - 3 a[3,6,7,7] >= 0
 R2599: alpha[7,7] + beta[6,7] + gamma[4,6] - 3 a[4,6,7,7] >= 0
 R2600: alpha[7,7] + beta[6,7] + gamma[5,6] - 3 a[5,6,7,7] >= 0
 R2601: alpha[7,7] + beta[6,7] + gamma[6,6] - 3 a[6,6,7,7] >= 0
 R2602: alpha[7,7] + beta[7,7] + gamma[0,7] - 3 a[0,7,7,7] >= 0
 R2603: alpha[7,7] + beta[7,7] + gamma[1,7] - 3 a[1,7,7,7] >= 0
 R2604: alpha[7,7] + beta[7,7] + gamma[2,7] - 3 a[2,7,7,7] >= 0
 R2605: alpha[7,7] + beta[7,7] + gamma[3,7] - 3 a[3,7,7,7] >= 0
 R2606: alpha[7,7] + beta[7,7] + gamma[4,7] - 3 a[4,7,7,7] >= 0
 R2607: alpha[7,7] + beta[7,7] + gamma[5,7] - 3 a[5,7,7,7] >= 0
 R2608: alpha[7,7] + beta[7,7] + gamma[6,7] - 3 a[6,7,7,7] >= 0
 R2609: alpha[7,7] + beta[7,7] + gamma[7,7] - 3 a[7,7,7,7] >= 0
 R2610: alpha[7,7] + beta[8,7] + gamma[0,8] - 3 a[0,8,7,7] >= 0
 R2611: alpha[7,7] + beta[8,7] + gamma[1,8] - 3 a[1,8,7,7] >= 0
 R2612: alpha[7,7] + beta[8,7] + gamma[2,8] - 3 a[2,8,7,7] >= 0
 R2613: alpha[7,7] + beta[8,7] + gamma[3,8] - 3 a[3,8,7,7] >= 0
 R2614: alpha[7,7] + beta[8,7] + gamma[4,8] - 3 a[4,8,7,7] >= 0
 R2615: alpha[7,7] + beta[8,7] + gamma[5,8] - 3 a[5,8,7,7] >= 0
 R2616: alpha[7,7] + beta[8,7] + gamma[6,8] - 3 a[6,8,7,7] >= 0
 R2617: alpha[7,7] + beta[8,7] + gamma[7,8] - 3 a[7,8,7,7] >= 0
 R2618: alpha[7,7] + beta[8,7] + gamma[8,8] - 3 a[8,8,7,7] >= 0
 R2619: alpha[7,7] + beta[9,7] + gamma[0,9] - 3 a[0,9,7,7] >= 0
 R2620: alpha[7,7] + beta[9,7] + gamma[1,9] - 3 a[1,9,7,7] >= 0
 R2621: alpha[7,7] + beta[9,7] + gamma[2,9] - 3 a[2,9,7,7] >= 0
 R2622: alpha[7,7] + beta[9,7] + gamma[3,9] - 3 a[3,9,7,7] >= 0
 R2623: alpha[7,7] + beta[9,7] + gamma[4,9] - 3 a[4,9,7,7] >= 0
 R2624: alpha[7,7] + beta[9,7] + gamma[5,9] - 3 a[5,9,7,7] >= 0
 R2625: alpha[7,7] + beta[9,7] + gamma[6,9] - 3 a[6,9,7,7] >= 0
 R2626: alpha[7,7] + beta[9,7] + gamma[7,9] - 3 a[7,9,7,7] >= 0
 R2627: alpha[7,7] + beta[9,7] + gamma[8,9] - 3 a[8,9,7,7] >= 0
 R2628: alpha[7,7] + beta[9,7] + gamma[9,9] - 3 a[9,9,7,7] >= 0
 R2629: alpha[8,8] + beta[0,8] + gamma[0,0] - 3 a[0,0,8,8] >= 0
 R2630: alpha[8,8] + beta[1,8] + gamma[0,1] - 3 a[0,1,8,8] >= 0
 R2631: alpha[8,8] + beta[1,8] + gamma[1,1] - 3 a[1,1,8,8] >= 0
 R2632: alpha[8,8] + beta[2,8] + gamma[0,2] - 3 a[0,2,8,8] >= 0
 R2633: alpha[8,8] + beta[2,8] + gamma[1,2] - 3 a[1,2,8,8] >= 0
 R2634: alpha[8,8] + beta[2,8] + gamma[2,2] - 3 a[2,2,8,8] >= 0
 R2635: alpha[8,8] + beta[3,8] + gamma[0,3] - 3 a[0,3,8,8] >= 0
 R2636: alpha[8,8] + beta[3,8] + gamma[1,3] - 3 a[1,3,8,8] >= 0
 R2637: alpha[8,8] + beta[3,8] + gamma[2,3] - 3 a[2,3,8,8] >= 0
 R2638: alpha[8,8] + beta[3,8] + gamma[3,3] - 3 a[3,3,8,8] >= 0
 R2639: alpha[8,8] + beta[4,8] + gamma[0,4] - 3 a[0,4,8,8] >= 0
 R2640: alpha[8,8] + beta[4,8] + gamma[1,4] - 3 a[1,4,8,8] >= 0
 R2641: alpha[8,8] + beta[4,8] + gamma[2,4] - 3 a[2,4,8,8] >= 0
 R2642: alpha[8,8] + beta[4,8] + gamma[3,4] - 3 a[3,4,8,8] >= 0
 R2643: alpha[8,8] + beta[4,8] + gamma[4,4] - 3 a[4,4,8,8] >= 0
 R2644: alpha[8,8] + beta[5,8] + gamma[0,5] - 3 a[0,5,8,8] >= 0
 R2645: alpha[8,8] + beta[5,8] + gamma[1,5] - 3 a[1,5,8,8] >= 0
 R2646: alpha[8,8] + beta[5,8] + gamma[2,5] - 3 a[2,5,8,8] >= 0
 R2647: alpha[8,8] + beta[5,8] + gamma[3,5] - 3 a[3,5,8,8] >= 0
 R2648: alpha[8,8] + beta[5,8] + gamma[4,5] - 3 a[4,5,8,8] >= 0
 R2649: alpha[8,8] + beta[5,8] + gamma[5,5] - 3 a[5,5,8,8] >= 0
 R2650: alpha[8,8] + beta[6,8] + gamma[0,6] - 3 a[0,6,8,8] >= 0
 R2651: alpha[8,8] + beta[6,8] + gamma[1,6] - 3 a[1,6,8,8] >= 0
 R2652: alpha[8,8] + beta[6,8] + gamma[2,6] - 3 a[2,6,8,8] >= 0
 R2653: alpha[8,8] + beta[6,8] + gamma[3,6] - 3 a[3,6,8,8] >= 0
 R2654: alpha[8,8] + beta[6,8] + gamma[4,6] - 3 a[4,6,8,8] >= 0
 R2655: alpha[8,8] + beta[6,8] + gamma[5,6] - 3 a[5,6,8,8] >= 0
 R2656: alpha[8,8] + beta[6,8] + gamma[6,6] - 3 a[6,6,8,8] >= 0
 R2657: alpha[8,8] + beta[7,8] + gamma[0,7] - 3 a[0,7,8,8] >= 0
 R2658: alpha[8,8] + beta[7,8] + gamma[1,7] - 3 a[1,7,8,8] >= 0
 R2659: alpha[8,8] + beta[7,8] + gamma[2,7] - 3 a[2,7,8,8] >= 0
 R2660: alpha[8,8] + beta[7,8] + gamma[3,7] - 3 a[3,7,8,8] >= 0
 R2661: alpha[8,8] + beta[7,8] + gamma[4,7] - 3 a[4,7,8,8] >= 0
 R2662: alpha[8,8] + beta[7,8] + gamma[5,7] - 3 a[5,7,8,8] >= 0
 R2663: alpha[8,8] + beta[7,8] + gamma[6,7] - 3 a[6,7,8,8] >= 0
 R2664: alpha[8,8] + beta[7,8] + gamma[7,7] - 3 a[7,7,8,8] >= 0
 R2665: alpha[8,8] + beta[8,8] + gamma[0,8] - 3 a[0,8,8,8] >= 0
 R2666: alpha[8,8] + beta[8,8] + gamma[1,8] - 3 a[1,8,8,8] >= 0
 R2667: alpha[8,8] + beta[8,8] + gamma[2,8] - 3 a[2,8,8,8] >= 0
 R2668: alpha[8,8] + beta[8,8] + gamma[3,8] - 3 a[3,8,8,8] >= 0
 R2669: alpha[8,8] + beta[8,8] + gamma[4,8] - 3 a[4,8,8,8] >= 0
 R2670: alpha[8,8] + beta[8,8] + gamma[5,8] - 3 a[5,8,8,8] >= 0
 R2671: alpha[8,8] + beta[8,8] + gamma[6,8] - 3 a[6,8,8,8] >= 0
 R2672: alpha[8,8] + beta[8,8] + gamma[7,8] - 3 a[7,8,8,8] >= 0
 R2673: alpha[8,8] + beta[8,8] + gamma[8,8] - 3 a[8,8,8,8] >= 0
 R2674: alpha[8,8] + beta[9,8] + gamma[0,9] - 3 a[0,9,8,8] >= 0
 R2675: alpha[8,8] + beta[9,8] + gamma[1,9] - 3 a[1,9,8,8] >= 0
 R2676: alpha[8,8] + beta[9,8] + gamma[2,9] - 3 a[2,9,8,8] >= 0
 R2677: alpha[8,8] + beta[9,8] + gamma[3,9] - 3 a[3,9,8,8] >= 0
 R2678: alpha[8,8] + beta[9,8] + gamma[4,9] - 3 a[4,9,8,8] >= 0
 R2679: alpha[8,8] + beta[9,8] + gamma[5,9] - 3 a[5,9,8,8] >= 0
 R2680: alpha[8,8] + beta[9,8] + gamma[6,9] - 3 a[6,9,8,8] >= 0
 R2681: alpha[8,8] + beta[9,8] + gamma[7,9] - 3 a[7,9,8,8] >= 0
 R2682: alpha[8,8] + beta[9,8] + gamma[8,9] - 3 a[8,9,8,8] >= 0
 R2683: alpha[8,8] + beta[9,8] + gamma[9,9] - 3 a[9,9,8,8] >= 0
 R2684: alpha[8,9] + beta[0,8] + gamma[0,0] - 3 a[0,0,8,9] >= 0
 R2685: alpha[8,9] + beta[1,8] + gamma[0,1] - 3 a[0,1,8,9] >= 0
 R2686: alpha[8,9] + beta[1,8] + gamma[1,1] - 3 a[1,1,8,9] >= 0
 R2687: alpha[8,9] + beta[2,8] + gamma[0,2] - 3 a[0,2,8,9] >= 0
 R2688: alpha[8,9] + beta[2,8] + gamma[1,2] - 3 a[1,2,8,9] >= 0
 R2689: alpha[8,9] + beta[2,8] + gamma[2,2] - 3 a[2,2,8,9] >= 0
 R2690: alpha[8,9] + beta[3,8] + gamma[0,3] - 3 a[0,3,8,9] >= 0
 R2691: alpha[8,9] + beta[3,8] + gamma[1,3] - 3 a[1,3,8,9] >= 0
 R2692: alpha[8,9] + beta[3,8] + gamma[2,3] - 3 a[2,3,8,9] >= 0
 R2693: alpha[8,9] + beta[3,8] + gamma[3,3] - 3 a[3,3,8,9] >= 0
 R2694: alpha[8,9] + beta[4,8] + gamma[0,4] - 3 a[0,4,8,9] >= 0
 R2695: alpha[8,9] + beta[4,8] + gamma[1,4] - 3 a[1,4,8,9] >= 0
 R2696: alpha[8,9] + beta[4,8] + gamma[2,4] - 3 a[2,4,8,9] >= 0
 R2697: alpha[8,9] + beta[4,8] + gamma[3,4] - 3 a[3,4,8,9] >= 0
 R2698: alpha[8,9] + beta[4,8] + gamma[4,4] - 3 a[4,4,8,9] >= 0
 R2699: alpha[8,9] + beta[5,8] + gamma[0,5] - 3 a[0,5,8,9] >= 0
 R2700: alpha[8,9] + beta[5,8] + gamma[1,5] - 3 a[1,5,8,9] >= 0
 R2701: alpha[8,9] + beta[5,8] + gamma[2,5] - 3 a[2,5,8,9] >= 0
 R2702: alpha[8,9] + beta[5,8] + gamma[3,5] - 3 a[3,5,8,9] >= 0
 R2703: alpha[8,9] + beta[5,8] + gamma[4,5] - 3 a[4,5,8,9] >= 0
 R2704: alpha[8,9] + beta[5,8] + gamma[5,5] - 3 a[5,5,8,9] >= 0
 R2705: alpha[8,9] + beta[6,8] + gamma[0,6] - 3 a[0,6,8,9] >= 0
 R2706: alpha[8,9] + beta[6,8] + gamma[1,6] - 3 a[1,6,8,9] >= 0
 R2707: alpha[8,9] + beta[6,8] + gamma[2,6] - 3 a[2,6,8,9] >= 0
 R2708: alpha[8,9] + beta[6,8] + gamma[3,6] - 3 a[3,6,8,9] >= 0
 R2709: alpha[8,9] + beta[6,8] + gamma[4,6] - 3 a[4,6,8,9] >= 0
 R2710: alpha[8,9] + beta[6,8] + gamma[5,6] - 3 a[5,6,8,9] >= 0
 R2711: alpha[8,9] + beta[6,8] + gamma[6,6] - 3 a[6,6,8,9] >= 0
 R2712: alpha[8,9] + beta[7,8] + gamma[0,7] - 3 a[0,7,8,9] >= 0
 R2713: alpha[8,9] + beta[7,8] + gamma[1,7] - 3 a[1,7,8,9] >= 0
 R2714: alpha[8,9] + beta[7,8] + gamma[2,7] - 3 a[2,7,8,9] >= 0
 R2715: alpha[8,9] + beta[7,8] + gamma[3,7] - 3 a[3,7,8,9] >= 0
 R2716: alpha[8,9] + beta[7,8] + gamma[4,7] - 3 a[4,7,8,9] >= 0
 R2717: alpha[8,9] + beta[7,8] + gamma[5,7] - 3 a[5,7,8,9] >= 0
 R2718: alpha[8,9] + beta[7,8] + gamma[6,7] - 3 a[6,7,8,9] >= 0
 R2719: alpha[8,9] + beta[7,8] + gamma[7,7] - 3 a[7,7,8,9] >= 0
 R2720: alpha[8,9] + beta[8,8] + gamma[0,8] - 3 a[0,8,8,9] >= 0
 R2721: alpha[8,9] + beta[8,8] + gamma[1,8] - 3 a[1,8,8,9] >= 0
 R2722: alpha[8,9] + beta[8,8] + gamma[2,8] - 3 a[2,8,8,9] >= 0
 R2723: alpha[8,9] + beta[8,8] + gamma[3,8] - 3 a[3,8,8,9] >= 0
 R2724: alpha[8,9] + beta[8,8] + gamma[4,8] - 3 a[4,8,8,9] >= 0
 R2725: alpha[8,9] + beta[8,8] + gamma[5,8] - 3 a[5,8,8,9] >= 0
 R2726: alpha[8,9] + beta[8,8] + gamma[6,8] - 3 a[6,8,8,9] >= 0
 R2727: alpha[8,9] + beta[8,8] + gamma[7,8] - 3 a[7,8,8,9] >= 0
 R2728: alpha[8,9] + beta[8,8] + gamma[8,8] - 3 a[8,8,8,9] >= 0
 R2729: alpha[8,9] + beta[9,8] + gamma[0,9] - 3 a[0,9,8,9] >= 0
 R2730: alpha[8,9] + beta[9,8] + gamma[1,9] - 3 a[1,9,8,9] >= 0
 R2731: alpha[8,9] + beta[9,8] + gamma[2,9] - 3 a[2,9,8,9] >= 0
 R2732: alpha[8,9] + beta[9,8] + gamma[3,9] - 3 a[3,9,8,9] >= 0
 R2733: alpha[8,9] + beta[9,8] + gamma[4,9] - 3 a[4,9,8,9] >= 0
 R2734: alpha[8,9] + beta[9,8] + gamma[5,9] - 3 a[5,9,8,9] >= 0
 R2735: alpha[8,9] + beta[9,8] + gamma[6,9] - 3 a[6,9,8,9] >= 0
 R2736: alpha[8,9] + beta[9,8] + gamma[7,9] - 3 a[7,9,8,9] >= 0
 R2737: alpha[8,9] + beta[9,8] + gamma[8,9] - 3 a[8,9,8,9] >= 0
 R2738: alpha[8,9] + beta[9,8] + gamma[9,9] - 3 a[9,9,8,9] >= 0
 R2739: alpha[9,9] + beta[0,9] + gamma[0,0] - 3 a[0,0,9,9] >= 0
 R2740: alpha[9,9] + beta[1,9] + gamma[0,1] - 3 a[0,1,9,9] >= 0
 R2741: alpha[9,9] + beta[1,9] + gamma[1,1] - 3 a[1,1,9,9] >= 0
 R2742: alpha[9,9] + beta[2,9] + gamma[0,2] - 3 a[0,2,9,9] >= 0
 R2743: alpha[9,9] + beta[2,9] + gamma[1,2] - 3 a[1,2,9,9] >= 0
 R2744: alpha[9,9] + beta[2,9] + gamma[2,2] - 3 a[2,2,9,9] >= 0
 R2745: alpha[9,9] + beta[3,9] + gamma[0,3] - 3 a[0,3,9,9] >= 0
 R2746: alpha[9,9] + beta[3,9] + gamma[1,3] - 3 a[1,3,9,9] >= 0
 R2747: alpha[9,9] + beta[3,9] + gamma[2,3] - 3 a[2,3,9,9] >= 0
 R2748: alpha[9,9] + beta[3,9] + gamma[3,3] - 3 a[3,3,9,9] >= 0
 R2749: alpha[9,9] + beta[4,9] + gamma[0,4] - 3 a[0,4,9,9] >= 0
 R2750: alpha[9,9] + beta[4,9] + gamma[1,4] - 3 a[1,4,9,9] >= 0
 R2751: alpha[9,9] + beta[4,9] + gamma[2,4] - 3 a[2,4,9,9] >= 0
 R2752: alpha[9,9] + beta[4,9] + gamma[3,4] - 3 a[3,4,9,9] >= 0
 R2753: alpha[9,9] + beta[4,9] + gamma[4,4] - 3 a[4,4,9,9] >= 0
 R2754: alpha[9,9] + beta[5,9] + gamma[0,5] - 3 a[0,5,9,9] >= 0
 R2755: alpha[9,9] + beta[5,9] + gamma[1,5] - 3 a[1,5,9,9] >= 0
 R2756: alpha[9,9] + beta[5,9] + gamma[2,5] - 3 a[2,5,9,9] >= 0
 R2757: alpha[9,9] + beta[5,9] + gamma[3,5] - 3 a[3,5,9,9] >= 0
 R2758: alpha[9,9] + beta[5,9] + gamma[4,5] - 3 a[4,5,9,9] >= 0
 R2759: alpha[9,9] + beta[5,9] + gamma[5,5] - 3 a[5,5,9,9] >= 0
 R2760: alpha[9,9] + beta[6,9] + gamma[0,6] - 3 a[0,6,9,9] >= 0
 R2761: alpha[9,9] + beta[6,9] + gamma[1,6] - 3 a[1,6,9,9] >= 0
 R2762: alpha[9,9] + beta[6,9] + gamma[2,6] - 3 a[2,6,9,9] >= 0
 R2763: alpha[9,9] + beta[6,9] + gamma[3,6] - 3 a[3,6,9,9] >= 0
 R2764: alpha[9,9] + beta[6,9] + gamma[4,6] - 3 a[4,6,9,9] >= 0
 R2765: alpha[9,9] + beta[6,9] + gamma[5,6] - 3 a[5,6,9,9] >= 0
 R2766: alpha[9,9] + beta[6,9] + gamma[6,6] - 3 a[6,6,9,9] >= 0
 R2767: alpha[9,9] + beta[7,9] + gamma[0,7] - 3 a[0,7,9,9] >= 0
 R2768: alpha[9,9] + beta[7,9] + gamma[1,7] - 3 a[1,7,9,9] >= 0
 R2769: alpha[9,9] + beta[7,9] + gamma[2,7] - 3 a[2,7,9,9] >= 0
 R2770: alpha[9,9] + beta[7,9] + gamma[3,7] - 3 a[3,7,9,9] >= 0
 R2771: alpha[9,9] + beta[7,9] + gamma[4,7] - 3 a[4,7,9,9] >= 0
 R2772: alpha[9,9] + beta[7,9] + gamma[5,7] - 3 a[5,7,9,9] >= 0
 R2773: alpha[9,9] + beta[7,9] + gamma[6,7] - 3 a[6,7,9,9] >= 0
 R2774: alpha[9,9] + beta[7,9] + gamma[7,7] - 3 a[7,7,9,9] >= 0
 R2775: alpha[9,9] + beta[8,9] + gamma[0,8] - 3 a[0,8,9,9] >= 0
 R2776: alpha[9,9] + beta[8,9] + gamma[1,8] - 3 a[1,8,9,9] >= 0
 R2777: alpha[9,9] + beta[8,9] + gamma[2,8] - 3 a[2,8,9,9] >= 0
 R2778: alpha[9,9] + beta[8,9] + gamma[3,8] - 3 a[3,8,9,9] >= 0
 R2779: alpha[9,9] + beta[8,9] + gamma[4,8] - 3 a[4,8,9,9] >= 0
 R2780: alpha[9,9] + beta[8,9] + gamma[5,8] - 3 a[5,8,9,9] >= 0
 R2781: alpha[9,9] + beta[8,9] + gamma[6,8] - 3 a[6,8,9,9] >= 0
 R2782: alpha[9,9] + beta[8,9] + gamma[7,8] - 3 a[7,8,9,9] >= 0
 R2783: alpha[9,9] + beta[8,9] + gamma[8,8] - 3 a[8,8,9,9] >= 0
 R2784: alpha[9,9] + beta[9,9] + gamma[0,9] - 3 a[0,9,9,9] >= 0
 R2785: alpha[9,9] + beta[9,9] + gamma[1,9] - 3 a[1,9,9,9] >= 0
 R2786: alpha[9,9] + beta[9,9] + gamma[2,9] - 3 a[2,9,9,9] >= 0
 R2787: alpha[9,9] + beta[9,9] + gamma[3,9] - 3 a[3,9,9,9] >= 0
 R2788: alpha[9,9] + beta[9,9] + gamma[4,9] - 3 a[4,9,9,9] >= 0
 R2789: alpha[9,9] + beta[9,9] + gamma[5,9] - 3 a[5,9,9,9] >= 0
 R2790: alpha[9,9] + beta[9,9] + gamma[6,9] - 3 a[6,9,9,9] >= 0
 R2791: alpha[9,9] + beta[9,9] + gamma[7,9] - 3 a[7,9,9,9] >= 0
 R2792: alpha[9,9] + beta[9,9] + gamma[8,9] - 3 a[8,9,9,9] >= 0
 R2793: alpha[9,9] + beta[9,9] + gamma[9,9] - 3 a[9,9,9,9] >= 0
 R2794: 5000 a[0,0,0,0] + 5000 a[0,1,0,0] + 5000 a[0,2,0,0]
   + 5000 a[0,3,0,0] + 5000 a[0,4,0,0] + 5000 a[0,5,0,0] + 5000 a[0,6,0,0]
   + 5000 a[0,7,0,0] + 5000 a[0,8,0,0] + 5000 a[0,9,0,0] - Q[0] <= 0
 R2795: 5000 a[1,1,0,0] + 5000 a[1,2,0,0] + 5000 a[1,3,0,0]
   + 5000 a[1,4,0,0] + 5000 a[1,5,0,0] + 5000 a[1,6,0,0] + 5000 a[1,7,0,0]
   + 5000 a[1,8,0,0] + 5000 a[1,9,0,0] - Q[1] <= 0
 R2796: 5000 a[2,2,0,0] + 5000 a[2,3,0,0] + 5000 a[2,4,0,0]
   + 5000 a[2,5,0,0] + 5000 a[2,6,0,0] + 5000 a[2,7,0,0] + 5000 a[2,8,0,0]
   + 5000 a[2,9,0,0] - Q[2] <= 0
 R2797: 5000 a[3,3,0,0] + 5000 a[3,4,0,0] + 5000 a[3,5,0,0]
   + 5000 a[3,6,0,0] + 5000 a[3,7,0,0] + 5000 a[3,8,0,0] + 5000 a[3,9,0,0]
   - Q[3] <= 0
 R2798: 5000 a[4,4,0,0] + 5000 a[4,5,0,0] + 5000 a[4,6,0,0]
   + 5000 a[4,7,0,0] + 5000 a[4,8,0,0] + 5000 a[4,9,0,0] - Q[4] <= 0
 R2799: 5000 a[5,5,0,0] + 5000 a[5,6,0,0] + 5000 a[5,7,0,0]
   + 5000 a[5,8,0,0] + 5000 a[5,9,0,0] - Q[5] <= 0
 R2800: 5000 a[6,6,0,0] + 5000 a[6,7,0,0] + 5000 a[6,8,0,0]
   + 5000 a[6,9,0,0] - Q[6] <= 0
 R2801: 5000 a[7,7,0,0] + 5000 a[7,8,0,0] + 5000 a[7,9,0,0] - Q[7] <= 0
 R2802: 5000 a[8,8,0,0] + 5000 a[8,9,0,0] - Q[8] <= 0
 R2803: 5000 a[9,9,0,0] - Q[9] <= 0
 R2804: 5000 a[0,0,1,1] + 5000 a[0,1,1,1] + 5000 a[0,2,1,1]
   + 5000 a[0,3,1,1] + 5000 a[0,4,1,1] + 5000 a[0,5,1,1] + 5000 a[0,6,1,1]
   + 5000 a[0,7,1,1] + 5000 a[0,8,1,1] + 5000 a[0,9,1,1] - Q[0] <= 0
 R2805: 5000 a[1,1,1,1] + 5000 a[1,2,1,1] + 5000 a[1,3,1,1]
   + 5000 a[1,4,1,1] + 5000 a[1,5,1,1] + 5000 a[1,6,1,1] + 5000 a[1,7,1,1]
   + 5000 a[1,8,1,1] + 5000 a[1,9,1,1] - Q[1] <= 0
 R2806: 5000 a[2,2,1,1] + 5000 a[2,3,1,1] + 5000 a[2,4,1,1]
   + 5000 a[2,5,1,1] + 5000 a[2,6,1,1] + 5000 a[2,7,1,1] + 5000 a[2,8,1,1]
   + 5000 a[2,9,1,1] - Q[2] <= 0
 R2807: 5000 a[3,3,1,1] + 5000 a[3,4,1,1] + 5000 a[3,5,1,1]
   + 5000 a[3,6,1,1] + 5000 a[3,7,1,1] + 5000 a[3,8,1,1] + 5000 a[3,9,1,1]
   - Q[3] <= 0
 R2808: 5000 a[4,4,1,1] + 5000 a[4,5,1,1] + 5000 a[4,6,1,1]
   + 5000 a[4,7,1,1] + 5000 a[4,8,1,1] + 5000 a[4,9,1,1] - Q[4] <= 0
 R2809: 5000 a[5,5,1,1] + 5000 a[5,6,1,1] + 5000 a[5,7,1,1]
   + 5000 a[5,8,1,1] + 5000 a[5,9,1,1] - Q[5] <= 0
 R2810: 5000 a[6,6,1,1] + 5000 a[6,7,1,1] + 5000 a[6,8,1,1]
   + 5000 a[6,9,1,1] - Q[6] <= 0
 R2811: 5000 a[7,7,1,1] + 5000 a[7,8,1,1] + 5000 a[7,9,1,1] - Q[7] <= 0
 R2812: 5000 a[8,8,1,1] + 5000 a[8,9,1,1] - Q[8] <= 0
 R2813: 5000 a[9,9,1,1] - Q[9] <= 0
 R2814: 5000 a[0,0,0,2] + 5000 a[0,1,0,2] + 5000 a[0,2,0,2]
   + 5000 a[0,3,0,2] + 5000 a[0,4,0,2] + 5000 a[0,5,0,2] + 5000 a[0,6,0,2]
   + 5000 a[0,7,0,2] + 5000 a[0,8,0,2] + 5000 a[0,9,0,2] + 5000 a[0,0,2,2]
   + 5000 a[0,1,2,2] + 5000 a[0,2,2,2] + 5000 a[0,3,2,2] + 5000 a[0,4,2,2]
   + 5000 a[0,5,2,2] + 5000 a[0,6,2,2] + 5000 a[0,7,2,2] + 5000 a[0,8,2,2]
   + 5000 a[0,9,2,2] - Q[0] <= 0
 R2815: 5000 a[1,1,0,2] + 5000 a[1,2,0,2] + 5000 a[1,3,0,2]
   + 5000 a[1,4,0,2] + 5000 a[1,5,0,2] + 5000 a[1,6,0,2] + 5000 a[1,7,0,2]
   + 5000 a[1,8,0,2] + 5000 a[1,9,0,2] + 5000 a[1,1,2,2] + 5000 a[1,2,2,2]
   + 5000 a[1,3,2,2] + 5000 a[1,4,2,2] + 5000 a[1,5,2,2] + 5000 a[1,6,2,2]
   + 5000 a[1,7,2,2] + 5000 a[1,8,2,2] + 5000 a[1,9,2,2] - Q[1] <= 0
 R2816: 5000 a[2,2,0,2] + 5000 a[2,3,0,2] + 5000 a[2,4,0,2]
   + 5000 a[2,5,0,2] + 5000 a[2,6,0,2] + 5000 a[2,7,0,2] + 5000 a[2,8,0,2]
   + 5000 a[2,9,0,2] + 5000 a[2,2,2,2] + 5000 a[2,3,2,2] + 5000 a[2,4,2,2]
   + 5000 a[2,5,2,2] + 5000 a[2,6,2,2] + 5000 a[2,7,2,2] + 5000 a[2,8,2,2]
   + 5000 a[2,9,2,2] - Q[2] <= 0
 R2817: 5000 a[3,3,0,2] + 5000 a[3,4,0,2] + 5000 a[3,5,0,2]
   + 5000 a[3,6,0,2] + 5000 a[3,7,0,2] + 5000 a[3,8,0,2] + 5000 a[3,9,0,2]
   + 5000 a[3,3,2,2] + 5000 a[3,4,2,2] + 5000 a[3,5,2,2] + 5000 a[3,6,2,2]
   + 5000 a[3,7,2,2] + 5000 a[3,8,2,2] + 5000 a[3,9,2,2] - Q[3] <= 0
 R2818: 5000 a[4,4,0,2] + 5000 a[4,5,0,2] + 5000 a[4,6,0,2]
   + 5000 a[4,7,0,2] + 5000 a[4,8,0,2] + 5000 a[4,9,0,2] + 5000 a[4,4,2,2]
   + 5000 a[4,5,2,2] + 5000 a[4,6,2,2] + 5000 a[4,7,2,2] + 5000 a[4,8,2,2]
   + 5000 a[4,9,2,2] - Q[4] <= 0
 R2819: 5000 a[5,5,0,2] + 5000 a[5,6,0,2] + 5000 a[5,7,0,2]
   + 5000 a[5,8,0,2] + 5000 a[5,9,0,2] + 5000 a[5,5,2,2] + 5000 a[5,6,2,2]
   + 5000 a[5,7,2,2] + 5000 a[5,8,2,2] + 5000 a[5,9,2,2] - Q[5] <= 0
 R2820: 5000 a[6,6,0,2] + 5000 a[6,7,0,2] + 5000 a[6,8,0,2]
   + 5000 a[6,9,0,2] + 5000 a[6,6,2,2] + 5000 a[6,7,2,2] + 5000 a[6,8,2,2]
   + 5000 a[6,9,2,2] - Q[6] <= 0
 R2821: 5000 a[7,7,0,2] + 5000 a[7,8,0,2] + 5000 a[7,9,0,2]
   + 5000 a[7,7,2,2] + 5000 a[7,8,2,2] + 5000 a[7,9,2,2] - Q[7] <= 0
 R2822: 5000 a[8,8,0,2] + 5000 a[8,9,0,2] + 5000 a[8,8,2,2]
   + 5000 a[8,9,2,2] - Q[8] <= 0
 R2823: 5000 a[9,9,0,2] + 5000 a[9,9,2,2] - Q[9] <= 0
 R2824: 5000 a[0,0,1,3] + 5000 a[0,1,1,3] + 5000 a[0,2,1,3]
   + 5000 a[0,3,1,3] + 5000 a[0,4,1,3] + 5000 a[0,5,1,3] + 5000 a[0,6,1,3]
   + 5000 a[0,7,1,3] + 5000 a[0,8,1,3] + 5000 a[0,9,1,3] + 5000 a[0,0,3,3]
   + 5000 a[0,1,3,3] + 5000 a[0,2,3,3] + 5000 a[0,3,3,3] + 5000 a[0,4,3,3]
   + 5000 a[0,5,3,3] + 5000 a[0,6,3,3] + 5000 a[0,7,3,3] + 5000 a[0,8,3,3]
   + 5000 a[0,9,3,3] - Q[0] <= 0
 R2825: 5000 a[1,1,1,3] + 5000 a[1,2,1,3] + 5000 a[1,3,1,3]
   + 5000 a[1,4,1,3] + 5000 a[1,5,1,3] + 5000 a[1,6,1,3] + 5000 a[1,7,1,3]
   + 5000 a[1,8,1,3] + 5000 a[1,9,1,3] + 5000 a[1,1,3,3] + 5000 a[1,2,3,3]
   + 5000 a[1,3,3,3] + 5000 a[1,4,3,3] + 5000 a[1,5,3,3] + 5000 a[1,6,3,3]
   + 5000 a[1,7,3,3] + 5000 a[1,8,3,3] + 5000 a[1,9,3,3] - Q[1] <= 0
 R2826: 5000 a[2,2,1,3] + 5000 a[2,3,1,3] + 5000 a[2,4,1,3]
   + 5000 a[2,5,1,3] + 5000 a[2,6,1,3] + 5000 a[2,7,1,3] + 5000 a[2,8,1,3]
   + 5000 a[2,9,1,3] + 5000 a[2,2,3,3] + 5000 a[2,3,3,3] + 5000 a[2,4,3,3]
   + 5000 a[2,5,3,3] + 5000 a[2,6,3,3] + 5000 a[2,7,3,3] + 5000 a[2,8,3,3]
   + 5000 a[2,9,3,3] - Q[2] <= 0
 R2827: 5000 a[3,3,1,3] + 5000 a[3,4,1,3] + 5000 a[3,5,1,3]
   + 5000 a[3,6,1,3] + 5000 a[3,7,1,3] + 5000 a[3,8,1,3] + 5000 a[3,9,1,3]
   + 5000 a[3,3,3,3] + 5000 a[3,4,3,3] + 5000 a[3,5,3,3] + 5000 a[3,6,3,3]
   + 5000 a[3,7,3,3] + 5000 a[3,8,3,3] + 5000 a[3,9,3,3] - Q[3] <= 0
 R2828: 5000 a[4,4,1,3] + 5000 a[4,5,1,3] + 5000 a[4,6,1,3]
   + 5000 a[4,7,1,3] + 5000 a[4,8,1,3] + 5000 a[4,9,1,3] + 5000 a[4,4,3,3]
   + 5000 a[4,5,3,3] + 5000 a[4,6,3,3] + 5000 a[4,7,3,3] + 5000 a[4,8,3,3]
   + 5000 a[4,9,3,3] - Q[4] <= 0
 R2829: 5000 a[5,5,1,3] + 5000 a[5,6,1,3] + 5000 a[5,7,1,3]
   + 5000 a[5,8,1,3] + 5000 a[5,9,1,3] + 5000 a[5,5,3,3] + 5000 a[5,6,3,3]
   + 5000 a[5,7,3,3] + 5000 a[5,8,3,3] + 5000 a[5,9,3,3] - Q[5] <= 0
 R2830: 5000 a[6,6,1,3] + 5000 a[6,7,1,3] + 5000 a[6,8,1,3]
   + 5000 a[6,9,1,3] + 5000 a[6,6,3,3] + 5000 a[6,7,3,3] + 5000 a[6,8,3,3]
   + 5000 a[6,9,3,3] - Q[6] <= 0
 R2831: 5000 a[7,7,1,3] + 5000 a[7,8,1,3] + 5000 a[7,9,1,3]
   + 5000 a[7,7,3,3] + 5000 a[7,8,3,3] + 5000 a[7,9,3,3] - Q[7] <= 0
 R2832: 5000 a[8,8,1,3] + 5000 a[8,9,1,3] + 5000 a[8,8,3,3]
   + 5000 a[8,9,3,3] - Q[8] <= 0
 R2833: 5000 a[9,9,1,3] + 5000 a[9,9,3,3] - Q[9] <= 0
 R2834: 1500 a[0,0,4,4] + 1500 a[0,1,4,4] + 1500 a[0,2,4,4]
   + 1500 a[0,3,4,4] + 1500 a[0,4,4,4] + 1500 a[0,5,4,4] + 1500 a[0,6,4,4]
   + 1500 a[0,7,4,4] + 1500 a[0,8,4,4] + 1500 a[0,9,4,4] - Q[0] <= 0
 R2835: 1500 a[1,1,4,4] + 1500 a[1,2,4,4] + 1500 a[1,3,4,4]
   + 1500 a[1,4,4,4] + 1500 a[1,5,4,4] + 1500 a[1,6,4,4] + 1500 a[1,7,4,4]
   + 1500 a[1,8,4,4] + 1500 a[1,9,4,4] - Q[1] <= 0
 R2836: 1500 a[2,2,4,4] + 1500 a[2,3,4,4] + 1500 a[2,4,4,4]
   + 1500 a[2,5,4,4] + 1500 a[2,6,4,4] + 1500 a[2,7,4,4] + 1500 a[2,8,4,4]
   + 1500 a[2,9,4,4] - Q[2] <= 0
 R2837: 1500 a[3,3,4,4] + 1500 a[3,4,4,4] + 1500 a[3,5,4,4]
   + 1500 a[3,6,4,4] + 1500 a[3,7,4,4] + 1500 a[3,8,4,4] + 1500 a[3,9,4,4]
   - Q[3] <= 0
 R2838: 1500 a[4,4,4,4] + 1500 a[4,5,4,4] + 1500 a[4,6,4,4]
   + 1500 a[4,7,4,4] + 1500 a[4,8,4,4] + 1500 a[4,9,4,4] - Q[4] <= 0
 R2839: 1500 a[5,5,4,4] + 1500 a[5,6,4,4] + 1500 a[5,7,4,4]
   + 1500 a[5,8,4,4] + 1500 a[5,9,4,4] - Q[5] <= 0
 R2840: 1500 a[6,6,4,4] + 1500 a[6,7,4,4] + 1500 a[6,8,4,4]
   + 1500 a[6,9,4,4] - Q[6] <= 0
 R2841: 1500 a[7,7,4,4] + 1500 a[7,8,4,4] + 1500 a[7,9,4,4] - Q[7] <= 0
 R2842: 1500 a[8,8,4,4] + 1500 a[8,9,4,4] - Q[8] <= 0
 R2843: 1500 a[9,9,4,4] - Q[9] <= 0
 R2844: 5000 a[0,0,5,5] + 5000 a[0,1,5,5] + 5000 a[0,2,5,5]
   + 5000 a[0,3,5,5] + 5000 a[0,4,5,5] + 5000 a[0,5,5,5] + 5000 a[0,6,5,5]
   + 5000 a[0,7,5,5] + 5000 a[0,8,5,5] + 5000 a[0,9,5,5] - Q[0] <= 0
 R2845: 5000 a[1,1,5,5] + 5000 a[1,2,5,5] + 5000 a[1,3,5,5]
   + 5000 a[1,4,5,5] + 5000 a[1,5,5,5] + 5000 a[1,6,5,5] + 5000 a[1,7,5,5]
   + 5000 a[1,8,5,5] + 5000 a[1,9,5,5] - Q[1] <= 0
 R2846: 5000 a[2,2,5,5] + 5000 a[2,3,5,5] + 5000 a[2,4,5,5]
   + 5000 a[2,5,5,5] + 5000 a[2,6,5,5] + 5000 a[2,7,5,5] + 5000 a[2,8,5,5]
   + 5000 a[2,9,5,5] - Q[2] <= 0
 R2847: 5000 a[3,3,5,5] + 5000 a[3,4,5,5] + 5000 a[3,5,5,5]
   + 5000 a[3,6,5,5] + 5000 a[3,7,5,5] + 5000 a[3,8,5,5] + 5000 a[3,9,5,5]
   - Q[3] <= 0
 R2848: 5000 a[4,4,5,5] + 5000 a[4,5,5,5] + 5000 a[4,6,5,5]
   + 5000 a[4,7,5,5] + 5000 a[4,8,5,5] + 5000 a[4,9,5,5] - Q[4] <= 0
 R2849: 5000 a[5,5,5,5] + 5000 a[5,6,5,5] + 5000 a[5,7,5,5]
   + 5000 a[5,8,5,5] + 5000 a[5,9,5,5] - Q[5] <= 0
 R2850: 5000 a[6,6,5,5] + 5000 a[6,7,5,5] + 5000 a[6,8,5,5]
   + 5000 a[6,9,5,5] - Q[6] <= 0
 R2851: 5000 a[7,7,5,5] + 5000 a[7,8,5,5] + 5000 a[7,9,5,5] - Q[7] <= 0
 R2852: 5000 a[8,8,5,5] + 5000 a[8,9,5,5] - Q[8] <= 0
 R2853: 5000 a[9,9,5,5] - Q[9] <= 0
 R2854: 10000 a[0,0,0,6] + 10000 a[0,1,0,6] + 10000 a[0,2,0,6]
   + 10000 a[0,3,0,6] + 10000 a[0,4,0,6] + 10000 a[0,5,0,6]
   + 10000 a[0,6,0,6] + 10000 a[0,7,0,6] + 10000 a[0,8,0,6]
   + 10000 a[0,9,0,6] + 10000 a[0,0,2,6] + 10000 a[0,1,2,6]
   + 10000 a[0,2,2,6] + 10000 a[0,3,2,6] + 10000 a[0,4,2,6]
   + 10000 a[0,5,2,6] + 10000 a[0,6,2,6] + 10000 a[0,7,2,6]
   + 10000 a[0,8,2,6] + 10000 a[0,9,2,6] + 10000 a[0,0,6,6]
   + 10000 a[0,1,6,6] + 10000 a[0,2,6,6] + 10000 a[0,3,6,6]
   + 10000 a[0,4,6,6] + 10000 a[0,5,6,6] + 10000 a[0,6,6,6]
   + 10000 a[0,7,6,6] + 10000 a[0,8,6,6] + 10000 a[0,9,6,6] - Q[0] <= 0
 R2855: 10000 a[1,1,0,6] + 10000 a[1,2,0,6] + 10000 a[1,3,0,6]
   + 10000 a[1,4,0,6] + 10000 a[1,5,0,6] + 10000 a[1,6,0,6]
   + 10000 a[1,7,0,6] + 10000 a[1,8,0,6] + 10000 a[1,9,0,6]
   + 10000 a[1,1,2,6] + 10000 a[1,2,2,6] + 10000 a[1,3,2,6]
   + 10000 a[1,4,2,6] + 10000 a[1,5,2,6] + 10000 a[1,6,2,6]
   + 10000 a[1,7,2,6] + 10000 a[1,8,2,6] + 10000 a[1,9,2,6]
   + 10000 a[1,1,6,6] + 10000 a[1,2,6,6] + 10000 a[1,3,6,6]
   + 10000 a[1,4,6,6] + 10000 a[1,5,6,6] + 10000 a[1,6,6,6]
   + 10000 a[1,7,6,6] + 10000 a[1,8,6,6] + 10000 a[1,9,6,6] - Q[1] <= 0
 R2856: 10000 a[2,2,0,6] + 10000 a[2,3,0,6] + 10000 a[2,4,0,6]
   + 10000 a[2,5,0,6] + 10000 a[2,6,0,6] + 10000 a[2,7,0,6]
   + 10000 a[2,8,0,6] + 10000 a[2,9,0,6] + 10000 a[2,2,2,6]
   + 10000 a[2,3,2,6] + 10000 a[2,4,2,6] + 10000 a[2,5,2,6]
   + 10000 a[2,6,2,6] + 10000 a[2,7,2,6] + 10000 a[2,8,2,6]
   + 10000 a[2,9,2,6] + 10000 a[2,2,6,6] + 10000 a[2,3,6,6]
   + 10000 a[2,4,6,6] + 10000 a[2,5,6,6] + 10000 a[2,6,6,6]
   + 10000 a[2,7,6,6] + 10000 a[2,8,6,6] + 10000 a[2,9,6,6] - Q[2] <= 0
 R2857: 10000 a[3,3,0,6] + 10000 a[3,4,0,6] + 10000 a[3,5,0,6]
   + 10000 a[3,6,0,6] + 10000 a[3,7,0,6] + 10000 a[3,8,0,6]
   + 10000 a[3,9,0,6] + 10000 a[3,3,2,6] + 10000 a[3,4,2,6]
   + 10000 a[3,5,2,6] + 10000 a[3,6,2,6] + 10000 a[3,7,2,6]
   + 10000 a[3,8,2,6] + 10000 a[3,9,2,6] + 10000 a[3,3,6,6]
   + 10000 a[3,4,6,6] + 10000 a[3,5,6,6] + 10000 a[3,6,6,6]
   + 10000 a[3,7,6,6] + 10000 a[3,8,6,6] + 10000 a[3,9,6,6] - Q[3] <= 0
 R2858: 10000 a[4,4,0,6] + 10000 a[4,5,0,6] + 10000 a[4,6,0,6]
   + 10000 a[4,7,0,6] + 10000 a[4,8,0,6] + 10000 a[4,9,0,6]
   + 10000 a[4,4,2,6] + 10000 a[4,5,2,6] + 10000 a[4,6,2,6]
   + 10000 a[4,7,2,6] + 10000 a[4,8,2,6] + 10000 a[4,9,2,6]
   + 10000 a[4,4,6,6] + 10000 a[4,5,6,6] + 10000 a[4,6,6,6]
   + 10000 a[4,7,6,6] + 10000 a[4,8,6,6] + 10000 a[4,9,6,6] - Q[4] <= 0
 R2859: 10000 a[5,5,0,6] + 10000 a[5,6,0,6] + 10000 a[5,7,0,6]
   + 10000 a[5,8,0,6] + 10000 a[5,9,0,6] + 10000 a[5,5,2,6]
   + 10000 a[5,6,2,6] + 10000 a[5,7,2,6] + 10000 a[5,8,2,6]
   + 10000 a[5,9,2,6] + 10000 a[5,5,6,6] + 10000 a[5,6,6,6]
   + 10000 a[5,7,6,6] + 10000 a[5,8,6,6] + 10000 a[5,9,6,6] - Q[5] <= 0
 R2860: 10000 a[6,6,0,6] + 10000 a[6,7,0,6] + 10000 a[6,8,0,6]
   + 10000 a[6,9,0,6] + 10000 a[6,6,2,6] + 10000 a[6,7,2,6]
   + 10000 a[6,8,2,6] + 10000 a[6,9,2,6] + 10000 a[6,6,6,6]
   + 10000 a[6,7,6,6] + 10000 a[6,8,6,6] + 10000 a[6,9,6,6] - Q[6] <= 0
 R2861: 10000 a[7,7,0,6] + 10000 a[7,8,0,6] + 10000 a[7,9,0,6]
   + 10000 a[7,7,2,6] + 10000 a[7,8,2,6] + 10000 a[7,9,2,6]
   + 10000 a[7,7,6,6] + 10000 a[7,8,6,6] + 10000 a[7,9,6,6] - Q[7] <= 0
 R2862: 10000 a[8,8,0,6] + 10000 a[8,9,0,6] + 10000 a[8,8,2,6]
   + 10000 a[8,9,2,6] + 10000 a[8,8,6,6] + 10000 a[8,9,6,6] - Q[8] <= 0
 R2863: 10000 a[9,9,0,6] + 10000 a[9,9,2,6] + 10000 a[9,9,6,6] - Q[9] <= 0
 R2864: 20000 a[0,0,1,7] + 20000 a[0,1,1,7] + 20000 a[0,2,1,7]
   + 20000 a[0,3,1,7] + 20000 a[0,4,1,7] + 20000 a[0,5,1,7]
   + 20000 a[0,6,1,7] + 20000 a[0,7,1,7] + 20000 a[0,8,1,7]
   + 20000 a[0,9,1,7] + 20000 a[0,0,3,7] + 20000 a[0,1,3,7]
   + 20000 a[0,2,3,7] + 20000 a[0,3,3,7] + 20000 a[0,4,3,7]
   + 20000 a[0,5,3,7] + 20000 a[0,6,3,7] + 20000 a[0,7,3,7]
   + 20000 a[0,8,3,7] + 20000 a[0,9,3,7] + 20000 a[0,0,7,7]
   + 20000 a[0,1,7,7] + 20000 a[0,2,7,7] + 20000 a[0,3,7,7]
   + 20000 a[0,4,7,7] + 20000 a[0,5,7,7] + 20000 a[0,6,7,7]
   + 20000 a[0,7,7,7] + 20000 a[0,8,7,7] + 20000 a[0,9,7,7] - Q[0] <= 0
 R2865: 20000 a[1,1,1,7] + 20000 a[1,2,1,7] + 20000 a[1,3,1,7]
   + 20000 a[1,4,1,7] + 20000 a[1,5,1,7] + 20000 a[1,6,1,7]
   + 20000 a[1,7,1,7] + 20000 a[1,8,1,7] + 20000 a[1,9,1,7]
   + 20000 a[1,1,3,7] + 20000 a[1,2,3,7] + 20000 a[1,3,3,7]
   + 20000 a[1,4,3,7] + 20000 a[1,5,3,7] + 20000 a[1,6,3,7]
   + 20000 a[1,7,3,7] + 20000 a[1,8,3,7] + 20000 a[1,9,3,7]
   + 20000 a[1,1,7,7] + 20000 a[1,2,7,7] + 20000 a[1,3,7,7]
   + 20000 a[1,4,7,7] + 20000 a[1,5,7,7] + 20000 a[1,6,7,7]
   + 20000 a[1,7,7,7] + 20000 a[1,8,7,7] + 20000 a[1,9,7,7] - Q[1] <= 0
 R2866: 20000 a[2,2,1,7] + 20000 a[2,3,1,7] + 20000 a[2,4,1,7]
   + 20000 a[2,5,1,7] + 20000 a[2,6,1,7] + 20000 a[2,7,1,7]
   + 20000 a[2,8,1,7] + 20000 a[2,9,1,7] + 20000 a[2,2,3,7]
   + 20000 a[2,3,3,7] + 20000 a[2,4,3,7] + 20000 a[2,5,3,7]
   + 20000 a[2,6,3,7] + 20000 a[2,7,3,7] + 20000 a[2,8,3,7]
   + 20000 a[2,9,3,7] + 20000 a[2,2,7,7] + 20000 a[2,3,7,7]
   + 20000 a[2,4,7,7] + 20000 a[2,5,7,7] + 20000 a[2,6,7,7]
   + 20000 a[2,7,7,7] + 20000 a[2,8,7,7] + 20000 a[2,9,7,7] - Q[2] <= 0
 R2867: 20000 a[3,3,1,7] + 20000 a[3,4,1,7] + 20000 a[3,5,1,7]
   + 20000 a[3,6,1,7] + 20000 a[3,7,1,7] + 20000 a[3,8,1,7]
   + 20000 a[3,9,1,7] + 20000 a[3,3,3,7] + 20000 a[3,4,3,7]
   + 20000 a[3,5,3,7] + 20000 a[3,6,3,7] + 20000 a[3,7,3,7]
   + 20000 a[3,8,3,7] + 20000 a[3,9,3,7] + 20000 a[3,3,7,7]
   + 20000 a[3,4,7,7] + 20000 a[3,5,7,7] + 20000 a[3,6,7,7]
   + 20000 a[3,7,7,7] + 20000 a[3,8,7,7] + 20000 a[3,9,7,7] - Q[3] <= 0
 R2868: 20000 a[4,4,1,7] + 20000 a[4,5,1,7] + 20000 a[4,6,1,7]
   + 20000 a[4,7,1,7] + 20000 a[4,8,1,7] + 20000 a[4,9,1,7]
   + 20000 a[4,4,3,7] + 20000 a[4,5,3,7] + 20000 a[4,6,3,7]
   + 20000 a[4,7,3,7] + 20000 a[4,8,3,7] + 20000 a[4,9,3,7]
   + 20000 a[4,4,7,7] + 20000 a[4,5,7,7] + 20000 a[4,6,7,7]
   + 20000 a[4,7,7,7] + 20000 a[4,8,7,7] + 20000 a[4,9,7,7] - Q[4] <= 0
 R2869: 20000 a[5,5,1,7] + 20000 a[5,6,1,7] + 20000 a[5,7,1,7]
   + 20000 a[5,8,1,7] + 20000 a[5,9,1,7] + 20000 a[5,5,3,7]
   + 20000 a[5,6,3,7] + 20000 a[5,7,3,7] + 20000 a[5,8,3,7]
   + 20000 a[5,9,3,7] + 20000 a[5,5,7,7] + 20000 a[5,6,7,7]
   + 20000 a[5,7,7,7] + 20000 a[5,8,7,7] + 20000 a[5,9,7,7] - Q[5] <= 0
 R2870: 20000 a[6,6,1,7] + 20000 a[6,7,1,7] + 20000 a[6,8,1,7]
   + 20000 a[6,9,1,7] + 20000 a[6,6,3,7] + 20000 a[6,7,3,7]
   + 20000 a[6,8,3,7] + 20000 a[6,9,3,7] + 20000 a[6,6,7,7]
   + 20000 a[6,7,7,7] + 20000 a[6,8,7,7] + 20000 a[6,9,7,7] - Q[6] <= 0
 R2871: 20000 a[7,7,1,7] + 20000 a[7,8,1,7] + 20000 a[7,9,1,7]
   + 20000 a[7,7,3,7] + 20000 a[7,8,3,7] + 20000 a[7,9,3,7]
   + 20000 a[7,7,7,7] + 20000 a[7,8,7,7] + 20000 a[7,9,7,7] - Q[7] <= 0
 R2872: 20000 a[8,8,1,7] + 20000 a[8,9,1,7] + 20000 a[8,8,3,7]
   + 20000 a[8,9,3,7] + 20000 a[8,8,7,7] + 20000 a[8,9,7,7] - Q[8] <= 0
 R2873: 20000 a[9,9,1,7] + 20000 a[9,9,3,7] + 20000 a[9,9,7,7] - Q[9] <= 0
 R2874: 10000 a[0,0,0,8] + 10000 a[0,1,0,8] + 10000 a[0,2,0,8]
   + 10000 a[0,3,0,8] + 10000 a[0,4,0,8] + 10000 a[0,5,0,8]
   + 10000 a[0,6,0,8] + 10000 a[0,7,0,8] + 10000 a[0,8,0,8]
   + 10000 a[0,9,0,8] + 10000 a[0,0,2,8] + 10000 a[0,1,2,8]
   + 10000 a[0,2,2,8] + 10000 a[0,3,2,8] + 10000 a[0,4,2,8]
   + 10000 a[0,5,2,8] + 10000 a[0,6,2,8] + 10000 a[0,7,2,8]
   + 10000 a[0,8,2,8] + 10000 a[0,9,2,8] + 10000 a[0,0,6,8]
   + 10000 a[0,1,6,8] + 10000 a[0,2,6,8] + 10000 a[0,3,6,8]
   + 10000 a[0,4,6,8] + 10000 a[0,5,6,8] + 10000 a[0,6,6,8]
   + 10000 a[0,7,6,8] + 10000 a[0,8,6,8] + 10000 a[0,9,6,8]
   + 10000 a[0,0,8,8] + 10000 a[0,1,8,8] + 10000 a[0,2,8,8]
   + 10000 a[0,3,8,8] + 10000 a[0,4,8,8] + 10000 a[0,5,8,8]
   + 10000 a[0,6,8,8] + 10000 a[0,7,8,8] + 10000 a[0,8,8,8]
   + 10000 a[0,9,8,8] - Q[0] <= 0
 R2875: 10000 a[1,1,0,8] + 10000 a[1,2,0,8] + 10000 a[1,3,0,8]
   + 10000 a[1,4,0,8] + 10000 a[1,5,0,8] + 10000 a[1,6,0,8]
   + 10000 a[1,7,0,8] + 10000 a[1,8,0,8] + 10000 a[1,9,0,8]
   + 10000 a[1,1,2,8] + 10000 a[1,2,2,8] + 10000 a[1,3,2,8]
   + 10000 a[1,4,2,8] + 10000 a[1,5,2,8] + 10000 a[1,6,2,8]
   + 10000 a[1,7,2,8] + 10000 a[1,8,2,8] + 10000 a[1,9,2,8]
   + 10000 a[1,1,6,8] + 10000 a[1,2,6,8] + 10000 a[1,3,6,8]
   + 10000 a[1,4,6,8] + 10000 a[1,5,6,8] + 10000 a[1,6,6,8]
   + 10000 a[1,7,6,8] + 10000 a[1,8,6,8] + 10000 a[1,9,6,8]
   + 10000 a[1,1,8,8] + 10000 a[1,2,8,8] + 10000 a[1,3,8,8]
   + 10000 a[1,4,8,8] + 10000 a[1,5,8,8] + 10000 a[1,6,8,8]
   + 10000 a[1,7,8,8] + 10000 a[1,8,8,8] + 10000 a[1,9,8,8] - Q[1] <= 0
 R2876: 10000 a[2,2,0,8] + 10000 a[2,3,0,8] + 10000 a[2,4,0,8]
   + 10000 a[2,5,0,8] + 10000 a[2,6,0,8] + 10000 a[2,7,0,8]
   + 10000 a[2,8,0,8] + 10000 a[2,9,0,8] + 10000 a[2,2,2,8]
   + 10000 a[2,3,2,8] + 10000 a[2,4,2,8] + 10000 a[2,5,2,8]
   + 10000 a[2,6,2,8] + 10000 a[2,7,2,8] + 10000 a[2,8,2,8]
   + 10000 a[2,9,2,8] + 10000 a[2,2,6,8] + 10000 a[2,3,6,8]
   + 10000 a[2,4,6,8] + 10000 a[2,5,6,8] + 10000 a[2,6,6,8]
   + 10000 a[2,7,6,8] + 10000 a[2,8,6,8] + 10000 a[2,9,6,8]
   + 10000 a[2,2,8,8] + 10000 a[2,3,8,8] + 10000 a[2,4,8,8]
   + 10000 a[2,5,8,8] + 10000 a[2,6,8,8] + 10000 a[2,7,8,8]
   + 10000 a[2,8,8,8] + 10000 a[2,9,8,8] - Q[2] <= 0
 R2877: 10000 a[3,3,0,8] + 10000 a[3,4,0,8] + 10000 a[3,5,0,8]
   + 10000 a[3,6,0,8] + 10000 a[3,7,0,8] + 10000 a[3,8,0,8]
   + 10000 a[3,9,0,8] + 10000 a[3,3,2,8] + 10000 a[3,4,2,8]
   + 10000 a[3,5,2,8] + 10000 a[3,6,2,8] + 10000 a[3,7,2,8]
   + 10000 a[3,8,2,8] + 10000 a[3,9,2,8] + 10000 a[3,3,6,8]
   + 10000 a[3,4,6,8] + 10000 a[3,5,6,8] + 10000 a[3,6,6,8]
   + 10000 a[3,7,6,8] + 10000 a[3,8,6,8] + 10000 a[3,9,6,8]
   + 10000 a[3,3,8,8] + 10000 a[3,4,8,8] + 10000 a[3,5,8,8]
   + 10000 a[3,6,8,8] + 10000 a[3,7,8,8] + 10000 a[3,8,8,8]
   + 10000 a[3,9,8,8] - Q[3] <= 0
 R2878: 10000 a[4,4,0,8] + 10000 a[4,5,0,8] + 10000 a[4,6,0,8]
   + 10000 a[4,7,0,8] + 10000 a[4,8,0,8] + 10000 a[4,9,0,8]
   + 10000 a[4,4,2,8] + 10000 a[4,5,2,8] + 10000 a[4,6,2,8]
   + 10000 a[4,7,2,8] + 10000 a[4,8,2,8] + 10000 a[4,9,2,8]
   + 10000 a[4,4,6,8] + 10000 a[4,5,6,8] + 10000 a[4,6,6,8]
   + 10000 a[4,7,6,8] + 10000 a[4,8,6,8] + 10000 a[4,9,6,8]
   + 10000 a[4,4,8,8] + 10000 a[4,5,8,8] + 10000 a[4,6,8,8]
   + 10000 a[4,7,8,8] + 10000 a[4,8,8,8] + 10000 a[4,9,8,8] - Q[4] <= 0
 R2879: 10000 a[5,5,0,8] + 10000 a[5,6,0,8] + 10000 a[5,7,0,8]
   + 10000 a[5,8,0,8] + 10000 a[5,9,0,8] + 10000 a[5,5,2,8]
   + 10000 a[5,6,2,8] + 10000 a[5,7,2,8] + 10000 a[5,8,2,8]
   + 10000 a[5,9,2,8] + 10000 a[5,5,6,8] + 10000 a[5,6,6,8]
   + 10000 a[5,7,6,8] + 10000 a[5,8,6,8] + 10000 a[5,9,6,8]
   + 10000 a[5,5,8,8] + 10000 a[5,6,8,8] + 10000 a[5,7,8,8]
   + 10000 a[5,8,8,8] + 10000 a[5,9,8,8] - Q[5] <= 0
 R2880: 10000 a[6,6,0,8] + 10000 a[6,7,0,8] + 10000 a[6,8,0,8]
   + 10000 a[6,9,0,8] + 10000 a[6,6,2,8] + 10000 a[6,7,2,8]
   + 10000 a[6,8,2,8] + 10000 a[6,9,2,8] + 10000 a[6,6,6,8]
   + 10000 a[6,7,6,8] + 10000 a[6,8,6,8] + 10000 a[6,9,6,8]
   + 10000 a[6,6,8,8] + 10000 a[6,7,8,8] + 10000 a[6,8,8,8]
   + 10000 a[6,9,8,8] - Q[6] <= 0
 R2881: 10000 a[7,7,0,8] + 10000 a[7,8,0,8] + 10000 a[7,9,0,8]
   + 10000 a[7,7,2,8] + 10000 a[7,8,2,8] + 10000 a[7,9,2,8]
   + 10000 a[7,7,6,8] + 10000 a[7,8,6,8] + 10000 a[7,9,6,8]
   + 10000 a[7,7,8,8] + 10000 a[7,8,8,8] + 10000 a[7,9,8,8] - Q[7] <= 0
 R2882: 10000 a[8,8,0,8] + 10000 a[8,9,0,8] + 10000 a[8,8,2,8]
   + 10000 a[8,9,2,8] + 10000 a[8,8,6,8] + 10000 a[8,9,6,8]
   + 10000 a[8,8,8,8] + 10000 a[8,9,8,8] - Q[8] <= 0
 R2883: 10000 a[9,9,0,8] + 10000 a[9,9,2,8] + 10000 a[9,9,6,8]
   + 10000 a[9,9,8,8] - Q[9] <= 0
 R2884: 10000 a[0,0,0,9] + 10000 a[0,1,0,9] + 10000 a[0,2,0,9]
   + 10000 a[0,3,0,9] + 10000 a[0,4,0,9] + 10000 a[0,5,0,9]
   + 10000 a[0,6,0,9] + 10000 a[0,7,0,9] + 10000 a[0,8,0,9]
   + 10000 a[0,9,0,9] + 10000 a[0,0,2,9] + 10000 a[0,1,2,9]
   + 10000 a[0,2,2,9] + 10000 a[0,3,2,9] + 10000 a[0,4,2,9]
   + 10000 a[0,5,2,9] + 10000 a[0,6,2,9] + 10000 a[0,7,2,9]
   + 10000 a[0,8,2,9] + 10000 a[0,9,2,9] + 10000 a[0,0,6,9]
   + 10000 a[0,1,6,9] + 10000 a[0,2,6,9] + 10000 a[0,3,6,9]
   + 10000 a[0,4,6,9] + 10000 a[0,5,6,9] + 10000 a[0,6,6,9]
   + 10000 a[0,7,6,9] + 10000 a[0,8,6,9] + 10000 a[0,9,6,9]
   + 10000 a[0,0,8,9] + 10000 a[0,1,8,9] + 10000 a[0,2,8,9]
   + 10000 a[0,3,8,9] + 10000 a[0,4,8,9] + 10000 a[0,5,8,9]
   + 10000 a[0,6,8,9] + 10000 a[0,7,8,9] + 10000 a[0,8,8,9]
   + 10000 a[0,9,8,9] + 10000 a[0,0,9,9] + 10000 a[0,1,9,9]
   + 10000 a[0,2,9,9] + 10000 a[0,3,9,9] + 10000 a[0,4,9,9]
   + 10000 a[0,5,9,9] + 10000 a[0,6,9,9] + 10000 a[0,7,9,9]
   + 10000 a[0,8,9,9] + 10000 a[0,9,9,9] - Q[0] <= 0
 R2885: 10000 a[1,1,0,9] + 10000 a[1,2,0,9] + 10000 a[1,3,0,9]
   + 10000 a[1,4,0,9] + 10000 a[1,5,0,9] + 10000 a[1,6,0,9]
   + 10000 a[1,7,0,9] + 10000 a[1,8,0,9] + 10000 a[1,9,0,9]
   + 10000 a[1,1,2,9] + 10000 a[1,2,2,9] + 10000 a[1,3,2,9]
   + 10000 a[1,4,2,9] + 10000 a[1,5,2,9] + 10000 a[1,6,2,9]
   + 10000 a[1,7,2,9] + 10000 a[1,8,2,9] + 10000 a[1,9,2,9]
   + 10000 a[1,1,6,9] + 10000 a[1,2,6,9] + 10000 a[1,3,6,9]
   + 10000 a[1,4,6,9] + 10000 a[1,5,6,9] + 10000 a[1,6,6,9]
   + 10000 a[1,7,6,9] + 10000 a[1,8,6,9] + 10000 a[1,9,6,9]
   + 10000 a[1,1,8,9] + 10000 a[1,2,8,9] + 10000 a[1,3,8,9]
   + 10000 a[1,4,8,9] + 10000 a[1,5,8,9] + 10000 a[1,6,8,9]
   + 10000 a[1,7,8,9] + 10000 a[1,8,8,9] + 10000 a[1,9,8,9]
   + 10000 a[1,1,9,9] + 10000 a[1,2,9,9] + 10000 a[1,3,9,9]
   + 10000 a[1,4,9,9] + 10000 a[1,5,9,9] + 10000 a[1,6,9,9]
   + 10000 a[1,7,9,9] + 10000 a[1,8,9,9] + 10000 a[1,9,9,9] - Q[1] <= 0
 R2886: 10000 a[2,2,0,9] + 10000 a[2,3,0,9] + 10000 a[2,4,0,9]
   + 10000 a[2,5,0,9] + 10000 a[2,6,0,9] + 10000 a[2,7,0,9]
   + 10000 a[2,8,0,9] + 10000 a[2,9,0,9] + 10000 a[2,2,2,9]
   + 10000 a[2,3,2,9] + 10000 a[2,4,2,9] + 10000 a[2,5,2,9]
   + 10000 a[2,6,2,9] + 10000 a[2,7,2,9] + 10000 a[2,8,2,9]
   + 10000 a[2,9,2,9] + 10000 a[2,2,6,9] + 10000 a[2,3,6,9]
   + 10000 a[2,4,6,9] + 10000 a[2,5,6,9] + 10000 a[2,6,6,9]
   + 10000 a[2,7,6,9] + 10000 a[2,8,6,9] + 10000 a[2,9,6,9]
   + 10000 a[2,2,8,9] + 10000 a[2,3,8,9] + 10000 a[2,4,8,9]
   + 10000 a[2,5,8,9] + 10000 a[2,6,8,9] + 10000 a[2,7,8,9]
   + 10000 a[2,8,8,9] + 10000 a[2,9,8,9] + 10000 a[2,2,9,9]
   + 10000 a[2,3,9,9] + 10000 a[2,4,9,9] + 10000 a[2,5,9,9]
   + 10000 a[2,6,9,9] + 10000 a[2,7,9,9] + 10000 a[2,8,9,9]
   + 10000 a[2,9,9,9] - Q[2] <= 0
 R2887: 10000 a[3,3,0,9] + 10000 a[3,4,0,9] + 10000 a[3,5,0,9]
   + 10000 a[3,6,0,9] + 10000 a[3,7,0,9] + 10000 a[3,8,0,9]
   + 10000 a[3,9,0,9] + 10000 a[3,3,2,9] + 10000 a[3,4,2,9]
   + 10000 a[3,5,2,9] + 10000 a[3,6,2,9] + 10000 a[3,7,2,9]
   + 10000 a[3,8,2,9] + 10000 a[3,9,2,9] + 10000 a[3,3,6,9]
   + 10000 a[3,4,6,9] + 10000 a[3,5,6,9] + 10000 a[3,6,6,9]
   + 10000 a[3,7,6,9] + 10000 a[3,8,6,9] + 10000 a[3,9,6,9]
   + 10000 a[3,3,8,9] + 10000 a[3,4,8,9] + 10000 a[3,5,8,9]
   + 10000 a[3,6,8,9] + 10000 a[3,7,8,9] + 10000 a[3,8,8,9]
   + 10000 a[3,9,8,9] + 10000 a[3,3,9,9] + 10000 a[3,4,9,9]
   + 10000 a[3,5,9,9] + 10000 a[3,6,9,9] + 10000 a[3,7,9,9]
   + 10000 a[3,8,9,9] + 10000 a[3,9,9,9] - Q[3] <= 0
 R2888: 10000 a[4,4,0,9] + 10000 a[4,5,0,9] + 10000 a[4,6,0,9]
   + 10000 a[4,7,0,9] + 10000 a[4,8,0,9] + 10000 a[4,9,0,9]
   + 10000 a[4,4,2,9] + 10000 a[4,5,2,9] + 10000 a[4,6,2,9]
   + 10000 a[4,7,2,9] + 10000 a[4,8,2,9] + 10000 a[4,9,2,9]
   + 10000 a[4,4,6,9] + 10000 a[4,5,6,9] + 10000 a[4,6,6,9]
   + 10000 a[4,7,6,9] + 10000 a[4,8,6,9] + 10000 a[4,9,6,9]
   + 10000 a[4,4,8,9] + 10000 a[4,5,8,9] + 10000 a[4,6,8,9]
   + 10000 a[4,7,8,9] + 10000 a[4,8,8,9] + 10000 a[4,9,8,9]
   + 10000 a[4,4,9,9] + 10000 a[4,5,9,9] + 10000 a[4,6,9,9]
   + 10000 a[4,7,9,9] + 10000 a[4,8,9,9] + 10000 a[4,9,9,9] - Q[4] <= 0
 R2889: 10000 a[5,5,0,9] + 10000 a[5,6,0,9] + 10000 a[5,7,0,9]
   + 10000 a[5,8,0,9] + 10000 a[5,9,0,9] + 10000 a[5,5,2,9]
   + 10000 a[5,6,2,9] + 10000 a[5,7,2,9] + 10000 a[5,8,2,9]
   + 10000 a[5,9,2,9] + 10000 a[5,5,6,9] + 10000 a[5,6,6,9]
   + 10000 a[5,7,6,9] + 10000 a[5,8,6,9] + 10000 a[5,9,6,9]
   + 10000 a[5,5,8,9] + 10000 a[5,6,8,9] + 10000 a[5,7,8,9]
   + 10000 a[5,8,8,9] + 10000 a[5,9,8,9] + 10000 a[5,5,9,9]
   + 10000 a[5,6,9,9] + 10000 a[5,7,9,9] + 10000 a[5,8,9,9]
   + 10000 a[5,9,9,9] - Q[5] <= 0
 R2890: 10000 a[6,6,0,9] + 10000 a[6,7,0,9] + 10000 a[6,8,0,9]
   + 10000 a[6,9,0,9] + 10000 a[6,6,2,9] + 10000 a[6,7,2,9]
   + 10000 a[6,8,2,9] + 10000 a[6,9,2,9] + 10000 a[6,6,6,9]
   + 10000 a[6,7,6,9] + 10000 a[6,8,6,9] + 10000 a[6,9,6,9]
   + 10000 a[6,6,8,9] + 10000 a[6,7,8,9] + 10000 a[6,8,8,9]
   + 10000 a[6,9,8,9] + 10000 a[6,6,9,9] + 10000 a[6,7,9,9]
   + 10000 a[6,8,9,9] + 10000 a[6,9,9,9] - Q[6] <= 0
 R2891: 10000 a[7,7,0,9] + 10000 a[7,8,0,9] + 10000 a[7,9,0,9]
   + 10000 a[7,7,2,9] + 10000 a[7,8,2,9] + 10000 a[7,9,2,9]
   + 10000 a[7,7,6,9] + 10000 a[7,8,6,9] + 10000 a[7,9,6,9]
   + 10000 a[7,7,8,9] + 10000 a[7,8,8,9] + 10000 a[7,9,8,9]
   + 10000 a[7,7,9,9] + 10000 a[7,8,9,9] + 10000 a[7,9,9,9] - Q[7] <= 0
 R2892: 10000 a[8,8,0,9] + 10000 a[8,9,0,9] + 10000 a[8,8,2,9]
   + 10000 a[8,9,2,9] + 10000 a[8,8,6,9] + 10000 a[8,9,6,9]
   + 10000 a[8,8,8,9] + 10000 a[8,9,8,9] + 10000 a[8,8,9,9]
   + 10000 a[8,9,9,9] - Q[8] <= 0
 R2893: 10000 a[9,9,0,9] + 10000 a[9,9,2,9] + 10000 a[9,9,6,9]
   + 10000 a[9,9,8,9] + 10000 a[9,9,9,9] - Q[9] <= 0
 R2894: - a[0,0,0,0] - a[0,1,0,0] - a[0,2,0,0] - a[0,3,0,0] - a[0,4,0,0]
   - a[0,5,0,0] - a[0,6,0,0] - a[0,7,0,0] - a[0,8,0,0] - a[0,9,0,0]
   + 2 b[0,0,0] - x[0,0] <= 0
 R2895: - a[0,0,0,0] - a[0,1,0,0] - a[0,2,0,0] - a[0,3,0,0] - a[0,4,0,0]
   - a[0,5,0,0] - a[0,6,0,0] - a[0,7,0,0] - a[0,8,0,0] - a[0,9,0,0]
   + 2 b[0,0,1] - x[0,1] <= 0
 R2896: - a[0,0,0,0] - a[0,1,0,0] - a[0,2,0,0] - a[0,3,0,0] - a[0,4,0,0]
   - a[0,5,0,0] - a[0,6,0,0] - a[0,7,0,0] - a[0,8,0,0] - a[0,9,0,0]
   + 2 b[0,0,2] - x[0,2] <= 0
 R2897: - a[1,1,0,0] - a[1,2,0,0] - a[1,3,0,0] - a[1,4,0,0] - a[1,5,0,0]
   - a[1,6,0,0] - a[1,7,0,0] - a[1,8,0,0] - a[1,9,0,0] + 2 b[0,1,0]
   - x[1,0] <= 0
 R2898: - a[1,1,0,0] - a[1,2,0,0] - a[1,3,0,0] - a[1,4,0,0] - a[1,5,0,0]
   - a[1,6,0,0] - a[1,7,0,0] - a[1,8,0,0] - a[1,9,0,0] + 2 b[0,1,1]
   - x[1,1] <= 0
 R2899: - a[1,1,0,0] - a[1,2,0,0] - a[1,3,0,0] - a[1,4,0,0] - a[1,5,0,0]
   - a[1,6,0,0] - a[1,7,0,0] - a[1,8,0,0] - a[1,9,0,0] + 2 b[0,1,2]
   - x[1,2] <= 0
 R2900: - a[2,2,0,0] - a[2,3,0,0] - a[2,4,0,0] - a[2,5,0,0] - a[2,6,0,0]
   - a[2,7,0,0] - a[2,8,0,0] - a[2,9,0,0] + 2 b[0,2,0] - x[2,0] <= 0
 R2901: - a[2,2,0,0] - a[2,3,0,0] - a[2,4,0,0] - a[2,5,0,0] - a[2,6,0,0]
   - a[2,7,0,0] - a[2,8,0,0] - a[2,9,0,0] + 2 b[0,2,1] - x[2,1] <= 0
 R2902: - a[2,2,0,0] - a[2,3,0,0] - a[2,4,0,0] - a[2,5,0,0] - a[2,6,0,0]
   - a[2,7,0,0] - a[2,8,0,0] - a[2,9,0,0] + 2 b[0,2,2] - x[2,2] <= 0
 R2903: - a[3,3,0,0] - a[3,4,0,0] - a[3,5,0,0] - a[3,6,0,0] - a[3,7,0,0]
   - a[3,8,0,0] - a[3,9,0,0] + 2 b[0,3,0] - x[3,0] <= 0
 R2904: - a[3,3,0,0] - a[3,4,0,0] - a[3,5,0,0] - a[3,6,0,0] - a[3,7,0,0]
   - a[3,8,0,0] - a[3,9,0,0] + 2 b[0,3,1] - x[3,1] <= 0
 R2905: - a[3,3,0,0] - a[3,4,0,0] - a[3,5,0,0] - a[3,6,0,0] - a[3,7,0,0]
   - a[3,8,0,0] - a[3,9,0,0] + 2 b[0,3,2] - x[3,2] <= 0
 R2906: - a[4,4,0,0] - a[4,5,0,0] - a[4,6,0,0] - a[4,7,0,0] - a[4,8,0,0]
   - a[4,9,0,0] + 2 b[0,4,0] - x[4,0] <= 0
 R2907: - a[4,4,0,0] - a[4,5,0,0] - a[4,6,0,0] - a[4,7,0,0] - a[4,8,0,0]
   - a[4,9,0,0] + 2 b[0,4,1] - x[4,1] <= 0
 R2908: - a[4,4,0,0] - a[4,5,0,0] - a[4,6,0,0] - a[4,7,0,0] - a[4,8,0,0]
   - a[4,9,0,0] + 2 b[0,4,2] - x[4,2] <= 0
 R2909: - a[5,5,0,0] - a[5,6,0,0] - a[5,7,0,0] - a[5,8,0,0] - a[5,9,0,0]
   + 2 b[0,5,0] - x[5,0] <= 0
 R2910: - a[5,5,0,0] - a[5,6,0,0] - a[5,7,0,0] - a[5,8,0,0] - a[5,9,0,0]
   + 2 b[0,5,1] - x[5,1] <= 0
 R2911: - a[5,5,0,0] - a[5,6,0,0] - a[5,7,0,0] - a[5,8,0,0] - a[5,9,0,0]
   + 2 b[0,5,2] - x[5,2] <= 0
 R2912: - a[6,6,0,0] - a[6,7,0,0] - a[6,8,0,0] - a[6,9,0,0] + 2 b[0,6,0]
   - x[6,0] <= 0
 R2913: - a[6,6,0,0] - a[6,7,0,0] - a[6,8,0,0] - a[6,9,0,0] + 2 b[0,6,1]
   - x[6,1] <= 0
 R2914: - a[6,6,0,0] - a[6,7,0,0] - a[6,8,0,0] - a[6,9,0,0] + 2 b[0,6,2]
   - x[6,2] <= 0
 R2915: - a[7,7,0,0] - a[7,8,0,0] - a[7,9,0,0] + 2 b[0,7,0] - x[7,0] <= 0
 R2916: - a[7,7,0,0] - a[7,8,0,0] - a[7,9,0,0] + 2 b[0,7,1] - x[7,1] <= 0
 R2917: - a[7,7,0,0] - a[7,8,0,0] - a[7,9,0,0] + 2 b[0,7,2] - x[7,2] <= 0
 R2918: - a[8,8,0,0] - a[8,9,0,0] + 2 b[0,8,0] - x[8,0] <= 0
 R2919: - a[8,8,0,0] - a[8,9,0,0] + 2 b[0,8,1] - x[8,1] <= 0
 R2920: - a[8,8,0,0] - a[8,9,0,0] + 2 b[0,8,2] - x[8,2] <= 0
 R2921: - a[9,9,0,0] + 2 b[0,9,0] - x[9,0] <= 0
 R2922: - a[9,9,0,0] + 2 b[0,9,1] - x[9,1] <= 0
 R2923: - a[9,9,0,0] + 2 b[0,9,2] - x[9,2] <= 0
 R2924: - a[0,0,0,2] - a[0,1,0,2] - a[0,2,0,2] - a[0,3,0,2] - a[0,4,0,2]
   - a[0,5,0,2] - a[0,6,0,2] - a[0,7,0,2] - a[0,8,0,2] - a[0,9,0,2]
   - a[0,0,2,2] - a[0,1,2,2] - a[0,2,2,2] - a[0,3,2,2] - a[0,4,2,2]
   - a[0,5,2,2] - a[0,6,2,2] - a[0,7,2,2] - a[0,8,2,2] - a[0,9,2,2]
   + 2 b[2,0,0] - x[0,0] <= 0
 R2925: - a[0,0,0,2] - a[0,1,0,2] - a[0,2,0,2] - a[0,3,0,2] - a[0,4,0,2]
   - a[0,5,0,2] - a[0,6,0,2] - a[0,7,0,2] - a[0,8,0,2] - a[0,9,0,2]
   - a[0,0,2,2] - a[0,1,2,2] - a[0,2,2,2] - a[0,3,2,2] - a[0,4,2,2]
   - a[0,5,2,2] - a[0,6,2,2] - a[0,7,2,2] - a[0,8,2,2] - a[0,9,2,2]
   + 2 b[2,0,1] - x[0,1] <= 0
 R2926: - a[0,0,0,2] - a[0,1,0,2] - a[0,2,0,2] - a[0,3,0,2] - a[0,4,0,2]
   - a[0,5,0,2] - a[0,6,0,2] - a[0,7,0,2] - a[0,8,0,2] - a[0,9,0,2]
   - a[0,0,2,2] - a[0,1,2,2] - a[0,2,2,2] - a[0,3,2,2] - a[0,4,2,2]
   - a[0,5,2,2] - a[0,6,2,2] - a[0,7,2,2] - a[0,8,2,2] - a[0,9,2,2]
   + 2 b[2,0,2] - x[0,2] <= 0
 R2927: - a[1,1,0,2] - a[1,2,0,2] - a[1,3,0,2] - a[1,4,0,2] - a[1,5,0,2]
   - a[1,6,0,2] - a[1,7,0,2] - a[1,8,0,2] - a[1,9,0,2] - a[1,1,2,2]
   - a[1,2,2,2] - a[1,3,2,2] - a[1,4,2,2] - a[1,5,2,2] - a[1,6,2,2]
   - a[1,7,2,2] - a[1,8,2,2] - a[1,9,2,2] + 2 b[2,1,0] - x[1,0] <= 0
 R2928: - a[1,1,0,2] - a[1,2,0,2] - a[1,3,0,2] - a[1,4,0,2] - a[1,5,0,2]
   - a[1,6,0,2] - a[1,7,0,2] - a[1,8,0,2] - a[1,9,0,2] - a[1,1,2,2]
   - a[1,2,2,2] - a[1,3,2,2] - a[1,4,2,2] - a[1,5,2,2] - a[1,6,2,2]
   - a[1,7,2,2] - a[1,8,2,2] - a[1,9,2,2] + 2 b[2,1,1] - x[1,1] <= 0
 R2929: - a[1,1,0,2] - a[1,2,0,2] - a[1,3,0,2] - a[1,4,0,2] - a[1,5,0,2]
   - a[1,6,0,2] - a[1,7,0,2] - a[1,8,0,2] - a[1,9,0,2] - a[1,1,2,2]
   - a[1,2,2,2] - a[1,3,2,2] - a[1,4,2,2] - a[1,5,2,2] - a[1,6,2,2]
   - a[1,7,2,2] - a[1,8,2,2] - a[1,9,2,2] + 2 b[2,1,2] - x[1,2] <= 0
 R2930: - a[2,2,0,2] - a[2,3,0,2] - a[2,4,0,2] - a[2,5,0,2] - a[2,6,0,2]
   - a[2,7,0,2] - a[2,8,0,2] - a[2,9,0,2] - a[2,2,2,2] - a[2,3,2,2]
   - a[2,4,2,2] - a[2,5,2,2] - a[2,6,2,2] - a[2,7,2,2] - a[2,8,2,2]
   - a[2,9,2,2] + 2 b[2,2,0] - x[2,0] <= 0
 R2931: - a[2,2,0,2] - a[2,3,0,2] - a[2,4,0,2] - a[2,5,0,2] - a[2,6,0,2]
   - a[2,7,0,2] - a[2,8,0,2] - a[2,9,0,2] - a[2,2,2,2] - a[2,3,2,2]
   - a[2,4,2,2] - a[2,5,2,2] - a[2,6,2,2] - a[2,7,2,2] - a[2,8,2,2]
   - a[2,9,2,2] + 2 b[2,2,1] - x[2,1] <= 0
 R2932: - a[2,2,0,2] - a[2,3,0,2] - a[2,4,0,2] - a[2,5,0,2] - a[2,6,0,2]
   - a[2,7,0,2] - a[2,8,0,2] - a[2,9,0,2] - a[2,2,2,2] - a[2,3,2,2]
   - a[2,4,2,2] - a[2,5,2,2] - a[2,6,2,2] - a[2,7,2,2] - a[2,8,2,2]
   - a[2,9,2,2] + 2 b[2,2,2] - x[2,2] <= 0
 R2933: - a[3,3,0,2] - a[3,4,0,2] - a[3,5,0,2] - a[3,6,0,2] - a[3,7,0,2]
   - a[3,8,0,2] - a[3,9,0,2] - a[3,3,2,2] - a[3,4,2,2] - a[3,5,2,2]
   - a[3,6,2,2] - a[3,7,2,2] - a[3,8,2,2] - a[3,9,2,2] + 2 b[2,3,0]
   - x[3,0] <= 0
 R2934: - a[3,3,0,2] - a[3,4,0,2] - a[3,5,0,2] - a[3,6,0,2] - a[3,7,0,2]
   - a[3,8,0,2] - a[3,9,0,2] - a[3,3,2,2] - a[3,4,2,2] - a[3,5,2,2]
   - a[3,6,2,2] - a[3,7,2,2] - a[3,8,2,2] - a[3,9,2,2] + 2 b[2,3,1]
   - x[3,1] <= 0
 R2935: - a[3,3,0,2] - a[3,4,0,2] - a[3,5,0,2] - a[3,6,0,2] - a[3,7,0,2]
   - a[3,8,0,2] - a[3,9,0,2] - a[3,3,2,2] - a[3,4,2,2] - a[3,5,2,2]
   - a[3,6,2,2] - a[3,7,2,2] - a[3,8,2,2] - a[3,9,2,2] + 2 b[2,3,2]
   - x[3,2] <= 0
 R2936: - a[4,4,0,2] - a[4,5,0,2] - a[4,6,0,2] - a[4,7,0,2] - a[4,8,0,2]
   - a[4,9,0,2] - a[4,4,2,2] - a[4,5,2,2] - a[4,6,2,2] - a[4,7,2,2]
   - a[4,8,2,2] - a[4,9,2,2] + 2 b[2,4,0] - x[4,0] <= 0
 R2937: - a[4,4,0,2] - a[4,5,0,2] - a[4,6,0,2] - a[4,7,0,2] - a[4,8,0,2]
   - a[4,9,0,2] - a[4,4,2,2] - a[4,5,2,2] - a[4,6,2,2] - a[4,7,2,2]
   - a[4,8,2,2] - a[4,9,2,2] + 2 b[2,4,1] - x[4,1] <= 0
 R2938: - a[4,4,0,2] - a[4,5,0,2] - a[4,6,0,2] - a[4,7,0,2] - a[4,8,0,2]
   - a[4,9,0,2] - a[4,4,2,2] - a[4,5,2,2] - a[4,6,2,2] - a[4,7,2,2]
   - a[4,8,2,2] - a[4,9,2,2] + 2 b[2,4,2] - x[4,2] <= 0
 R2939: - a[5,5,0,2] - a[5,6,0,2] - a[5,7,0,2] - a[5,8,0,2] - a[5,9,0,2]
   - a[5,5,2,2] - a[5,6,2,2] - a[5,7,2,2] - a[5,8,2,2] - a[5,9,2,2]
   + 2 b[2,5,0] - x[5,0] <= 0
 R2940: - a[5,5,0,2] - a[5,6,0,2] - a[5,7,0,2] - a[5,8,0,2] - a[5,9,0,2]
   - a[5,5,2,2] - a[5,6,2,2] - a[5,7,2,2] - a[5,8,2,2] - a[5,9,2,2]
   + 2 b[2,5,1] - x[5,1] <= 0
 R2941: - a[5,5,0,2] - a[5,6,0,2] - a[5,7,0,2] - a[5,8,0,2] - a[5,9,0,2]
   - a[5,5,2,2] - a[5,6,2,2] - a[5,7,2,2] - a[5,8,2,2] - a[5,9,2,2]
   + 2 b[2,5,2] - x[5,2] <= 0
 R2942: - a[6,6,0,2] - a[6,7,0,2] - a[6,8,0,2] - a[6,9,0,2] - a[6,6,2,2]
   - a[6,7,2,2] - a[6,8,2,2] - a[6,9,2,2] + 2 b[2,6,0] - x[6,0] <= 0
 R2943: - a[6,6,0,2] - a[6,7,0,2] - a[6,8,0,2] - a[6,9,0,2] - a[6,6,2,2]
   - a[6,7,2,2] - a[6,8,2,2] - a[6,9,2,2] + 2 b[2,6,1] - x[6,1] <= 0
 R2944: - a[6,6,0,2] - a[6,7,0,2] - a[6,8,0,2] - a[6,9,0,2] - a[6,6,2,2]
   - a[6,7,2,2] - a[6,8,2,2] - a[6,9,2,2] + 2 b[2,6,2] - x[6,2] <= 0
 R2945: - a[7,7,0,2] - a[7,8,0,2] - a[7,9,0,2] - a[7,7,2,2] - a[7,8,2,2]
   - a[7,9,2,2] + 2 b[2,7,0] - x[7,0] <= 0
 R2946: - a[7,7,0,2] - a[7,8,0,2] - a[7,9,0,2] - a[7,7,2,2] - a[7,8,2,2]
   - a[7,9,2,2] + 2 b[2,7,1] - x[7,1] <= 0
 R2947: - a[7,7,0,2] - a[7,8,0,2] - a[7,9,0,2] - a[7,7,2,2] - a[7,8,2,2]
   - a[7,9,2,2] + 2 b[2,7,2] - x[7,2] <= 0
 R2948: - a[8,8,0,2] - a[8,9,0,2] - a[8,8,2,2] - a[8,9,2,2] + 2 b[2,8,0]
   - x[8,0] <= 0
 R2949: - a[8,8,0,2] - a[8,9,0,2] - a[8,8,2,2] - a[8,9,2,2] + 2 b[2,8,1]
   - x[8,1] <= 0
 R2950: - a[8,8,0,2] - a[8,9,0,2] - a[8,8,2,2] - a[8,9,2,2] + 2 b[2,8,2]
   - x[8,2] <= 0
 R2951: - a[9,9,0,2] - a[9,9,2,2] + 2 b[2,9,0] - x[9,0] <= 0
 R2952: - a[9,9,0,2] - a[9,9,2,2] + 2 b[2,9,1] - x[9,1] <= 0
 R2953: - a[9,9,0,2] - a[9,9,2,2] + 2 b[2,9,2] - x[9,2] <= 0
 R2954: - a[0,0,0,6] - a[0,1,0,6] - a[0,2,0,6] - a[0,3,0,6] - a[0,4,0,6]
   - a[0,5,0,6] - a[0,6,0,6] - a[0,7,0,6] - a[0,8,0,6] - a[0,9,0,6]
   - a[0,0,2,6] - a[0,1,2,6] - a[0,2,2,6] - a[0,3,2,6] - a[0,4,2,6]
   - a[0,5,2,6] - a[0,6,2,6] - a[0,7,2,6] - a[0,8,2,6] - a[0,9,2,6]
   - a[0,0,6,6] - a[0,1,6,6] - a[0,2,6,6] - a[0,3,6,6] - a[0,4,6,6]
   - a[0,5,6,6] - a[0,6,6,6] - a[0,7,6,6] - a[0,8,6,6] - a[0,9,6,6]
   + 2 b[6,0,0] - x[0,0] <= 0
 R2955: - a[0,0,0,6] - a[0,1,0,6] - a[0,2,0,6] - a[0,3,0,6] - a[0,4,0,6]
   - a[0,5,0,6] - a[0,6,0,6] - a[0,7,0,6] - a[0,8,0,6] - a[0,9,0,6]
   - a[0,0,2,6] - a[0,1,2,6] - a[0,2,2,6] - a[0,3,2,6] - a[0,4,2,6]
   - a[0,5,2,6] - a[0,6,2,6] - a[0,7,2,6] - a[0,8,2,6] - a[0,9,2,6]
   - a[0,0,6,6] - a[0,1,6,6] - a[0,2,6,6] - a[0,3,6,6] - a[0,4,6,6]
   - a[0,5,6,6] - a[0,6,6,6] - a[0,7,6,6] - a[0,8,6,6] - a[0,9,6,6]
   + 2 b[6,0,1] - x[0,1] <= 0
 R2956: - a[0,0,0,6] - a[0,1,0,6] - a[0,2,0,6] - a[0,3,0,6] - a[0,4,0,6]
   - a[0,5,0,6] - a[0,6,0,6] - a[0,7,0,6] - a[0,8,0,6] - a[0,9,0,6]
   - a[0,0,2,6] - a[0,1,2,6] - a[0,2,2,6] - a[0,3,2,6] - a[0,4,2,6]
   - a[0,5,2,6] - a[0,6,2,6] - a[0,7,2,6] - a[0,8,2,6] - a[0,9,2,6]
   - a[0,0,6,6] - a[0,1,6,6] - a[0,2,6,6] - a[0,3,6,6] - a[0,4,6,6]
   - a[0,5,6,6] - a[0,6,6,6] - a[0,7,6,6] - a[0,8,6,6] - a[0,9,6,6]
   + 2 b[6,0,2] - x[0,2] <= 0
 R2957: - a[1,1,0,6] - a[1,2,0,6] - a[1,3,0,6] - a[1,4,0,6] - a[1,5,0,6]
   - a[1,6,0,6] - a[1,7,0,6] - a[1,8,0,6] - a[1,9,0,6] - a[1,1,2,6]
   - a[1,2,2,6] - a[1,3,2,6] - a[1,4,2,6] - a[1,5,2,6] - a[1,6,2,6]
   - a[1,7,2,6] - a[1,8,2,6] - a[1,9,2,6] - a[1,1,6,6] - a[1,2,6,6]
   - a[1,3,6,6] - a[1,4,6,6] - a[1,5,6,6] - a[1,6,6,6] - a[1,7,6,6]
   - a[1,8,6,6] - a[1,9,6,6] + 2 b[6,1,0] - x[1,0] <= 0
 R2958: - a[1,1,0,6] - a[1,2,0,6] - a[1,3,0,6] - a[1,4,0,6] - a[1,5,0,6]
   - a[1,6,0,6] - a[1,7,0,6] - a[1,8,0,6] - a[1,9,0,6] - a[1,1,2,6]
   - a[1,2,2,6] - a[1,3,2,6] - a[1,4,2,6] - a[1,5,2,6] - a[1,6,2,6]
   - a[1,7,2,6] - a[1,8,2,6] - a[1,9,2,6] - a[1,1,6,6] - a[1,2,6,6]
   - a[1,3,6,6] - a[1,4,6,6] - a[1,5,6,6] - a[1,6,6,6] - a[1,7,6,6]
   - a[1,8,6,6] - a[1,9,6,6] + 2 b[6,1,1] - x[1,1] <= 0
 R2959: - a[1,1,0,6] - a[1,2,0,6] - a[1,3,0,6] - a[1,4,0,6] - a[1,5,0,6]
   - a[1,6,0,6] - a[1,7,0,6] - a[1,8,0,6] - a[1,9,0,6] - a[1,1,2,6]
   - a[1,2,2,6] - a[1,3,2,6] - a[1,4,2,6] - a[1,5,2,6] - a[1,6,2,6]
   - a[1,7,2,6] - a[1,8,2,6] - a[1,9,2,6] - a[1,1,6,6] - a[1,2,6,6]
   - a[1,3,6,6] - a[1,4,6,6] - a[1,5,6,6] - a[1,6,6,6] - a[1,7,6,6]
   - a[1,8,6,6] - a[1,9,6,6] + 2 b[6,1,2] - x[1,2] <= 0
 R2960: - a[2,2,0,6] - a[2,3,0,6] - a[2,4,0,6] - a[2,5,0,6] - a[2,6,0,6]
   - a[2,7,0,6] - a[2,8,0,6] - a[2,9,0,6] - a[2,2,2,6] - a[2,3,2,6]
   - a[2,4,2,6] - a[2,5,2,6] - a[2,6,2,6] - a[2,7,2,6] - a[2,8,2,6]
   - a[2,9,2,6] - a[2,2,6,6] - a[2,3,6,6] - a[2,4,6,6] - a[2,5,6,6]
   - a[2,6,6,6] - a[2,7,6,6] - a[2,8,6,6] - a[2,9,6,6] + 2 b[6,2,0]
   - x[2,0] <= 0
 R2961: - a[2,2,0,6] - a[2,3,0,6] - a[2,4,0,6] - a[2,5,0,6] - a[2,6,0,6]
   - a[2,7,0,6] - a[2,8,0,6] - a[2,9,0,6] - a[2,2,2,6] - a[2,3,2,6]
   - a[2,4,2,6] - a[2,5,2,6] - a[2,6,2,6] - a[2,7,2,6] - a[2,8,2,6]
   - a[2,9,2,6] - a[2,2,6,6] - a[2,3,6,6] - a[2,4,6,6] - a[2,5,6,6]
   - a[2,6,6,6] - a[2,7,6,6] - a[2,8,6,6] - a[2,9,6,6] + 2 b[6,2,1]
   - x[2,1] <= 0
 R2962: - a[2,2,0,6] - a[2,3,0,6] - a[2,4,0,6] - a[2,5,0,6] - a[2,6,0,6]
   - a[2,7,0,6] - a[2,8,0,6] - a[2,9,0,6] - a[2,2,2,6] - a[2,3,2,6]
   - a[2,4,2,6] - a[2,5,2,6] - a[2,6,2,6] - a[2,7,2,6] - a[2,8,2,6]
   - a[2,9,2,6] - a[2,2,6,6] - a[2,3,6,6] - a[2,4,6,6] - a[2,5,6,6]
   - a[2,6,6,6] - a[2,7,6,6] - a[2,8,6,6] - a[2,9,6,6] + 2 b[6,2,2]
   - x[2,2] <= 0
 R2963: - a[3,3,0,6] - a[3,4,0,6] - a[3,5,0,6] - a[3,6,0,6] - a[3,7,0,6]
   - a[3,8,0,6] - a[3,9,0,6] - a[3,3,2,6] - a[3,4,2,6] - a[3,5,2,6]
   - a[3,6,2,6] - a[3,7,2,6] - a[3,8,2,6] - a[3,9,2,6] - a[3,3,6,6]
   - a[3,4,6,6] - a[3,5,6,6] - a[3,6,6,6] - a[3,7,6,6] - a[3,8,6,6]
   - a[3,9,6,6] + 2 b[6,3,0] - x[3,0] <= 0
 R2964: - a[3,3,0,6] - a[3,4,0,6] - a[3,5,0,6] - a[3,6,0,6] - a[3,7,0,6]
   - a[3,8,0,6] - a[3,9,0,6] - a[3,3,2,6] - a[3,4,2,6] - a[3,5,2,6]
   - a[3,6,2,6] - a[3,7,2,6] - a[3,8,2,6] - a[3,9,2,6] - a[3,3,6,6]
   - a[3,4,6,6] - a[3,5,6,6] - a[3,6,6,6] - a[3,7,6,6] - a[3,8,6,6]
   - a[3,9,6,6] + 2 b[6,3,1] - x[3,1] <= 0
 R2965: - a[3,3,0,6] - a[3,4,0,6] - a[3,5,0,6] - a[3,6,0,6] - a[3,7,0,6]
   - a[3,8,0,6] - a[3,9,0,6] - a[3,3,2,6] - a[3,4,2,6] - a[3,5,2,6]
   - a[3,6,2,6] - a[3,7,2,6] - a[3,8,2,6] - a[3,9,2,6] - a[3,3,6,6]
   - a[3,4,6,6] - a[3,5,6,6] - a[3,6,6,6] - a[3,7,6,6] - a[3,8,6,6]
   - a[3,9,6,6] + 2 b[6,3,2] - x[3,2] <= 0
 R2966: - a[4,4,0,6] - a[4,5,0,6] - a[4,6,0,6] - a[4,7,0,6] - a[4,8,0,6]
   - a[4,9,0,6] - a[4,4,2,6] - a[4,5,2,6] - a[4,6,2,6] - a[4,7,2,6]
   - a[4,8,2,6] - a[4,9,2,6] - a[4,4,6,6] - a[4,5,6,6] - a[4,6,6,6]
   - a[4,7,6,6] - a[4,8,6,6] - a[4,9,6,6] + 2 b[6,4,0] - x[4,0] <= 0
 R2967: - a[4,4,0,6] - a[4,5,0,6] - a[4,6,0,6] - a[4,7,0,6] - a[4,8,0,6]
   - a[4,9,0,6] - a[4,4,2,6] - a[4,5,2,6] - a[4,6,2,6] - a[4,7,2,6]
   - a[4,8,2,6] - a[4,9,2,6] - a[4,4,6,6] - a[4,5,6,6] - a[4,6,6,6]
   - a[4,7,6,6] - a[4,8,6,6] - a[4,9,6,6] + 2 b[6,4,1] - x[4,1] <= 0
 R2968: - a[4,4,0,6] - a[4,5,0,6] - a[4,6,0,6] - a[4,7,0,6] - a[4,8,0,6]
   - a[4,9,0,6] - a[4,4,2,6] - a[4,5,2,6] - a[4,6,2,6] - a[4,7,2,6]
   - a[4,8,2,6] - a[4,9,2,6] - a[4,4,6,6] - a[4,5,6,6] - a[4,6,6,6]
   - a[4,7,6,6] - a[4,8,6,6] - a[4,9,6,6] + 2 b[6,4,2] - x[4,2] <= 0
 R2969: - a[5,5,0,6] - a[5,6,0,6] - a[5,7,0,6] - a[5,8,0,6] - a[5,9,0,6]
   - a[5,5,2,6] - a[5,6,2,6] - a[5,7,2,6] - a[5,8,2,6] - a[5,9,2,6]
   - a[5,5,6,6] - a[5,6,6,6] - a[5,7,6,6] - a[5,8,6,6] - a[5,9,6,6]
   + 2 b[6,5,0] - x[5,0] <= 0
 R2970: - a[5,5,0,6] - a[5,6,0,6] - a[5,7,0,6] - a[5,8,0,6] - a[5,9,0,6]
   - a[5,5,2,6] - a[5,6,2,6] - a[5,7,2,6] - a[5,8,2,6] - a[5,9,2,6]
   - a[5,5,6,6] - a[5,6,6,6] - a[5,7,6,6] - a[5,8,6,6] - a[5,9,6,6]
   + 2 b[6,5,1] - x[5,1] <= 0
 R2971: - a[5,5,0,6] - a[5,6,0,6] - a[5,7,0,6] - a[5,8,0,6] - a[5,9,0,6]
   - a[5,5,2,6] - a[5,6,2,6] - a[5,7,2,6] - a[5,8,2,6] - a[5,9,2,6]
   - a[5,5,6,6] - a[5,6,6,6] - a[5,7,6,6] - a[5,8,6,6] - a[5,9,6,6]
   + 2 b[6,5,2] - x[5,2] <= 0
 R2972: - a[6,6,0,6] - a[6,7,0,6] - a[6,8,0,6] - a[6,9,0,6] - a[6,6,2,6]
   - a[6,7,2,6] - a[6,8,2,6] - a[6,9,2,6] - a[6,6,6,6] - a[6,7,6,6]
   - a[6,8,6,6] - a[6,9,6,6] + 2 b[6,6,0] - x[6,0] <= 0
 R2973: - a[6,6,0,6] - a[6,7,0,6] - a[6,8,0,6] - a[6,9,0,6] - a[6,6,2,6]
   - a[6,7,2,6] - a[6,8,2,6] - a[6,9,2,6] - a[6,6,6,6] - a[6,7,6,6]
   - a[6,8,6,6] - a[6,9,6,6] + 2 b[6,6,1] - x[6,1] <= 0
 R2974: - a[6,6,0,6] - a[6,7,0,6] - a[6,8,0,6] - a[6,9,0,6] - a[6,6,2,6]
   - a[6,7,2,6] - a[6,8,2,6] - a[6,9,2,6] - a[6,6,6,6] - a[6,7,6,6]
   - a[6,8,6,6] - a[6,9,6,6] + 2 b[6,6,2] - x[6,2] <= 0
 R2975: - a[7,7,0,6] - a[7,8,0,6] - a[7,9,0,6] - a[7,7,2,6] - a[7,8,2,6]
   - a[7,9,2,6] - a[7,7,6,6] - a[7,8,6,6] - a[7,9,6,6] + 2 b[6,7,0]
   - x[7,0] <= 0
 R2976: - a[7,7,0,6] - a[7,8,0,6] - a[7,9,0,6] - a[7,7,2,6] - a[7,8,2,6]
   - a[7,9,2,6] - a[7,7,6,6] - a[7,8,6,6] - a[7,9,6,6] + 2 b[6,7,1]
   - x[7,1] <= 0
 R2977: - a[7,7,0,6] - a[7,8,0,6] - a[7,9,0,6] - a[7,7,2,6] - a[7,8,2,6]
   - a[7,9,2,6] - a[7,7,6,6] - a[7,8,6,6] - a[7,9,6,6] + 2 b[6,7,2]
   - x[7,2] <= 0
 R2978: - a[8,8,0,6] - a[8,9,0,6] - a[8,8,2,6] - a[8,9,2,6] - a[8,8,6,6]
   - a[8,9,6,6] + 2 b[6,8,0] - x[8,0] <= 0
 R2979: - a[8,8,0,6] - a[8,9,0,6] - a[8,8,2,6] - a[8,9,2,6] - a[8,8,6,6]
   - a[8,9,6,6] + 2 b[6,8,1] - x[8,1] <= 0
 R2980: - a[8,8,0,6] - a[8,9,0,6] - a[8,8,2,6] - a[8,9,2,6] - a[8,8,6,6]
   - a[8,9,6,6] + 2 b[6,8,2] - x[8,2] <= 0
 R2981: - a[9,9,0,6] - a[9,9,2,6] - a[9,9,6,6] + 2 b[6,9,0] - x[9,0] <= 0
 R2982: - a[9,9,0,6] - a[9,9,2,6] - a[9,9,6,6] + 2 b[6,9,1] - x[9,1] <= 0
 R2983: - a[9,9,0,6] - a[9,9,2,6] - a[9,9,6,6] + 2 b[6,9,2] - x[9,2] <= 0
 R2984: - a[0,0,0,8] - a[0,1,0,8] - a[0,2,0,8] - a[0,3,0,8] - a[0,4,0,8]
   - a[0,5,0,8] - a[0,6,0,8] - a[0,7,0,8] - a[0,8,0,8] - a[0,9,0,8]
   - a[0,0,2,8] - a[0,1,2,8] - a[0,2,2,8] - a[0,3,2,8] - a[0,4,2,8]
   - a[0,5,2,8] - a[0,6,2,8] - a[0,7,2,8] - a[0,8,2,8] - a[0,9,2,8]
   - a[0,0,6,8] - a[0,1,6,8] - a[0,2,6,8] - a[0,3,6,8] - a[0,4,6,8]
   - a[0,5,6,8] - a[0,6,6,8] - a[0,7,6,8] - a[0,8,6,8] - a[0,9,6,8]
   - a[0,0,8,8] - a[0,1,8,8] - a[0,2,8,8] - a[0,3,8,8] - a[0,4,8,8]
   - a[0,5,8,8] - a[0,6,8,8] - a[0,7,8,8] - a[0,8,8,8] - a[0,9,8,8]
   + 2 b[8,0,0] - x[0,0] <= 0
 R2985: - a[0,0,0,8] - a[0,1,0,8] - a[0,2,0,8] - a[0,3,0,8] - a[0,4,0,8]
   - a[0,5,0,8] - a[0,6,0,8] - a[0,7,0,8] - a[0,8,0,8] - a[0,9,0,8]
   - a[0,0,2,8] - a[0,1,2,8] - a[0,2,2,8] - a[0,3,2,8] - a[0,4,2,8]
   - a[0,5,2,8] - a[0,6,2,8] - a[0,7,2,8] - a[0,8,2,8] - a[0,9,2,8]
   - a[0,0,6,8] - a[0,1,6,8] - a[0,2,6,8] - a[0,3,6,8] - a[0,4,6,8]
   - a[0,5,6,8] - a[0,6,6,8] - a[0,7,6,8] - a[0,8,6,8] - a[0,9,6,8]
   - a[0,0,8,8] - a[0,1,8,8] - a[0,2,8,8] - a[0,3,8,8] - a[0,4,8,8]
   - a[0,5,8,8] - a[0,6,8,8] - a[0,7,8,8] - a[0,8,8,8] - a[0,9,8,8]
   + 2 b[8,0,1] - x[0,1] <= 0
 R2986: - a[0,0,0,8] - a[0,1,0,8] - a[0,2,0,8] - a[0,3,0,8] - a[0,4,0,8]
   - a[0,5,0,8] - a[0,6,0,8] - a[0,7,0,8] - a[0,8,0,8] - a[0,9,0,8]
   - a[0,0,2,8] - a[0,1,2,8] - a[0,2,2,8] - a[0,3,2,8] - a[0,4,2,8]
   - a[0,5,2,8] - a[0,6,2,8] - a[0,7,2,8] - a[0,8,2,8] - a[0,9,2,8]
   - a[0,0,6,8] - a[0,1,6,8] - a[0,2,6,8] - a[0,3,6,8] - a[0,4,6,8]
   - a[0,5,6,8] - a[0,6,6,8] - a[0,7,6,8] - a[0,8,6,8] - a[0,9,6,8]
   - a[0,0,8,8] - a[0,1,8,8] - a[0,2,8,8] - a[0,3,8,8] - a[0,4,8,8]
   - a[0,5,8,8] - a[0,6,8,8] - a[0,7,8,8] - a[0,8,8,8] - a[0,9,8,8]
   + 2 b[8,0,2] - x[0,2] <= 0
 R2987: - a[1,1,0,8] - a[1,2,0,8] - a[1,3,0,8] - a[1,4,0,8] - a[1,5,0,8]
   - a[1,6,0,8] - a[1,7,0,8] - a[1,8,0,8] - a[1,9,0,8] - a[1,1,2,8]
   - a[1,2,2,8] - a[1,3,2,8] - a[1,4,2,8] - a[1,5,2,8] - a[1,6,2,8]
   - a[1,7,2,8] - a[1,8,2,8] - a[1,9,2,8] - a[1,1,6,8] - a[1,2,6,8]
   - a[1,3,6,8] - a[1,4,6,8] - a[1,5,6,8] - a[1,6,6,8] - a[1,7,6,8]
   - a[1,8,6,8] - a[1,9,6,8] - a[1,1,8,8] - a[1,2,8,8] - a[1,3,8,8]
   - a[1,4,8,8] - a[1,5,8,8] - a[1,6,8,8] - a[1,7,8,8] - a[1,8,8,8]
   - a[1,9,8,8] + 2 b[8,1,0] - x[1,0] <= 0
 R2988: - a[1,1,0,8] - a[1,2,0,8] - a[1,3,0,8] - a[1,4,0,8] - a[1,5,0,8]
   - a[1,6,0,8] - a[1,7,0,8] - a[1,8,0,8] - a[1,9,0,8] - a[1,1,2,8]
   - a[1,2,2,8] - a[1,3,2,8] - a[1,4,2,8] - a[1,5,2,8] - a[1,6,2,8]
   - a[1,7,2,8] - a[1,8,2,8] - a[1,9,2,8] - a[1,1,6,8] - a[1,2,6,8]
   - a[1,3,6,8] - a[1,4,6,8] - a[1,5,6,8] - a[1,6,6,8] - a[1,7,6,8]
   - a[1,8,6,8] - a[1,9,6,8] - a[1,1,8,8] - a[1,2,8,8] - a[1,3,8,8]
   - a[1,4,8,8] - a[1,5,8,8] - a[1,6,8,8] - a[1,7,8,8] - a[1,8,8,8]
   - a[1,9,8,8] + 2 b[8,1,1] - x[1,1] <= 0
 R2989: - a[1,1,0,8] - a[1,2,0,8] - a[1,3,0,8] - a[1,4,0,8] - a[1,5,0,8]
   - a[1,6,0,8] - a[1,7,0,8] - a[1,8,0,8] - a[1,9,0,8] - a[1,1,2,8]
   - a[1,2,2,8] - a[1,3,2,8] - a[1,4,2,8] - a[1,5,2,8] - a[1,6,2,8]
   - a[1,7,2,8] - a[1,8,2,8] - a[1,9,2,8] - a[1,1,6,8] - a[1,2,6,8]
   - a[1,3,6,8] - a[1,4,6,8] - a[1,5,6,8] - a[1,6,6,8] - a[1,7,6,8]
   - a[1,8,6,8] - a[1,9,6,8] - a[1,1,8,8] - a[1,2,8,8] - a[1,3,8,8]
   - a[1,4,8,8] - a[1,5,8,8] - a[1,6,8,8] - a[1,7,8,8] - a[1,8,8,8]
   - a[1,9,8,8] + 2 b[8,1,2] - x[1,2] <= 0
 R2990: - a[2,2,0,8] - a[2,3,0,8] - a[2,4,0,8] - a[2,5,0,8] - a[2,6,0,8]
   - a[2,7,0,8] - a[2,8,0,8] - a[2,9,0,8] - a[2,2,2,8] - a[2,3,2,8]
   - a[2,4,2,8] - a[2,5,2,8] - a[2,6,2,8] - a[2,7,2,8] - a[2,8,2,8]
   - a[2,9,2,8] - a[2,2,6,8] - a[2,3,6,8] - a[2,4,6,8] - a[2,5,6,8]
   - a[2,6,6,8] - a[2,7,6,8] - a[2,8,6,8] - a[2,9,6,8] - a[2,2,8,8]
   - a[2,3,8,8] - a[2,4,8,8] - a[2,5,8,8] - a[2,6,8,8] - a[2,7,8,8]
   - a[2,8,8,8] - a[2,9,8,8] + 2 b[8,2,0] - x[2,0] <= 0
 R2991: - a[2,2,0,8] - a[2,3,0,8] - a[2,4,0,8] - a[2,5,0,8] - a[2,6,0,8]
   - a[2,7,0,8] - a[2,8,0,8] - a[2,9,0,8] - a[2,2,2,8] - a[2,3,2,8]
   - a[2,4,2,8] - a[2,5,2,8] - a[2,6,2,8] - a[2,7,2,8] - a[2,8,2,8]
   - a[2,9,2,8] - a[2,2,6,8] - a[2,3,6,8] - a[2,4,6,8] - a[2,5,6,8]
   - a[2,6,6,8] - a[2,7,6,8] - a[2,8,6,8] - a[2,9,6,8] - a[2,2,8,8]
   - a[2,3,8,8] - a[2,4,8,8] - a[2,5,8,8] - a[2,6,8,8] - a[2,7,8,8]
   - a[2,8,8,8] - a[2,9,8,8] + 2 b[8,2,1] - x[2,1] <= 0
 R2992: - a[2,2,0,8] - a[2,3,0,8] - a[2,4,0,8] - a[2,5,0,8] - a[2,6,0,8]
   - a[2,7,0,8] - a[2,8,0,8] - a[2,9,0,8] - a[2,2,2,8] - a[2,3,2,8]
   - a[2,4,2,8] - a[2,5,2,8] - a[2,6,2,8] - a[2,7,2,8] - a[2,8,2,8]
   - a[2,9,2,8] - a[2,2,6,8] - a[2,3,6,8] - a[2,4,6,8] - a[2,5,6,8]
   - a[2,6,6,8] - a[2,7,6,8] - a[2,8,6,8] - a[2,9,6,8] - a[2,2,8,8]
   - a[2,3,8,8] - a[2,4,8,8] - a[2,5,8,8] - a[2,6,8,8] - a[2,7,8,8]
   - a[2,8,8,8] - a[2,9,8,8] + 2 b[8,2,2] - x[2,2] <= 0
 R2993: - a[3,3,0,8] - a[3,4,0,8] - a[3,5,0,8] - a[3,6,0,8] - a[3,7,0,8]
   - a[3,8,0,8] - a[3,9,0,8] - a[3,3,2,8] - a[3,4,2,8] - a[3,5,2,8]
   - a[3,6,2,8] - a[3,7,2,8] - a[3,8,2,8] - a[3,9,2,8] - a[3,3,6,8]
   - a[3,4,6,8] - a[3,5,6,8] - a[3,6,6,8] - a[3,7,6,8] - a[3,8,6,8]
   - a[3,9,6,8] - a[3,3,8,8] - a[3,4,8,8] - a[3,5,8,8] - a[3,6,8,8]
   - a[3,7,8,8] - a[3,8,8,8] - a[3,9,8,8] + 2 b[8,3,0] - x[3,0] <= 0
 R2994: - a[3,3,0,8] - a[3,4,0,8] - a[3,5,0,8] - a[3,6,0,8] - a[3,7,0,8]
   - a[3,8,0,8] - a[3,9,0,8] - a[3,3,2,8] - a[3,4,2,8] - a[3,5,2,8]
   - a[3,6,2,8] - a[3,7,2,8] - a[3,8,2,8] - a[3,9,2,8] - a[3,3,6,8]
   - a[3,4,6,8] - a[3,5,6,8] - a[3,6,6,8] - a[3,7,6,8] - a[3,8,6,8]
   - a[3,9,6,8] - a[3,3,8,8] - a[3,4,8,8] - a[3,5,8,8] - a[3,6,8,8]
   - a[3,7,8,8] - a[3,8,8,8] - a[3,9,8,8] + 2 b[8,3,1] - x[3,1] <= 0
 R2995: - a[3,3,0,8] - a[3,4,0,8] - a[3,5,0,8] - a[3,6,0,8] - a[3,7,0,8]
   - a[3,8,0,8] - a[3,9,0,8] - a[3,3,2,8] - a[3,4,2,8] - a[3,5,2,8]
   - a[3,6,2,8] - a[3,7,2,8] - a[3,8,2,8] - a[3,9,2,8] - a[3,3,6,8]
   - a[3,4,6,8] - a[3,5,6,8] - a[3,6,6,8] - a[3,7,6,8] - a[3,8,6,8]
   - a[3,9,6,8] - a[3,3,8,8] - a[3,4,8,8] - a[3,5,8,8] - a[3,6,8,8]
   - a[3,7,8,8] - a[3,8,8,8] - a[3,9,8,8] + 2 b[8,3,2] - x[3,2] <= 0
 R2996: - a[4,4,0,8] - a[4,5,0,8] - a[4,6,0,8] - a[4,7,0,8] - a[4,8,0,8]
   - a[4,9,0,8] - a[4,4,2,8] - a[4,5,2,8] - a[4,6,2,8] - a[4,7,2,8]
   - a[4,8,2,8] - a[4,9,2,8] - a[4,4,6,8] - a[4,5,6,8] - a[4,6,6,8]
   - a[4,7,6,8] - a[4,8,6,8] - a[4,9,6,8] - a[4,4,8,8] - a[4,5,8,8]
   - a[4,6,8,8] - a[4,7,8,8] - a[4,8,8,8] - a[4,9,8,8] + 2 b[8,4,0]
   - x[4,0] <= 0
 R2997: - a[4,4,0,8] - a[4,5,0,8] - a[4,6,0,8] - a[4,7,0,8] - a[4,8,0,8]
   - a[4,9,0,8] - a[4,4,2,8] - a[4,5,2,8] - a[4,6,2,8] - a[4,7,2,8]
   - a[4,8,2,8] - a[4,9,2,8] - a[4,4,6,8] - a[4,5,6,8] - a[4,6,6,8]
   - a[4,7,6,8] - a[4,8,6,8] - a[4,9,6,8] - a[4,4,8,8] - a[4,5,8,8]
   - a[4,6,8,8] - a[4,7,8,8] - a[4,8,8,8] - a[4,9,8,8] + 2 b[8,4,1]
   - x[4,1] <= 0
 R2998: - a[4,4,0,8] - a[4,5,0,8] - a[4,6,0,8] - a[4,7,0,8] - a[4,8,0,8]
   - a[4,9,0,8] - a[4,4,2,8] - a[4,5,2,8] - a[4,6,2,8] - a[4,7,2,8]
   - a[4,8,2,8] - a[4,9,2,8] - a[4,4,6,8] - a[4,5,6,8] - a[4,6,6,8]
   - a[4,7,6,8] - a[4,8,6,8] - a[4,9,6,8] - a[4,4,8,8] - a[4,5,8,8]
   - a[4,6,8,8] - a[4,7,8,8] - a[4,8,8,8] - a[4,9,8,8] + 2 b[8,4,2]
   - x[4,2] <= 0
 R2999: - a[5,5,0,8] - a[5,6,0,8] - a[5,7,0,8] - a[5,8,0,8] - a[5,9,0,8]
   - a[5,5,2,8] - a[5,6,2,8] - a[5,7,2,8] - a[5,8,2,8] - a[5,9,2,8]
   - a[5,5,6,8] - a[5,6,6,8] - a[5,7,6,8] - a[5,8,6,8] - a[5,9,6,8]
   - a[5,5,8,8] - a[5,6,8,8] - a[5,7,8,8] - a[5,8,8,8] - a[5,9,8,8]
   + 2 b[8,5,0] - x[5,0] <= 0
 R3000: - a[5,5,0,8] - a[5,6,0,8] - a[5,7,0,8] - a[5,8,0,8] - a[5,9,0,8]
   - a[5,5,2,8] - a[5,6,2,8] - a[5,7,2,8] - a[5,8,2,8] - a[5,9,2,8]
   - a[5,5,6,8] - a[5,6,6,8] - a[5,7,6,8] - a[5,8,6,8] - a[5,9,6,8]
   - a[5,5,8,8] - a[5,6,8,8] - a[5,7,8,8] - a[5,8,8,8] - a[5,9,8,8]
   + 2 b[8,5,1] - x[5,1] <= 0
 R3001: - a[5,5,0,8] - a[5,6,0,8] - a[5,7,0,8] - a[5,8,0,8] - a[5,9,0,8]
   - a[5,5,2,8] - a[5,6,2,8] - a[5,7,2,8] - a[5,8,2,8] - a[5,9,2,8]
   - a[5,5,6,8] - a[5,6,6,8] - a[5,7,6,8] - a[5,8,6,8] - a[5,9,6,8]
   - a[5,5,8,8] - a[5,6,8,8] - a[5,7,8,8] - a[5,8,8,8] - a[5,9,8,8]
   + 2 b[8,5,2] - x[5,2] <= 0
 R3002: - a[6,6,0,8] - a[6,7,0,8] - a[6,8,0,8] - a[6,9,0,8] - a[6,6,2,8]
   - a[6,7,2,8] - a[6,8,2,8] - a[6,9,2,8] - a[6,6,6,8] - a[6,7,6,8]
   - a[6,8,6,8] - a[6,9,6,8] - a[6,6,8,8] - a[6,7,8,8] - a[6,8,8,8]
   - a[6,9,8,8] + 2 b[8,6,0] - x[6,0] <= 0
 R3003: - a[6,6,0,8] - a[6,7,0,8] - a[6,8,0,8] - a[6,9,0,8] - a[6,6,2,8]
   - a[6,7,2,8] - a[6,8,2,8] - a[6,9,2,8] - a[6,6,6,8] - a[6,7,6,8]
   - a[6,8,6,8] - a[6,9,6,8] - a[6,6,8,8] - a[6,7,8,8] - a[6,8,8,8]
   - a[6,9,8,8] + 2 b[8,6,1] - x[6,1] <= 0
 R3004: - a[6,6,0,8] - a[6,7,0,8] - a[6,8,0,8] - a[6,9,0,8] - a[6,6,2,8]
   - a[6,7,2,8] - a[6,8,2,8] - a[6,9,2,8] - a[6,6,6,8] - a[6,7,6,8]
   - a[6,8,6,8] - a[6,9,6,8] - a[6,6,8,8] - a[6,7,8,8] - a[6,8,8,8]
   - a[6,9,8,8] + 2 b[8,6,2] - x[6,2] <= 0
 R3005: - a[7,7,0,8] - a[7,8,0,8] - a[7,9,0,8] - a[7,7,2,8] - a[7,8,2,8]
   - a[7,9,2,8] - a[7,7,6,8] - a[7,8,6,8] - a[7,9,6,8] - a[7,7,8,8]
   - a[7,8,8,8] - a[7,9,8,8] + 2 b[8,7,0] - x[7,0] <= 0
 R3006: - a[7,7,0,8] - a[7,8,0,8] - a[7,9,0,8] - a[7,7,2,8] - a[7,8,2,8]
   - a[7,9,2,8] - a[7,7,6,8] - a[7,8,6,8] - a[7,9,6,8] - a[7,7,8,8]
   - a[7,8,8,8] - a[7,9,8,8] + 2 b[8,7,1] - x[7,1] <= 0
 R3007: - a[7,7,0,8] - a[7,8,0,8] - a[7,9,0,8] - a[7,7,2,8] - a[7,8,2,8]
   - a[7,9,2,8] - a[7,7,6,8] - a[7,8,6,8] - a[7,9,6,8] - a[7,7,8,8]
   - a[7,8,8,8] - a[7,9,8,8] + 2 b[8,7,2] - x[7,2] <= 0
 R3008: - a[8,8,0,8] - a[8,9,0,8] - a[8,8,2,8] - a[8,9,2,8] - a[8,8,6,8]
   - a[8,9,6,8] - a[8,8,8,8] - a[8,9,8,8] + 2 b[8,8,0] - x[8,0] <= 0
 R3009: - a[8,8,0,8] - a[8,9,0,8] - a[8,8,2,8] - a[8,9,2,8] - a[8,8,6,8]
   - a[8,9,6,8] - a[8,8,8,8] - a[8,9,8,8] + 2 b[8,8,1] - x[8,1] <= 0
 R3010: - a[8,8,0,8] - a[8,9,0,8] - a[8,8,2,8] - a[8,9,2,8] - a[8,8,6,8]
   - a[8,9,6,8] - a[8,8,8,8] - a[8,9,8,8] + 2 b[8,8,2] - x[8,2] <= 0
 R3011: - a[9,9,0,8] - a[9,9,2,8] - a[9,9,6,8] - a[9,9,8,8] + 2 b[8,9,0]
   - x[9,0] <= 0
 R3012: - a[9,9,0,8] - a[9,9,2,8] - a[9,9,6,8] - a[9,9,8,8] + 2 b[8,9,1]
   - x[9,1] <= 0
 R3013: - a[9,9,0,8] - a[9,9,2,8] - a[9,9,6,8] - a[9,9,8,8] + 2 b[8,9,2]
   - x[9,2] <= 0
 R3014: - a[0,0,0,9] - a[0,1,0,9] - a[0,2,0,9] - a[0,3,0,9] - a[0,4,0,9]
   - a[0,5,0,9] - a[0,6,0,9] - a[0,7,0,9] - a[0,8,0,9] - a[0,9,0,9]
   - a[0,0,2,9] - a[0,1,2,9] - a[0,2,2,9] - a[0,3,2,9] - a[0,4,2,9]
   - a[0,5,2,9] - a[0,6,2,9] - a[0,7,2,9] - a[0,8,2,9] - a[0,9,2,9]
   - a[0,0,6,9] - a[0,1,6,9] - a[0,2,6,9] - a[0,3,6,9] - a[0,4,6,9]
   - a[0,5,6,9] - a[0,6,6,9] - a[0,7,6,9] - a[0,8,6,9] - a[0,9,6,9]
   - a[0,0,8,9] - a[0,1,8,9] - a[0,2,8,9] - a[0,3,8,9] - a[0,4,8,9]
   - a[0,5,8,9] - a[0,6,8,9] - a[0,7,8,9] - a[0,8,8,9] - a[0,9,8,9]
   - a[0,0,9,9] - a[0,1,9,9] - a[0,2,9,9] - a[0,3,9,9] - a[0,4,9,9]
   - a[0,5,9,9] - a[0,6,9,9] - a[0,7,9,9] - a[0,8,9,9] - a[0,9,9,9]
   + 2 b[9,0,0] - x[0,0] <= 0
 R3015: - a[0,0,0,9] - a[0,1,0,9] - a[0,2,0,9] - a[0,3,0,9] - a[0,4,0,9]
   - a[0,5,0,9] - a[0,6,0,9] - a[0,7,0,9] - a[0,8,0,9] - a[0,9,0,9]
   - a[0,0,2,9] - a[0,1,2,9] - a[0,2,2,9] - a[0,3,2,9] - a[0,4,2,9]
   - a[0,5,2,9] - a[0,6,2,9] - a[0,7,2,9] - a[0,8,2,9] - a[0,9,2,9]
   - a[0,0,6,9] - a[0,1,6,9] - a[0,2,6,9] - a[0,3,6,9] - a[0,4,6,9]
   - a[0,5,6,9] - a[0,6,6,9] - a[0,7,6,9] - a[0,8,6,9] - a[0,9,6,9]
   - a[0,0,8,9] - a[0,1,8,9] - a[0,2,8,9] - a[0,3,8,9] - a[0,4,8,9]
   - a[0,5,8,9] - a[0,6,8,9] - a[0,7,8,9] - a[0,8,8,9] - a[0,9,8,9]
   - a[0,0,9,9] - a[0,1,9,9] - a[0,2,9,9] - a[0,3,9,9] - a[0,4,9,9]
   - a[0,5,9,9] - a[0,6,9,9] - a[0,7,9,9] - a[0,8,9,9] - a[0,9,9,9]
   + 2 b[9,0,1] - x[0,1] <= 0
 R3016: - a[0,0,0,9] - a[0,1,0,9] - a[0,2,0,9] - a[0,3,0,9] - a[0,4,0,9]
   - a[0,5,0,9] - a[0,6,0,9] - a[0,7,0,9] - a[0,8,0,9] - a[0,9,0,9]
   - a[0,0,2,9] - a[0,1,2,9] - a[0,2,2,9] - a[0,3,2,9] - a[0,4,2,9]
   - a[0,5,2,9] - a[0,6,2,9] - a[0,7,2,9] - a[0,8,2,9] - a[0,9,2,9]
   - a[0,0,6,9] - a[0,1,6,9] - a[0,2,6,9] - a[0,3,6,9] - a[0,4,6,9]
   - a[0,5,6,9] - a[0,6,6,9] - a[0,7,6,9] - a[0,8,6,9] - a[0,9,6,9]
   - a[0,0,8,9] - a[0,1,8,9] - a[0,2,8,9] - a[0,3,8,9] - a[0,4,8,9]
   - a[0,5,8,9] - a[0,6,8,9] - a[0,7,8,9] - a[0,8,8,9] - a[0,9,8,9]
   - a[0,0,9,9] - a[0,1,9,9] - a[0,2,9,9] - a[0,3,9,9] - a[0,4,9,9]
   - a[0,5,9,9] - a[0,6,9,9] - a[0,7,9,9] - a[0,8,9,9] - a[0,9,9,9]
   + 2 b[9,0,2] - x[0,2] <= 0
 R3017: - a[1,1,0,9] - a[1,2,0,9] - a[1,3,0,9] - a[1,4,0,9] - a[1,5,0,9]
   - a[1,6,0,9] - a[1,7,0,9] - a[1,8,0,9] - a[1,9,0,9] - a[1,1,2,9]
   - a[1,2,2,9] - a[1,3,2,9] - a[1,4,2,9] - a[1,5,2,9] - a[1,6,2,9]
   - a[1,7,2,9] - a[1,8,2,9] - a[1,9,2,9] - a[1,1,6,9] - a[1,2,6,9]
   - a[1,3,6,9] - a[1,4,6,9] - a[1,5,6,9] - a[1,6,6,9] - a[1,7,6,9]
   - a[1,8,6,9] - a[1,9,6,9] - a[1,1,8,9] - a[1,2,8,9] - a[1,3,8,9]
   - a[1,4,8,9] - a[1,5,8,9] - a[1,6,8,9] - a[1,7,8,9] - a[1,8,8,9]
   - a[1,9,8,9] - a[1,1,9,9] - a[1,2,9,9] - a[1,3,9,9] - a[1,4,9,9]
   - a[1,5,9,9] - a[1,6,9,9] - a[1,7,9,9] - a[1,8,9,9] - a[1,9,9,9]
   + 2 b[9,1,0] - x[1,0] <= 0
 R3018: - a[1,1,0,9] - a[1,2,0,9] - a[1,3,0,9] - a[1,4,0,9] - a[1,5,0,9]
   - a[1,6,0,9] - a[1,7,0,9] - a[1,8,0,9] - a[1,9,0,9] - a[1,1,2,9]
   - a[1,2,2,9] - a[1,3,2,9] - a[1,4,2,9] - a[1,5,2,9] - a[1,6,2,9]
   - a[1,7,2,9] - a[1,8,2,9] - a[1,9,2,9] - a[1,1,6,9] - a[1,2,6,9]
   - a[1,3,6,9] - a[1,4,6,9] - a[1,5,6,9] - a[1,6,6,9] - a[1,7,6,9]
   - a[1,8,6,9] - a[1,9,6,9] - a[1,1,8,9] - a[1,2,8,9] - a[1,3,8,9]
   - a[1,4,8,9] - a[1,5,8,9] - a[1,6,8,9] - a[1,7,8,9] - a[1,8,8,9]
   - a[1,9,8,9] - a[1,1,9,9] - a[1,2,9,9] - a[1,3,9,9] - a[1,4,9,9]
   - a[1,5,9,9] - a[1,6,9,9] - a[1,7,9,9] - a[1,8,9,9] - a[1,9,9,9]
   + 2 b[9,1,1] - x[1,1] <= 0
 R3019: - a[1,1,0,9] - a[1,2,0,9] - a[1,3,0,9] - a[1,4,0,9] - a[1,5,0,9]
   - a[1,6,0,9] - a[1,7,0,9] - a[1,8,0,9] - a[1,9,0,9] - a[1,1,2,9]
   - a[1,2,2,9] - a[1,3,2,9] - a[1,4,2,9] - a[1,5,2,9] - a[1,6,2,9]
   - a[1,7,2,9] - a[1,8,2,9] - a[1,9,2,9] - a[1,1,6,9] - a[1,2,6,9]
   - a[1,3,6,9] - a[1,4,6,9] - a[1,5,6,9] - a[1,6,6,9] - a[1,7,6,9]
   - a[1,8,6,9] - a[1,9,6,9] - a[1,1,8,9] - a[1,2,8,9] - a[1,3,8,9]
   - a[1,4,8,9] - a[1,5,8,9] - a[1,6,8,9] - a[1,7,8,9] - a[1,8,8,9]
   - a[1,9,8,9] - a[1,1,9,9] - a[1,2,9,9] - a[1,3,9,9] - a[1,4,9,9]
   - a[1,5,9,9] - a[1,6,9,9] - a[1,7,9,9] - a[1,8,9,9] - a[1,9,9,9]
   + 2 b[9,1,2] - x[1,2] <= 0
 R3020: - a[2,2,0,9] - a[2,3,0,9] - a[2,4,0,9] - a[2,5,0,9] - a[2,6,0,9]
   - a[2,7,0,9] - a[2,8,0,9] - a[2,9,0,9] - a[2,2,2,9] - a[2,3,2,9]
   - a[2,4,2,9] - a[2,5,2,9] - a[2,6,2,9] - a[2,7,2,9] - a[2,8,2,9]
   - a[2,9,2,9] - a[2,2,6,9] - a[2,3,6,9] - a[2,4,6,9] - a[2,5,6,9]
   - a[2,6,6,9] - a[2,7,6,9] - a[2,8,6,9] - a[2,9,6,9] - a[2,2,8,9]
   - a[2,3,8,9] - a[2,4,8,9] - a[2,5,8,9] - a[2,6,8,9] - a[2,7,8,9]
   - a[2,8,8,9] - a[2,9,8,9] - a[2,2,9,9] - a[2,3,9,9] - a[2,4,9,9]
   - a[2,5,9,9] - a[2,6,9,9] - a[2,7,9,9] - a[2,8,9,9] - a[2,9,9,9]
   + 2 b[9,2,0] - x[2,0] <= 0
 R3021: - a[2,2,0,9] - a[2,3,0,9] - a[2,4,0,9] - a[2,5,0,9] - a[2,6,0,9]
   - a[2,7,0,9] - a[2,8,0,9] - a[2,9,0,9] - a[2,2,2,9] - a[2,3,2,9]
   - a[2,4,2,9] - a[2,5,2,9] - a[2,6,2,9] - a[2,7,2,9] - a[2,8,2,9]
   - a[2,9,2,9] - a[2,2,6,9] - a[2,3,6,9] - a[2,4,6,9] - a[2,5,6,9]
   - a[2,6,6,9] - a[2,7,6,9] - a[2,8,6,9] - a[2,9,6,9] - a[2,2,8,9]
   - a[2,3,8,9] - a[2,4,8,9] - a[2,5,8,9] - a[2,6,8,9] - a[2,7,8,9]
   - a[2,8,8,9] - a[2,9,8,9] - a[2,2,9,9] - a[2,3,9,9] - a[2,4,9,9]
   - a[2,5,9,9] - a[2,6,9,9] - a[2,7,9,9] - a[2,8,9,9] - a[2,9,9,9]
   + 2 b[9,2,1] - x[2,1] <= 0
 R3022: - a[2,2,0,9] - a[2,3,0,9] - a[2,4,0,9] - a[2,5,0,9] - a[2,6,0,9]
   - a[2,7,0,9] - a[2,8,0,9] - a[2,9,0,9] - a[2,2,2,9] - a[2,3,2,9]
   - a[2,4,2,9] - a[2,5,2,9] - a[2,6,2,9] - a[2,7,2,9] - a[2,8,2,9]
   - a[2,9,2,9] - a[2,2,6,9] - a[2,3,6,9] - a[2,4,6,9] - a[2,5,6,9]
   - a[2,6,6,9] - a[2,7,6,9] - a[2,8,6,9] - a[2,9,6,9] - a[2,2,8,9]
   - a[2,3,8,9] - a[2,4,8,9] - a[2,5,8,9] - a[2,6,8,9] - a[2,7,8,9]
   - a[2,8,8,9] - a[2,9,8,9] - a[2,2,9,9] - a[2,3,9,9] - a[2,4,9,9]
   - a[2,5,9,9] - a[2,6,9,9] - a[2,7,9,9] - a[2,8,9,9] - a[2,9,9,9]
   + 2 b[9,2,2] - x[2,2] <= 0
 R3023: - a[3,3,0,9] - a[3,4,0,9] - a[3,5,0,9] - a[3,6,0,9] - a[3,7,0,9]
   - a[3,8,0,9] - a[3,9,0,9] - a[3,3,2,9] - a[3,4,2,9] - a[3,5,2,9]
   - a[3,6,2,9] - a[3,7,2,9] - a[3,8,2,9] - a[3,9,2,9] - a[3,3,6,9]
   - a[3,4,6,9] - a[3,5,6,9] - a[3,6,6,9] - a[3,7,6,9] - a[3,8,6,9]
   - a[3,9,6,9] - a[3,3,8,9] - a[3,4,8,9] - a[3,5,8,9] - a[3,6,8,9]
   - a[3,7,8,9] - a[3,8,8,9] - a[3,9,8,9] - a[3,3,9,9] - a[3,4,9,9]
   - a[3,5,9,9] - a[3,6,9,9] - a[3,7,9,9] - a[3,8,9,9] - a[3,9,9,9]
   + 2 b[9,3,0] - x[3,0] <= 0
 R3024: - a[3,3,0,9] - a[3,4,0,9] - a[3,5,0,9] - a[3,6,0,9] - a[3,7,0,9]
   - a[3,8,0,9] - a[3,9,0,9] - a[3,3,2,9] - a[3,4,2,9] - a[3,5,2,9]
   - a[3,6,2,9] - a[3,7,2,9] - a[3,8,2,9] - a[3,9,2,9] - a[3,3,6,9]
   - a[3,4,6,9] - a[3,5,6,9] - a[3,6,6,9] - a[3,7,6,9] - a[3,8,6,9]
   - a[3,9,6,9] - a[3,3,8,9] - a[3,4,8,9] - a[3,5,8,9] - a[3,6,8,9]
   - a[3,7,8,9] - a[3,8,8,9] - a[3,9,8,9] - a[3,3,9,9] - a[3,4,9,9]
   - a[3,5,9,9] - a[3,6,9,9] - a[3,7,9,9] - a[3,8,9,9] - a[3,9,9,9]
   + 2 b[9,3,1] - x[3,1] <= 0
 R3025: - a[3,3,0,9] - a[3,4,0,9] - a[3,5,0,9] - a[3,6,0,9] - a[3,7,0,9]
   - a[3,8,0,9] - a[3,9,0,9] - a[3,3,2,9] - a[3,4,2,9] - a[3,5,2,9]
   - a[3,6,2,9] - a[3,7,2,9] - a[3,8,2,9] - a[3,9,2,9] - a[3,3,6,9]
   - a[3,4,6,9] - a[3,5,6,9] - a[3,6,6,9] - a[3,7,6,9] - a[3,8,6,9]
   - a[3,9,6,9] - a[3,3,8,9] - a[3,4,8,9] - a[3,5,8,9] - a[3,6,8,9]
   - a[3,7,8,9] - a[3,8,8,9] - a[3,9,8,9] - a[3,3,9,9] - a[3,4,9,9]
   - a[3,5,9,9] - a[3,6,9,9] - a[3,7,9,9] - a[3,8,9,9] - a[3,9,9,9]
   + 2 b[9,3,2] - x[3,2] <= 0
 R3026: - a[4,4,0,9] - a[4,5,0,9] - a[4,6,0,9] - a[4,7,0,9] - a[4,8,0,9]
   - a[4,9,0,9] - a[4,4,2,9] - a[4,5,2,9] - a[4,6,2,9] - a[4,7,2,9]
   - a[4,8,2,9] - a[4,9,2,9] - a[4,4,6,9] - a[4,5,6,9] - a[4,6,6,9]
   - a[4,7,6,9] - a[4,8,6,9] - a[4,9,6,9] - a[4,4,8,9] - a[4,5,8,9]
   - a[4,6,8,9] - a[4,7,8,9] - a[4,8,8,9] - a[4,9,8,9] - a[4,4,9,9]
   - a[4,5,9,9] - a[4,6,9,9] - a[4,7,9,9] - a[4,8,9,9] - a[4,9,9,9]
   + 2 b[9,4,0] - x[4,0] <= 0
 R3027: - a[4,4,0,9] - a[4,5,0,9] - a[4,6,0,9] - a[4,7,0,9] - a[4,8,0,9]
   - a[4,9,0,9] - a[4,4,2,9] - a[4,5,2,9] - a[4,6,2,9] - a[4,7,2,9]
   - a[4,8,2,9] - a[4,9,2,9] - a[4,4,6,9] - a[4,5,6,9] - a[4,6,6,9]
   - a[4,7,6,9] - a[4,8,6,9] - a[4,9,6,9] - a[4,4,8,9] - a[4,5,8,9]
   - a[4,6,8,9] - a[4,7,8,9] - a[4,8,8,9] - a[4,9,8,9] - a[4,4,9,9]
   - a[4,5,9,9] - a[4,6,9,9] - a[4,7,9,9] - a[4,8,9,9] - a[4,9,9,9]
   + 2 b[9,4,1] - x[4,1] <= 0
 R3028: - a[4,4,0,9] - a[4,5,0,9] - a[4,6,0,9] - a[4,7,0,9] - a[4,8,0,9]
   - a[4,9,0,9] - a[4,4,2,9] - a[4,5,2,9] - a[4,6,2,9] - a[4,7,2,9]
   - a[4,8,2,9] - a[4,9,2,9] - a[4,4,6,9] - a[4,5,6,9] - a[4,6,6,9]
   - a[4,7,6,9] - a[4,8,6,9] - a[4,9,6,9] - a[4,4,8,9] - a[4,5,8,9]
   - a[4,6,8,9] - a[4,7,8,9] - a[4,8,8,9] - a[4,9,8,9] - a[4,4,9,9]
   - a[4,5,9,9] - a[4,6,9,9] - a[4,7,9,9] - a[4,8,9,9] - a[4,9,9,9]
   + 2 b[9,4,2] - x[4,2] <= 0
 R3029: - a[5,5,0,9] - a[5,6,0,9] - a[5,7,0,9] - a[5,8,0,9] - a[5,9,0,9]
   - a[5,5,2,9] - a[5,6,2,9] - a[5,7,2,9] - a[5,8,2,9] - a[5,9,2,9]
   - a[5,5,6,9] - a[5,6,6,9] - a[5,7,6,9] - a[5,8,6,9] - a[5,9,6,9]
   - a[5,5,8,9] - a[5,6,8,9] - a[5,7,8,9] - a[5,8,8,9] - a[5,9,8,9]
   - a[5,5,9,9] - a[5,6,9,9] - a[5,7,9,9] - a[5,8,9,9] - a[5,9,9,9]
   + 2 b[9,5,0] - x[5,0] <= 0
 R3030: - a[5,5,0,9] - a[5,6,0,9] - a[5,7,0,9] - a[5,8,0,9] - a[5,9,0,9]
   - a[5,5,2,9] - a[5,6,2,9] - a[5,7,2,9] - a[5,8,2,9] - a[5,9,2,9]
   - a[5,5,6,9] - a[5,6,6,9] - a[5,7,6,9] - a[5,8,6,9] - a[5,9,6,9]
   - a[5,5,8,9] - a[5,6,8,9] - a[5,7,8,9] - a[5,8,8,9] - a[5,9,8,9]
   - a[5,5,9,9] - a[5,6,9,9] - a[5,7,9,9] - a[5,8,9,9] - a[5,9,9,9]
   + 2 b[9,5,1] - x[5,1] <= 0
 R3031: - a[5,5,0,9] - a[5,6,0,9] - a[5,7,0,9] - a[5,8,0,9] - a[5,9,0,9]
   - a[5,5,2,9] - a[5,6,2,9] - a[5,7,2,9] - a[5,8,2,9] - a[5,9,2,9]
   - a[5,5,6,9] - a[5,6,6,9] - a[5,7,6,9] - a[5,8,6,9] - a[5,9,6,9]
   - a[5,5,8,9] - a[5,6,8,9] - a[5,7,8,9] - a[5,8,8,9] - a[5,9,8,9]
   - a[5,5,9,9] - a[5,6,9,9] - a[5,7,9,9] - a[5,8,9,9] - a[5,9,9,9]
   + 2 b[9,5,2] - x[5,2] <= 0
 R3032: - a[6,6,0,9] - a[6,7,0,9] - a[6,8,0,9] - a[6,9,0,9] - a[6,6,2,9]
   - a[6,7,2,9] - a[6,8,2,9] - a[6,9,2,9] - a[6,6,6,9] - a[6,7,6,9]
   - a[6,8,6,9] - a[6,9,6,9] - a[6,6,8,9] - a[6,7,8,9] - a[6,8,8,9]
   - a[6,9,8,9] - a[6,6,9,9] - a[6,7,9,9] - a[6,8,9,9] - a[6,9,9,9]
   + 2 b[9,6,0] - x[6,0] <= 0
 R3033: - a[6,6,0,9] - a[6,7,0,9] - a[6,8,0,9] - a[6,9,0,9] - a[6,6,2,9]
   - a[6,7,2,9] - a[6,8,2,9] - a[6,9,2,9] - a[6,6,6,9] - a[6,7,6,9]
   - a[6,8,6,9] - a[6,9,6,9] - a[6,6,8,9] - a[6,7,8,9] - a[6,8,8,9]
   - a[6,9,8,9] - a[6,6,9,9] - a[6,7,9,9] - a[6,8,9,9] - a[6,9,9,9]
   + 2 b[9,6,1] - x[6,1] <= 0
 R3034: - a[6,6,0,9] - a[6,7,0,9] - a[6,8,0,9] - a[6,9,0,9] - a[6,6,2,9]
   - a[6,7,2,9] - a[6,8,2,9] - a[6,9,2,9] - a[6,6,6,9] - a[6,7,6,9]
   - a[6,8,6,9] - a[6,9,6,9] - a[6,6,8,9] - a[6,7,8,9] - a[6,8,8,9]
   - a[6,9,8,9] - a[6,6,9,9] - a[6,7,9,9] - a[6,8,9,9] - a[6,9,9,9]
   + 2 b[9,6,2] - x[6,2] <= 0
 R3035: - a[7,7,0,9] - a[7,8,0,9] - a[7,9,0,9] - a[7,7,2,9] - a[7,8,2,9]
   - a[7,9,2,9] - a[7,7,6,9] - a[7,8,6,9] - a[7,9,6,9] - a[7,7,8,9]
   - a[7,8,8,9] - a[7,9,8,9] - a[7,7,9,9] - a[7,8,9,9] - a[7,9,9,9]
   + 2 b[9,7,0] - x[7,0] <= 0
 R3036: - a[7,7,0,9] - a[7,8,0,9] - a[7,9,0,9] - a[7,7,2,9] - a[7,8,2,9]
   - a[7,9,2,9] - a[7,7,6,9] - a[7,8,6,9] - a[7,9,6,9] - a[7,7,8,9]
   - a[7,8,8,9] - a[7,9,8,9] - a[7,7,9,9] - a[7,8,9,9] - a[7,9,9,9]
   + 2 b[9,7,1] - x[7,1] <= 0
 R3037: - a[7,7,0,9] - a[7,8,0,9] - a[7,9,0,9] - a[7,7,2,9] - a[7,8,2,9]
   - a[7,9,2,9] - a[7,7,6,9] - a[7,8,6,9] - a[7,9,6,9] - a[7,7,8,9]
   - a[7,8,8,9] - a[7,9,8,9] - a[7,7,9,9] - a[7,8,9,9] - a[7,9,9,9]
   + 2 b[9,7,2] - x[7,2] <= 0
 R3038: - a[8,8,0,9] - a[8,9,0,9] - a[8,8,2,9] - a[8,9,2,9] - a[8,8,6,9]
   - a[8,9,6,9] - a[8,8,8,9] - a[8,9,8,9] - a[8,8,9,9] - a[8,9,9,9]
   + 2 b[9,8,0] - x[8,0] <= 0
 R3039: - a[8,8,0,9] - a[8,9,0,9] - a[8,8,2,9] - a[8,9,2,9] - a[8,8,6,9]
   - a[8,9,6,9] - a[8,8,8,9] - a[8,9,8,9] - a[8,8,9,9] - a[8,9,9,9]
   + 2 b[9,8,1] - x[8,1] <= 0
 R3040: - a[8,8,0,9] - a[8,9,0,9] - a[8,8,2,9] - a[8,9,2,9] - a[8,8,6,9]
   - a[8,9,6,9] - a[8,8,8,9] - a[8,9,8,9] - a[8,8,9,9] - a[8,9,9,9]
   + 2 b[9,8,2] - x[8,2] <= 0
 R3041: - a[9,9,0,9] - a[9,9,2,9] - a[9,9,6,9] - a[9,9,8,9] - a[9,9,9,9]
   + 2 b[9,9,0] - x[9,0] <= 0
 R3042: - a[9,9,0,9] - a[9,9,2,9] - a[9,9,6,9] - a[9,9,8,9] - a[9,9,9,9]
   + 2 b[9,9,1] - x[9,1] <= 0
 R3043: - a[9,9,0,9] - a[9,9,2,9] - a[9,9,6,9] - a[9,9,8,9] - a[9,9,9,9]
   + 2 b[9,9,2] - x[9,2] <= 0
 R3044: - a[0,0,1,1] - a[0,1,1,1] - a[0,2,1,1] - a[0,3,1,1] - a[0,4,1,1]
   - a[0,5,1,1] - a[0,6,1,1] - a[0,7,1,1] - a[0,8,1,1] - a[0,9,1,1]
   + 2 b[1,0,0] - x[0,0] <= 0
 R3045: - a[0,0,1,1] - a[0,1,1,1] - a[0,2,1,1] - a[0,3,1,1] - a[0,4,1,1]
   - a[0,5,1,1] - a[0,6,1,1] - a[0,7,1,1] - a[0,8,1,1] - a[0,9,1,1]
   + 2 b[1,0,1] - x[0,1] <= 0
 R3046: - a[0,0,1,1] - a[0,1,1,1] - a[0,2,1,1] - a[0,3,1,1] - a[0,4,1,1]
   - a[0,5,1,1] - a[0,6,1,1] - a[0,7,1,1] - a[0,8,1,1] - a[0,9,1,1]
   + 2 b[1,0,2] - x[0,2] <= 0
 R3047: - a[1,1,1,1] - a[1,2,1,1] - a[1,3,1,1] - a[1,4,1,1] - a[1,5,1,1]
   - a[1,6,1,1] - a[1,7,1,1] - a[1,8,1,1] - a[1,9,1,1] + 2 b[1,1,0]
   - x[1,0] <= 0
 R3048: - a[1,1,1,1] - a[1,2,1,1] - a[1,3,1,1] - a[1,4,1,1] - a[1,5,1,1]
   - a[1,6,1,1] - a[1,7,1,1] - a[1,8,1,1] - a[1,9,1,1] + 2 b[1,1,1]
   - x[1,1] <= 0
 R3049: - a[1,1,1,1] - a[1,2,1,1] - a[1,3,1,1] - a[1,4,1,1] - a[1,5,1,1]
   - a[1,6,1,1] - a[1,7,1,1] - a[1,8,1,1] - a[1,9,1,1] + 2 b[1,1,2]
   - x[1,2] <= 0
 R3050: - a[2,2,1,1] - a[2,3,1,1] - a[2,4,1,1] - a[2,5,1,1] - a[2,6,1,1]
   - a[2,7,1,1] - a[2,8,1,1] - a[2,9,1,1] + 2 b[1,2,0] - x[2,0] <= 0
 R3051: - a[2,2,1,1] - a[2,3,1,1] - a[2,4,1,1] - a[2,5,1,1] - a[2,6,1,1]
   - a[2,7,1,1] - a[2,8,1,1] - a[2,9,1,1] + 2 b[1,2,1] - x[2,1] <= 0
 R3052: - a[2,2,1,1] - a[2,3,1,1] - a[2,4,1,1] - a[2,5,1,1] - a[2,6,1,1]
   - a[2,7,1,1] - a[2,8,1,1] - a[2,9,1,1] + 2 b[1,2,2] - x[2,2] <= 0
 R3053: - a[3,3,1,1] - a[3,4,1,1] - a[3,5,1,1] - a[3,6,1,1] - a[3,7,1,1]
   - a[3,8,1,1] - a[3,9,1,1] + 2 b[1,3,0] - x[3,0] <= 0
 R3054: - a[3,3,1,1] - a[3,4,1,1] - a[3,5,1,1] - a[3,6,1,1] - a[3,7,1,1]
   - a[3,8,1,1] - a[3,9,1,1] + 2 b[1,3,1] - x[3,1] <= 0
 R3055: - a[3,3,1,1] - a[3,4,1,1] - a[3,5,1,1] - a[3,6,1,1] - a[3,7,1,1]
   - a[3,8,1,1] - a[3,9,1,1] + 2 b[1,3,2] - x[3,2] <= 0
 R3056: - a[4,4,1,1] - a[4,5,1,1] - a[4,6,1,1] - a[4,7,1,1] - a[4,8,1,1]
   - a[4,9,1,1] + 2 b[1,4,0] - x[4,0] <= 0
 R3057: - a[4,4,1,1] - a[4,5,1,1] - a[4,6,1,1] - a[4,7,1,1] - a[4,8,1,1]
   - a[4,9,1,1] + 2 b[1,4,1] - x[4,1] <= 0
 R3058: - a[4,4,1,1] - a[4,5,1,1] - a[4,6,1,1] - a[4,7,1,1] - a[4,8,1,1]
   - a[4,9,1,1] + 2 b[1,4,2] - x[4,2] <= 0
 R3059: - a[5,5,1,1] - a[5,6,1,1] - a[5,7,1,1] - a[5,8,1,1] - a[5,9,1,1]
   + 2 b[1,5,0] - x[5,0] <= 0
 R3060: - a[5,5,1,1] - a[5,6,1,1] - a[5,7,1,1] - a[5,8,1,1] - a[5,9,1,1]
   + 2 b[1,5,1] - x[5,1] <= 0
 R3061: - a[5,5,1,1] - a[5,6,1,1] - a[5,7,1,1] - a[5,8,1,1] - a[5,9,1,1]
   + 2 b[1,5,2] - x[5,2] <= 0
 R3062: - a[6,6,1,1] - a[6,7,1,1] - a[6,8,1,1] - a[6,9,1,1] + 2 b[1,6,0]
   - x[6,0] <= 0
 R3063: - a[6,6,1,1] - a[6,7,1,1] - a[6,8,1,1] - a[6,9,1,1] + 2 b[1,6,1]
   - x[6,1] <= 0
 R3064: - a[6,6,1,1] - a[6,7,1,1] - a[6,8,1,1] - a[6,9,1,1] + 2 b[1,6,2]
   - x[6,2] <= 0
 R3065: - a[7,7,1,1] - a[7,8,1,1] - a[7,9,1,1] + 2 b[1,7,0] - x[7,0] <= 0
 R3066: - a[7,7,1,1] - a[7,8,1,1] - a[7,9,1,1] + 2 b[1,7,1] - x[7,1] <= 0
 R3067: - a[7,7,1,1] - a[7,8,1,1] - a[7,9,1,1] + 2 b[1,7,2] - x[7,2] <= 0
 R3068: - a[8,8,1,1] - a[8,9,1,1] + 2 b[1,8,0] - x[8,0] <= 0
 R3069: - a[8,8,1,1] - a[8,9,1,1] + 2 b[1,8,1] - x[8,1] <= 0
 R3070: - a[8,8,1,1] - a[8,9,1,1] + 2 b[1,8,2] - x[8,2] <= 0
 R3071: - a[9,9,1,1] + 2 b[1,9,0] - x[9,0] <= 0
 R3072: - a[9,9,1,1] + 2 b[1,9,1] - x[9,1] <= 0
 R3073: - a[9,9,1,1] + 2 b[1,9,2] - x[9,2] <= 0
 R3074: - a[0,0,1,3] - a[0,1,1,3] - a[0,2,1,3] - a[0,3,1,3] - a[0,4,1,3]
   - a[0,5,1,3] - a[0,6,1,3] - a[0,7,1,3] - a[0,8,1,3] - a[0,9,1,3]
   - a[0,0,3,3] - a[0,1,3,3] - a[0,2,3,3] - a[0,3,3,3] - a[0,4,3,3]
   - a[0,5,3,3] - a[0,6,3,3] - a[0,7,3,3] - a[0,8,3,3] - a[0,9,3,3]
   + 2 b[3,0,0] - x[0,0] <= 0
 R3075: - a[0,0,1,3] - a[0,1,1,3] - a[0,2,1,3] - a[0,3,1,3] - a[0,4,1,3]
   - a[0,5,1,3] - a[0,6,1,3] - a[0,7,1,3] - a[0,8,1,3] - a[0,9,1,3]
   - a[0,0,3,3] - a[0,1,3,3] - a[0,2,3,3] - a[0,3,3,3] - a[0,4,3,3]
   - a[0,5,3,3] - a[0,6,3,3] - a[0,7,3,3] - a[0,8,3,3] - a[0,9,3,3]
   + 2 b[3,0,1] - x[0,1] <= 0
 R3076: - a[0,0,1,3] - a[0,1,1,3] - a[0,2,1,3] - a[0,3,1,3] - a[0,4,1,3]
   - a[0,5,1,3] - a[0,6,1,3] - a[0,7,1,3] - a[0,8,1,3] - a[0,9,1,3]
   - a[0,0,3,3] - a[0,1,3,3] - a[0,2,3,3] - a[0,3,3,3] - a[0,4,3,3]
   - a[0,5,3,3] - a[0,6,3,3] - a[0,7,3,3] - a[0,8,3,3] - a[0,9,3,3]
   + 2 b[3,0,2] - x[0,2] <= 0
 R3077: - a[1,1,1,3] - a[1,2,1,3] - a[1,3,1,3] - a[1,4,1,3] - a[1,5,1,3]
   - a[1,6,1,3] - a[1,7,1,3] - a[1,8,1,3] - a[1,9,1,3] - a[1,1,3,3]
   - a[1,2,3,3] - a[1,3,3,3] - a[1,4,3,3] - a[1,5,3,3] - a[1,6,3,3]
   - a[1,7,3,3] - a[1,8,3,3] - a[1,9,3,3] + 2 b[3,1,0] - x[1,0] <= 0
 R3078: - a[1,1,1,3] - a[1,2,1,3] - a[1,3,1,3] - a[1,4,1,3] - a[1,5,1,3]
   - a[1,6,1,3] - a[1,7,1,3] - a[1,8,1,3] - a[1,9,1,3] - a[1,1,3,3]
   - a[1,2,3,3] - a[1,3,3,3] - a[1,4,3,3] - a[1,5,3,3] - a[1,6,3,3]
   - a[1,7,3,3] - a[1,8,3,3] - a[1,9,3,3] + 2 b[3,1,1] - x[1,1] <= 0
 R3079: - a[1,1,1,3] - a[1,2,1,3] - a[1,3,1,3] - a[1,4,1,3] - a[1,5,1,3]
   - a[1,6,1,3] - a[1,7,1,3] - a[1,8,1,3] - a[1,9,1,3] - a[1,1,3,3]
   - a[1,2,3,3] - a[1,3,3,3] - a[1,4,3,3] - a[1,5,3,3] - a[1,6,3,3]
   - a[1,7,3,3] - a[1,8,3,3] - a[1,9,3,3] + 2 b[3,1,2] - x[1,2] <= 0
 R3080: - a[2,2,1,3] - a[2,3,1,3] - a[2,4,1,3] - a[2,5,1,3] - a[2,6,1,3]
   - a[2,7,1,3] - a[2,8,1,3] - a[2,9,1,3] - a[2,2,3,3] - a[2,3,3,3]
   - a[2,4,3,3] - a[2,5,3,3] - a[2,6,3,3] - a[2,7,3,3] - a[2,8,3,3]
   - a[2,9,3,3] + 2 b[3,2,0] - x[2,0] <= 0
 R3081: - a[2,2,1,3] - a[2,3,1,3] - a[2,4,1,3] - a[2,5,1,3] - a[2,6,1,3]
   - a[2,7,1,3] - a[2,8,1,3] - a[2,9,1,3] - a[2,2,3,3] - a[2,3,3,3]
   - a[2,4,3,3] - a[2,5,3,3] - a[2,6,3,3] - a[2,7,3,3] - a[2,8,3,3]
   - a[2,9,3,3] + 2 b[3,2,1] - x[2,1] <= 0
 R3082: - a[2,2,1,3] - a[2,3,1,3] - a[2,4,1,3] - a[2,5,1,3] - a[2,6,1,3]
   - a[2,7,1,3] - a[2,8,1,3] - a[2,9,1,3] - a[2,2,3,3] - a[2,3,3,3]
   - a[2,4,3,3] - a[2,5,3,3] - a[2,6,3,3] - a[2,7,3,3] - a[2,8,3,3]
   - a[2,9,3,3] + 2 b[3,2,2] - x[2,2] <= 0
 R3083: - a[3,3,1,3] - a[3,4,1,3] - a[3,5,1,3] - a[3,6,1,3] - a[3,7,1,3]
   - a[3,8,1,3] - a[3,9,1,3] - a[3,3,3,3] - a[3,4,3,3] - a[3,5,3,3]
   - a[3,6,3,3] - a[3,7,3,3] - a[3,8,3,3] - a[3,9,3,3] + 2 b[3,3,0]
   - x[3,0] <= 0
 R3084: - a[3,3,1,3] - a[3,4,1,3] - a[3,5,1,3] - a[3,6,1,3] - a[3,7,1,3]
   - a[3,8,1,3] - a[3,9,1,3] - a[3,3,3,3] - a[3,4,3,3] - a[3,5,3,3]
   - a[3,6,3,3] - a[3,7,3,3] - a[3,8,3,3] - a[3,9,3,3] + 2 b[3,3,1]
   - x[3,1] <= 0
 R3085: - a[3,3,1,3] - a[3,4,1,3] - a[3,5,1,3] - a[3,6,1,3] - a[3,7,1,3]
   - a[3,8,1,3] - a[3,9,1,3] - a[3,3,3,3] - a[3,4,3,3] - a[3,5,3,3]
   - a[3,6,3,3] - a[3,7,3,3] - a[3,8,3,3] - a[3,9,3,3] + 2 b[3,3,2]
   - x[3,2] <= 0
 R3086: - a[4,4,1,3] - a[4,5,1,3] - a[4,6,1,3] - a[4,7,1,3] - a[4,8,1,3]
   - a[4,9,1,3] - a[4,4,3,3] - a[4,5,3,3] - a[4,6,3,3] - a[4,7,3,3]
   - a[4,8,3,3] - a[4,9,3,3] + 2 b[3,4,0] - x[4,0] <= 0
 R3087: - a[4,4,1,3] - a[4,5,1,3] - a[4,6,1,3] - a[4,7,1,3] - a[4,8,1,3]
   - a[4,9,1,3] - a[4,4,3,3] - a[4,5,3,3] - a[4,6,3,3] - a[4,7,3,3]
   - a[4,8,3,3] - a[4,9,3,3] + 2 b[3,4,1] - x[4,1] <= 0
 R3088: - a[4,4,1,3] - a[4,5,1,3] - a[4,6,1,3] - a[4,7,1,3] - a[4,8,1,3]
   - a[4,9,1,3] - a[4,4,3,3] - a[4,5,3,3] - a[4,6,3,3] - a[4,7,3,3]
   - a[4,8,3,3] - a[4,9,3,3] + 2 b[3,4,2] - x[4,2] <= 0
 R3089: - a[5,5,1,3] - a[5,6,1,3] - a[5,7,1,3] - a[5,8,1,3] - a[5,9,1,3]
   - a[5,5,3,3] - a[5,6,3,3] - a[5,7,3,3] - a[5,8,3,3] - a[5,9,3,3]
   + 2 b[3,5,0] - x[5,0] <= 0
 R3090: - a[5,5,1,3] - a[5,6,1,3] - a[5,7,1,3] - a[5,8,1,3] - a[5,9,1,3]
   - a[5,5,3,3] - a[5,6,3,3] - a[5,7,3,3] - a[5,8,3,3] - a[5,9,3,3]
   + 2 b[3,5,1] - x[5,1] <= 0
 R3091: - a[5,5,1,3] - a[5,6,1,3] - a[5,7,1,3] - a[5,8,1,3] - a[5,9,1,3]
   - a[5,5,3,3] - a[5,6,3,3] - a[5,7,3,3] - a[5,8,3,3] - a[5,9,3,3]
   + 2 b[3,5,2] - x[5,2] <= 0
 R3092: - a[6,6,1,3] - a[6,7,1,3] - a[6,8,1,3] - a[6,9,1,3] - a[6,6,3,3]
   - a[6,7,3,3] - a[6,8,3,3] - a[6,9,3,3] + 2 b[3,6,0] - x[6,0] <= 0
 R3093: - a[6,6,1,3] - a[6,7,1,3] - a[6,8,1,3] - a[6,9,1,3] - a[6,6,3,3]
   - a[6,7,3,3] - a[6,8,3,3] - a[6,9,3,3] + 2 b[3,6,1] - x[6,1] <= 0
 R3094: - a[6,6,1,3] - a[6,7,1,3] - a[6,8,1,3] - a[6,9,1,3] - a[6,6,3,3]
   - a[6,7,3,3] - a[6,8,3,3] - a[6,9,3,3] + 2 b[3,6,2] - x[6,2] <= 0
 R3095: - a[7,7,1,3] - a[7,8,1,3] - a[7,9,1,3] - a[7,7,3,3] - a[7,8,3,3]
   - a[7,9,3,3] + 2 b[3,7,0] - x[7,0] <= 0
 R3096: - a[7,7,1,3] - a[7,8,1,3] - a[7,9,1,3] - a[7,7,3,3] - a[7,8,3,3]
   - a[7,9,3,3] + 2 b[3,7,1] - x[7,1] <= 0
 R3097: - a[7,7,1,3] - a[7,8,1,3] - a[7,9,1,3] - a[7,7,3,3] - a[7,8,3,3]
   - a[7,9,3,3] + 2 b[3,7,2] - x[7,2] <= 0
 R3098: - a[8,8,1,3] - a[8,9,1,3] - a[8,8,3,3] - a[8,9,3,3] + 2 b[3,8,0]
   - x[8,0] <= 0
 R3099: - a[8,8,1,3] - a[8,9,1,3] - a[8,8,3,3] - a[8,9,3,3] + 2 b[3,8,1]
   - x[8,1] <= 0
 R3100: - a[8,8,1,3] - a[8,9,1,3] - a[8,8,3,3] - a[8,9,3,3] + 2 b[3,8,2]
   - x[8,2] <= 0
 R3101: - a[9,9,1,3] - a[9,9,3,3] + 2 b[3,9,0] - x[9,0] <= 0
 R3102: - a[9,9,1,3] - a[9,9,3,3] + 2 b[3,9,1] - x[9,1] <= 0
 R3103: - a[9,9,1,3] - a[9,9,3,3] + 2 b[3,9,2] - x[9,2] <= 0
 R3104: - a[0,0,1,7] - a[0,1,1,7] - a[0,2,1,7] - a[0,3,1,7] - a[0,4,1,7]
   - a[0,5,1,7] - a[0,6,1,7] - a[0,7,1,7] - a[0,8,1,7] - a[0,9,1,7]
   - a[0,0,3,7] - a[0,1,3,7] - a[0,2,3,7] - a[0,3,3,7] - a[0,4,3,7]
   - a[0,5,3,7] - a[0,6,3,7] - a[0,7,3,7] - a[0,8,3,7] - a[0,9,3,7]
   - a[0,0,7,7] - a[0,1,7,7] - a[0,2,7,7] - a[0,3,7,7] - a[0,4,7,7]
   - a[0,5,7,7] - a[0,6,7,7] - a[0,7,7,7] - a[0,8,7,7] - a[0,9,7,7]
   + 2 b[7,0,0] - x[0,0] <= 0
 R3105: - a[0,0,1,7] - a[0,1,1,7] - a[0,2,1,7] - a[0,3,1,7] - a[0,4,1,7]
   - a[0,5,1,7] - a[0,6,1,7] - a[0,7,1,7] - a[0,8,1,7] - a[0,9,1,7]
   - a[0,0,3,7] - a[0,1,3,7] - a[0,2,3,7] - a[0,3,3,7] - a[0,4,3,7]
   - a[0,5,3,7] - a[0,6,3,7] - a[0,7,3,7] - a[0,8,3,7] - a[0,9,3,7]
   - a[0,0,7,7] - a[0,1,7,7] - a[0,2,7,7] - a[0,3,7,7] - a[0,4,7,7]
   - a[0,5,7,7] - a[0,6,7,7] - a[0,7,7,7] - a[0,8,7,7] - a[0,9,7,7]
   + 2 b[7,0,1] - x[0,1] <= 0
 R3106: - a[0,0,1,7] - a[0,1,1,7] - a[0,2,1,7] - a[0,3,1,7] - a[0,4,1,7]
   - a[0,5,1,7] - a[0,6,1,7] - a[0,7,1,7] - a[0,8,1,7] - a[0,9,1,7]
   - a[0,0,3,7] - a[0,1,3,7] - a[0,2,3,7] - a[0,3,3,7] - a[0,4,3,7]
   - a[0,5,3,7] - a[0,6,3,7] - a[0,7,3,7] - a[0,8,3,7] - a[0,9,3,7]
   - a[0,0,7,7] - a[0,1,7,7] - a[0,2,7,7] - a[0,3,7,7] - a[0,4,7,7]
   - a[0,5,7,7] - a[0,6,7,7] - a[0,7,7,7] - a[0,8,7,7] - a[0,9,7,7]
   + 2 b[7,0,2] - x[0,2] <= 0
 R3107: - a[1,1,1,7] - a[1,2,1,7] - a[1,3,1,7] - a[1,4,1,7] - a[1,5,1,7]
   - a[1,6,1,7] - a[1,7,1,7] - a[1,8,1,7] - a[1,9,1,7] - a[1,1,3,7]
   - a[1,2,3,7] - a[1,3,3,7] - a[1,4,3,7] - a[1,5,3,7] - a[1,6,3,7]
   - a[1,7,3,7] - a[1,8,3,7] - a[1,9,3,7] - a[1,1,7,7] - a[1,2,7,7]
   - a[1,3,7,7] - a[1,4,7,7] - a[1,5,7,7] - a[1,6,7,7] - a[1,7,7,7]
   - a[1,8,7,7] - a[1,9,7,7] + 2 b[7,1,0] - x[1,0] <= 0
 R3108: - a[1,1,1,7] - a[1,2,1,7] - a[1,3,1,7] - a[1,4,1,7] - a[1,5,1,7]
   - a[1,6,1,7] - a[1,7,1,7] - a[1,8,1,7] - a[1,9,1,7] - a[1,1,3,7]
   - a[1,2,3,7] - a[1,3,3,7] - a[1,4,3,7] - a[1,5,3,7] - a[1,6,3,7]
   - a[1,7,3,7] - a[1,8,3,7] - a[1,9,3,7] - a[1,1,7,7] - a[1,2,7,7]
   - a[1,3,7,7] - a[1,4,7,7] - a[1,5,7,7] - a[1,6,7,7] - a[1,7,7,7]
   - a[1,8,7,7] - a[1,9,7,7] + 2 b[7,1,1] - x[1,1] <= 0
 R3109: - a[1,1,1,7] - a[1,2,1,7] - a[1,3,1,7] - a[1,4,1,7] - a[1,5,1,7]
   - a[1,6,1,7] - a[1,7,1,7] - a[1,8,1,7] - a[1,9,1,7] - a[1,1,3,7]
   - a[1,2,3,7] - a[1,3,3,7] - a[1,4,3,7] - a[1,5,3,7] - a[1,6,3,7]
   - a[1,7,3,7] - a[1,8,3,7] - a[1,9,3,7] - a[1,1,7,7] - a[1,2,7,7]
   - a[1,3,7,7] - a[1,4,7,7] - a[1,5,7,7] - a[1,6,7,7] - a[1,7,7,7]
   - a[1,8,7,7] - a[1,9,7,7] + 2 b[7,1,2] - x[1,2] <= 0
 R3110: - a[2,2,1,7] - a[2,3,1,7] - a[2,4,1,7] - a[2,5,1,7] - a[2,6,1,7]
   - a[2,7,1,7] - a[2,8,1,7] - a[2,9,1,7] - a[2,2,3,7] - a[2,3,3,7]
   - a[2,4,3,7] - a[2,5,3,7] - a[2,6,3,7] - a[2,7,3,7] - a[2,8,3,7]
   - a[2,9,3,7] - a[2,2,7,7] - a[2,3,7,7] - a[2,4,7,7] - a[2,5,7,7]
   - a[2,6,7,7] - a[2,7,7,7] - a[2,8,7,7] - a[2,9,7,7] + 2 b[7,2,0]
   - x[2,0] <= 0
 R3111: - a[2,2,1,7] - a[2,3,1,7] - a[2,4,1,7] - a[2,5,1,7] - a[2,6,1,7]
   - a[2,7,1,7] - a[2,8,1,7] - a[2,9,1,7] - a[2,2,3,7] - a[2,3,3,7]
   - a[2,4,3,7] - a[2,5,3,7] - a[2,6,3,7] - a[2,7,3,7] - a[2,8,3,7]
   - a[2,9,3,7] - a[2,2,7,7] - a[2,3,7,7] - a[2,4,7,7] - a[2,5,7,7]
   - a[2,6,7,7] - a[2,7,7,7] - a[2,8,7,7] - a[2,9,7,7] + 2 b[7,2,1]
   - x[2,1] <= 0
 R3112: - a[2,2,1,7] - a[2,3,1,7] - a[2,4,1,7] - a[2,5,1,7] - a[2,6,1,7]
   - a[2,7,1,7] - a[2,8,1,7] - a[2,9,1,7] - a[2,2,3,7] - a[2,3,3,7]
   - a[2,4,3,7] - a[2,5,3,7] - a[2,6,3,7] - a[2,7,3,7] - a[2,8,3,7]
   - a[2,9,3,7] - a[2,2,7,7] - a[2,3,7,7] - a[2,4,7,7] - a[2,5,7,7]
   - a[2,6,7,7] - a[2,7,7,7] - a[2,8,7,7] - a[2,9,7,7] + 2 b[7,2,2]
   - x[2,2] <= 0
 R3113: - a[3,3,1,7] - a[3,4,1,7] - a[3,5,1,7] - a[3,6,1,7] - a[3,7,1,7]
   - a[3,8,1,7] - a[3,9,1,7] - a[3,3,3,7] - a[3,4,3,7] - a[3,5,3,7]
   - a[3,6,3,7] - a[3,7,3,7] - a[3,8,3,7] - a[3,9,3,7] - a[3,3,7,7]
   - a[3,4,7,7] - a[3,5,7,7] - a[3,6,7,7] - a[3,7,7,7] - a[3,8,7,7]
   - a[3,9,7,7] + 2 b[7,3,0] - x[3,0] <= 0
 R3114: - a[3,3,1,7] - a[3,4,1,7] - a[3,5,1,7] - a[3,6,1,7] - a[3,7,1,7]
   - a[3,8,1,7] - a[3,9,1,7] - a[3,3,3,7] - a[3,4,3,7] - a[3,5,3,7]
   - a[3,6,3,7] - a[3,7,3,7] - a[3,8,3,7] - a[3,9,3,7] - a[3,3,7,7]
   - a[3,4,7,7] - a[3,5,7,7] - a[3,6,7,7] - a[3,7,7,7] - a[3,8,7,7]
   - a[3,9,7,7] + 2 b[7,3,1] - x[3,1] <= 0
 R3115: - a[3,3,1,7] - a[3,4,1,7] - a[3,5,1,7] - a[3,6,1,7] - a[3,7,1,7]
   - a[3,8,1,7] - a[3,9,1,7] - a[3,3,3,7] - a[3,4,3,7] - a[3,5,3,7]
   - a[3,6,3,7] - a[3,7,3,7] - a[3,8,3,7] - a[3,9,3,7] - a[3,3,7,7]
   - a[3,4,7,7] - a[3,5,7,7] - a[3,6,7,7] - a[3,7,7,7] - a[3,8,7,7]
   - a[3,9,7,7] + 2 b[7,3,2] - x[3,2] <= 0
 R3116: - a[4,4,1,7] - a[4,5,1,7] - a[4,6,1,7] - a[4,7,1,7] - a[4,8,1,7]
   - a[4,9,1,7] - a[4,4,3,7] - a[4,5,3,7] - a[4,6,3,7] - a[4,7,3,7]
   - a[4,8,3,7] - a[4,9,3,7] - a[4,4,7,7] - a[4,5,7,7] - a[4,6,7,7]
   - a[4,7,7,7] - a[4,8,7,7] - a[4,9,7,7] + 2 b[7,4,0] - x[4,0] <= 0
 R3117: - a[4,4,1,7] - a[4,5,1,7] - a[4,6,1,7] - a[4,7,1,7] - a[4,8,1,7]
   - a[4,9,1,7] - a[4,4,3,7] - a[4,5,3,7] - a[4,6,3,7] - a[4,7,3,7]
   - a[4,8,3,7] - a[4,9,3,7] - a[4,4,7,7] - a[4,5,7,7] - a[4,6,7,7]
   - a[4,7,7,7] - a[4,8,7,7] - a[4,9,7,7] + 2 b[7,4,1] - x[4,1] <= 0
 R3118: - a[4,4,1,7] - a[4,5,1,7] - a[4,6,1,7] - a[4,7,1,7] - a[4,8,1,7]
   - a[4,9,1,7] - a[4,4,3,7] - a[4,5,3,7] - a[4,6,3,7] - a[4,7,3,7]
   - a[4,8,3,7] - a[4,9,3,7] - a[4,4,7,7] - a[4,5,7,7] - a[4,6,7,7]
   - a[4,7,7,7] - a[4,8,7,7] - a[4,9,7,7] + 2 b[7,4,2] - x[4,2] <= 0
 R3119: - a[5,5,1,7] - a[5,6,1,7] - a[5,7,1,7] - a[5,8,1,7] - a[5,9,1,7]
   - a[5,5,3,7] - a[5,6,3,7] - a[5,7,3,7] - a[5,8,3,7] - a[5,9,3,7]
   - a[5,5,7,7] - a[5,6,7,7] - a[5,7,7,7] - a[5,8,7,7] - a[5,9,7,7]
   + 2 b[7,5,0] - x[5,0] <= 0
 R3120: - a[5,5,1,7] - a[5,6,1,7] - a[5,7,1,7] - a[5,8,1,7] - a[5,9,1,7]
   - a[5,5,3,7] - a[5,6,3,7] - a[5,7,3,7] - a[5,8,3,7] - a[5,9,3,7]
   - a[5,5,7,7] - a[5,6,7,7] - a[5,7,7,7] - a[5,8,7,7] - a[5,9,7,7]
   + 2 b[7,5,1] - x[5,1] <= 0
 R3121: - a[5,5,1,7] - a[5,6,1,7] - a[5,7,1,7] - a[5,8,1,7] - a[5,9,1,7]
   - a[5,5,3,7] - a[5,6,3,7] - a[5,7,3,7] - a[5,8,3,7] - a[5,9,3,7]
   - a[5,5,7,7] - a[5,6,7,7] - a[5,7,7,7] - a[5,8,7,7] - a[5,9,7,7]
   + 2 b[7,5,2] - x[5,2] <= 0
 R3122: - a[6,6,1,7] - a[6,7,1,7] - a[6,8,1,7] - a[6,9,1,7] - a[6,6,3,7]
   - a[6,7,3,7] - a[6,8,3,7] - a[6,9,3,7] - a[6,6,7,7] - a[6,7,7,7]
   - a[6,8,7,7] - a[6,9,7,7] + 2 b[7,6,0] - x[6,0] <= 0
 R3123: - a[6,6,1,7] - a[6,7,1,7] - a[6,8,1,7] - a[6,9,1,7] - a[6,6,3,7]
   - a[6,7,3,7] - a[6,8,3,7] - a[6,9,3,7] - a[6,6,7,7] - a[6,7,7,7]
   - a[6,8,7,7] - a[6,9,7,7] + 2 b[7,6,1] - x[6,1] <= 0
 R3124: - a[6,6,1,7] - a[6,7,1,7] - a[6,8,1,7] - a[6,9,1,7] - a[6,6,3,7]
   - a[6,7,3,7] - a[6,8,3,7] - a[6,9,3,7] - a[6,6,7,7] - a[6,7,7,7]
   - a[6,8,7,7] - a[6,9,7,7] + 2 b[7,6,2] - x[6,2] <= 0
 R3125: - a[7,7,1,7] - a[7,8,1,7] - a[7,9,1,7] - a[7,7,3,7] - a[7,8,3,7]
   - a[7,9,3,7] - a[7,7,7,7] - a[7,8,7,7] - a[7,9,7,7] + 2 b[7,7,0]
   - x[7,0] <= 0
 R3126: - a[7,7,1,7] - a[7,8,1,7] - a[7,9,1,7] - a[7,7,3,7] - a[7,8,3,7]
   - a[7,9,3,7] - a[7,7,7,7] - a[7,8,7,7] - a[7,9,7,7] + 2 b[7,7,1]
   - x[7,1] <= 0
 R3127: - a[7,7,1,7] - a[7,8,1,7] - a[7,9,1,7] - a[7,7,3,7] - a[7,8,3,7]
   - a[7,9,3,7] - a[7,7,7,7] - a[7,8,7,7] - a[7,9,7,7] + 2 b[7,7,2]
   - x[7,2] <= 0
 R3128: - a[8,8,1,7] - a[8,9,1,7] - a[8,8,3,7] - a[8,9,3,7] - a[8,8,7,7]
   - a[8,9,7,7] + 2 b[7,8,0] - x[8,0] <= 0
 R3129: - a[8,8,1,7] - a[8,9,1,7] - a[8,8,3,7] - a[8,9,3,7] - a[8,8,7,7]
   - a[8,9,7,7] + 2 b[7,8,1] - x[8,1] <= 0
 R3130: - a[8,8,1,7] - a[8,9,1,7] - a[8,8,3,7] - a[8,9,3,7] - a[8,8,7,7]
   - a[8,9,7,7] + 2 b[7,8,2] - x[8,2] <= 0
 R3131: - a[9,9,1,7] - a[9,9,3,7] - a[9,9,7,7] + 2 b[7,9,0] - x[9,0] <= 0
 R3132: - a[9,9,1,7] - a[9,9,3,7] - a[9,9,7,7] + 2 b[7,9,1] - x[9,1] <= 0
 R3133: - a[9,9,1,7] - a[9,9,3,7] - a[9,9,7,7] + 2 b[7,9,2] - x[9,2] <= 0
 R3134: - a[0,0,4,4] - a[0,1,4,4] - a[0,2,4,4] - a[0,3,4,4] - a[0,4,4,4]
   - a[0,5,4,4] - a[0,6,4,4] - a[0,7,4,4] - a[0,8,4,4] - a[0,9,4,4]
   + 2 b[4,0,0] - x[0,0] <= 0
 R3135: - a[0,0,4,4] - a[0,1,4,4] - a[0,2,4,4] - a[0,3,4,4] - a[0,4,4,4]
   - a[0,5,4,4] - a[0,6,4,4] - a[0,7,4,4] - a[0,8,4,4] - a[0,9,4,4]
   + 2 b[4,0,1] - x[0,1] <= 0
 R3136: - a[0,0,4,4] - a[0,1,4,4] - a[0,2,4,4] - a[0,3,4,4] - a[0,4,4,4]
   - a[0,5,4,4] - a[0,6,4,4] - a[0,7,4,4] - a[0,8,4,4] - a[0,9,4,4]
   + 2 b[4,0,2] - x[0,2] <= 0
 R3137: - a[1,1,4,4] - a[1,2,4,4] - a[1,3,4,4] - a[1,4,4,4] - a[1,5,4,4]
   - a[1,6,4,4] - a[1,7,4,4] - a[1,8,4,4] - a[1,9,4,4] + 2 b[4,1,0]
   - x[1,0] <= 0
 R3138: - a[1,1,4,4] - a[1,2,4,4] - a[1,3,4,4] - a[1,4,4,4] - a[1,5,4,4]
   - a[1,6,4,4] - a[1,7,4,4] - a[1,8,4,4] - a[1,9,4,4] + 2 b[4,1,1]
   - x[1,1] <= 0
 R3139: - a[1,1,4,4] - a[1,2,4,4] - a[1,3,4,4] - a[1,4,4,4] - a[1,5,4,4]
   - a[1,6,4,4] - a[1,7,4,4] - a[1,8,4,4] - a[1,9,4,4] + 2 b[4,1,2]
   - x[1,2] <= 0
 R3140: - a[2,2,4,4] - a[2,3,4,4] - a[2,4,4,4] - a[2,5,4,4] - a[2,6,4,4]
   - a[2,7,4,4] - a[2,8,4,4] - a[2,9,4,4] + 2 b[4,2,0] - x[2,0] <= 0
 R3141: - a[2,2,4,4] - a[2,3,4,4] - a[2,4,4,4] - a[2,5,4,4] - a[2,6,4,4]
   - a[2,7,4,4] - a[2,8,4,4] - a[2,9,4,4] + 2 b[4,2,1] - x[2,1] <= 0
 R3142: - a[2,2,4,4] - a[2,3,4,4] - a[2,4,4,4] - a[2,5,4,4] - a[2,6,4,4]
   - a[2,7,4,4] - a[2,8,4,4] - a[2,9,4,4] + 2 b[4,2,2] - x[2,2] <= 0
 R3143: - a[3,3,4,4] - a[3,4,4,4] - a[3,5,4,4] - a[3,6,4,4] - a[3,7,4,4]
   - a[3,8,4,4] - a[3,9,4,4] + 2 b[4,3,0] - x[3,0] <= 0
 R3144: - a[3,3,4,4] - a[3,4,4,4] - a[3,5,4,4] - a[3,6,4,4] - a[3,7,4,4]
   - a[3,8,4,4] - a[3,9,4,4] + 2 b[4,3,1] - x[3,1] <= 0
 R3145: - a[3,3,4,4] - a[3,4,4,4] - a[3,5,4,4] - a[3,6,4,4] - a[3,7,4,4]
   - a[3,8,4,4] - a[3,9,4,4] + 2 b[4,3,2] - x[3,2] <= 0
 R3146: - a[4,4,4,4] - a[4,5,4,4] - a[4,6,4,4] - a[4,7,4,4] - a[4,8,4,4]
   - a[4,9,4,4] + 2 b[4,4,0] - x[4,0] <= 0
 R3147: - a[4,4,4,4] - a[4,5,4,4] - a[4,6,4,4] - a[4,7,4,4] - a[4,8,4,4]
   - a[4,9,4,4] + 2 b[4,4,1] - x[4,1] <= 0
 R3148: - a[4,4,4,4] - a[4,5,4,4] - a[4,6,4,4] - a[4,7,4,4] - a[4,8,4,4]
   - a[4,9,4,4] + 2 b[4,4,2] - x[4,2] <= 0
 R3149: - a[5,5,4,4] - a[5,6,4,4] - a[5,7,4,4] - a[5,8,4,4] - a[5,9,4,4]
   + 2 b[4,5,0] - x[5,0] <= 0
 R3150: - a[5,5,4,4] - a[5,6,4,4] - a[5,7,4,4] - a[5,8,4,4] - a[5,9,4,4]
   + 2 b[4,5,1] - x[5,1] <= 0
 R3151: - a[5,5,4,4] - a[5,6,4,4] - a[5,7,4,4] - a[5,8,4,4] - a[5,9,4,4]
   + 2 b[4,5,2] - x[5,2] <= 0
 R3152: - a[6,6,4,4] - a[6,7,4,4] - a[6,8,4,4] - a[6,9,4,4] + 2 b[4,6,0]
   - x[6,0] <= 0
 R3153: - a[6,6,4,4] - a[6,7,4,4] - a[6,8,4,4] - a[6,9,4,4] + 2 b[4,6,1]
   - x[6,1] <= 0
 R3154: - a[6,6,4,4] - a[6,7,4,4] - a[6,8,4,4] - a[6,9,4,4] + 2 b[4,6,2]
   - x[6,2] <= 0
 R3155: - a[7,7,4,4] - a[7,8,4,4] - a[7,9,4,4] + 2 b[4,7,0] - x[7,0] <= 0
 R3156: - a[7,7,4,4] - a[7,8,4,4] - a[7,9,4,4] + 2 b[4,7,1] - x[7,1] <= 0
 R3157: - a[7,7,4,4] - a[7,8,4,4] - a[7,9,4,4] + 2 b[4,7,2] - x[7,2] <= 0
 R3158: - a[8,8,4,4] - a[8,9,4,4] + 2 b[4,8,0] - x[8,0] <= 0
 R3159: - a[8,8,4,4] - a[8,9,4,4] + 2 b[4,8,1] - x[8,1] <= 0
 R3160: - a[8,8,4,4] - a[8,9,4,4] + 2 b[4,8,2] - x[8,2] <= 0
 R3161: - a[9,9,4,4] + 2 b[4,9,0] - x[9,0] <= 0
 R3162: - a[9,9,4,4] + 2 b[4,9,1] - x[9,1] <= 0
 R3163: - a[9,9,4,4] + 2 b[4,9,2] - x[9,2] <= 0
 R3164: - a[0,0,5,5] - a[0,1,5,5] - a[0,2,5,5] - a[0,3,5,5] - a[0,4,5,5]
   - a[0,5,5,5] - a[0,6,5,5] - a[0,7,5,5] - a[0,8,5,5] - a[0,9,5,5]
   + 2 b[5,0,0] - x[0,0] <= 0
 R3165: - a[0,0,5,5] - a[0,1,5,5] - a[0,2,5,5] - a[0,3,5,5] - a[0,4,5,5]
   - a[0,5,5,5] - a[0,6,5,5] - a[0,7,5,5] - a[0,8,5,5] - a[0,9,5,5]
   + 2 b[5,0,1] - x[0,1] <= 0
 R3166: - a[0,0,5,5] - a[0,1,5,5] - a[0,2,5,5] - a[0,3,5,5] - a[0,4,5,5]
   - a[0,5,5,5] - a[0,6,5,5] - a[0,7,5,5] - a[0,8,5,5] - a[0,9,5,5]
   + 2 b[5,0,2] - x[0,2] <= 0
 R3167: - a[1,1,5,5] - a[1,2,5,5] - a[1,3,5,5] - a[1,4,5,5] - a[1,5,5,5]
   - a[1,6,5,5] - a[1,7,5,5] - a[1,8,5,5] - a[1,9,5,5] + 2 b[5,1,0]
   - x[1,0] <= 0
 R3168: - a[1,1,5,5] - a[1,2,5,5] - a[1,3,5,5] - a[1,4,5,5] - a[1,5,5,5]
   - a[1,6,5,5] - a[1,7,5,5] - a[1,8,5,5] - a[1,9,5,5] + 2 b[5,1,1]
   - x[1,1] <= 0
 R3169: - a[1,1,5,5] - a[1,2,5,5] - a[1,3,5,5] - a[1,4,5,5] - a[1,5,5,5]
   - a[1,6,5,5] - a[1,7,5,5] - a[1,8,5,5] - a[1,9,5,5] + 2 b[5,1,2]
   - x[1,2] <= 0
 R3170: - a[2,2,5,5] - a[2,3,5,5] - a[2,4,5,5] - a[2,5,5,5] - a[2,6,5,5]
   - a[2,7,5,5] - a[2,8,5,5] - a[2,9,5,5] + 2 b[5,2,0] - x[2,0] <= 0
 R3171: - a[2,2,5,5] - a[2,3,5,5] - a[2,4,5,5] - a[2,5,5,5] - a[2,6,5,5]
   - a[2,7,5,5] - a[2,8,5,5] - a[2,9,5,5] + 2 b[5,2,1] - x[2,1] <= 0
 R3172: - a[2,2,5,5] - a[2,3,5,5] - a[2,4,5,5] - a[2,5,5,5] - a[2,6,5,5]
   - a[2,7,5,5] - a[2,8,5,5] - a[2,9,5,5] + 2 b[5,2,2] - x[2,2] <= 0
 R3173: - a[3,3,5,5] - a[3,4,5,5] - a[3,5,5,5] - a[3,6,5,5] - a[3,7,5,5]
   - a[3,8,5,5] - a[3,9,5,5] + 2 b[5,3,0] - x[3,0] <= 0
 R3174: - a[3,3,5,5] - a[3,4,5,5] - a[3,5,5,5] - a[3,6,5,5] - a[3,7,5,5]
   - a[3,8,5,5] - a[3,9,5,5] + 2 b[5,3,1] - x[3,1] <= 0
 R3175: - a[3,3,5,5] - a[3,4,5,5] - a[3,5,5,5] - a[3,6,5,5] - a[3,7,5,5]
   - a[3,8,5,5] - a[3,9,5,5] + 2 b[5,3,2] - x[3,2] <= 0
 R3176: - a[4,4,5,5] - a[4,5,5,5] - a[4,6,5,5] - a[4,7,5,5] - a[4,8,5,5]
   - a[4,9,5,5] + 2 b[5,4,0] - x[4,0] <= 0
 R3177: - a[4,4,5,5] - a[4,5,5,5] - a[4,6,5,5] - a[4,7,5,5] - a[4,8,5,5]
   - a[4,9,5,5] + 2 b[5,4,1] - x[4,1] <= 0
 R3178: - a[4,4,5,5] - a[4,5,5,5] - a[4,6,5,5] - a[4,7,5,5] - a[4,8,5,5]
   - a[4,9,5,5] + 2 b[5,4,2] - x[4,2] <= 0
 R3179: - a[5,5,5,5] - a[5,6,5,5] - a[5,7,5,5] - a[5,8,5,5] - a[5,9,5,5]
   + 2 b[5,5,0] - x[5,0] <= 0
 R3180: - a[5,5,5,5] - a[5,6,5,5] - a[5,7,5,5] - a[5,8,5,5] - a[5,9,5,5]
   + 2 b[5,5,1] - x[5,1] <= 0
 R3181: - a[5,5,5,5] - a[5,6,5,5] - a[5,7,5,5] - a[5,8,5,5] - a[5,9,5,5]
   + 2 b[5,5,2] - x[5,2] <= 0
 R3182: - a[6,6,5,5] - a[6,7,5,5] - a[6,8,5,5] - a[6,9,5,5] + 2 b[5,6,0]
   - x[6,0] <= 0
 R3183: - a[6,6,5,5] - a[6,7,5,5] - a[6,8,5,5] - a[6,9,5,5] + 2 b[5,6,1]
   - x[6,1] <= 0
 R3184: - a[6,6,5,5] - a[6,7,5,5] - a[6,8,5,5] - a[6,9,5,5] + 2 b[5,6,2]
   - x[6,2] <= 0
 R3185: - a[7,7,5,5] - a[7,8,5,5] - a[7,9,5,5] + 2 b[5,7,0] - x[7,0] <= 0
 R3186: - a[7,7,5,5] - a[7,8,5,5] - a[7,9,5,5] + 2 b[5,7,1] - x[7,1] <= 0
 R3187: - a[7,7,5,5] - a[7,8,5,5] - a[7,9,5,5] + 2 b[5,7,2] - x[7,2] <= 0
 R3188: - a[8,8,5,5] - a[8,9,5,5] + 2 b[5,8,0] - x[8,0] <= 0
 R3189: - a[8,8,5,5] - a[8,9,5,5] + 2 b[5,8,1] - x[8,1] <= 0
 R3190: - a[8,8,5,5] - a[8,9,5,5] + 2 b[5,8,2] - x[8,2] <= 0
 R3191: - a[9,9,5,5] + 2 b[5,9,0] - x[9,0] <= 0
 R3192: - a[9,9,5,5] + 2 b[5,9,1] - x[9,1] <= 0
 R3193: - a[9,9,5,5] + 2 b[5,9,2] - x[9,2] <= 0
 R3194: - a[0,0,0,0] - a[0,1,0,0] - a[0,2,0,0] - a[0,3,0,0] - a[0,4,0,0]
   - a[0,5,0,0] - a[0,6,0,0] - a[0,7,0,0] - a[0,8,0,0] - a[0,9,0,0]
   + b[0,0,0] - x[0,0] >= -1
 R3195: - a[0,0,0,0] - a[0,1,0,0] - a[0,2,0,0] - a[0,3,0,0] - a[0,4,0,0]
   - a[0,5,0,0] - a[0,6,0,0] - a[0,7,0,0] - a[0,8,0,0] - a[0,9,0,0]
   + b[0,0,1] - x[0,1] >= -1
 R3196: - a[0,0,0,0] - a[0,1,0,0] - a[0,2,0,0] - a[0,3,0,0] - a[0,4,0,0]
   - a[0,5,0,0] - a[0,6,0,0] - a[0,7,0,0] - a[0,8,0,0] - a[0,9,0,0]
   + b[0,0,2] - x[0,2] >= -1
 R3197: - a[1,1,0,0] - a[1,2,0,0] - a[1,3,0,0] - a[1,4,0,0] - a[1,5,0,0]
   - a[1,6,0,0] - a[1,7,0,0] - a[1,8,0,0] - a[1,9,0,0] + b[0,1,0] - x[1,0]
   >= -1
 R3198: - a[1,1,0,0] - a[1,2,0,0] - a[1,3,0,0] - a[1,4,0,0] - a[1,5,0,0]
   - a[1,6,0,0] - a[1,7,0,0] - a[1,8,0,0] - a[1,9,0,0] + b[0,1,1] - x[1,1]
   >= -1
 R3199: - a[1,1,0,0] - a[1,2,0,0] - a[1,3,0,0] - a[1,4,0,0] - a[1,5,0,0]
   - a[1,6,0,0] - a[1,7,0,0] - a[1,8,0,0] - a[1,9,0,0] + b[0,1,2] - x[1,2]
   >= -1
 R3200: - a[2,2,0,0] - a[2,3,0,0] - a[2,4,0,0] - a[2,5,0,0] - a[2,6,0,0]
   - a[2,7,0,0] - a[2,8,0,0] - a[2,9,0,0] + b[0,2,0] - x[2,0] >= -1
 R3201: - a[2,2,0,0] - a[2,3,0,0] - a[2,4,0,0] - a[2,5,0,0] - a[2,6,0,0]
   - a[2,7,0,0] - a[2,8,0,0] - a[2,9,0,0] + b[0,2,1] - x[2,1] >= -1
 R3202: - a[2,2,0,0] - a[2,3,0,0] - a[2,4,0,0] - a[2,5,0,0] - a[2,6,0,0]
   - a[2,7,0,0] - a[2,8,0,0] - a[2,9,0,0] + b[0,2,2] - x[2,2] >= -1
 R3203: - a[3,3,0,0] - a[3,4,0,0] - a[3,5,0,0] - a[3,6,0,0] - a[3,7,0,0]
   - a[3,8,0,0] - a[3,9,0,0] + b[0,3,0] - x[3,0] >= -1
 R3204: - a[3,3,0,0] - a[3,4,0,0] - a[3,5,0,0] - a[3,6,0,0] - a[3,7,0,0]
   - a[3,8,0,0] - a[3,9,0,0] + b[0,3,1] - x[3,1] >= -1
 R3205: - a[3,3,0,0] - a[3,4,0,0] - a[3,5,0,0] - a[3,6,0,0] - a[3,7,0,0]
   - a[3,8,0,0] - a[3,9,0,0] + b[0,3,2] - x[3,2] >= -1
 R3206: - a[4,4,0,0] - a[4,5,0,0] - a[4,6,0,0] - a[4,7,0,0] - a[4,8,0,0]
   - a[4,9,0,0] + b[0,4,0] - x[4,0] >= -1
 R3207: - a[4,4,0,0] - a[4,5,0,0] - a[4,6,0,0] - a[4,7,0,0] - a[4,8,0,0]
   - a[4,9,0,0] + b[0,4,1] - x[4,1] >= -1
 R3208: - a[4,4,0,0] - a[4,5,0,0] - a[4,6,0,0] - a[4,7,0,0] - a[4,8,0,0]
   - a[4,9,0,0] + b[0,4,2] - x[4,2] >= -1
 R3209: - a[5,5,0,0] - a[5,6,0,0] - a[5,7,0,0] - a[5,8,0,0] - a[5,9,0,0]
   + b[0,5,0] - x[5,0] >= -1
 R3210: - a[5,5,0,0] - a[5,6,0,0] - a[5,7,0,0] - a[5,8,0,0] - a[5,9,0,0]
   + b[0,5,1] - x[5,1] >= -1
 R3211: - a[5,5,0,0] - a[5,6,0,0] - a[5,7,0,0] - a[5,8,0,0] - a[5,9,0,0]
   + b[0,5,2] - x[5,2] >= -1
 R3212: - a[6,6,0,0] - a[6,7,0,0] - a[6,8,0,0] - a[6,9,0,0] + b[0,6,0]
   - x[6,0] >= -1
 R3213: - a[6,6,0,0] - a[6,7,0,0] - a[6,8,0,0] - a[6,9,0,0] + b[0,6,1]
   - x[6,1] >= -1
 R3214: - a[6,6,0,0] - a[6,7,0,0] - a[6,8,0,0] - a[6,9,0,0] + b[0,6,2]
   - x[6,2] >= -1
 R3215: - a[7,7,0,0] - a[7,8,0,0] - a[7,9,0,0] + b[0,7,0] - x[7,0] >= -1
 R3216: - a[7,7,0,0] - a[7,8,0,0] - a[7,9,0,0] + b[0,7,1] - x[7,1] >= -1
 R3217: - a[7,7,0,0] - a[7,8,0,0] - a[7,9,0,0] + b[0,7,2] - x[7,2] >= -1
 R3218: - a[8,8,0,0] - a[8,9,0,0] + b[0,8,0] - x[8,0] >= -1
 R3219: - a[8,8,0,0] - a[8,9,0,0] + b[0,8,1] - x[8,1] >= -1
 R3220: - a[8,8,0,0] - a[8,9,0,0] + b[0,8,2] - x[8,2] >= -1
 R3221: - a[9,9,0,0] + b[0,9,0] - x[9,0] >= -1
 R3222: - a[9,9,0,0] + b[0,9,1] - x[9,1] >= -1
 R3223: - a[9,9,0,0] + b[0,9,2] - x[9,2] >= -1
 R3224: - a[0,0,0,2] - a[0,1,0,2] - a[0,2,0,2] - a[0,3,0,2] - a[0,4,0,2]
   - a[0,5,0,2] - a[0,6,0,2] - a[0,7,0,2] - a[0,8,0,2] - a[0,9,0,2]
   - a[0,0,2,2] - a[0,1,2,2] - a[0,2,2,2] - a[0,3,2,2] - a[0,4,2,2]
   - a[0,5,2,2] - a[0,6,2,2] - a[0,7,2,2] - a[0,8,2,2] - a[0,9,2,2]
   + b[2,0,0] - x[0,0] >= -1
 R3225: - a[0,0,0,2] - a[0,1,0,2] - a[0,2,0,2] - a[0,3,0,2] - a[0,4,0,2]
   - a[0,5,0,2] - a[0,6,0,2] - a[0,7,0,2] - a[0,8,0,2] - a[0,9,0,2]
   - a[0,0,2,2] - a[0,1,2,2] - a[0,2,2,2] - a[0,3,2,2] - a[0,4,2,2]
   - a[0,5,2,2] - a[0,6,2,2] - a[0,7,2,2] - a[0,8,2,2] - a[0,9,2,2]
   + b[2,0,1] - x[0,1] >= -1
 R3226: - a[0,0,0,2] - a[0,1,0,2] - a[0,2,0,2] - a[0,3,0,2] - a[0,4,0,2]
   - a[0,5,0,2] - a[0,6,0,2] - a[0,7,0,2] - a[0,8,0,2] - a[0,9,0,2]
   - a[0,0,2,2] - a[0,1,2,2] - a[0,2,2,2] - a[0,3,2,2] - a[0,4,2,2]
   - a[0,5,2,2] - a[0,6,2,2] - a[0,7,2,2] - a[0,8,2,2] - a[0,9,2,2]
   + b[2,0,2] - x[0,2] >= -1
 R3227: - a[1,1,0,2] - a[1,2,0,2] - a[1,3,0,2] - a[1,4,0,2] - a[1,5,0,2]
   - a[1,6,0,2] - a[1,7,0,2] - a[1,8,0,2] - a[1,9,0,2] - a[1,1,2,2]
   - a[1,2,2,2] - a[1,3,2,2] - a[1,4,2,2] - a[1,5,2,2] - a[1,6,2,2]
   - a[1,7,2,2] - a[1,8,2,2] - a[1,9,2,2] + b[2,1,0] - x[1,0] >= -1
 R3228: - a[1,1,0,2] - a[1,2,0,2] - a[1,3,0,2] - a[1,4,0,2] - a[1,5,0,2]
   - a[1,6,0,2] - a[1,7,0,2] - a[1,8,0,2] - a[1,9,0,2] - a[1,1,2,2]
   - a[1,2,2,2] - a[1,3,2,2] - a[1,4,2,2] - a[1,5,2,2] - a[1,6,2,2]
   - a[1,7,2,2] - a[1,8,2,2] - a[1,9,2,2] + b[2,1,1] - x[1,1] >= -1
 R3229: - a[1,1,0,2] - a[1,2,0,2] - a[1,3,0,2] - a[1,4,0,2] - a[1,5,0,2]
   - a[1,6,0,2] - a[1,7,0,2] - a[1,8,0,2] - a[1,9,0,2] - a[1,1,2,2]
   - a[1,2,2,2] - a[1,3,2,2] - a[1,4,2,2] - a[1,5,2,2] - a[1,6,2,2]
   - a[1,7,2,2] - a[1,8,2,2] - a[1,9,2,2] + b[2,1,2] - x[1,2] >= -1
 R3230: - a[2,2,0,2] - a[2,3,0,2] - a[2,4,0,2] - a[2,5,0,2] - a[2,6,0,2]
   - a[2,7,0,2] - a[2,8,0,2] - a[2,9,0,2] - a[2,2,2,2] - a[2,3,2,2]
   - a[2,4,2,2] - a[2,5,2,2] - a[2,6,2,2] - a[2,7,2,2] - a[2,8,2,2]
   - a[2,9,2,2] + b[2,2,0] - x[2,0] >= -1
 R3231: - a[2,2,0,2] - a[2,3,0,2] - a[2,4,0,2] - a[2,5,0,2] - a[2,6,0,2]
   - a[2,7,0,2] - a[2,8,0,2] - a[2,9,0,2] - a[2,2,2,2] - a[2,3,2,2]
   - a[2,4,2,2] - a[2,5,2,2] - a[2,6,2,2] - a[2,7,2,2] - a[2,8,2,2]
   - a[2,9,2,2] + b[2,2,1] - x[2,1] >= -1
 R3232: - a[2,2,0,2] - a[2,3,0,2] - a[2,4,0,2] - a[2,5,0,2] - a[2,6,0,2]
   - a[2,7,0,2] - a[2,8,0,2] - a[2,9,0,2] - a[2,2,2,2] - a[2,3,2,2]
   - a[2,4,2,2] - a[2,5,2,2] - a[2,6,2,2] - a[2,7,2,2] - a[2,8,2,2]
   - a[2,9,2,2] + b[2,2,2] - x[2,2] >= -1
 R3233: - a[3,3,0,2] - a[3,4,0,2] - a[3,5,0,2] - a[3,6,0,2] - a[3,7,0,2]
   - a[3,8,0,2] - a[3,9,0,2] - a[3,3,2,2] - a[3,4,2,2] - a[3,5,2,2]
   - a[3,6,2,2] - a[3,7,2,2] - a[3,8,2,2] - a[3,9,2,2] + b[2,3,0] - x[3,0]
   >= -1
 R3234: - a[3,3,0,2] - a[3,4,0,2] - a[3,5,0,2] - a[3,6,0,2] - a[3,7,0,2]
   - a[3,8,0,2] - a[3,9,0,2] - a[3,3,2,2] - a[3,4,2,2] - a[3,5,2,2]
   - a[3,6,2,2] - a[3,7,2,2] - a[3,8,2,2] - a[3,9,2,2] + b[2,3,1] - x[3,1]
   >= -1
 R3235: - a[3,3,0,2] - a[3,4,0,2] - a[3,5,0,2] - a[3,6,0,2] - a[3,7,0,2]
   - a[3,8,0,2] - a[3,9,0,2] - a[3,3,2,2] - a[3,4,2,2] - a[3,5,2,2]
   - a[3,6,2,2] - a[3,7,2,2] - a[3,8,2,2] - a[3,9,2,2] + b[2,3,2] - x[3,2]
   >= -1
 R3236: - a[4,4,0,2] - a[4,5,0,2] - a[4,6,0,2] - a[4,7,0,2] - a[4,8,0,2]
   - a[4,9,0,2] - a[4,4,2,2] - a[4,5,2,2] - a[4,6,2,2] - a[4,7,2,2]
   - a[4,8,2,2] - a[4,9,2,2] + b[2,4,0] - x[4,0] >= -1
 R3237: - a[4,4,0,2] - a[4,5,0,2] - a[4,6,0,2] - a[4,7,0,2] - a[4,8,0,2]
   - a[4,9,0,2] - a[4,4,2,2] - a[4,5,2,2] - a[4,6,2,2] - a[4,7,2,2]
   - a[4,8,2,2] - a[4,9,2,2] + b[2,4,1] - x[4,1] >= -1
 R3238: - a[4,4,0,2] - a[4,5,0,2] - a[4,6,0,2] - a[4,7,0,2] - a[4,8,0,2]
   - a[4,9,0,2] - a[4,4,2,2] - a[4,5,2,2] - a[4,6,2,2] - a[4,7,2,2]
   - a[4,8,2,2] - a[4,9,2,2] + b[2,4,2] - x[4,2] >= -1
 R3239: - a[5,5,0,2] - a[5,6,0,2] - a[5,7,0,2] - a[5,8,0,2] - a[5,9,0,2]
   - a[5,5,2,2] - a[5,6,2,2] - a[5,7,2,2] - a[5,8,2,2] - a[5,9,2,2]
   + b[2,5,0] - x[5,0] >= -1
 R3240: - a[5,5,0,2] - a[5,6,0,2] - a[5,7,0,2] - a[5,8,0,2] - a[5,9,0,2]
   - a[5,5,2,2] - a[5,6,2,2] - a[5,7,2,2] - a[5,8,2,2] - a[5,9,2,2]
   + b[2,5,1] - x[5,1] >= -1
 R3241: - a[5,5,0,2] - a[5,6,0,2] - a[5,7,0,2] - a[5,8,0,2] - a[5,9,0,2]
   - a[5,5,2,2] - a[5,6,2,2] - a[5,7,2,2] - a[5,8,2,2] - a[5,9,2,2]
   + b[2,5,2] - x[5,2] >= -1
 R3242: - a[6,6,0,2] - a[6,7,0,2] - a[6,8,0,2] - a[6,9,0,2] - a[6,6,2,2]
   - a[6,7,2,2] - a[6,8,2,2] - a[6,9,2,2] + b[2,6,0] - x[6,0] >= -1
 R3243: - a[6,6,0,2] - a[6,7,0,2] - a[6,8,0,2] - a[6,9,0,2] - a[6,6,2,2]
   - a[6,7,2,2] - a[6,8,2,2] - a[6,9,2,2] + b[2,6,1] - x[6,1] >= -1
 R3244: - a[6,6,0,2] - a[6,7,0,2] - a[6,8,0,2] - a[6,9,0,2] - a[6,6,2,2]
   - a[6,7,2,2] - a[6,8,2,2] - a[6,9,2,2] + b[2,6,2] - x[6,2] >= -1
 R3245: - a[7,7,0,2] - a[7,8,0,2] - a[7,9,0,2] - a[7,7,2,2] - a[7,8,2,2]
   - a[7,9,2,2] + b[2,7,0] - x[7,0] >= -1
 R3246: - a[7,7,0,2] - a[7,8,0,2] - a[7,9,0,2] - a[7,7,2,2] - a[7,8,2,2]
   - a[7,9,2,2] + b[2,7,1] - x[7,1] >= -1
 R3247: - a[7,7,0,2] - a[7,8,0,2] - a[7,9,0,2] - a[7,7,2,2] - a[7,8,2,2]
   - a[7,9,2,2] + b[2,7,2] - x[7,2] >= -1
 R3248: - a[8,8,0,2] - a[8,9,0,2] - a[8,8,2,2] - a[8,9,2,2] + b[2,8,0]
   - x[8,0] >= -1
 R3249: - a[8,8,0,2] - a[8,9,0,2] - a[8,8,2,2] - a[8,9,2,2] + b[2,8,1]
   - x[8,1] >= -1
 R3250: - a[8,8,0,2] - a[8,9,0,2] - a[8,8,2,2] - a[8,9,2,2] + b[2,8,2]
   - x[8,2] >= -1
 R3251: - a[9,9,0,2] - a[9,9,2,2] + b[2,9,0] - x[9,0] >= -1
 R3252: - a[9,9,0,2] - a[9,9,2,2] + b[2,9,1] - x[9,1] >= -1
 R3253: - a[9,9,0,2] - a[9,9,2,2] + b[2,9,2] - x[9,2] >= -1
 R3254: - a[0,0,0,6] - a[0,1,0,6] - a[0,2,0,6] - a[0,3,0,6] - a[0,4,0,6]
   - a[0,5,0,6] - a[0,6,0,6] - a[0,7,0,6] - a[0,8,0,6] - a[0,9,0,6]
   - a[0,0,2,6] - a[0,1,2,6] - a[0,2,2,6] - a[0,3,2,6] - a[0,4,2,6]
   - a[0,5,2,6] - a[0,6,2,6] - a[0,7,2,6] - a[0,8,2,6] - a[0,9,2,6]
   - a[0,0,6,6] - a[0,1,6,6] - a[0,2,6,6] - a[0,3,6,6] - a[0,4,6,6]
   - a[0,5,6,6] - a[0,6,6,6] - a[0,7,6,6] - a[0,8,6,6] - a[0,9,6,6]
   + b[6,0,0] - x[0,0] >= -1
 R3255: - a[0,0,0,6] - a[0,1,0,6] - a[0,2,0,6] - a[0,3,0,6] - a[0,4,0,6]
   - a[0,5,0,6] - a[0,6,0,6] - a[0,7,0,6] - a[0,8,0,6] - a[0,9,0,6]
   - a[0,0,2,6] - a[0,1,2,6] - a[0,2,2,6] - a[0,3,2,6] - a[0,4,2,6]
   - a[0,5,2,6] - a[0,6,2,6] - a[0,7,2,6] - a[0,8,2,6] - a[0,9,2,6]
   - a[0,0,6,6] - a[0,1,6,6] - a[0,2,6,6] - a[0,3,6,6] - a[0,4,6,6]
   - a[0,5,6,6] - a[0,6,6,6] - a[0,7,6,6] - a[0,8,6,6] - a[0,9,6,6]
   + b[6,0,1] - x[0,1] >= -1
 R3256: - a[0,0,0,6] - a[0,1,0,6] - a[0,2,0,6] - a[0,3,0,6] - a[0,4,0,6]
   - a[0,5,0,6] - a[0,6,0,6] - a[0,7,0,6] - a[0,8,0,6] - a[0,9,0,6]
   - a[0,0,2,6] - a[0,1,2,6] - a[0,2,2,6] - a[0,3,2,6] - a[0,4,2,6]
   - a[0,5,2,6] - a[0,6,2,6] - a[0,7,2,6] - a[0,8,2,6] - a[0,9,2,6]
   - a[0,0,6,6] - a[0,1,6,6] - a[0,2,6,6] - a[0,3,6,6] - a[0,4,6,6]
   - a[0,5,6,6] - a[0,6,6,6] - a[0,7,6,6] - a[0,8,6,6] - a[0,9,6,6]
   + b[6,0,2] - x[0,2] >= -1
 R3257: - a[1,1,0,6] - a[1,2,0,6] - a[1,3,0,6] - a[1,4,0,6] - a[1,5,0,6]
   - a[1,6,0,6] - a[1,7,0,6] - a[1,8,0,6] - a[1,9,0,6] - a[1,1,2,6]
   - a[1,2,2,6] - a[1,3,2,6] - a[1,4,2,6] - a[1,5,2,6] - a[1,6,2,6]
   - a[1,7,2,6] - a[1,8,2,6] - a[1,9,2,6] - a[1,1,6,6] - a[1,2,6,6]
   - a[1,3,6,6] - a[1,4,6,6] - a[1,5,6,6] - a[1,6,6,6] - a[1,7,6,6]
   - a[1,8,6,6] - a[1,9,6,6] + b[6,1,0] - x[1,0] >= -1
 R3258: - a[1,1,0,6] - a[1,2,0,6] - a[1,3,0,6] - a[1,4,0,6] - a[1,5,0,6]
   - a[1,6,0,6] - a[1,7,0,6] - a[1,8,0,6] - a[1,9,0,6] - a[1,1,2,6]
   - a[1,2,2,6] - a[1,3,2,6] - a[1,4,2,6] - a[1,5,2,6] - a[1,6,2,6]
   - a[1,7,2,6] - a[1,8,2,6] - a[1,9,2,6] - a[1,1,6,6] - a[1,2,6,6]
   - a[1,3,6,6] - a[1,4,6,6] - a[1,5,6,6] - a[1,6,6,6] - a[1,7,6,6]
   - a[1,8,6,6] - a[1,9,6,6] + b[6,1,1] - x[1,1] >= -1
 R3259: - a[1,1,0,6] - a[1,2,0,6] - a[1,3,0,6] - a[1,4,0,6] - a[1,5,0,6]
   - a[1,6,0,6] - a[1,7,0,6] - a[1,8,0,6] - a[1,9,0,6] - a[1,1,2,6]
   - a[1,2,2,6] - a[1,3,2,6] - a[1,4,2,6] - a[1,5,2,6] - a[1,6,2,6]
   - a[1,7,2,6] - a[1,8,2,6] - a[1,9,2,6] - a[1,1,6,6] - a[1,2,6,6]
   - a[1,3,6,6] - a[1,4,6,6] - a[1,5,6,6] - a[1,6,6,6] - a[1,7,6,6]
   - a[1,8,6,6] - a[1,9,6,6] + b[6,1,2] - x[1,2] >= -1
 R3260: - a[2,2,0,6] - a[2,3,0,6] - a[2,4,0,6] - a[2,5,0,6] - a[2,6,0,6]
   - a[2,7,0,6] - a[2,8,0,6] - a[2,9,0,6] - a[2,2,2,6] - a[2,3,2,6]
   - a[2,4,2,6] - a[2,5,2,6] - a[2,6,2,6] - a[2,7,2,6] - a[2,8,2,6]
   - a[2,9,2,6] - a[2,2,6,6] - a[2,3,6,6] - a[2,4,6,6] - a[2,5,6,6]
   - a[2,6,6,6] - a[2,7,6,6] - a[2,8,6,6] - a[2,9,6,6] + b[6,2,0] - x[2,0]
   >= -1
 R3261: - a[2,2,0,6] - a[2,3,0,6] - a[2,4,0,6] - a[2,5,0,6] - a[2,6,0,6]
   - a[2,7,0,6] - a[2,8,0,6] - a[2,9,0,6] - a[2,2,2,6] - a[2,3,2,6]
   - a[2,4,2,6] - a[2,5,2,6] - a[2,6,2,6] - a[2,7,2,6] - a[2,8,2,6]
   - a[2,9,2,6] - a[2,2,6,6] - a[2,3,6,6] - a[2,4,6,6] - a[2,5,6,6]
   - a[2,6,6,6] - a[2,7,6,6] - a[2,8,6,6] - a[2,9,6,6] + b[6,2,1] - x[2,1]
   >= -1
 R3262: - a[2,2,0,6] - a[2,3,0,6] - a[2,4,0,6] - a[2,5,0,6] - a[2,6,0,6]
   - a[2,7,0,6] - a[2,8,0,6] - a[2,9,0,6] - a[2,2,2,6] - a[2,3,2,6]
   - a[2,4,2,6] - a[2,5,2,6] - a[2,6,2,6] - a[2,7,2,6] - a[2,8,2,6]
   - a[2,9,2,6] - a[2,2,6,6] - a[2,3,6,6] - a[2,4,6,6] - a[2,5,6,6]
   - a[2,6,6,6] - a[2,7,6,6] - a[2,8,6,6] - a[2,9,6,6] + b[6,2,2] - x[2,2]
   >= -1
 R3263: - a[3,3,0,6] - a[3,4,0,6] - a[3,5,0,6] - a[3,6,0,6] - a[3,7,0,6]
   - a[3,8,0,6] - a[3,9,0,6] - a[3,3,2,6] - a[3,4,2,6] - a[3,5,2,6]
   - a[3,6,2,6] - a[3,7,2,6] - a[3,8,2,6] - a[3,9,2,6] - a[3,3,6,6]
   - a[3,4,6,6] - a[3,5,6,6] - a[3,6,6,6] - a[3,7,6,6] - a[3,8,6,6]
   - a[3,9,6,6] + b[6,3,0] - x[3,0] >= -1
 R3264: - a[3,3,0,6] - a[3,4,0,6] - a[3,5,0,6] - a[3,6,0,6] - a[3,7,0,6]
   - a[3,8,0,6] - a[3,9,0,6] - a[3,3,2,6] - a[3,4,2,6] - a[3,5,2,6]
   - a[3,6,2,6] - a[3,7,2,6] - a[3,8,2,6] - a[3,9,2,6] - a[3,3,6,6]
   - a[3,4,6,6] - a[3,5,6,6] - a[3,6,6,6] - a[3,7,6,6] - a[3,8,6,6]
   - a[3,9,6,6] + b[6,3,1] - x[3,1] >= -1
 R3265: - a[3,3,0,6] - a[3,4,0,6] - a[3,5,0,6] - a[3,6,0,6] - a[3,7,0,6]
   - a[3,8,0,6] - a[3,9,0,6] - a[3,3,2,6] - a[3,4,2,6] - a[3,5,2,6]
   - a[3,6,2,6] - a[3,7,2,6] - a[3,8,2,6] - a[3,9,2,6] - a[3,3,6,6]
   - a[3,4,6,6] - a[3,5,6,6] - a[3,6,6,6] - a[3,7,6,6] - a[3,8,6,6]
   - a[3,9,6,6] + b[6,3,2] - x[3,2] >= -1
 R3266: - a[4,4,0,6] - a[4,5,0,6] - a[4,6,0,6] - a[4,7,0,6] - a[4,8,0,6]
   - a[4,9,0,6] - a[4,4,2,6] - a[4,5,2,6] - a[4,6,2,6] - a[4,7,2,6]
   - a[4,8,2,6] - a[4,9,2,6] - a[4,4,6,6] - a[4,5,6,6] - a[4,6,6,6]
   - a[4,7,6,6] - a[4,8,6,6] - a[4,9,6,6] + b[6,4,0] - x[4,0] >= -1
 R3267: - a[4,4,0,6] - a[4,5,0,6] - a[4,6,0,6] - a[4,7,0,6] - a[4,8,0,6]
   - a[4,9,0,6] - a[4,4,2,6] - a[4,5,2,6] - a[4,6,2,6] - a[4,7,2,6]
   - a[4,8,2,6] - a[4,9,2,6] - a[4,4,6,6] - a[4,5,6,6] - a[4,6,6,6]
   - a[4,7,6,6] - a[4,8,6,6] - a[4,9,6,6] + b[6,4,1] - x[4,1] >= -1
 R3268: - a[4,4,0,6] - a[4,5,0,6] - a[4,6,0,6] - a[4,7,0,6] - a[4,8,0,6]
   - a[4,9,0,6] - a[4,4,2,6] - a[4,5,2,6] - a[4,6,2,6] - a[4,7,2,6]
   - a[4,8,2,6] - a[4,9,2,6] - a[4,4,6,6] - a[4,5,6,6] - a[4,6,6,6]
   - a[4,7,6,6] - a[4,8,6,6] - a[4,9,6,6] + b[6,4,2] - x[4,2] >= -1
 R3269: - a[5,5,0,6] - a[5,6,0,6] - a[5,7,0,6] - a[5,8,0,6] - a[5,9,0,6]
   - a[5,5,2,6] - a[5,6,2,6] - a[5,7,2,6] - a[5,8,2,6] - a[5,9,2,6]
   - a[5,5,6,6] - a[5,6,6,6] - a[5,7,6,6] - a[5,8,6,6] - a[5,9,6,6]
   + b[6,5,0] - x[5,0] >= -1
 R3270: - a[5,5,0,6] - a[5,6,0,6] - a[5,7,0,6] - a[5,8,0,6] - a[5,9,0,6]
   - a[5,5,2,6] - a[5,6,2,6] - a[5,7,2,6] - a[5,8,2,6] - a[5,9,2,6]
   - a[5,5,6,6] - a[5,6,6,6] - a[5,7,6,6] - a[5,8,6,6] - a[5,9,6,6]
   + b[6,5,1] - x[5,1] >= -1
 R3271: - a[5,5,0,6] - a[5,6,0,6] - a[5,7,0,6] - a[5,8,0,6] - a[5,9,0,6]
   - a[5,5,2,6] - a[5,6,2,6] - a[5,7,2,6] - a[5,8,2,6] - a[5,9,2,6]
   - a[5,5,6,6] - a[5,6,6,6] - a[5,7,6,6] - a[5,8,6,6] - a[5,9,6,6]
   + b[6,5,2] - x[5,2] >= -1
 R3272: - a[6,6,0,6] - a[6,7,0,6] - a[6,8,0,6] - a[6,9,0,6] - a[6,6,2,6]
   - a[6,7,2,6] - a[6,8,2,6] - a[6,9,2,6] - a[6,6,6,6] - a[6,7,6,6]
   - a[6,8,6,6] - a[6,9,6,6] + b[6,6,0] - x[6,0] >= -1
 R3273: - a[6,6,0,6] - a[6,7,0,6] - a[6,8,0,6] - a[6,9,0,6] - a[6,6,2,6]
   - a[6,7,2,6] - a[6,8,2,6] - a[6,9,2,6] - a[6,6,6,6] - a[6,7,6,6]
   - a[6,8,6,6] - a[6,9,6,6] + b[6,6,1] - x[6,1] >= -1
 R3274: - a[6,6,0,6] - a[6,7,0,6] - a[6,8,0,6] - a[6,9,0,6] - a[6,6,2,6]
   - a[6,7,2,6] - a[6,8,2,6] - a[6,9,2,6] - a[6,6,6,6] - a[6,7,6,6]
   - a[6,8,6,6] - a[6,9,6,6] + b[6,6,2] - x[6,2] >= -1
 R3275: - a[7,7,0,6] - a[7,8,0,6] - a[7,9,0,6] - a[7,7,2,6] - a[7,8,2,6]
   - a[7,9,2,6] - a[7,7,6,6] - a[7,8,6,6] - a[7,9,6,6] + b[6,7,0] - x[7,0]
   >= -1
 R3276: - a[7,7,0,6] - a[7,8,0,6] - a[7,9,0,6] - a[7,7,2,6] - a[7,8,2,6]
   - a[7,9,2,6] - a[7,7,6,6] - a[7,8,6,6] - a[7,9,6,6] + b[6,7,1] - x[7,1]
   >= -1
 R3277: - a[7,7,0,6] - a[7,8,0,6] - a[7,9,0,6] - a[7,7,2,6] - a[7,8,2,6]
   - a[7,9,2,6] - a[7,7,6,6] - a[7,8,6,6] - a[7,9,6,6] + b[6,7,2] - x[7,2]
   >= -1
 R3278: - a[8,8,0,6] - a[8,9,0,6] - a[8,8,2,6] - a[8,9,2,6] - a[8,8,6,6]
   - a[8,9,6,6] + b[6,8,0] - x[8,0] >= -1
 R3279: - a[8,8,0,6] - a[8,9,0,6] - a[8,8,2,6] - a[8,9,2,6] - a[8,8,6,6]
   - a[8,9,6,6] + b[6,8,1] - x[8,1] >= -1
 R3280: - a[8,8,0,6] - a[8,9,0,6] - a[8,8,2,6] - a[8,9,2,6] - a[8,8,6,6]
   - a[8,9,6,6] + b[6,8,2] - x[8,2] >= -1
 R3281: - a[9,9,0,6] - a[9,9,2,6] - a[9,9,6,6] + b[6,9,0] - x[9,0] >= -1
 R3282: - a[9,9,0,6] - a[9,9,2,6] - a[9,9,6,6] + b[6,9,1] - x[9,1] >= -1
 R3283: - a[9,9,0,6] - a[9,9,2,6] - a[9,9,6,6] + b[6,9,2] - x[9,2] >= -1
 R3284: - a[0,0,0,8] - a[0,1,0,8] - a[0,2,0,8] - a[0,3,0,8] - a[0,4,0,8]
   - a[0,5,0,8] - a[0,6,0,8] - a[0,7,0,8] - a[0,8,0,8] - a[0,9,0,8]
   - a[0,0,2,8] - a[0,1,2,8] - a[0,2,2,8] - a[0,3,2,8] - a[0,4,2,8]
   - a[0,5,2,8] - a[0,6,2,8] - a[0,7,2,8] - a[0,8,2,8] - a[0,9,2,8]
   - a[0,0,6,8] - a[0,1,6,8] - a[0,2,6,8] - a[0,3,6,8] - a[0,4,6,8]
   - a[0,5,6,8] - a[0,6,6,8] - a[0,7,6,8] - a[0,8,6,8] - a[0,9,6,8]
   - a[0,0,8,8] - a[0,1,8,8] - a[0,2,8,8] - a[0,3,8,8] - a[0,4,8,8]
   - a[0,5,8,8] - a[0,6,8,8] - a[0,7,8,8] - a[0,8,8,8] - a[0,9,8,8]
   + b[8,0,0] - x[0,0] >= -1
 R3285: - a[0,0,0,8] - a[0,1,0,8] - a[0,2,0,8] - a[0,3,0,8] - a[0,4,0,8]
   - a[0,5,0,8] - a[0,6,0,8] - a[0,7,0,8] - a[0,8,0,8] - a[0,9,0,8]
   - a[0,0,2,8] - a[0,1,2,8] - a[0,2,2,8] - a[0,3,2,8] - a[0,4,2,8]
   - a[0,5,2,8] - a[0,6,2,8] - a[0,7,2,8] - a[0,8,2,8] - a[0,9,2,8]
   - a[0,0,6,8] - a[0,1,6,8] - a[0,2,6,8] - a[0,3,6,8] - a[0,4,6,8]
   - a[0,5,6,8] - a[0,6,6,8] - a[0,7,6,8] - a[0,8,6,8] - a[0,9,6,8]
   - a[0,0,8,8] - a[0,1,8,8] - a[0,2,8,8] - a[0,3,8,8] - a[0,4,8,8]
   - a[0,5,8,8] - a[0,6,8,8] - a[0,7,8,8] - a[0,8,8,8] - a[0,9,8,8]
   + b[8,0,1] - x[0,1] >= -1
 R3286: - a[0,0,0,8] - a[0,1,0,8] - a[0,2,0,8] - a[0,3,0,8] - a[0,4,0,8]
   - a[0,5,0,8] - a[0,6,0,8] - a[0,7,0,8] - a[0,8,0,8] - a[0,9,0,8]
   - a[0,0,2,8] - a[0,1,2,8] - a[0,2,2,8] - a[0,3,2,8] - a[0,4,2,8]
   - a[0,5,2,8] - a[0,6,2,8] - a[0,7,2,8] - a[0,8,2,8] - a[0,9,2,8]
   - a[0,0,6,8] - a[0,1,6,8] - a[0,2,6,8] - a[0,3,6,8] - a[0,4,6,8]
   - a[0,5,6,8] - a[0,6,6,8] - a[0,7,6,8] - a[0,8,6,8] - a[0,9,6,8]
   - a[0,0,8,8] - a[0,1,8,8] - a[0,2,8,8] - a[0,3,8,8] - a[0,4,8,8]
   - a[0,5,8,8] - a[0,6,8,8] - a[0,7,8,8] - a[0,8,8,8] - a[0,9,8,8]
   + b[8,0,2] - x[0,2] >= -1
 R3287: - a[1,1,0,8] - a[1,2,0,8] - a[1,3,0,8] - a[1,4,0,8] - a[1,5,0,8]
   - a[1,6,0,8] - a[1,7,0,8] - a[1,8,0,8] - a[1,9,0,8] - a[1,1,2,8]
   - a[1,2,2,8] - a[1,3,2,8] - a[1,4,2,8] - a[1,5,2,8] - a[1,6,2,8]
   - a[1,7,2,8] - a[1,8,2,8] - a[1,9,2,8] - a[1,1,6,8] - a[1,2,6,8]
   - a[1,3,6,8] - a[1,4,6,8] - a[1,5,6,8] - a[1,6,6,8] - a[1,7,6,8]
   - a[1,8,6,8] - a[1,9,6,8] - a[1,1,8,8] - a[1,2,8,8] - a[1,3,8,8]
   - a[1,4,8,8] - a[1,5,8,8] - a[1,6,8,8] - a[1,7,8,8] - a[1,8,8,8]
   - a[1,9,8,8] + b[8,1,0] - x[1,0] >= -1
 R3288: - a[1,1,0,8] - a[1,2,0,8] - a[1,3,0,8] - a[1,4,0,8] - a[1,5,0,8]
   - a[1,6,0,8] - a[1,7,0,8] - a[1,8,0,8] - a[1,9,0,8] - a[1,1,2,8]
   - a[1,2,2,8] - a[1,3,2,8] - a[1,4,2,8] - a[1,5,2,8] - a[1,6,2,8]
   - a[1,7,2,8] - a[1,8,2,8] - a[1,9,2,8] - a[1,1,6,8] - a[1,2,6,8]
   - a[1,3,6,8] - a[1,4,6,8] - a[1,5,6,8] - a[1,6,6,8] - a[1,7,6,8]
   - a[1,8,6,8] - a[1,9,6,8] - a[1,1,8,8] - a[1,2,8,8] - a[1,3,8,8]
   - a[1,4,8,8] - a[1,5,8,8] - a[1,6,8,8] - a[1,7,8,8] - a[1,8,8,8]
   - a[1,9,8,8] + b[8,1,1] - x[1,1] >= -1
 R3289: - a[1,1,0,8] - a[1,2,0,8] - a[1,3,0,8] - a[1,4,0,8] - a[1,5,0,8]
   - a[1,6,0,8] - a[1,7,0,8] - a[1,8,0,8] - a[1,9,0,8] - a[1,1,2,8]
   - a[1,2,2,8] - a[1,3,2,8] - a[1,4,2,8] - a[1,5,2,8] - a[1,6,2,8]
   - a[1,7,2,8] - a[1,8,2,8] - a[1,9,2,8] - a[1,1,6,8] - a[1,2,6,8]
   - a[1,3,6,8] - a[1,4,6,8] - a[1,5,6,8] - a[1,6,6,8] - a[1,7,6,8]
   - a[1,8,6,8] - a[1,9,6,8] - a[1,1,8,8] - a[1,2,8,8] - a[1,3,8,8]
   - a[1,4,8,8] - a[1,5,8,8] - a[1,6,8,8] - a[1,7,8,8] - a[1,8,8,8]
   - a[1,9,8,8] + b[8,1,2] - x[1,2] >= -1
 R3290: - a[2,2,0,8] - a[2,3,0,8] - a[2,4,0,8] - a[2,5,0,8] - a[2,6,0,8]
   - a[2,7,0,8] - a[2,8,0,8] - a[2,9,0,8] - a[2,2,2,8] - a[2,3,2,8]
   - a[2,4,2,8] - a[2,5,2,8] - a[2,6,2,8] - a[2,7,2,8] - a[2,8,2,8]
   - a[2,9,2,8] - a[2,2,6,8] - a[2,3,6,8] - a[2,4,6,8] - a[2,5,6,8]
   - a[2,6,6,8] - a[2,7,6,8] - a[2,8,6,8] - a[2,9,6,8] - a[2,2,8,8]
   - a[2,3,8,8] - a[2,4,8,8] - a[2,5,8,8] - a[2,6,8,8] - a[2,7,8,8]
   - a[2,8,8,8] - a[2,9,8,8] + b[8,2,0] - x[2,0] >= -1
 R3291: - a[2,2,0,8] - a[2,3,0,8] - a[2,4,0,8] - a[2,5,0,8] - a[2,6,0,8]
   - a[2,7,0,8] - a[2,8,0,8] - a[2,9,0,8] - a[2,2,2,8] - a[2,3,2,8]
   - a[2,4,2,8] - a[2,5,2,8] - a[2,6,2,8] - a[2,7,2,8] - a[2,8,2,8]
   - a[2,9,2,8] - a[2,2,6,8] - a[2,3,6,8] - a[2,4,6,8] - a[2,5,6,8]
   - a[2,6,6,8] - a[2,7,6,8] - a[2,8,6,8] - a[2,9,6,8] - a[2,2,8,8]
   - a[2,3,8,8] - a[2,4,8,8] - a[2,5,8,8] - a[2,6,8,8] - a[2,7,8,8]
   - a[2,8,8,8] - a[2,9,8,8] + b[8,2,1] - x[2,1] >= -1
 R3292: - a[2,2,0,8] - a[2,3,0,8] - a[2,4,0,8] - a[2,5,0,8] - a[2,6,0,8]
   - a[2,7,0,8] - a[2,8,0,8] - a[2,9,0,8] - a[2,2,2,8] - a[2,3,2,8]
   - a[2,4,2,8] - a[2,5,2,8] - a[2,6,2,8] - a[2,7,2,8] - a[2,8,2,8]
   - a[2,9,2,8] - a[2,2,6,8] - a[2,3,6,8] - a[2,4,6,8] - a[2,5,6,8]
   - a[2,6,6,8] - a[2,7,6,8] - a[2,8,6,8] - a[2,9,6,8] - a[2,2,8,8]
   - a[2,3,8,8] - a[2,4,8,8] - a[2,5,8,8] - a[2,6,8,8] - a[2,7,8,8]
   - a[2,8,8,8] - a[2,9,8,8] + b[8,2,2] - x[2,2] >= -1
 R3293: - a[3,3,0,8] - a[3,4,0,8] - a[3,5,0,8] - a[3,6,0,8] - a[3,7,0,8]
   - a[3,8,0,8] - a[3,9,0,8] - a[3,3,2,8] - a[3,4,2,8] - a[3,5,2,8]
   - a[3,6,2,8] - a[3,7,2,8] - a[3,8,2,8] - a[3,9,2,8] - a[3,3,6,8]
   - a[3,4,6,8] - a[3,5,6,8] - a[3,6,6,8] - a[3,7,6,8] - a[3,8,6,8]
   - a[3,9,6,8] - a[3,3,8,8] - a[3,4,8,8] - a[3,5,8,8] - a[3,6,8,8]
   - a[3,7,8,8] - a[3,8,8,8] - a[3,9,8,8] + b[8,3,0] - x[3,0] >= -1
 R3294: - a[3,3,0,8] - a[3,4,0,8] - a[3,5,0,8] - a[3,6,0,8] - a[3,7,0,8]
   - a[3,8,0,8] - a[3,9,0,8] - a[3,3,2,8] - a[3,4,2,8] - a[3,5,2,8]
   - a[3,6,2,8] - a[3,7,2,8] - a[3,8,2,8] - a[3,9,2,8] - a[3,3,6,8]
   - a[3,4,6,8] - a[3,5,6,8] - a[3,6,6,8] - a[3,7,6,8] - a[3,8,6,8]
   - a[3,9,6,8] - a[3,3,8,8] - a[3,4,8,8] - a[3,5,8,8] - a[3,6,8,8]
   - a[3,7,8,8] - a[3,8,8,8] - a[3,9,8,8] + b[8,3,1] - x[3,1] >= -1
 R3295: - a[3,3,0,8] - a[3,4,0,8] - a[3,5,0,8] - a[3,6,0,8] - a[3,7,0,8]
   - a[3,8,0,8] - a[3,9,0,8] - a[3,3,2,8] - a[3,4,2,8] - a[3,5,2,8]
   - a[3,6,2,8] - a[3,7,2,8] - a[3,8,2,8] - a[3,9,2,8] - a[3,3,6,8]
   - a[3,4,6,8] - a[3,5,6,8] - a[3,6,6,8] - a[3,7,6,8] - a[3,8,6,8]
   - a[3,9,6,8] - a[3,3,8,8] - a[3,4,8,8] - a[3,5,8,8] - a[3,6,8,8]
   - a[3,7,8,8] - a[3,8,8,8] - a[3,9,8,8] + b[8,3,2] - x[3,2] >= -1
 R3296: - a[4,4,0,8] - a[4,5,0,8] - a[4,6,0,8] - a[4,7,0,8] - a[4,8,0,8]
   - a[4,9,0,8] - a[4,4,2,8] - a[4,5,2,8] - a[4,6,2,8] - a[4,7,2,8]
   - a[4,8,2,8] - a[4,9,2,8] - a[4,4,6,8] - a[4,5,6,8] - a[4,6,6,8]
   - a[4,7,6,8] - a[4,8,6,8] - a[4,9,6,8] - a[4,4,8,8] - a[4,5,8,8]
   - a[4,6,8,8] - a[4,7,8,8] - a[4,8,8,8] - a[4,9,8,8] + b[8,4,0] - x[4,0]
   >= -1
 R3297: - a[4,4,0,8] - a[4,5,0,8] - a[4,6,0,8] - a[4,7,0,8] - a[4,8,0,8]
   - a[4,9,0,8] - a[4,4,2,8] - a[4,5,2,8] - a[4,6,2,8] - a[4,7,2,8]
   - a[4,8,2,8] - a[4,9,2,8] - a[4,4,6,8] - a[4,5,6,8] - a[4,6,6,8]
   - a[4,7,6,8] - a[4,8,6,8] - a[4,9,6,8] - a[4,4,8,8] - a[4,5,8,8]
   - a[4,6,8,8] - a[4,7,8,8] - a[4,8,8,8] - a[4,9,8,8] + b[8,4,1] - x[4,1]
   >= -1
 R3298: - a[4,4,0,8] - a[4,5,0,8] - a[4,6,0,8] - a[4,7,0,8] - a[4,8,0,8]
   - a[4,9,0,8] - a[4,4,2,8] - a[4,5,2,8] - a[4,6,2,8] - a[4,7,2,8]
   - a[4,8,2,8] - a[4,9,2,8] - a[4,4,6,8] - a[4,5,6,8] - a[4,6,6,8]
   - a[4,7,6,8] - a[4,8,6,8] - a[4,9,6,8] - a[4,4,8,8] - a[4,5,8,8]
   - a[4,6,8,8] - a[4,7,8,8] - a[4,8,8,8] - a[4,9,8,8] + b[8,4,2] - x[4,2]
   >= -1
 R3299: - a[5,5,0,8] - a[5,6,0,8] - a[5,7,0,8] - a[5,8,0,8] - a[5,9,0,8]
   - a[5,5,2,8] - a[5,6,2,8] - a[5,7,2,8] - a[5,8,2,8] - a[5,9,2,8]
   - a[5,5,6,8] - a[5,6,6,8] - a[5,7,6,8] - a[5,8,6,8] - a[5,9,6,8]
   - a[5,5,8,8] - a[5,6,8,8] - a[5,7,8,8] - a[5,8,8,8] - a[5,9,8,8]
   + b[8,5,0] - x[5,0] >= -1
 R3300: - a[5,5,0,8] - a[5,6,0,8] - a[5,7,0,8] - a[5,8,0,8] - a[5,9,0,8]
   - a[5,5,2,8] - a[5,6,2,8] - a[5,7,2,8] - a[5,8,2,8] - a[5,9,2,8]
   - a[5,5,6,8] - a[5,6,6,8] - a[5,7,6,8] - a[5,8,6,8] - a[5,9,6,8]
   - a[5,5,8,8] - a[5,6,8,8] - a[5,7,8,8] - a[5,8,8,8] - a[5,9,8,8]
   + b[8,5,1] - x[5,1] >= -1
 R3301: - a[5,5,0,8] - a[5,6,0,8] - a[5,7,0,8] - a[5,8,0,8] - a[5,9,0,8]
   - a[5,5,2,8] - a[5,6,2,8] - a[5,7,2,8] - a[5,8,2,8] - a[5,9,2,8]
   - a[5,5,6,8] - a[5,6,6,8] - a[5,7,6,8] - a[5,8,6,8] - a[5,9,6,8]
   - a[5,5,8,8] - a[5,6,8,8] - a[5,7,8,8] - a[5,8,8,8] - a[5,9,8,8]
   + b[8,5,2] - x[5,2] >= -1
 R3302: - a[6,6,0,8] - a[6,7,0,8] - a[6,8,0,8] - a[6,9,0,8] - a[6,6,2,8]
   - a[6,7,2,8] - a[6,8,2,8] - a[6,9,2,8] - a[6,6,6,8] - a[6,7,6,8]
   - a[6,8,6,8] - a[6,9,6,8] - a[6,6,8,8] - a[6,7,8,8] - a[6,8,8,8]
   - a[6,9,8,8] + b[8,6,0] - x[6,0] >= -1
 R3303: - a[6,6,0,8] - a[6,7,0,8] - a[6,8,0,8] - a[6,9,0,8] - a[6,6,2,8]
   - a[6,7,2,8] - a[6,8,2,8] - a[6,9,2,8] - a[6,6,6,8] - a[6,7,6,8]
   - a[6,8,6,8] - a[6,9,6,8] - a[6,6,8,8] - a[6,7,8,8] - a[6,8,8,8]
   - a[6,9,8,8] + b[8,6,1] - x[6,1] >= -1
 R3304: - a[6,6,0,8] - a[6,7,0,8] - a[6,8,0,8] - a[6,9,0,8] - a[6,6,2,8]
   - a[6,7,2,8] - a[6,8,2,8] - a[6,9,2,8] - a[6,6,6,8] - a[6,7,6,8]
   - a[6,8,6,8] - a[6,9,6,8] - a[6,6,8,8] - a[6,7,8,8] - a[6,8,8,8]
   - a[6,9,8,8] + b[8,6,2] - x[6,2] >= -1
 R3305: - a[7,7,0,8] - a[7,8,0,8] - a[7,9,0,8] - a[7,7,2,8] - a[7,8,2,8]
   - a[7,9,2,8] - a[7,7,6,8] - a[7,8,6,8] - a[7,9,6,8] - a[7,7,8,8]
   - a[7,8,8,8] - a[7,9,8,8] + b[8,7,0] - x[7,0] >= -1
 R3306: - a[7,7,0,8] - a[7,8,0,8] - a[7,9,0,8] - a[7,7,2,8] - a[7,8,2,8]
   - a[7,9,2,8] - a[7,7,6,8] - a[7,8,6,8] - a[7,9,6,8] - a[7,7,8,8]
   - a[7,8,8,8] - a[7,9,8,8] + b[8,7,1] - x[7,1] >= -1
 R3307: - a[7,7,0,8] - a[7,8,0,8] - a[7,9,0,8] - a[7,7,2,8] - a[7,8,2,8]
   - a[7,9,2,8] - a[7,7,6,8] - a[7,8,6,8] - a[7,9,6,8] - a[7,7,8,8]
   - a[7,8,8,8] - a[7,9,8,8] + b[8,7,2] - x[7,2] >= -1
 R3308: - a[8,8,0,8] - a[8,9,0,8] - a[8,8,2,8] - a[8,9,2,8] - a[8,8,6,8]
   - a[8,9,6,8] - a[8,8,8,8] - a[8,9,8,8] + b[8,8,0] - x[8,0] >= -1
 R3309: - a[8,8,0,8] - a[8,9,0,8] - a[8,8,2,8] - a[8,9,2,8] - a[8,8,6,8]
   - a[8,9,6,8] - a[8,8,8,8] - a[8,9,8,8] + b[8,8,1] - x[8,1] >= -1
 R3310: - a[8,8,0,8] - a[8,9,0,8] - a[8,8,2,8] - a[8,9,2,8] - a[8,8,6,8]
   - a[8,9,6,8] - a[8,8,8,8] - a[8,9,8,8] + b[8,8,2] - x[8,2] >= -1
 R3311: - a[9,9,0,8] - a[9,9,2,8] - a[9,9,6,8] - a[9,9,8,8] + b[8,9,0]
   - x[9,0] >= -1
 R3312: - a[9,9,0,8] - a[9,9,2,8] - a[9,9,6,8] - a[9,9,8,8] + b[8,9,1]
   - x[9,1] >= -1
 R3313: - a[9,9,0,8] - a[9,9,2,8] - a[9,9,6,8] - a[9,9,8,8] + b[8,9,2]
   - x[9,2] >= -1
 R3314: - a[0,0,0,9] - a[0,1,0,9] - a[0,2,0,9] - a[0,3,0,9] - a[0,4,0,9]
   - a[0,5,0,9] - a[0,6,0,9] - a[0,7,0,9] - a[0,8,0,9] - a[0,9,0,9]
   - a[0,0,2,9] - a[0,1,2,9] - a[0,2,2,9] - a[0,3,2,9] - a[0,4,2,9]
   - a[0,5,2,9] - a[0,6,2,9] - a[0,7,2,9] - a[0,8,2,9] - a[0,9,2,9]
   - a[0,0,6,9] - a[0,1,6,9] - a[0,2,6,9] - a[0,3,6,9] - a[0,4,6,9]
   - a[0,5,6,9] - a[0,6,6,9] - a[0,7,6,9] - a[0,8,6,9] - a[0,9,6,9]
   - a[0,0,8,9] - a[0,1,8,9] - a[0,2,8,9] - a[0,3,8,9] - a[0,4,8,9]
   - a[0,5,8,9] - a[0,6,8,9] - a[0,7,8,9] - a[0,8,8,9] - a[0,9,8,9]
   - a[0,0,9,9] - a[0,1,9,9] - a[0,2,9,9] - a[0,3,9,9] - a[0,4,9,9]
   - a[0,5,9,9] - a[0,6,9,9] - a[0,7,9,9] - a[0,8,9,9] - a[0,9,9,9]
   + b[9,0,0] - x[0,0] >= -1
 R3315: - a[0,0,0,9] - a[0,1,0,9] - a[0,2,0,9] - a[0,3,0,9] - a[0,4,0,9]
   - a[0,5,0,9] - a[0,6,0,9] - a[0,7,0,9] - a[0,8,0,9] - a[0,9,0,9]
   - a[0,0,2,9] - a[0,1,2,9] - a[0,2,2,9] - a[0,3,2,9] - a[0,4,2,9]
   - a[0,5,2,9] - a[0,6,2,9] - a[0,7,2,9] - a[0,8,2,9] - a[0,9,2,9]
   - a[0,0,6,9] - a[0,1,6,9] - a[0,2,6,9] - a[0,3,6,9] - a[0,4,6,9]
   - a[0,5,6,9] - a[0,6,6,9] - a[0,7,6,9] - a[0,8,6,9] - a[0,9,6,9]
   - a[0,0,8,9] - a[0,1,8,9] - a[0,2,8,9] - a[0,3,8,9] - a[0,4,8,9]
   - a[0,5,8,9] - a[0,6,8,9] - a[0,7,8,9] - a[0,8,8,9] - a[0,9,8,9]
   - a[0,0,9,9] - a[0,1,9,9] - a[0,2,9,9] - a[0,3,9,9] - a[0,4,9,9]
   - a[0,5,9,9] - a[0,6,9,9] - a[0,7,9,9] - a[0,8,9,9] - a[0,9,9,9]
   + b[9,0,1] - x[0,1] >= -1
 R3316: - a[0,0,0,9] - a[0,1,0,9] - a[0,2,0,9] - a[0,3,0,9] - a[0,4,0,9]
   - a[0,5,0,9] - a[0,6,0,9] - a[0,7,0,9] - a[0,8,0,9] - a[0,9,0,9]
   - a[0,0,2,9] - a[0,1,2,9] - a[0,2,2,9] - a[0,3,2,9] - a[0,4,2,9]
   - a[0,5,2,9] - a[0,6,2,9] - a[0,7,2,9] - a[0,8,2,9] - a[0,9,2,9]
   - a[0,0,6,9] - a[0,1,6,9] - a[0,2,6,9] - a[0,3,6,9] - a[0,4,6,9]
   - a[0,5,6,9] - a[0,6,6,9] - a[0,7,6,9] - a[0,8,6,9] - a[0,9,6,9]
   - a[0,0,8,9] - a[0,1,8,9] - a[0,2,8,9] - a[0,3,8,9] - a[0,4,8,9]
   - a[0,5,8,9] - a[0,6,8,9] - a[0,7,8,9] - a[0,8,8,9] - a[0,9,8,9]
   - a[0,0,9,9] - a[0,1,9,9] - a[0,2,9,9] - a[0,3,9,9] - a[0,4,9,9]
   - a[0,5,9,9] - a[0,6,9,9] - a[0,7,9,9] - a[0,8,9,9] - a[0,9,9,9]
   + b[9,0,2] - x[0,2] >= -1
 R3317: - a[1,1,0,9] - a[1,2,0,9] - a[1,3,0,9] - a[1,4,0,9] - a[1,5,0,9]
   - a[1,6,0,9] - a[1,7,0,9] - a[1,8,0,9] - a[1,9,0,9] - a[1,1,2,9]
   - a[1,2,2,9] - a[1,3,2,9] - a[1,4,2,9] - a[1,5,2,9] - a[1,6,2,9]
   - a[1,7,2,9] - a[1,8,2,9] - a[1,9,2,9] - a[1,1,6,9] - a[1,2,6,9]
   - a[1,3,6,9] - a[1,4,6,9] - a[1,5,6,9] - a[1,6,6,9] - a[1,7,6,9]
   - a[1,8,6,9] - a[1,9,6,9] - a[1,1,8,9] - a[1,2,8,9] - a[1,3,8,9]
   - a[1,4,8,9] - a[1,5,8,9] - a[1,6,8,9] - a[1,7,8,9] - a[1,8,8,9]
   - a[1,9,8,9] - a[1,1,9,9] - a[1,2,9,9] - a[1,3,9,9] - a[1,4,9,9]
   - a[1,5,9,9] - a[1,6,9,9] - a[1,7,9,9] - a[1,8,9,9] - a[1,9,9,9]
   + b[9,1,0] - x[1,0] >= -1
 R3318: - a[1,1,0,9] - a[1,2,0,9] - a[1,3,0,9] - a[1,4,0,9] - a[1,5,0,9]
   - a[1,6,0,9] - a[1,7,0,9] - a[1,8,0,9] - a[1,9,0,9] - a[1,1,2,9]
   - a[1,2,2,9] - a[1,3,2,9] - a[1,4,2,9] - a[1,5,2,9] - a[1,6,2,9]
   - a[1,7,2,9] - a[1,8,2,9] - a[1,9,2,9] - a[1,1,6,9] - a[1,2,6,9]
   - a[1,3,6,9] - a[1,4,6,9] - a[1,5,6,9] - a[1,6,6,9] - a[1,7,6,9]
   - a[1,8,6,9] - a[1,9,6,9] - a[1,1,8,9] - a[1,2,8,9] - a[1,3,8,9]
   - a[1,4,8,9] - a[1,5,8,9] - a[1,6,8,9] - a[1,7,8,9] - a[1,8,8,9]
   - a[1,9,8,9] - a[1,1,9,9] - a[1,2,9,9] - a[1,3,9,9] - a[1,4,9,9]
   - a[1,5,9,9] - a[1,6,9,9] - a[1,7,9,9] - a[1,8,9,9] - a[1,9,9,9]
   + b[9,1,1] - x[1,1] >= -1
 R3319: - a[1,1,0,9] - a[1,2,0,9] - a[1,3,0,9] - a[1,4,0,9] - a[1,5,0,9]
   - a[1,6,0,9] - a[1,7,0,9] - a[1,8,0,9] - a[1,9,0,9] - a[1,1,2,9]
   - a[1,2,2,9] - a[1,3,2,9] - a[1,4,2,9] - a[1,5,2,9] - a[1,6,2,9]
   - a[1,7,2,9] - a[1,8,2,9] - a[1,9,2,9] - a[1,1,6,9] - a[1,2,6,9]
   - a[1,3,6,9] - a[1,4,6,9] - a[1,5,6,9] - a[1,6,6,9] - a[1,7,6,9]
   - a[1,8,6,9] - a[1,9,6,9] - a[1,1,8,9] - a[1,2,8,9] - a[1,3,8,9]
   - a[1,4,8,9] - a[1,5,8,9] - a[1,6,8,9] - a[1,7,8,9] - a[1,8,8,9]
   - a[1,9,8,9] - a[1,1,9,9] - a[1,2,9,9] - a[1,3,9,9] - a[1,4,9,9]
   - a[1,5,9,9] - a[1,6,9,9] - a[1,7,9,9] - a[1,8,9,9] - a[1,9,9,9]
   + b[9,1,2] - x[1,2] >= -1
 R3320: - a[2,2,0,9] - a[2,3,0,9] - a[2,4,0,9] - a[2,5,0,9] - a[2,6,0,9]
   - a[2,7,0,9] - a[2,8,0,9] - a[2,9,0,9] - a[2,2,2,9] - a[2,3,2,9]
   - a[2,4,2,9] - a[2,5,2,9] - a[2,6,2,9] - a[2,7,2,9] - a[2,8,2,9]
   - a[2,9,2,9] - a[2,2,6,9] - a[2,3,6,9] - a[2,4,6,9] - a[2,5,6,9]
   - a[2,6,6,9] - a[2,7,6,9] - a[2,8,6,9] - a[2,9,6,9] - a[2,2,8,9]
   - a[2,3,8,9] - a[2,4,8,9] - a[2,5,8,9] - a[2,6,8,9] - a[2,7,8,9]
   - a[2,8,8,9] - a[2,9,8,9] - a[2,2,9,9] - a[2,3,9,9] - a[2,4,9,9]
   - a[2,5,9,9] - a[2,6,9,9] - a[2,7,9,9] - a[2,8,9,9] - a[2,9,9,9]
   + b[9,2,0] - x[2,0] >= -1
 R3321: - a[2,2,0,9] - a[2,3,0,9] - a[2,4,0,9] - a[2,5,0,9] - a[2,6,0,9]
   - a[2,7,0,9] - a[2,8,0,9] - a[2,9,0,9] - a[2,2,2,9] - a[2,3,2,9]
   - a[2,4,2,9] - a[2,5,2,9] - a[2,6,2,9] - a[2,7,2,9] - a[2,8,2,9]
   - a[2,9,2,9] - a[2,2,6,9] - a[2,3,6,9] - a[2,4,6,9] - a[2,5,6,9]
   - a[2,6,6,9] - a[2,7,6,9] - a[2,8,6,9] - a[2,9,6,9] - a[2,2,8,9]
   - a[2,3,8,9] - a[2,4,8,9] - a[2,5,8,9] - a[2,6,8,9] - a[2,7,8,9]
   - a[2,8,8,9] - a[2,9,8,9] - a[2,2,9,9] - a[2,3,9,9] - a[2,4,9,9]
   - a[2,5,9,9] - a[2,6,9,9] - a[2,7,9,9] - a[2,8,9,9] - a[2,9,9,9]
   + b[9,2,1] - x[2,1] >= -1
 R3322: - a[2,2,0,9] - a[2,3,0,9] - a[2,4,0,9] - a[2,5,0,9] - a[2,6,0,9]
   - a[2,7,0,9] - a[2,8,0,9] - a[2,9,0,9] - a[2,2,2,9] - a[2,3,2,9]
   - a[2,4,2,9] - a[2,5,2,9] - a[2,6,2,9] - a[2,7,2,9] - a[2,8,2,9]
   - a[2,9,2,9] - a[2,2,6,9] - a[2,3,6,9] - a[2,4,6,9] - a[2,5,6,9]
   - a[2,6,6,9] - a[2,7,6,9] - a[2,8,6,9] - a[2,9,6,9] - a[2,2,8,9]
   - a[2,3,8,9] - a[2,4,8,9] - a[2,5,8,9] - a[2,6,8,9] - a[2,7,8,9]
   - a[2,8,8,9] - a[2,9,8,9] - a[2,2,9,9] - a[2,3,9,9] - a[2,4,9,9]
   - a[2,5,9,9] - a[2,6,9,9] - a[2,7,9,9] - a[2,8,9,9] - a[2,9,9,9]
   + b[9,2,2] - x[2,2] >= -1
 R3323: - a[3,3,0,9] - a[3,4,0,9] - a[3,5,0,9] - a[3,6,0,9] - a[3,7,0,9]
   - a[3,8,0,9] - a[3,9,0,9] - a[3,3,2,9] - a[3,4,2,9] - a[3,5,2,9]
   - a[3,6,2,9] - a[3,7,2,9] - a[3,8,2,9] - a[3,9,2,9] - a[3,3,6,9]
   - a[3,4,6,9] - a[3,5,6,9] - a[3,6,6,9] - a[3,7,6,9] - a[3,8,6,9]
   - a[3,9,6,9] - a[3,3,8,9] - a[3,4,8,9] - a[3,5,8,9] - a[3,6,8,9]
   - a[3,7,8,9] - a[3,8,8,9] - a[3,9,8,9] - a[3,3,9,9] - a[3,4,9,9]
   - a[3,5,9,9] - a[3,6,9,9] - a[3,7,9,9] - a[3,8,9,9] - a[3,9,9,9]
   + b[9,3,0] - x[3,0] >= -1
 R3324: - a[3,3,0,9] - a[3,4,0,9] - a[3,5,0,9] - a[3,6,0,9] - a[3,7,0,9]
   - a[3,8,0,9] - a[3,9,0,9] - a[3,3,2,9] - a[3,4,2,9] - a[3,5,2,9]
   - a[3,6,2,9] - a[3,7,2,9] - a[3,8,2,9] - a[3,9,2,9] - a[3,3,6,9]
   - a[3,4,6,9] - a[3,5,6,9] - a[3,6,6,9] - a[3,7,6,9] - a[3,8,6,9]
   - a[3,9,6,9] - a[3,3,8,9] - a[3,4,8,9] - a[3,5,8,9] - a[3,6,8,9]
   - a[3,7,8,9] - a[3,8,8,9] - a[3,9,8,9] - a[3,3,9,9] - a[3,4,9,9]
   - a[3,5,9,9] - a[3,6,9,9] - a[3,7,9,9] - a[3,8,9,9] - a[3,9,9,9]
   + b[9,3,1] - x[3,1] >= -1
 R3325: - a[3,3,0,9] - a[3,4,0,9] - a[3,5,0,9] - a[3,6,0,9] - a[3,7,0,9]
   - a[3,8,0,9] - a[3,9,0,9] - a[3,3,2,9] - a[3,4,2,9] - a[3,5,2,9]
   - a[3,6,2,9] - a[3,7,2,9] - a[3,8,2,9] - a[3,9,2,9] - a[3,3,6,9]
   - a[3,4,6,9] - a[3,5,6,9] - a[3,6,6,9] - a[3,7,6,9] - a[3,8,6,9]
   - a[3,9,6,9] - a[3,3,8,9] - a[3,4,8,9] - a[3,5,8,9] - a[3,6,8,9]
   - a[3,7,8,9] - a[3,8,8,9] - a[3,9,8,9] - a[3,3,9,9] - a[3,4,9,9]
   - a[3,5,9,9] - a[3,6,9,9] - a[3,7,9,9] - a[3,8,9,9] - a[3,9,9,9]
   + b[9,3,2] - x[3,2] >= -1
 R3326: - a[4,4,0,9] - a[4,5,0,9] - a[4,6,0,9] - a[4,7,0,9] - a[4,8,0,9]
   - a[4,9,0,9] - a[4,4,2,9] - a[4,5,2,9] - a[4,6,2,9] - a[4,7,2,9]
   - a[4,8,2,9] - a[4,9,2,9] - a[4,4,6,9] - a[4,5,6,9] - a[4,6,6,9]
   - a[4,7,6,9] - a[4,8,6,9] - a[4,9,6,9] - a[4,4,8,9] - a[4,5,8,9]
   - a[4,6,8,9] - a[4,7,8,9] - a[4,8,8,9] - a[4,9,8,9] - a[4,4,9,9]
   - a[4,5,9,9] - a[4,6,9,9] - a[4,7,9,9] - a[4,8,9,9] - a[4,9,9,9]
   + b[9,4,0] - x[4,0] >= -1
 R3327: - a[4,4,0,9] - a[4,5,0,9] - a[4,6,0,9] - a[4,7,0,9] - a[4,8,0,9]
   - a[4,9,0,9] - a[4,4,2,9] - a[4,5,2,9] - a[4,6,2,9] - a[4,7,2,9]
   - a[4,8,2,9] - a[4,9,2,9] - a[4,4,6,9] - a[4,5,6,9] - a[4,6,6,9]
   - a[4,7,6,9] - a[4,8,6,9] - a[4,9,6,9] - a[4,4,8,9] - a[4,5,8,9]
   - a[4,6,8,9] - a[4,7,8,9] - a[4,8,8,9] - a[4,9,8,9] - a[4,4,9,9]
   - a[4,5,9,9] - a[4,6,9,9] - a[4,7,9,9] - a[4,8,9,9] - a[4,9,9,9]
   + b[9,4,1] - x[4,1] >= -1
 R3328: - a[4,4,0,9] - a[4,5,0,9] - a[4,6,0,9] - a[4,7,0,9] - a[4,8,0,9]
   - a[4,9,0,9] - a[4,4,2,9] - a[4,5,2,9] - a[4,6,2,9] - a[4,7,2,9]
   - a[4,8,2,9] - a[4,9,2,9] - a[4,4,6,9] - a[4,5,6,9] - a[4,6,6,9]
   - a[4,7,6,9] - a[4,8,6,9] - a[4,9,6,9] - a[4,4,8,9] - a[4,5,8,9]
   - a[4,6,8,9] - a[4,7,8,9] - a[4,8,8,9] - a[4,9,8,9] - a[4,4,9,9]
   - a[4,5,9,9] - a[4,6,9,9] - a[4,7,9,9] - a[4,8,9,9] - a[4,9,9,9]
   + b[9,4,2] - x[4,2] >= -1
 R3329: - a[5,5,0,9] - a[5,6,0,9] - a[5,7,0,9] - a[5,8,0,9] - a[5,9,0,9]
   - a[5,5,2,9] - a[5,6,2,9] - a[5,7,2,9] - a[5,8,2,9] - a[5,9,2,9]
   - a[5,5,6,9] - a[5,6,6,9] - a[5,7,6,9] - a[5,8,6,9] - a[5,9,6,9]
   - a[5,5,8,9] - a[5,6,8,9] - a[5,7,8,9] - a[5,8,8,9] - a[5,9,8,9]
   - a[5,5,9,9] - a[5,6,9,9] - a[5,7,9,9] - a[5,8,9,9] - a[5,9,9,9]
   + b[9,5,0] - x[5,0] >= -1
 R3330: - a[5,5,0,9] - a[5,6,0,9] - a[5,7,0,9] - a[5,8,0,9] - a[5,9,0,9]
   - a[5,5,2,9] - a[5,6,2,9] - a[5,7,2,9] - a[5,8,2,9] - a[5,9,2,9]
   - a[5,5,6,9] - a[5,6,6,9] - a[5,7,6,9] - a[5,8,6,9] - a[5,9,6,9]
   - a[5,5,8,9] - a[5,6,8,9] - a[5,7,8,9] - a[5,8,8,9] - a[5,9,8,9]
   - a[5,5,9,9] - a[5,6,9,9] - a[5,7,9,9] - a[5,8,9,9] - a[5,9,9,9]
   + b[9,5,1] - x[5,1] >= -1
 R3331: - a[5,5,0,9] - a[5,6,0,9] - a[5,7,0,9] - a[5,8,0,9] - a[5,9,0,9]
   - a[5,5,2,9] - a[5,6,2,9] - a[5,7,2,9] - a[5,8,2,9] - a[5,9,2,9]
   - a[5,5,6,9] - a[5,6,6,9] - a[5,7,6,9] - a[5,8,6,9] - a[5,9,6,9]
   - a[5,5,8,9] - a[5,6,8,9] - a[5,7,8,9] - a[5,8,8,9] - a[5,9,8,9]
   - a[5,5,9,9] - a[5,6,9,9] - a[5,7,9,9] - a[5,8,9,9] - a[5,9,9,9]
   + b[9,5,2] - x[5,2] >= -1
 R3332: - a[6,6,0,9] - a[6,7,0,9] - a[6,8,0,9] - a[6,9,0,9] - a[6,6,2,9]
   - a[6,7,2,9] - a[6,8,2,9] - a[6,9,2,9] - a[6,6,6,9] - a[6,7,6,9]
   - a[6,8,6,9] - a[6,9,6,9] - a[6,6,8,9] - a[6,7,8,9] - a[6,8,8,9]
   - a[6,9,8,9] - a[6,6,9,9] - a[6,7,9,9] - a[6,8,9,9] - a[6,9,9,9]
   + b[9,6,0] - x[6,0] >= -1
 R3333: - a[6,6,0,9] - a[6,7,0,9] - a[6,8,0,9] - a[6,9,0,9] - a[6,6,2,9]
   - a[6,7,2,9] - a[6,8,2,9] - a[6,9,2,9] - a[6,6,6,9] - a[6,7,6,9]
   - a[6,8,6,9] - a[6,9,6,9] - a[6,6,8,9] - a[6,7,8,9] - a[6,8,8,9]
   - a[6,9,8,9] - a[6,6,9,9] - a[6,7,9,9] - a[6,8,9,9] - a[6,9,9,9]
   + b[9,6,1] - x[6,1] >= -1
 R3334: - a[6,6,0,9] - a[6,7,0,9] - a[6,8,0,9] - a[6,9,0,9] - a[6,6,2,9]
   - a[6,7,2,9] - a[6,8,2,9] - a[6,9,2,9] - a[6,6,6,9] - a[6,7,6,9]
   - a[6,8,6,9] - a[6,9,6,9] - a[6,6,8,9] - a[6,7,8,9] - a[6,8,8,9]
   - a[6,9,8,9] - a[6,6,9,9] - a[6,7,9,9] - a[6,8,9,9] - a[6,9,9,9]
   + b[9,6,2] - x[6,2] >= -1
 R3335: - a[7,7,0,9] - a[7,8,0,9] - a[7,9,0,9] - a[7,7,2,9] - a[7,8,2,9]
   - a[7,9,2,9] - a[7,7,6,9] - a[7,8,6,9] - a[7,9,6,9] - a[7,7,8,9]
   - a[7,8,8,9] - a[7,9,8,9] - a[7,7,9,9] - a[7,8,9,9] - a[7,9,9,9]
   + b[9,7,0] - x[7,0] >= -1
 R3336: - a[7,7,0,9] - a[7,8,0,9] - a[7,9,0,9] - a[7,7,2,9] - a[7,8,2,9]
   - a[7,9,2,9] - a[7,7,6,9] - a[7,8,6,9] - a[7,9,6,9] - a[7,7,8,9]
   - a[7,8,8,9] - a[7,9,8,9] - a[7,7,9,9] - a[7,8,9,9] - a[7,9,9,9]
   + b[9,7,1] - x[7,1] >= -1
 R3337: - a[7,7,0,9] - a[7,8,0,9] - a[7,9,0,9] - a[7,7,2,9] - a[7,8,2,9]
   - a[7,9,2,9] - a[7,7,6,9] - a[7,8,6,9] - a[7,9,6,9] - a[7,7,8,9]
   - a[7,8,8,9] - a[7,9,8,9] - a[7,7,9,9] - a[7,8,9,9] - a[7,9,9,9]
   + b[9,7,2] - x[7,2] >= -1
 R3338: - a[8,8,0,9] - a[8,9,0,9] - a[8,8,2,9] - a[8,9,2,9] - a[8,8,6,9]
   - a[8,9,6,9] - a[8,8,8,9] - a[8,9,8,9] - a[8,8,9,9] - a[8,9,9,9]
   + b[9,8,0] - x[8,0] >= -1
 R3339: - a[8,8,0,9] - a[8,9,0,9] - a[8,8,2,9] - a[8,9,2,9] - a[8,8,6,9]
   - a[8,9,6,9] - a[8,8,8,9] - a[8,9,8,9] - a[8,8,9,9] - a[8,9,9,9]
   + b[9,8,1] - x[8,1] >= -1
 R3340: - a[8,8,0,9] - a[8,9,0,9] - a[8,8,2,9] - a[8,9,2,9] - a[8,8,6,9]
   - a[8,9,6,9] - a[8,8,8,9] - a[8,9,8,9] - a[8,8,9,9] - a[8,9,9,9]
   + b[9,8,2] - x[8,2] >= -1
 R3341: - a[9,9,0,9] - a[9,9,2,9] - a[9,9,6,9] - a[9,9,8,9] - a[9,9,9,9]
   + b[9,9,0] - x[9,0] >= -1
 R3342: - a[9,9,0,9] - a[9,9,2,9] - a[9,9,6,9] - a[9,9,8,9] - a[9,9,9,9]
   + b[9,9,1] - x[9,1] >= -1
 R3343: - a[9,9,0,9] - a[9,9,2,9] - a[9,9,6,9] - a[9,9,8,9] - a[9,9,9,9]
   + b[9,9,2] - x[9,2] >= -1
 R3344: - a[0,0,1,1] - a[0,1,1,1] - a[0,2,1,1] - a[0,3,1,1] - a[0,4,1,1]
   - a[0,5,1,1] - a[0,6,1,1] - a[0,7,1,1] - a[0,8,1,1] - a[0,9,1,1]
   + b[1,0,0] - x[0,0] >= -1
 R3345: - a[0,0,1,1] - a[0,1,1,1] - a[0,2,1,1] - a[0,3,1,1] - a[0,4,1,1]
   - a[0,5,1,1] - a[0,6,1,1] - a[0,7,1,1] - a[0,8,1,1] - a[0,9,1,1]
   + b[1,0,1] - x[0,1] >= -1
 R3346: - a[0,0,1,1] - a[0,1,1,1] - a[0,2,1,1] - a[0,3,1,1] - a[0,4,1,1]
   - a[0,5,1,1] - a[0,6,1,1] - a[0,7,1,1] - a[0,8,1,1] - a[0,9,1,1]
   + b[1,0,2] - x[0,2] >= -1
 R3347: - a[1,1,1,1] - a[1,2,1,1] - a[1,3,1,1] - a[1,4,1,1] - a[1,5,1,1]
   - a[1,6,1,1] - a[1,7,1,1] - a[1,8,1,1] - a[1,9,1,1] + b[1,1,0] - x[1,0]
   >= -1
 R3348: - a[1,1,1,1] - a[1,2,1,1] - a[1,3,1,1] - a[1,4,1,1] - a[1,5,1,1]
   - a[1,6,1,1] - a[1,7,1,1] - a[1,8,1,1] - a[1,9,1,1] + b[1,1,1] - x[1,1]
   >= -1
 R3349: - a[1,1,1,1] - a[1,2,1,1] - a[1,3,1,1] - a[1,4,1,1] - a[1,5,1,1]
   - a[1,6,1,1] - a[1,7,1,1] - a[1,8,1,1] - a[1,9,1,1] + b[1,1,2] - x[1,2]
   >= -1
 R3350: - a[2,2,1,1] - a[2,3,1,1] - a[2,4,1,1] - a[2,5,1,1] - a[2,6,1,1]
   - a[2,7,1,1] - a[2,8,1,1] - a[2,9,1,1] + b[1,2,0] - x[2,0] >= -1
 R3351: - a[2,2,1,1] - a[2,3,1,1] - a[2,4,1,1] - a[2,5,1,1] - a[2,6,1,1]
   - a[2,7,1,1] - a[2,8,1,1] - a[2,9,1,1] + b[1,2,1] - x[2,1] >= -1
 R3352: - a[2,2,1,1] - a[2,3,1,1] - a[2,4,1,1] - a[2,5,1,1] - a[2,6,1,1]
   - a[2,7,1,1] - a[2,8,1,1] - a[2,9,1,1] + b[1,2,2] - x[2,2] >= -1
 R3353: - a[3,3,1,1] - a[3,4,1,1] - a[3,5,1,1] - a[3,6,1,1] - a[3,7,1,1]
   - a[3,8,1,1] - a[3,9,1,1] + b[1,3,0] - x[3,0] >= -1
 R3354: - a[3,3,1,1] - a[3,4,1,1] - a[3,5,1,1] - a[3,6,1,1] - a[3,7,1,1]
   - a[3,8,1,1] - a[3,9,1,1] + b[1,3,1] - x[3,1] >= -1
 R3355: - a[3,3,1,1] - a[3,4,1,1] - a[3,5,1,1] - a[3,6,1,1] - a[3,7,1,1]
   - a[3,8,1,1] - a[3,9,1,1] + b[1,3,2] - x[3,2] >= -1
 R3356: - a[4,4,1,1] - a[4,5,1,1] - a[4,6,1,1] - a[4,7,1,1] - a[4,8,1,1]
   - a[4,9,1,1] + b[1,4,0] - x[4,0] >= -1
 R3357: - a[4,4,1,1] - a[4,5,1,1] - a[4,6,1,1] - a[4,7,1,1] - a[4,8,1,1]
   - a[4,9,1,1] + b[1,4,1] - x[4,1] >= -1
 R3358: - a[4,4,1,1] - a[4,5,1,1] - a[4,6,1,1] - a[4,7,1,1] - a[4,8,1,1]
   - a[4,9,1,1] + b[1,4,2] - x[4,2] >= -1
 R3359: - a[5,5,1,1] - a[5,6,1,1] - a[5,7,1,1] - a[5,8,1,1] - a[5,9,1,1]
   + b[1,5,0] - x[5,0] >= -1
 R3360: - a[5,5,1,1] - a[5,6,1,1] - a[5,7,1,1] - a[5,8,1,1] - a[5,9,1,1]
   + b[1,5,1] - x[5,1] >= -1
 R3361: - a[5,5,1,1] - a[5,6,1,1] - a[5,7,1,1] - a[5,8,1,1] - a[5,9,1,1]
   + b[1,5,2] - x[5,2] >= -1
 R3362: - a[6,6,1,1] - a[6,7,1,1] - a[6,8,1,1] - a[6,9,1,1] + b[1,6,0]
   - x[6,0] >= -1
 R3363: - a[6,6,1,1] - a[6,7,1,1] - a[6,8,1,1] - a[6,9,1,1] + b[1,6,1]
   - x[6,1] >= -1
 R3364: - a[6,6,1,1] - a[6,7,1,1] - a[6,8,1,1] - a[6,9,1,1] + b[1,6,2]
   - x[6,2] >= -1
 R3365: - a[7,7,1,1] - a[7,8,1,1] - a[7,9,1,1] + b[1,7,0] - x[7,0] >= -1
 R3366: - a[7,7,1,1] - a[7,8,1,1] - a[7,9,1,1] + b[1,7,1] - x[7,1] >= -1
 R3367: - a[7,7,1,1] - a[7,8,1,1] - a[7,9,1,1] + b[1,7,2] - x[7,2] >= -1
 R3368: - a[8,8,1,1] - a[8,9,1,1] + b[1,8,0] - x[8,0] >= -1
 R3369: - a[8,8,1,1] - a[8,9,1,1] + b[1,8,1] - x[8,1] >= -1
 R3370: - a[8,8,1,1] - a[8,9,1,1] + b[1,8,2] - x[8,2] >= -1
 R3371: - a[9,9,1,1] + b[1,9,0] - x[9,0] >= -1
 R3372: - a[9,9,1,1] + b[1,9,1] - x[9,1] >= -1
 R3373: - a[9,9,1,1] + b[1,9,2] - x[9,2] >= -1
 R3374: - a[0,0,1,3] - a[0,1,1,3] - a[0,2,1,3] - a[0,3,1,3] - a[0,4,1,3]
   - a[0,5,1,3] - a[0,6,1,3] - a[0,7,1,3] - a[0,8,1,3] - a[0,9,1,3]
   - a[0,0,3,3] - a[0,1,3,3] - a[0,2,3,3] - a[0,3,3,3] - a[0,4,3,3]
   - a[0,5,3,3] - a[0,6,3,3] - a[0,7,3,3] - a[0,8,3,3] - a[0,9,3,3]
   + b[3,0,0] - x[0,0] >= -1
 R3375: - a[0,0,1,3] - a[0,1,1,3] - a[0,2,1,3] - a[0,3,1,3] - a[0,4,1,3]
   - a[0,5,1,3] - a[0,6,1,3] - a[0,7,1,3] - a[0,8,1,3] - a[0,9,1,3]
   - a[0,0,3,3] - a[0,1,3,3] - a[0,2,3,3] - a[0,3,3,3] - a[0,4,3,3]
   - a[0,5,3,3] - a[0,6,3,3] - a[0,7,3,3] - a[0,8,3,3] - a[0,9,3,3]
   + b[3,0,1] - x[0,1] >= -1
 R3376: - a[0,0,1,3] - a[0,1,1,3] - a[0,2,1,3] - a[0,3,1,3] - a[0,4,1,3]
   - a[0,5,1,3] - a[0,6,1,3] - a[0,7,1,3] - a[0,8,1,3] - a[0,9,1,3]
   - a[0,0,3,3] - a[0,1,3,3] - a[0,2,3,3] - a[0,3,3,3] - a[0,4,3,3]
   - a[0,5,3,3] - a[0,6,3,3] - a[0,7,3,3] - a[0,8,3,3] - a[0,9,3,3]
   + b[3,0,2] - x[0,2] >= -1
 R3377: - a[1,1,1,3] - a[1,2,1,3] - a[1,3,1,3] - a[1,4,1,3] - a[1,5,1,3]
   - a[1,6,1,3] - a[1,7,1,3] - a[1,8,1,3] - a[1,9,1,3] - a[1,1,3,3]
   - a[1,2,3,3] - a[1,3,3,3] - a[1,4,3,3] - a[1,5,3,3] - a[1,6,3,3]
   - a[1,7,3,3] - a[1,8,3,3] - a[1,9,3,3] + b[3,1,0] - x[1,0] >= -1
 R3378: - a[1,1,1,3] - a[1,2,1,3] - a[1,3,1,3] - a[1,4,1,3] - a[1,5,1,3]
   - a[1,6,1,3] - a[1,7,1,3] - a[1,8,1,3] - a[1,9,1,3] - a[1,1,3,3]
   - a[1,2,3,3] - a[1,3,3,3] - a[1,4,3,3] - a[1,5,3,3] - a[1,6,3,3]
   - a[1,7,3,3] - a[1,8,3,3] - a[1,9,3,3] + b[3,1,1] - x[1,1] >= -1
 R3379: - a[1,1,1,3] - a[1,2,1,3] - a[1,3,1,3] - a[1,4,1,3] - a[1,5,1,3]
   - a[1,6,1,3] - a[1,7,1,3] - a[1,8,1,3] - a[1,9,1,3] - a[1,1,3,3]
   - a[1,2,3,3] - a[1,3,3,3] - a[1,4,3,3] - a[1,5,3,3] - a[1,6,3,3]
   - a[1,7,3,3] - a[1,8,3,3] - a[1,9,3,3] + b[3,1,2] - x[1,2] >= -1
 R3380: - a[2,2,1,3] - a[2,3,1,3] - a[2,4,1,3] - a[2,5,1,3] - a[2,6,1,3]
   - a[2,7,1,3] - a[2,8,1,3] - a[2,9,1,3] - a[2,2,3,3] - a[2,3,3,3]
   - a[2,4,3,3] - a[2,5,3,3] - a[2,6,3,3] - a[2,7,3,3] - a[2,8,3,3]
   - a[2,9,3,3] + b[3,2,0] - x[2,0] >= -1
 R3381: - a[2,2,1,3] - a[2,3,1,3] - a[2,4,1,3] - a[2,5,1,3] - a[2,6,1,3]
   - a[2,7,1,3] - a[2,8,1,3] - a[2,9,1,3] - a[2,2,3,3] - a[2,3,3,3]
   - a[2,4,3,3] - a[2,5,3,3] - a[2,6,3,3] - a[2,7,3,3] - a[2,8,3,3]
   - a[2,9,3,3] + b[3,2,1] - x[2,1] >= -1
 R3382: - a[2,2,1,3] - a[2,3,1,3] - a[2,4,1,3] - a[2,5,1,3] - a[2,6,1,3]
   - a[2,7,1,3] - a[2,8,1,3] - a[2,9,1,3] - a[2,2,3,3] - a[2,3,3,3]
   - a[2,4,3,3] - a[2,5,3,3] - a[2,6,3,3] - a[2,7,3,3] - a[2,8,3,3]
   - a[2,9,3,3] + b[3,2,2] - x[2,2] >= -1
 R3383: - a[3,3,1,3] - a[3,4,1,3] - a[3,5,1,3] - a[3,6,1,3] - a[3,7,1,3]
   - a[3,8,1,3] - a[3,9,1,3] - a[3,3,3,3] - a[3,4,3,3] - a[3,5,3,3]
   - a[3,6,3,3] - a[3,7,3,3] - a[3,8,3,3] - a[3,9,3,3] + b[3,3,0] - x[3,0]
   >= -1
 R3384: - a[3,3,1,3] - a[3,4,1,3] - a[3,5,1,3] - a[3,6,1,3] - a[3,7,1,3]
   - a[3,8,1,3] - a[3,9,1,3] - a[3,3,3,3] - a[3,4,3,3] - a[3,5,3,3]
   - a[3,6,3,3] - a[3,7,3,3] - a[3,8,3,3] - a[3,9,3,3] + b[3,3,1] - x[3,1]
   >= -1
 R3385: - a[3,3,1,3] - a[3,4,1,3] - a[3,5,1,3] - a[3,6,1,3] - a[3,7,1,3]
   - a[3,8,1,3] - a[3,9,1,3] - a[3,3,3,3] - a[3,4,3,3] - a[3,5,3,3]
   - a[3,6,3,3] - a[3,7,3,3] - a[3,8,3,3] - a[3,9,3,3] + b[3,3,2] - x[3,2]
   >= -1
 R3386: - a[4,4,1,3] - a[4,5,1,3] - a[4,6,1,3] - a[4,7,1,3] - a[4,8,1,3]
   - a[4,9,1,3] - a[4,4,3,3] - a[4,5,3,3] - a[4,6,3,3] - a[4,7,3,3]
   - a[4,8,3,3] - a[4,9,3,3] + b[3,4,0] - x[4,0] >= -1
 R3387: - a[4,4,1,3] - a[4,5,1,3] - a[4,6,1,3] - a[4,7,1,3] - a[4,8,1,3]
   - a[4,9,1,3] - a[4,4,3,3] - a[4,5,3,3] - a[4,6,3,3] - a[4,7,3,3]
   - a[4,8,3,3] - a[4,9,3,3] + b[3,4,1] - x[4,1] >= -1
 R3388: - a[4,4,1,3] - a[4,5,1,3] - a[4,6,1,3] - a[4,7,1,3] - a[4,8,1,3]
   - a[4,9,1,3] - a[4,4,3,3] - a[4,5,3,3] - a[4,6,3,3] - a[4,7,3,3]
   - a[4,8,3,3] - a[4,9,3,3] + b[3,4,2] - x[4,2] >= -1
 R3389: - a[5,5,1,3] - a[5,6,1,3] - a[5,7,1,3] - a[5,8,1,3] - a[5,9,1,3]
   - a[5,5,3,3] - a[5,6,3,3] - a[5,7,3,3] - a[5,8,3,3] - a[5,9,3,3]
   + b[3,5,0] - x[5,0] >= -1
 R3390: - a[5,5,1,3] - a[5,6,1,3] - a[5,7,1,3] - a[5,8,1,3] - a[5,9,1,3]
   - a[5,5,3,3] - a[5,6,3,3] - a[5,7,3,3] - a[5,8,3,3] - a[5,9,3,3]
   + b[3,5,1] - x[5,1] >= -1
 R3391: - a[5,5,1,3] - a[5,6,1,3] - a[5,7,1,3] - a[5,8,1,3] - a[5,9,1,3]
   - a[5,5,3,3] - a[5,6,3,3] - a[5,7,3,3] - a[5,8,3,3] - a[5,9,3,3]
   + b[3,5,2] - x[5,2] >= -1
 R3392: - a[6,6,1,3] - a[6,7,1,3] - a[6,8,1,3] - a[6,9,1,3] - a[6,6,3,3]
   - a[6,7,3,3] - a[6,8,3,3] - a[6,9,3,3] + b[3,6,0] - x[6,0] >= -1
 R3393: - a[6,6,1,3] - a[6,7,1,3] - a[6,8,1,3] - a[6,9,1,3] - a[6,6,3,3]
   - a[6,7,3,3] - a[6,8,3,3] - a[6,9,3,3] + b[3,6,1] - x[6,1] >= -1
 R3394: - a[6,6,1,3] - a[6,7,1,3] - a[6,8,1,3] - a[6,9,1,3] - a[6,6,3,3]
   - a[6,7,3,3] - a[6,8,3,3] - a[6,9,3,3] + b[3,6,2] - x[6,2] >= -1
 R3395: - a[7,7,1,3] - a[7,8,1,3] - a[7,9,1,3] - a[7,7,3,3] - a[7,8,3,3]
   - a[7,9,3,3] + b[3,7,0] - x[7,0] >= -1
 R3396: - a[7,7,1,3] - a[7,8,1,3] - a[7,9,1,3] - a[7,7,3,3] - a[7,8,3,3]
   - a[7,9,3,3] + b[3,7,1] - x[7,1] >= -1
 R3397: - a[7,7,1,3] - a[7,8,1,3] - a[7,9,1,3] - a[7,7,3,3] - a[7,8,3,3]
   - a[7,9,3,3] + b[3,7,2] - x[7,2] >= -1
 R3398: - a[8,8,1,3] - a[8,9,1,3] - a[8,8,3,3] - a[8,9,3,3] + b[3,8,0]
   - x[8,0] >= -1
 R3399: - a[8,8,1,3] - a[8,9,1,3] - a[8,8,3,3] - a[8,9,3,3] + b[3,8,1]
   - x[8,1] >= -1
 R3400: - a[8,8,1,3] - a[8,9,1,3] - a[8,8,3,3] - a[8,9,3,3] + b[3,8,2]
   - x[8,2] >= -1
 R3401: - a[9,9,1,3] - a[9,9,3,3] + b[3,9,0] - x[9,0] >= -1
 R3402: - a[9,9,1,3] - a[9,9,3,3] + b[3,9,1] - x[9,1] >= -1
 R3403: - a[9,9,1,3] - a[9,9,3,3] + b[3,9,2] - x[9,2] >= -1
 R3404: - a[0,0,1,7] - a[0,1,1,7] - a[0,2,1,7] - a[0,3,1,7] - a[0,4,1,7]
   - a[0,5,1,7] - a[0,6,1,7] - a[0,7,1,7] - a[0,8,1,7] - a[0,9,1,7]
   - a[0,0,3,7] - a[0,1,3,7] - a[0,2,3,7] - a[0,3,3,7] - a[0,4,3,7]
   - a[0,5,3,7] - a[0,6,3,7] - a[0,7,3,7] - a[0,8,3,7] - a[0,9,3,7]
   - a[0,0,7,7] - a[0,1,7,7] - a[0,2,7,7] - a[0,3,7,7] - a[0,4,7,7]
   - a[0,5,7,7] - a[0,6,7,7] - a[0,7,7,7] - a[0,8,7,7] - a[0,9,7,7]
   + b[7,0,0] - x[0,0] >= -1
 R3405: - a[0,0,1,7] - a[0,1,1,7] - a[0,2,1,7] - a[0,3,1,7] - a[0,4,1,7]
   - a[0,5,1,7] - a[0,6,1,7] - a[0,7,1,7] - a[0,8,1,7] - a[0,9,1,7]
   - a[0,0,3,7] - a[0,1,3,7] - a[0,2,3,7] - a[0,3,3,7] - a[0,4,3,7]
   - a[0,5,3,7] - a[0,6,3,7] - a[0,7,3,7] - a[0,8,3,7] - a[0,9,3,7]
   - a[0,0,7,7] - a[0,1,7,7] - a[0,2,7,7] - a[0,3,7,7] - a[0,4,7,7]
   - a[0,5,7,7] - a[0,6,7,7] - a[0,7,7,7] - a[0,8,7,7] - a[0,9,7,7]
   + b[7,0,1] - x[0,1] >= -1
 R3406: - a[0,0,1,7] - a[0,1,1,7] - a[0,2,1,7] - a[0,3,1,7] - a[0,4,1,7]
   - a[0,5,1,7] - a[0,6,1,7] - a[0,7,1,7] - a[0,8,1,7] - a[0,9,1,7]
   - a[0,0,3,7] - a[0,1,3,7] - a[0,2,3,7] - a[0,3,3,7] - a[0,4,3,7]
   - a[0,5,3,7] - a[0,6,3,7] - a[0,7,3,7] - a[0,8,3,7] - a[0,9,3,7]
   - a[0,0,7,7] - a[0,1,7,7] - a[0,2,7,7] - a[0,3,7,7] - a[0,4,7,7]
   - a[0,5,7,7] - a[0,6,7,7] - a[0,7,7,7] - a[0,8,7,7] - a[0,9,7,7]
   + b[7,0,2] - x[0,2] >= -1
 R3407: - a[1,1,1,7] - a[1,2,1,7] - a[1,3,1,7] - a[1,4,1,7] - a[1,5,1,7]
   - a[1,6,1,7] - a[1,7,1,7] - a[1,8,1,7] - a[1,9,1,7] - a[1,1,3,7]
   - a[1,2,3,7] - a[1,3,3,7] - a[1,4,3,7] - a[1,5,3,7] - a[1,6,3,7]
   - a[1,7,3,7] - a[1,8,3,7] - a[1,9,3,7] - a[1,1,7,7] - a[1,2,7,7]
   - a[1,3,7,7] - a[1,4,7,7] - a[1,5,7,7] - a[1,6,7,7] - a[1,7,7,7]
   - a[1,8,7,7] - a[1,9,7,7] + b[7,1,0] - x[1,0] >= -1
 R3408: - a[1,1,1,7] - a[1,2,1,7] - a[1,3,1,7] - a[1,4,1,7] - a[1,5,1,7]
   - a[1,6,1,7] - a[1,7,1,7] - a[1,8,1,7] - a[1,9,1,7] - a[1,1,3,7]
   - a[1,2,3,7] - a[1,3,3,7] - a[1,4,3,7] - a[1,5,3,7] - a[1,6,3,7]
   - a[1,7,3,7] - a[1,8,3,7] - a[1,9,3,7] - a[1,1,7,7] - a[1,2,7,7]
   - a[1,3,7,7] - a[1,4,7,7] - a[1,5,7,7] - a[1,6,7,7] - a[1,7,7,7]
   - a[1,8,7,7] - a[1,9,7,7] + b[7,1,1] - x[1,1] >= -1
 R3409: - a[1,1,1,7] - a[1,2,1,7] - a[1,3,1,7] - a[1,4,1,7] - a[1,5,1,7]
   - a[1,6,1,7] - a[1,7,1,7] - a[1,8,1,7] - a[1,9,1,7] - a[1,1,3,7]
   - a[1,2,3,7] - a[1,3,3,7] - a[1,4,3,7] - a[1,5,3,7] - a[1,6,3,7]
   - a[1,7,3,7] - a[1,8,3,7] - a[1,9,3,7] - a[1,1,7,7] - a[1,2,7,7]
   - a[1,3,7,7] - a[1,4,7,7] - a[1,5,7,7] - a[1,6,7,7] - a[1,7,7,7]
   - a[1,8,7,7] - a[1,9,7,7] + b[7,1,2] - x[1,2] >= -1
 R3410: - a[2,2,1,7] - a[2,3,1,7] - a[2,4,1,7] - a[2,5,1,7] - a[2,6,1,7]
   - a[2,7,1,7] - a[2,8,1,7] - a[2,9,1,7] - a[2,2,3,7] - a[2,3,3,7]
   - a[2,4,3,7] - a[2,5,3,7] - a[2,6,3,7] - a[2,7,3,7] - a[2,8,3,7]
   - a[2,9,3,7] - a[2,2,7,7] - a[2,3,7,7] - a[2,4,7,7] - a[2,5,7,7]
   - a[2,6,7,7] - a[2,7,7,7] - a[2,8,7,7] - a[2,9,7,7] + b[7,2,0] - x[2,0]
   >= -1
 R3411: - a[2,2,1,7] - a[2,3,1,7] - a[2,4,1,7] - a[2,5,1,7] - a[2,6,1,7]
   - a[2,7,1,7] - a[2,8,1,7] - a[2,9,1,7] - a[2,2,3,7] - a[2,3,3,7]
   - a[2,4,3,7] - a[2,5,3,7] - a[2,6,3,7] - a[2,7,3,7] - a[2,8,3,7]
   - a[2,9,3,7] - a[2,2,7,7] - a[2,3,7,7] - a[2,4,7,7] - a[2,5,7,7]
   - a[2,6,7,7] - a[2,7,7,7] - a[2,8,7,7] - a[2,9,7,7] + b[7,2,1] - x[2,1]
   >= -1
 R3412: - a[2,2,1,7] - a[2,3,1,7] - a[2,4,1,7] - a[2,5,1,7] - a[2,6,1,7]
   - a[2,7,1,7] - a[2,8,1,7] - a[2,9,1,7] - a[2,2,3,7] - a[2,3,3,7]
   - a[2,4,3,7] - a[2,5,3,7] - a[2,6,3,7] - a[2,7,3,7] - a[2,8,3,7]
   - a[2,9,3,7] - a[2,2,7,7] - a[2,3,7,7] - a[2,4,7,7] - a[2,5,7,7]
   - a[2,6,7,7] - a[2,7,7,7] - a[2,8,7,7] - a[2,9,7,7] + b[7,2,2] - x[2,2]
   >= -1
 R3413: - a[3,3,1,7] - a[3,4,1,7] - a[3,5,1,7] - a[3,6,1,7] - a[3,7,1,7]
   - a[3,8,1,7] - a[3,9,1,7] - a[3,3,3,7] - a[3,4,3,7] - a[3,5,3,7]
   - a[3,6,3,7] - a[3,7,3,7] - a[3,8,3,7] - a[3,9,3,7] - a[3,3,7,7]
   - a[3,4,7,7] - a[3,5,7,7] - a[3,6,7,7] - a[3,7,7,7] - a[3,8,7,7]
   - a[3,9,7,7] + b[7,3,0] - x[3,0] >= -1
 R3414: - a[3,3,1,7] - a[3,4,1,7] - a[3,5,1,7] - a[3,6,1,7] - a[3,7,1,7]
   - a[3,8,1,7] - a[3,9,1,7] - a[3,3,3,7] - a[3,4,3,7] - a[3,5,3,7]
   - a[3,6,3,7] - a[3,7,3,7] - a[3,8,3,7] - a[3,9,3,7] - a[3,3,7,7]
   - a[3,4,7,7] - a[3,5,7,7] - a[3,6,7,7] - a[3,7,7,7] - a[3,8,7,7]
   - a[3,9,7,7] + b[7,3,1] - x[3,1] >= -1
 R3415: - a[3,3,1,7] - a[3,4,1,7] - a[3,5,1,7] - a[3,6,1,7] - a[3,7,1,7]
   - a[3,8,1,7] - a[3,9,1,7] - a[3,3,3,7] - a[3,4,3,7] - a[3,5,3,7]
   - a[3,6,3,7] - a[3,7,3,7] - a[3,8,3,7] - a[3,9,3,7] - a[3,3,7,7]
   - a[3,4,7,7] - a[3,5,7,7] - a[3,6,7,7] - a[3,7,7,7] - a[3,8,7,7]
   - a[3,9,7,7] + b[7,3,2] - x[3,2] >= -1
 R3416: - a[4,4,1,7] - a[4,5,1,7] - a[4,6,1,7] - a[4,7,1,7] - a[4,8,1,7]
   - a[4,9,1,7] - a[4,4,3,7] - a[4,5,3,7] - a[4,6,3,7] - a[4,7,3,7]
   - a[4,8,3,7] - a[4,9,3,7] - a[4,4,7,7] - a[4,5,7,7] - a[4,6,7,7]
   - a[4,7,7,7] - a[4,8,7,7] - a[4,9,7,7] + b[7,4,0] - x[4,0] >= -1
 R3417: - a[4,4,1,7] - a[4,5,1,7] - a[4,6,1,7] - a[4,7,1,7] - a[4,8,1,7]
   - a[4,9,1,7] - a[4,4,3,7] - a[4,5,3,7] - a[4,6,3,7] - a[4,7,3,7]
   - a[4,8,3,7] - a[4,9,3,7] - a[4,4,7,7] - a[4,5,7,7] - a[4,6,7,7]
   - a[4,7,7,7] - a[4,8,7,7] - a[4,9,7,7] + b[7,4,1] - x[4,1] >= -1
 R3418: - a[4,4,1,7] - a[4,5,1,7] - a[4,6,1,7] - a[4,7,1,7] - a[4,8,1,7]
   - a[4,9,1,7] - a[4,4,3,7] - a[4,5,3,7] - a[4,6,3,7] - a[4,7,3,7]
   - a[4,8,3,7] - a[4,9,3,7] - a[4,4,7,7] - a[4,5,7,7] - a[4,6,7,7]
   - a[4,7,7,7] - a[4,8,7,7] - a[4,9,7,7] + b[7,4,2] - x[4,2] >= -1
 R3419: - a[5,5,1,7] - a[5,6,1,7] - a[5,7,1,7] - a[5,8,1,7] - a[5,9,1,7]
   - a[5,5,3,7] - a[5,6,3,7] - a[5,7,3,7] - a[5,8,3,7] - a[5,9,3,7]
   - a[5,5,7,7] - a[5,6,7,7] - a[5,7,7,7] - a[5,8,7,7] - a[5,9,7,7]
   + b[7,5,0] - x[5,0] >= -1
 R3420: - a[5,5,1,7] - a[5,6,1,7] - a[5,7,1,7] - a[5,8,1,7] - a[5,9,1,7]
   - a[5,5,3,7] - a[5,6,3,7] - a[5,7,3,7] - a[5,8,3,7] - a[5,9,3,7]
   - a[5,5,7,7] - a[5,6,7,7] - a[5,7,7,7] - a[5,8,7,7] - a[5,9,7,7]
   + b[7,5,1] - x[5,1] >= -1
 R3421: - a[5,5,1,7] - a[5,6,1,7] - a[5,7,1,7] - a[5,8,1,7] - a[5,9,1,7]
   - a[5,5,3,7] - a[5,6,3,7] - a[5,7,3,7] - a[5,8,3,7] - a[5,9,3,7]
   - a[5,5,7,7] - a[5,6,7,7] - a[5,7,7,7] - a[5,8,7,7] - a[5,9,7,7]
   + b[7,5,2] - x[5,2] >= -1
 R3422: - a[6,6,1,7] - a[6,7,1,7] - a[6,8,1,7] - a[6,9,1,7] - a[6,6,3,7]
   - a[6,7,3,7] - a[6,8,3,7] - a[6,9,3,7] - a[6,6,7,7] - a[6,7,7,7]
   - a[6,8,7,7] - a[6,9,7,7] + b[7,6,0] - x[6,0] >= -1
 R3423: - a[6,6,1,7] - a[6,7,1,7] - a[6,8,1,7] - a[6,9,1,7] - a[6,6,3,7]
   - a[6,7,3,7] - a[6,8,3,7] - a[6,9,3,7] - a[6,6,7,7] - a[6,7,7,7]
   - a[6,8,7,7] - a[6,9,7,7] + b[7,6,1] - x[6,1] >= -1
 R3424: - a[6,6,1,7] - a[6,7,1,7] - a[6,8,1,7] - a[6,9,1,7] - a[6,6,3,7]
   - a[6,7,3,7] - a[6,8,3,7] - a[6,9,3,7] - a[6,6,7,7] - a[6,7,7,7]
   - a[6,8,7,7] - a[6,9,7,7] + b[7,6,2] - x[6,2] >= -1
 R3425: - a[7,7,1,7] - a[7,8,1,7] - a[7,9,1,7] - a[7,7,3,7] - a[7,8,3,7]
   - a[7,9,3,7] - a[7,7,7,7] - a[7,8,7,7] - a[7,9,7,7] + b[7,7,0] - x[7,0]
   >= -1
 R3426: - a[7,7,1,7] - a[7,8,1,7] - a[7,9,1,7] - a[7,7,3,7] - a[7,8,3,7]
   - a[7,9,3,7] - a[7,7,7,7] - a[7,8,7,7] - a[7,9,7,7] + b[7,7,1] - x[7,1]
   >= -1
 R3427: - a[7,7,1,7] - a[7,8,1,7] - a[7,9,1,7] - a[7,7,3,7] - a[7,8,3,7]
   - a[7,9,3,7] - a[7,7,7,7] - a[7,8,7,7] - a[7,9,7,7] + b[7,7,2] - x[7,2]
   >= -1
 R3428: - a[8,8,1,7] - a[8,9,1,7] - a[8,8,3,7] - a[8,9,3,7] - a[8,8,7,7]
   - a[8,9,7,7] + b[7,8,0] - x[8,0] >= -1
 R3429: - a[8,8,1,7] - a[8,9,1,7] - a[8,8,3,7] - a[8,9,3,7] - a[8,8,7,7]
   - a[8,9,7,7] + b[7,8,1] - x[8,1] >= -1
 R3430: - a[8,8,1,7] - a[8,9,1,7] - a[8,8,3,7] - a[8,9,3,7] - a[8,8,7,7]
   - a[8,9,7,7] + b[7,8,2] - x[8,2] >= -1
 R3431: - a[9,9,1,7] - a[9,9,3,7] - a[9,9,7,7] + b[7,9,0] - x[9,0] >= -1
 R3432: - a[9,9,1,7] - a[9,9,3,7] - a[9,9,7,7] + b[7,9,1] - x[9,1] >= -1
 R3433: - a[9,9,1,7] - a[9,9,3,7] - a[9,9,7,7] + b[7,9,2] - x[9,2] >= -1
 R3434: - a[0,0,4,4] - a[0,1,4,4] - a[0,2,4,4] - a[0,3,4,4] - a[0,4,4,4]
   - a[0,5,4,4] - a[0,6,4,4] - a[0,7,4,4] - a[0,8,4,4] - a[0,9,4,4]
   + b[4,0,0] - x[0,0] >= -1
 R3435: - a[0,0,4,4] - a[0,1,4,4] - a[0,2,4,4] - a[0,3,4,4] - a[0,4,4,4]
   - a[0,5,4,4] - a[0,6,4,4] - a[0,7,4,4] - a[0,8,4,4] - a[0,9,4,4]
   + b[4,0,1] - x[0,1] >= -1
 R3436: - a[0,0,4,4] - a[0,1,4,4] - a[0,2,4,4] - a[0,3,4,4] - a[0,4,4,4]
   - a[0,5,4,4] - a[0,6,4,4] - a[0,7,4,4] - a[0,8,4,4] - a[0,9,4,4]
   + b[4,0,2] - x[0,2] >= -1
 R3437: - a[1,1,4,4] - a[1,2,4,4] - a[1,3,4,4] - a[1,4,4,4] - a[1,5,4,4]
   - a[1,6,4,4] - a[1,7,4,4] - a[1,8,4,4] - a[1,9,4,4] + b[4,1,0] - x[1,0]
   >= -1
 R3438: - a[1,1,4,4] - a[1,2,4,4] - a[1,3,4,4] - a[1,4,4,4] - a[1,5,4,4]
   - a[1,6,4,4] - a[1,7,4,4] - a[1,8,4,4] - a[1,9,4,4] + b[4,1,1] - x[1,1]
   >= -1
 R3439: - a[1,1,4,4] - a[1,2,4,4] - a[1,3,4,4] - a[1,4,4,4] - a[1,5,4,4]
   - a[1,6,4,4] - a[1,7,4,4] - a[1,8,4,4] - a[1,9,4,4] + b[4,1,2] - x[1,2]
   >= -1
 R3440: - a[2,2,4,4] - a[2,3,4,4] - a[2,4,4,4] - a[2,5,4,4] - a[2,6,4,4]
   - a[2,7,4,4] - a[2,8,4,4] - a[2,9,4,4] + b[4,2,0] - x[2,0] >= -1
 R3441: - a[2,2,4,4] - a[2,3,4,4] - a[2,4,4,4] - a[2,5,4,4] - a[2,6,4,4]
   - a[2,7,4,4] - a[2,8,4,4] - a[2,9,4,4] + b[4,2,1] - x[2,1] >= -1
 R3442: - a[2,2,4,4] - a[2,3,4,4] - a[2,4,4,4] - a[2,5,4,4] - a[2,6,4,4]
   - a[2,7,4,4] - a[2,8,4,4] - a[2,9,4,4] + b[4,2,2] - x[2,2] >= -1
 R3443: - a[3,3,4,4] - a[3,4,4,4] - a[3,5,4,4] - a[3,6,4,4] - a[3,7,4,4]
   - a[3,8,4,4] - a[3,9,4,4] + b[4,3,0] - x[3,0] >= -1
 R3444: - a[3,3,4,4] - a[3,4,4,4] - a[3,5,4,4] - a[3,6,4,4] - a[3,7,4,4]
   - a[3,8,4,4] - a[3,9,4,4] + b[4,3,1] - x[3,1] >= -1
 R3445: - a[3,3,4,4] - a[3,4,4,4] - a[3,5,4,4] - a[3,6,4,4] - a[3,7,4,4]
   - a[3,8,4,4] - a[3,9,4,4] + b[4,3,2] - x[3,2] >= -1
 R3446: - a[4,4,4,4] - a[4,5,4,4] - a[4,6,4,4] - a[4,7,4,4] - a[4,8,4,4]
   - a[4,9,4,4] + b[4,4,0] - x[4,0] >= -1
 R3447: - a[4,4,4,4] - a[4,5,4,4] - a[4,6,4,4] - a[4,7,4,4] - a[4,8,4,4]
   - a[4,9,4,4] + b[4,4,1] - x[4,1] >= -1
 R3448: - a[4,4,4,4] - a[4,5,4,4] - a[4,6,4,4] - a[4,7,4,4] - a[4,8,4,4]
   - a[4,9,4,4] + b[4,4,2] - x[4,2] >= -1
 R3449: - a[5,5,4,4] - a[5,6,4,4] - a[5,7,4,4] - a[5,8,4,4] - a[5,9,4,4]
   + b[4,5,0] - x[5,0] >= -1
 R3450: - a[5,5,4,4] - a[5,6,4,4] - a[5,7,4,4] - a[5,8,4,4] - a[5,9,4,4]
   + b[4,5,1] - x[5,1] >= -1
 R3451: - a[5,5,4,4] - a[5,6,4,4] - a[5,7,4,4] - a[5,8,4,4] - a[5,9,4,4]
   + b[4,5,2] - x[5,2] >= -1
 R3452: - a[6,6,4,4] - a[6,7,4,4] - a[6,8,4,4] - a[6,9,4,4] + b[4,6,0]
   - x[6,0] >= -1
 R3453: - a[6,6,4,4] - a[6,7,4,4] - a[6,8,4,4] - a[6,9,4,4] + b[4,6,1]
   - x[6,1] >= -1
 R3454: - a[6,6,4,4] - a[6,7,4,4] - a[6,8,4,4] - a[6,9,4,4] + b[4,6,2]
   - x[6,2] >= -1
 R3455: - a[7,7,4,4] - a[7,8,4,4] - a[7,9,4,4] + b[4,7,0] - x[7,0] >= -1
 R3456: - a[7,7,4,4] - a[7,8,4,4] - a[7,9,4,4] + b[4,7,1] - x[7,1] >= -1
 R3457: - a[7,7,4,4] - a[7,8,4,4] - a[7,9,4,4] + b[4,7,2] - x[7,2] >= -1
 R3458: - a[8,8,4,4] - a[8,9,4,4] + b[4,8,0] - x[8,0] >= -1
 R3459: - a[8,8,4,4] - a[8,9,4,4] + b[4,8,1] - x[8,1] >= -1
 R3460: - a[8,8,4,4] - a[8,9,4,4] + b[4,8,2] - x[8,2] >= -1
 R3461: - a[9,9,4,4] + b[4,9,0] - x[9,0] >= -1
 R3462: - a[9,9,4,4] + b[4,9,1] - x[9,1] >= -1
 R3463: - a[9,9,4,4] + b[4,9,2] - x[9,2] >= -1
 R3464: - a[0,0,5,5] - a[0,1,5,5] - a[0,2,5,5] - a[0,3,5,5] - a[0,4,5,5]
   - a[0,5,5,5] - a[0,6,5,5] - a[0,7,5,5] - a[0,8,5,5] - a[0,9,5,5]
   + b[5,0,0] - x[0,0] >= -1
 R3465: - a[0,0,5,5] - a[0,1,5,5] - a[0,2,5,5] - a[0,3,5,5] - a[0,4,5,5]
   - a[0,5,5,5] - a[0,6,5,5] - a[0,7,5,5] - a[0,8,5,5] - a[0,9,5,5]
   + b[5,0,1] - x[0,1] >= -1
 R3466: - a[0,0,5,5] - a[0,1,5,5] - a[0,2,5,5] - a[0,3,5,5] - a[0,4,5,5]
   - a[0,5,5,5] - a[0,6,5,5] - a[0,7,5,5] - a[0,8,5,5] - a[0,9,5,5]
   + b[5,0,2] - x[0,2] >= -1
 R3467: - a[1,1,5,5] - a[1,2,5,5] - a[1,3,5,5] - a[1,4,5,5] - a[1,5,5,5]
   - a[1,6,5,5] - a[1,7,5,5] - a[1,8,5,5] - a[1,9,5,5] + b[5,1,0] - x[1,0]
   >= -1
 R3468: - a[1,1,5,5] - a[1,2,5,5] - a[1,3,5,5] - a[1,4,5,5] - a[1,5,5,5]
   - a[1,6,5,5] - a[1,7,5,5] - a[1,8,5,5] - a[1,9,5,5] + b[5,1,1] - x[1,1]
   >= -1
 R3469: - a[1,1,5,5] - a[1,2,5,5] - a[1,3,5,5] - a[1,4,5,5] - a[1,5,5,5]
   - a[1,6,5,5] - a[1,7,5,5] - a[1,8,5,5] - a[1,9,5,5] + b[5,1,2] - x[1,2]
   >= -1
 R3470: - a[2,2,5,5] - a[2,3,5,5] - a[2,4,5,5] - a[2,5,5,5] - a[2,6,5,5]
   - a[2,7,5,5] - a[2,8,5,5] - a[2,9,5,5] + b[5,2,0] - x[2,0] >= -1
 R3471: - a[2,2,5,5] - a[2,3,5,5] - a[2,4,5,5] - a[2,5,5,5] - a[2,6,5,5]
   - a[2,7,5,5] - a[2,8,5,5] - a[2,9,5,5] + b[5,2,1] - x[2,1] >= -1
 R3472: - a[2,2,5,5] - a[2,3,5,5] - a[2,4,5,5] - a[2,5,5,5] - a[2,6,5,5]
   - a[2,7,5,5] - a[2,8,5,5] - a[2,9,5,5] + b[5,2,2] - x[2,2] >= -1
 R3473: - a[3,3,5,5] - a[3,4,5,5] - a[3,5,5,5] - a[3,6,5,5] - a[3,7,5,5]
   - a[3,8,5,5] - a[3,9,5,5] + b[5,3,0] - x[3,0] >= -1
 R3474: - a[3,3,5,5] - a[3,4,5,5] - a[3,5,5,5] - a[3,6,5,5] - a[3,7,5,5]
   - a[3,8,5,5] - a[3,9,5,5] + b[5,3,1] - x[3,1] >= -1
 R3475: - a[3,3,5,5] - a[3,4,5,5] - a[3,5,5,5] - a[3,6,5,5] - a[3,7,5,5]
   - a[3,8,5,5] - a[3,9,5,5] + b[5,3,2] - x[3,2] >= -1
 R3476: - a[4,4,5,5] - a[4,5,5,5] - a[4,6,5,5] - a[4,7,5,5] - a[4,8,5,5]
   - a[4,9,5,5] + b[5,4,0] - x[4,0] >= -1
 R3477: - a[4,4,5,5] - a[4,5,5,5] - a[4,6,5,5] - a[4,7,5,5] - a[4,8,5,5]
   - a[4,9,5,5] + b[5,4,1] - x[4,1] >= -1
 R3478: - a[4,4,5,5] - a[4,5,5,5] - a[4,6,5,5] - a[4,7,5,5] - a[4,8,5,5]
   - a[4,9,5,5] + b[5,4,2] - x[4,2] >= -1
 R3479: - a[5,5,5,5] - a[5,6,5,5] - a[5,7,5,5] - a[5,8,5,5] - a[5,9,5,5]
   + b[5,5,0] - x[5,0] >= -1
 R3480: - a[5,5,5,5] - a[5,6,5,5] - a[5,7,5,5] - a[5,8,5,5] - a[5,9,5,5]
   + b[5,5,1] - x[5,1] >= -1
 R3481: - a[5,5,5,5] - a[5,6,5,5] - a[5,7,5,5] - a[5,8,5,5] - a[5,9,5,5]
   + b[5,5,2] - x[5,2] >= -1
 R3482: - a[6,6,5,5] - a[6,7,5,5] - a[6,8,5,5] - a[6,9,5,5] + b[5,6,0]
   - x[6,0] >= -1
 R3483: - a[6,6,5,5] - a[6,7,5,5] - a[6,8,5,5] - a[6,9,5,5] + b[5,6,1]
   - x[6,1] >= -1
 R3484: - a[6,6,5,5] - a[6,7,5,5] - a[6,8,5,5] - a[6,9,5,5] + b[5,6,2]
   - x[6,2] >= -1
 R3485: - a[7,7,5,5] - a[7,8,5,5] - a[7,9,5,5] + b[5,7,0] - x[7,0] >= -1
 R3486: - a[7,7,5,5] - a[7,8,5,5] - a[7,9,5,5] + b[5,7,1] - x[7,1] >= -1
 R3487: - a[7,7,5,5] - a[7,8,5,5] - a[7,9,5,5] + b[5,7,2] - x[7,2] >= -1
 R3488: - a[8,8,5,5] - a[8,9,5,5] + b[5,8,0] - x[8,0] >= -1
 R3489: - a[8,8,5,5] - a[8,9,5,5] + b[5,8,1] - x[8,1] >= -1
 R3490: - a[8,8,5,5] - a[8,9,5,5] + b[5,8,2] - x[8,2] >= -1
 R3491: - a[9,9,5,5] + b[5,9,0] - x[9,0] >= -1
 R3492: - a[9,9,5,5] + b[5,9,1] - x[9,1] >= -1
 R3493: - a[9,9,5,5] + b[5,9,2] - x[9,2] >= -1
 R3494: - Q[0] + z[0,0] <= 0
 R3495: - Q[1] + z[1,0] <= 0
 R3496: - Q[2] + z[2,0] <= 0
 R3497: - Q[3] + z[3,0] <= 0
 R3498: - Q[4] + z[4,0] <= 0
 R3499: - Q[5] + z[5,0] <= 0
 R3500: - Q[6] + z[6,0] <= 0
 R3501: - Q[7] + z[7,0] <= 0
 R3502: - Q[8] + z[8,0] <= 0
 R3503: - Q[9] + z[9,0] <= 0
 R3504: - Q[0] + z[0,1] <= 0
 R3505: - Q[1] + z[1,1] <= 0
 R3506: - Q[2] + z[2,1] <= 0
 R3507: - Q[3] + z[3,1] <= 0
 R3508: - Q[4] + z[4,1] <= 0
 R3509: - Q[5] + z[5,1] <= 0
 R3510: - Q[6] + z[6,1] <= 0
 R3511: - Q[7] + z[7,1] <= 0
 R3512: - Q[8] + z[8,1] <= 0
 R3513: - Q[9] + z[9,1] <= 0
 R3514: - Q[0] + z[0,2] <= 0
 R3515: - Q[1] + z[1,2] <= 0
 R3516: - Q[2] + z[2,2] <= 0
 R3517: - Q[3] + z[3,2] <= 0
 R3518: - Q[4] + z[4,2] <= 0
 R3519: - Q[5] + z[5,2] <= 0
 R3520: - Q[6] + z[6,2] <= 0
 R3521: - Q[7] + z[7,2] <= 0
 R3522: - Q[8] + z[8,2] <= 0
 R3523: - Q[9] + z[9,2] <= 0
 R3524: z[0,0] - 20000 x[0,0] <= 0
 R3525: z[1,0] - 20000 x[1,0] <= 0
 R3526: z[2,0] - 20000 x[2,0] <= 0
 R3527: z[3,0] - 20000 x[3,0] <= 0
 R3528: z[4,0] - 20000 x[4,0] <= 0
 R3529: z[5,0] - 20000 x[5,0] <= 0
 R3530: z[6,0] - 20000 x[6,0] <= 0
 R3531: z[7,0] - 20000 x[7,0] <= 0
 R3532: z[8,0] - 20000 x[8,0] <= 0
 R3533: z[9,0] - 20000 x[9,0] <= 0
 R3534: z[0,1] - 20000 x[0,1] <= 0
 R3535: z[1,1] - 20000 x[1,1] <= 0
 R3536: z[2,1] - 20000 x[2,1] <= 0
 R3537: z[3,1] - 20000 x[3,1] <= 0
 R3538: z[4,1] - 20000 x[4,1] <= 0
 R3539: z[5,1] - 20000 x[5,1] <= 0
 R3540: z[6,1] - 20000 x[6,1] <= 0
 R3541: z[7,1] - 20000 x[7,1] <= 0
 R3542: z[8,1] - 20000 x[8,1] <= 0
 R3543: z[9,1] - 20000 x[9,1] <= 0
 R3544: z[0,2] - 20000 x[0,2] <= 0
 R3545: z[1,2] - 20000 x[1,2] <= 0
 R3546: z[2,2] - 20000 x[2,2] <= 0
 R3547: z[3,2] - 20000 x[3,2] <= 0
 R3548: z[4,2] - 20000 x[4,2] <= 0
 R3549: z[5,2] - 20000 x[5,2] <= 0
 R3550: z[6,2] - 20000 x[6,2] <= 0
 R3551: z[7,2] - 20000 x[7,2] <= 0
 R3552: z[8,2] - 20000 x[8,2] <= 0
 R3553: z[9,2] - 20000 x[9,2] <= 0
 R3554: - Q[0] + z[0,0] - 20000 x[0,0] >= -20000
 R3555: - Q[1] + z[1,0] - 20000 x[1,0] >= -20000
 R3556: - Q[2] + z[2,0] - 20000 x[2,0] >= -20000
 R3557: - Q[3] + z[3,0] - 20000 x[3,0] >= -20000
 R3558: - Q[4] + z[4,0] - 20000 x[4,0] >= -20000
 R3559: - Q[5] + z[5,0] - 20000 x[5,0] >= -20000
 R3560: - Q[6] + z[6,0] - 20000 x[6,0] >= -20000
 R3561: - Q[7] + z[7,0] - 20000 x[7,0] >= -20000
 R3562: - Q[8] + z[8,0] - 20000 x[8,0] >= -20000
 R3563: - Q[9] + z[9,0] - 20000 x[9,0] >= -20000
 R3564: - Q[0] + z[0,1] - 20000 x[0,1] >= -20000
 R3565: - Q[1] + z[1,1] - 20000 x[1,1] >= -20000
 R3566: - Q[2] + z[2,1] - 20000 x[2,1] >= -20000
 R3567: - Q[3] + z[3,1] - 20000 x[3,1] >= -20000
 R3568: - Q[4] + z[4,1] - 20000 x[4,1] >= -20000
 R3569: - Q[5] + z[5,1] - 20000 x[5,1] >= -20000
 R3570: - Q[6] + z[6,1] - 20000 x[6,1] >= -20000
 R3571: - Q[7] + z[7,1] - 20000 x[7,1] >= -20000
 R3572: - Q[8] + z[8,1] - 20000 x[8,1] >= -20000
 R3573: - Q[9] + z[9,1] - 20000 x[9,1] >= -20000
 R3574: - Q[0] + z[0,2] - 20000 x[0,2] >= -20000
 R3575: - Q[1] + z[1,2] - 20000 x[1,2] >= -20000
 R3576: - Q[2] + z[2,2] - 20000 x[2,2] >= -20000
 R3577: - Q[3] + z[3,2] - 20000 x[3,2] >= -20000
 R3578: - Q[4] + z[4,2] - 20000 x[4,2] >= -20000
 R3579: - Q[5] + z[5,2] - 20000 x[5,2] >= -20000
 R3580: - Q[6] + z[6,2] - 20000 x[6,2] >= -20000
 R3581: - Q[7] + z[7,2] - 20000 x[7,2] >= -20000
 R3582: - Q[8] + z[8,2] - 20000 x[8,2] >= -20000
 R3583: - Q[9] + z[9,2] - 20000 x[9,2] >= -20000
 R3584: z[0,0] + z[1,0] + z[2,0] + z[3,0] + z[4,0] + z[5,0] + z[6,0]
   + z[7,0] + z[8,0] + z[9,0] <= 70000
 R3585: z[0,1] + z[1,1] + z[2,1] + z[3,1] + z[4,1] + z[5,1] + z[6,1]
   + z[7,1] + z[8,1] + z[9,1] <= 70000
 R3586: z[0,2] + z[1,2] + z[2,2] + z[3,2] + z[4,2] + z[5,2] + z[6,2]
   + z[7,2] + z[8,2] + z[9,2] <= 70000
 R3587: - b[0,0,1] - 2 b[0,0,2] - b[0,1,1] - 2 b[0,1,2] - b[0,2,1]
   - 2 b[0,2,2] - b[0,3,1] - 2 b[0,3,2] - b[0,4,1] - 2 b[0,4,2] - b[0,5,1]
   - 2 b[0,5,2] - b[0,6,1] - 2 b[0,6,2] - b[0,7,1] - 2 b[0,7,2] - b[0,8,1]
   - 2 b[0,8,2] - b[0,9,1] - 2 b[0,9,2] + tp[0] >= 0
 R3588: - b[1,0,1] - 2 b[1,0,2] - b[1,1,1] - 2 b[1,1,2] - b[1,2,1]
   - 2 b[1,2,2] - b[1,3,1] - 2 b[1,3,2] - b[1,4,1] - 2 b[1,4,2] - b[1,5,1]
   - 2 b[1,5,2] - b[1,6,1] - 2 b[1,6,2] - b[1,7,1] - 2 b[1,7,2] - b[1,8,1]
   - 2 b[1,8,2] - b[1,9,1] - 2 b[1,9,2] + tp[1] >= 0
 R3589: - b[2,0,1] - 2 b[2,0,2] - b[2,1,1] - 2 b[2,1,2] - b[2,2,1]
   - 2 b[2,2,2] - b[2,3,1] - 2 b[2,3,2] - b[2,4,1] - 2 b[2,4,2] - b[2,5,1]
   - 2 b[2,5,2] - b[2,6,1] - 2 b[2,6,2] - b[2,7,1] - 2 b[2,7,2] - b[2,8,1]
   - 2 b[2,8,2] - b[2,9,1] - 2 b[2,9,2] + tp[2] >= 0
 R3590: - b[3,0,1] - 2 b[3,0,2] - b[3,1,1] - 2 b[3,1,2] - b[3,2,1]
   - 2 b[3,2,2] - b[3,3,1] - 2 b[3,3,2] - b[3,4,1] - 2 b[3,4,2] - b[3,5,1]
   - 2 b[3,5,2] - b[3,6,1] - 2 b[3,6,2] - b[3,7,1] - 2 b[3,7,2] - b[3,8,1]
   - 2 b[3,8,2] - b[3,9,1] - 2 b[3,9,2] + tp[3] >= 0
 R3591: - b[4,0,1] - 2 b[4,0,2] - b[4,1,1] - 2 b[4,1,2] - b[4,2,1]
   - 2 b[4,2,2] - b[4,3,1] - 2 b[4,3,2] - b[4,4,1] - 2 b[4,4,2] - b[4,5,1]
   - 2 b[4,5,2] - b[4,6,1] - 2 b[4,6,2] - b[4,7,1] - 2 b[4,7,2] - b[4,8,1]
   - 2 b[4,8,2] - b[4,9,1] - 2 b[4,9,2] + tp[4] >= -1
 R3592: - b[5,0,1] - 2 b[5,0,2] - b[5,1,1] - 2 b[5,1,2] - b[5,2,1]
   - 2 b[5,2,2] - b[5,3,1] - 2 b[5,3,2] - b[5,4,1] - 2 b[5,4,2] - b[5,5,1]
   - 2 b[5,5,2] - b[5,6,1] - 2 b[5,6,2] - b[5,7,1] - 2 b[5,7,2] - b[5,8,1]
   - 2 b[5,8,2] - b[5,9,1] - 2 b[5,9,2] + tp[5] >= -1
 R3593: - b[6,0,1] - 2 b[6,0,2] - b[6,1,1] - 2 b[6,1,2] - b[6,2,1]
   - 2 b[6,2,2] - b[6,3,1] - 2 b[6,3,2] - b[6,4,1] - 2 b[6,4,2] - b[6,5,1]
   - 2 b[6,5,2] - b[6,6,1] - 2 b[6,6,2] - b[6,7,1] - 2 b[6,7,2] - b[6,8,1]
   - 2 b[6,8,2] - b[6,9,1] - 2 b[6,9,2] + tp[6] >= -1
 R3594: - b[7,0,1] - 2 b[7,0,2] - b[7,1,1] - 2 b[7,1,2] - b[7,2,1]
   - 2 b[7,2,2] - b[7,3,1] - 2 b[7,3,2] - b[7,4,1] - 2 b[7,4,2] - b[7,5,1]
   - 2 b[7,5,2] - b[7,6,1] - 2 b[7,6,2] - b[7,7,1] - 2 b[7,7,2] - b[7,8,1]
   - 2 b[7,8,2] - b[7,9,1] - 2 b[7,9,2] + tp[7] >= -1
 R3595: - b[8,0,1] - 2 b[8,0,2] - b[8,1,1] - 2 b[8,1,2] - b[8,2,1]
   - 2 b[8,2,2] - b[8,3,1] - 2 b[8,3,2] - b[8,4,1] - 2 b[8,4,2] - b[8,5,1]
   - 2 b[8,5,2] - b[8,6,1] - 2 b[8,6,2] - b[8,7,1] - 2 b[8,7,2] - b[8,8,1]
   - 2 b[8,8,2] - b[8,9,1] - 2 b[8,9,2] + tp[8] >= -1
 R3596: - b[9,0,1] - 2 b[9,0,2] - b[9,1,1] - 2 b[9,1,2] - b[9,2,1]
   - 2 b[9,2,2] - b[9,3,1] - 2 b[9,3,2] - b[9,4,1] - 2 b[9,4,2] - b[9,5,1]
   - 2 b[9,5,2] - b[9,6,1] - 2 b[9,6,2] - b[9,7,1] - 2 b[9,7,2] - b[9,8,1]
   - 2 b[9,8,2] - b[9,9,1] - 2 b[9,9,2] + tp[9] >= -1
 R3597: b[0,0,1] + 2 b[0,0,2] + b[0,1,1] + 2 b[0,1,2] + b[0,2,1]
   + 2 b[0,2,2] + b[0,3,1] + 2 b[0,3,2] + b[0,4,1] + 2 b[0,4,2] + b[0,5,1]
   + 2 b[0,5,2] + b[0,6,1] + 2 b[0,6,2] + b[0,7,1] + 2 b[0,7,2] + b[0,8,1]
   + 2 b[0,8,2] + b[0,9,1] + 2 b[0,9,2] + tn[0] >= 0
 R3598: b[1,0,1] + 2 b[1,0,2] + b[1,1,1] + 2 b[1,1,2] + b[1,2,1]
   + 2 b[1,2,2] + b[1,3,1] + 2 b[1,3,2] + b[1,4,1] + 2 b[1,4,2] + b[1,5,1]
   + 2 b[1,5,2] + b[1,6,1] + 2 b[1,6,2] + b[1,7,1] + 2 b[1,7,2] + b[1,8,1]
   + 2 b[1,8,2] + b[1,9,1] + 2 b[1,9,2] + tn[1] >= 0
 R3599: b[2,0,1] + 2 b[2,0,2] + b[2,1,1] + 2 b[2,1,2] + b[2,2,1]
   + 2 b[2,2,2] + b[2,3,1] + 2 b[2,3,2] + b[2,4,1] + 2 b[2,4,2] + b[2,5,1]
   + 2 b[2,5,2] + b[2,6,1] + 2 b[2,6,2] + b[2,7,1] + 2 b[2,7,2] + b[2,8,1]
   + 2 b[2,8,2] + b[2,9,1] + 2 b[2,9,2] + tn[2] >= 0
 R3600: b[3,0,1] + 2 b[3,0,2] + b[3,1,1] + 2 b[3,1,2] + b[3,2,1]
   + 2 b[3,2,2] + b[3,3,1] + 2 b[3,3,2] + b[3,4,1] + 2 b[3,4,2] + b[3,5,1]
   + 2 b[3,5,2] + b[3,6,1] + 2 b[3,6,2] + b[3,7,1] + 2 b[3,7,2] + b[3,8,1]
   + 2 b[3,8,2] + b[3,9,1] + 2 b[3,9,2] + tn[3] >= 0
 R3601: b[4,0,1] + 2 b[4,0,2] + b[4,1,1] + 2 b[4,1,2] + b[4,2,1]
   + 2 b[4,2,2] + b[4,3,1] + 2 b[4,3,2] + b[4,4,1] + 2 b[4,4,2] + b[4,5,1]
   + 2 b[4,5,2] + b[4,6,1] + 2 b[4,6,2] + b[4,7,1] + 2 b[4,7,2] + b[4,8,1]
   + 2 b[4,8,2] + b[4,9,1] + 2 b[4,9,2] + tn[4] >= 1
 R3602: b[5,0,1] + 2 b[5,0,2] + b[5,1,1] + 2 b[5,1,2] + b[5,2,1]
   + 2 b[5,2,2] + b[5,3,1] + 2 b[5,3,2] + b[5,4,1] + 2 b[5,4,2] + b[5,5,1]
   + 2 b[5,5,2] + b[5,6,1] + 2 b[5,6,2] + b[5,7,1] + 2 b[5,7,2] + b[5,8,1]
   + 2 b[5,8,2] + b[5,9,1] + 2 b[5,9,2] + tn[5] >= 1
 R3603: b[6,0,1] + 2 b[6,0,2] + b[6,1,1] + 2 b[6,1,2] + b[6,2,1]
   + 2 b[6,2,2] + b[6,3,1] + 2 b[6,3,2] + b[6,4,1] + 2 b[6,4,2] + b[6,5,1]
   + 2 b[6,5,2] + b[6,6,1] + 2 b[6,6,2] + b[6,7,1] + 2 b[6,7,2] + b[6,8,1]
   + 2 b[6,8,2] + b[6,9,1] + 2 b[6,9,2] + tn[6] >= 1
 R3604: b[7,0,1] + 2 b[7,0,2] + b[7,1,1] + 2 b[7,1,2] + b[7,2,1]
   + 2 b[7,2,2] + b[7,3,1] + 2 b[7,3,2] + b[7,4,1] + 2 b[7,4,2] + b[7,5,1]
   + 2 b[7,5,2] + b[7,6,1] + 2 b[7,6,2] + b[7,7,1] + 2 b[7,7,2] + b[7,8,1]
   + 2 b[7,8,2] + b[7,9,1] + 2 b[7,9,2] + tn[7] >= 1
 R3605: b[8,0,1] + 2 b[8,0,2] + b[8,1,1] + 2 b[8,1,2] + b[8,2,1]
   + 2 b[8,2,2] + b[8,3,1] + 2 b[8,3,2] + b[8,4,1] + 2 b[8,4,2] + b[8,5,1]
   + 2 b[8,5,2] + b[8,6,1] + 2 b[8,6,2] + b[8,7,1] + 2 b[8,7,2] + b[8,8,1]
   + 2 b[8,8,2] + b[8,9,1] + 2 b[8,9,2] + tn[8] >= 1
 R3606: b[9,0,1] + 2 b[9,0,2] + b[9,1,1] + 2 b[9,1,2] + b[9,2,1]
   + 2 b[9,2,2] + b[9,3,1] + 2 b[9,3,2] + b[9,4,1] + 2 b[9,4,2] + b[9,5,1]
   + 2 b[9,5,2] + b[9,6,1] + 2 b[9,6,2] + b[9,7,1] + 2 b[9,7,2] + b[9,8,1]
   + 2 b[9,8,2] + b[9,9,1] + 2 b[9,9,2] + tn[9] >= 1
 R3607: tp[0] <= 1
 R3608: tp[1] <= 1
 R3609: tp[2] <= 1
 R3610: tp[3] <= 1
 R3611: tp[4] <= 1
 R3612: tp[5] <= 1
 R3613: tp[6] <= 1
 R3614: tp[7] <= 1
 R3615: tp[8] <= 1
 R3616: tp[9] <= 1
 R3617: tn[0] <= 2
 R3618: tn[1] <= 2
 R3619: tn[2] <= 2
 R3620: tn[3] <= 2
 R3621: tn[4] <= 2
 R3622: tn[5] <= 2
 R3623: tn[6] <= 2
 R3624: tn[7] <= 2
 R3625: tn[8] <= 2
 R3626: tn[9] <= 2
Bounds
Binaries
 alpha[0,0] alpha[0,2] alpha[0,6] alpha[0,8] alpha[0,9] alpha[1,1]
 alpha[1,3] alpha[1,7] alpha[2,2] alpha[2,6] alpha[2,8] alpha[2,9]
 alpha[3,3] alpha[3,7] alpha[4,4] alpha[5,5] alpha[6,6] alpha[6,8]
 alpha[6,9] alpha[7,7] alpha[8,8] alpha[8,9] alpha[9,9] beta[0,0] beta[0,1]
 beta[0,2] beta[0,3] beta[0,4] beta[0,5] beta[0,6] beta[0,7] beta[0,8]
 beta[0,9] beta[1,0] beta[1,1] beta[1,2] beta[1,3] beta[1,4] beta[1,5]
 beta[1,6] beta[1,7] beta[1,8] beta[1,9] beta[2,0] beta[2,1] beta[2,2]
 beta[2,3] beta[2,4] beta[2,5] beta[2,6] beta[2,7] beta[2,8] beta[2,9]
 beta[3,0] beta[3,1] beta[3,2] beta[3,3] beta[3,4] beta[3,5] beta[3,6]
 beta[3,7] beta[3,8] beta[3,9] beta[4,0] beta[4,1] beta[4,2] beta[4,3]
 beta[4,4] beta[4,5] beta[4,6] beta[4,7] beta[4,8] beta[4,9] beta[5,0]
 beta[5,1] beta[5,2] beta[5,3] beta[5,4] beta[5,5] beta[5,6] beta[5,7]
 beta[5,8] beta[5,9] beta[6,0] beta[6,1] beta[6,2] beta[6,3] beta[6,4]
 beta[6,5] beta[6,6] beta[6,7] beta[6,8] beta[6,9] beta[7,0] beta[7,1]
 beta[7,2] beta[7,3] beta[7,4] beta[7,5] beta[7,6] beta[7,7] beta[7,8]
 beta[7,9] beta[8,0] beta[8,1] beta[8,2] beta[8,3] beta[8,4] beta[8,5]
 beta[8,6] beta[8,7] beta[8,8] beta[8,9] beta[9,0] beta[9,1] beta[9,2]
 beta[9,3] beta[9,4] beta[9,5] beta[9,6] beta[9,7] beta[9,8] beta[9,9]
 gamma[0,0] gamma[0,1] gamma[0,2] gamma[0,3] gamma[0,4] gamma[0,5]
 gamma[0,6] gamma[0,7] gamma[0,8] gamma[0,9] gamma[1,1] gamma[1,2]
 gamma[1,3] gamma[1,4] gamma[1,5] gamma[1,6] gamma[1,7] gamma[1,8]
 gamma[1,9] gamma[2,2] gamma[2,3] gamma[2,4] gamma[2,5] gamma[2,6]
 gamma[2,7] gamma[2,8] gamma[2,9] gamma[3,3] gamma[3,4] gamma[3,5]
 gamma[3,6] gamma[3,7] gamma[3,8] gamma[3,9] gamma[4,4] gamma[4,5]
 gamma[4,6] gamma[4,7] gamma[4,8] gamma[4,9] gamma[5,5] gamma[5,6]
 gamma[5,7] gamma[5,8] gamma[5,9] gamma[6,6] gamma[6,7] gamma[6,8]
 gamma[6,9] gamma[7,7] gamma[7,8] gamma[7,9] gamma[8,8] gamma[8,9]
 gamma[9,9] delta[0,1,0] delta[0,2,0] delta[0,2,1] delta[0,3,0]
 delta[0,3,1] delta[0,3,2] delta[0,4,0] delta[0,4,1] delta[0,4,2]
 delta[0,4,3] delta[0,5,0] delta[0,5,1] delta[0,5,2] delta[0,5,3]
 delta[0,5,4] delta[0,6,0] delta[0,6,1] delta[0,6,2] delta[0,6,3]
 delta[0,6,4] delta[0,6,5] delta[0,7,0] delta[0,7,1] delta[0,7,2]
 delta[0,7,3] delta[0,7,4] delta[0,7,5] delta[0,7,6] delta[0,8,0]
 delta[0,8,1] delta[0,8,2] delta[0,8,3] delta[0,8,4] delta[0,8,5]
 delta[0,8,6] delta[0,8,7] delta[0,9,0] delta[0,9,1] delta[0,9,2]
 delta[0,9,3] delta[0,9,4] delta[0,9,5] delta[0,9,6] delta[0,9,7]
 delta[0,9,8] delta[1,2,1] delta[1,3,1] delta[1,3,2] delta[1,4,1]
 delta[1,4,2] delta[1,4,3] delta[1,5,1] delta[1,5,2] delta[1,5,3]
 delta[1,5,4] delta[1,6,1] delta[1,6,2] delta[1,6,3] delta[1,6,4]
 delta[1,6,5] delta[1,7,1] delta[1,7,2] delta[1,7,3] delta[1,7,4]
 delta[1,7,5] delta[1,7,6] delta[1,8,1] delta[1,8,2] delta[1,8,3]
 delta[1,8,4] delta[1,8,5] delta[1,8,6] delta[1,8,7] delta[1,9,1]
 delta[1,9,2] delta[1,9,3] delta[1,9,4] delta[1,9,5] delta[1,9,6]
 delta[1,9,7] delta[1,9,8] delta[2,3,2] delta[2,4,2] delta[2,4,3]
 delta[2,5,2] delta[2,5,3] delta[2,5,4] delta[2,6,2] delta[2,6,3]
 delta[2,6,4] delta[2,6,5] delta[2,7,2] delta[2,7,3] delta[2,7,4]
 delta[2,7,5] delta[2,7,6] delta[2,8,2] delta[2,8,3] delta[2,8,4]
 delta[2,8,5] delta[2,8,6] delta[2,8,7] delta[2,9,2] delta[2,9,3]
 delta[2,9,4] delta[2,9,5] delta[2,9,6] delta[2,9,7] delta[2,9,8]
 delta[3,4,3] delta[3,5,3] delta[3,5,4] delta[3,6,3] delta[3,6,4]
 delta[3,6,5] delta[3,7,3] delta[3,7,4] delta[3,7,5] delta[3,7,6]
 delta[3,8,3] delta[3,8,4] delta[3,8,5] delta[3,8,6] delta[3,8,7]
 delta[3,9,3] delta[3,9,4] delta[3,9,5] delta[3,9,6] delta[3,9,7]
 delta[3,9,8] delta[4,5,4] delta[4,6,4] delta[4,6,5] delta[4,7,4]
 delta[4,7,5] delta[4,7,6] delta[4,8,4] delta[4,8,5] delta[4,8,6]
 delta[4,8,7] delta[4,9,4] delta[4,9,5] delta[4,9,6] delta[4,9,7]
 delta[4,9,8] delta[5,6,5] delta[5,7,5] delta[5,7,6] delta[5,8,5]
 delta[5,8,6] delta[5,8,7] delta[5,9,5] delta[5,9,6] delta[5,9,7]
 delta[5,9,8] delta[6,7,6] delta[6,8,6] delta[6,8,7] delta[6,9,6]
 delta[6,9,7] delta[6,9,8] delta[7,8,7] delta[7,9,7] delta[7,9,8]
 delta[8,9,8] a[0,0,0,0] a[0,1,0,0] a[1,1,0,0] a[0,2,0,0] a[1,2,0,0]
 a[2,2,0,0] a[0,3,0,0] a[1,3,0,0] a[2,3,0,0] a[3,3,0,0] a[0,4,0,0]
 a[1,4,0,0] a[2,4,0,0] a[3,4,0,0] a[4,4,0,0] a[0,5,0,0] a[1,5,0,0]
 a[2,5,0,0] a[3,5,0,0] a[4,5,0,0] a[5,5,0,0] a[0,6,0,0] a[1,6,0,0]
 a[2,6,0,0] a[3,6,0,0] a[4,6,0,0] a[5,6,0,0] a[6,6,0,0] a[0,7,0,0]
 a[1,7,0,0] a[2,7,0,0] a[3,7,0,0] a[4,7,0,0] a[5,7,0,0] a[6,7,0,0]
 a[7,7,0,0] a[0,8,0,0] a[1,8,0,0] a[2,8,0,0] a[3,8,0,0] a[4,8,0,0]
 a[5,8,0,0] a[6,8,0,0] a[7,8,0,0] a[8,8,0,0] a[0,9,0,0] a[1,9,0,0]
 a[2,9,0,0] a[3,9,0,0] a[4,9,0,0] a[5,9,0,0] a[6,9,0,0] a[7,9,0,0]
 a[8,9,0,0] a[9,9,0,0] a[0,0,0,2] a[0,1,0,2] a[1,1,0,2] a[0,2,0,2]
 a[1,2,0,2] a[2,2,0,2] a[0,3,0,2] a[1,3,0,2] a[2,3,0,2] a[3,3,0,2]
 a[0,4,0,2] a[1,4,0,2] a[2,4,0,2] a[3,4,0,2] a[4,4,0,2] a[0,5,0,2]
 a[1,5,0,2] a[2,5,0,2] a[3,5,0,2] a[4,5,0,2] a[5,5,0,2] a[0,6,0,2]
 a[1,6,0,2] a[2,6,0,2] a[3,6,0,2] a[4,6,0,2] a[5,6,0,2] a[6,6,0,2]
 a[0,7,0,2] a[1,7,0,2] a[2,7,0,2] a[3,7,0,2] a[4,7,0,2] a[5,7,0,2]
 a[6,7,0,2] a[7,7,0,2] a[0,8,0,2] a[1,8,0,2] a[2,8,0,2] a[3,8,0,2]
 a[4,8,0,2] a[5,8,0,2] a[6,8,0,2] a[7,8,0,2] a[8,8,0,2] a[0,9,0,2]
 a[1,9,0,2] a[2,9,0,2] a[3,9,0,2] a[4,9,0,2] a[5,9,0,2] a[6,9,0,2]
 a[7,9,0,2] a[8,9,0,2] a[9,9,0,2] a[0,0,0,6] a[0,1,0,6] a[1,1,0,6]
 a[0,2,0,6] a[1,2,0,6] a[2,2,0,6] a[0,3,0,6] a[1,3,0,6] a[2,3,0,6]
 a[3,3,0,6] a[0,4,0,6] a[1,4,0,6] a[2,4,0,6] a[3,4,0,6] a[4,4,0,6]
 a[0,5,0,6] a[1,5,0,6] a[2,5,0,6] a[3,5,0,6] a[4,5,0,6] a[5,5,0,6]
 a[0,6,0,6] a[1,6,0,6] a[2,6,0,6] a[3,6,0,6] a[4,6,0,6] a[5,6,0,6]
 a[6,6,0,6] a[0,7,0,6] a[1,7,0,6] a[2,7,0,6] a[3,7,0,6] a[4,7,0,6]
 a[5,7,0,6] a[6,7,0,6] a[7,7,0,6] a[0,8,0,6] a[1,8,0,6] a[2,8,0,6]
 a[3,8,0,6] a[4,8,0,6] a[5,8,0,6] a[6,8,0,6] a[7,8,0,6] a[8,8,0,6]
 a[0,9,0,6] a[1,9,0,6] a[2,9,0,6] a[3,9,0,6] a[4,9,0,6] a[5,9,0,6]
 a[6,9,0,6] a[7,9,0,6] a[8,9,0,6] a[9,9,0,6] a[0,0,0,8] a[0,1,0,8]
 a[1,1,0,8] a[0,2,0,8] a[1,2,0,8] a[2,2,0,8] a[0,3,0,8] a[1,3,0,8]
 a[2,3,0,8] a[3,3,0,8] a[0,4,0,8] a[1,4,0,8] a[2,4,0,8] a[3,4,0,8]
 a[4,4,0,8] a[0,5,0,8] a[1,5,0,8] a[2,5,0,8] a[3,5,0,8] a[4,5,0,8]
 a[5,5,0,8] a[0,6,0,8] a[1,6,0,8] a[2,6,0,8] a[3,6,0,8] a[4,6,0,8]
 a[5,6,0,8] a[6,6,0,8] a[0,7,0,8] a[1,7,0,8] a[2,7,0,8] a[3,7,0,8]
 a[4,7,0,8] a[5,7,0,8] a[6,7,0,8] a[7,7,0,8] a[0,8,0,8] a[1,8,0,8]
 a[2,8,0,8] a[3,8,0,8] a[4,8,0,8] a[5,8,0,8] a[6,8,0,8] a[7,8,0,8]
 a[8,8,0,8] a[0,9,0,8] a[1,9,0,8] a[2,9,0,8] a[3,9,0,8] a[4,9,0,8]
 a[5,9,0,8] a[6,9,0,8] a[7,9,0,8] a[8,9,0,8] a[9,9,0,8] a[0,0,0,9]
 a[0,1,0,9] a[1,1,0,9] a[0,2,0,9] a[1,2,0,9] a[2,2,0,9] a[0,3,0,9]
 a[1,3,0,9] a[2,3,0,9] a[3,3,0,9] a[0,4,0,9] a[1,4,0,9] a[2,4,0,9]
 a[3,4,0,9] a[4,4,0,9] a[0,5,0,9] a[1,5,0,9] a[2,5,0,9] a[3,5,0,9]
 a[4,5,0,9] a[5,5,0,9] a[0,6,0,9] a[1,6,0,9] a[2,6,0,9] a[3,6,0,9]
 a[4,6,0,9] a[5,6,0,9] a[6,6,0,9] a[0,7,0,9] a[1,7,0,9] a[2,7,0,9]
 a[3,7,0,9] a[4,7,0,9] a[5,7,0,9] a[6,7,0,9] a[7,7,0,9] a[0,8,0,9]
 a[1,8,0,9] a[2,8,0,9] a[3,8,0,9] a[4,8,0,9] a[5,8,0,9] a[6,8,0,9]
 a[7,8,0,9] a[8,8,0,9] a[0,9,0,9] a[1,9,0,9] a[2,9,0,9] a[3,9,0,9]
 a[4,9,0,9] a[5,9,0,9] a[6,9,0,9] a[7,9,0,9] a[8,9,0,9] a[9,9,0,9]
 a[0,0,1,1] a[0,1,1,1] a[1,1,1,1] a[0,2,1,1] a[1,2,1,1] a[2,2,1,1]
 a[0,3,1,1] a[1,3,1,1] a[2,3,1,1] a[3,3,1,1] a[0,4,1,1] a[1,4,1,1]
 a[2,4,1,1] a[3,4,1,1] a[4,4,1,1] a[0,5,1,1] a[1,5,1,1] a[2,5,1,1]
 a[3,5,1,1] a[4,5,1,1] a[5,5,1,1] a[0,6,1,1] a[1,6,1,1] a[2,6,1,1]
 a[3,6,1,1] a[4,6,1,1] a[5,6,1,1] a[6,6,1,1] a[0,7,1,1] a[1,7,1,1]
 a[2,7,1,1] a[3,7,1,1] a[4,7,1,1] a[5,7,1,1] a[6,7,1,1] a[7,7,1,1]
 a[0,8,1,1] a[1,8,1,1] a[2,8,1,1] a[3,8,1,1] a[4,8,1,1] a[5,8,1,1]
 a[6,8,1,1] a[7,8,1,1] a[8,8,1,1] a[0,9,1,1] a[1,9,1,1] a[2,9,1,1]
 a[3,9,1,1] a[4,9,1,1] a[5,9,1,1] a[6,9,1,1] a[7,9,1,1] a[8,9,1,1]
 a[9,9,1,1] a[0,0,1,3] a[0,1,1,3] a[1,1,1,3] a[0,2,1,3] a[1,2,1,3]
 a[2,2,1,3] a[0,3,1,3] a[1,3,1,3] a[2,3,1,3] a[3,3,1,3] a[0,4,1,3]
 a[1,4,1,3] a[2,4,1,3] a[3,4,1,3] a[4,4,1,3] a[0,5,1,3] a[1,5,1,3]
 a[2,5,1,3] a[3,5,1,3] a[4,5,1,3] a[5,5,1,3] a[0,6,1,3] a[1,6,1,3]
 a[2,6,1,3] a[3,6,1,3] a[4,6,1,3] a[5,6,1,3] a[6,6,1,3] a[0,7,1,3]
 a[1,7,1,3] a[2,7,1,3] a[3,7,1,3] a[4,7,1,3] a[5,7,1,3] a[6,7,1,3]
 a[7,7,1,3] a[0,8,1,3] a[1,8,1,3] a[2,8,1,3] a[3,8,1,3] a[4,8,1,3]
 a[5,8,1,3] a[6,8,1,3] a[7,8,1,3] a[8,8,1,3] a[0,9,1,3] a[1,9,1,3]
 a[2,9,1,3] a[3,9,1,3] a[4,9,1,3] a[5,9,1,3] a[6,9,1,3] a[7,9,1,3]
 a[8,9,1,3] a[9,9,1,3] a[0,0,1,7] a[0,1,1,7] a[1,1,1,7] a[0,2,1,7]
 a[1,2,1,7] a[2,2,1,7] a[0,3,1,7] a[1,3,1,7] a[2,3,1,7] a[3,3,1,7]
 a[0,4,1,7] a[1,4,1,7] a[2,4,1,7] a[3,4,1,7] a[4,4,1,7] a[0,5,1,7]
 a[1,5,1,7] a[2,5,1,7] a[3,5,1,7] a[4,5,1,7] a[5,5,1,7] a[0,6,1,7]
 a[1,6,1,7] a[2,6,1,7] a[3,6,1,7] a[4,6,1,7] a[5,6,1,7] a[6,6,1,7]
 a[0,7,1,7] a[1,7,1,7] a[2,7,1,7] a[3,7,1,7] a[4,7,1,7] a[5,7,1,7]
 a[6,7,1,7] a[7,7,1,7] a[0,8,1,7] a[1,8,1,7] a[2,8,1,7] a[3,8,1,7]
 a[4,8,1,7] a[5,8,1,7] a[6,8,1,7] a[7,8,1,7] a[8,8,1,7] a[0,9,1,7]
 a[1,9,1,7] a[2,9,1,7] a[3,9,1,7] a[4,9,1,7] a[5,9,1,7] a[6,9,1,7]
 a[7,9,1,7] a[8,9,1,7] a[9,9,1,7] a[0,0,2,2] a[0,1,2,2] a[1,1,2,2]
 a[0,2,2,2] a[1,2,2,2] a[2,2,2,2] a[0,3,2,2] a[1,3,2,2] a[2,3,2,2]
 a[3,3,2,2] a[0,4,2,2] a[1,4,2,2] a[2,4,2,2] a[3,4,2,2] a[4,4,2,2]
 a[0,5,2,2] a[1,5,2,2] a[2,5,2,2] a[3,5,2,2] a[4,5,2,2] a[5,5,2,2]
 a[0,6,2,2] a[1,6,2,2] a[2,6,2,2] a[3,6,2,2] a[4,6,2,2] a[5,6,2,2]
 a[6,6,2,2] a[0,7,2,2] a[1,7,2,2] a[2,7,2,2] a[3,7,2,2] a[4,7,2,2]
 a[5,7,2,2] a[6,7,2,2] a[7,7,2,2] a[0,8,2,2] a[1,8,2,2] a[2,8,2,2]
 a[3,8,2,2] a[4,8,2,2] a[5,8,2,2] a[6,8,2,2] a[7,8,2,2] a[8,8,2,2]
 a[0,9,2,2] a[1,9,2,2] a[2,9,2,2] a[3,9,2,2] a[4,9,2,2] a[5,9,2,2]
 a[6,9,2,2] a[7,9,2,2] a[8,9,2,2] a[9,9,2,2] a[0,0,2,6] a[0,1,2,6]
 a[1,1,2,6] a[0,2,2,6] a[1,2,2,6] a[2,2,2,6] a[0,3,2,6] a[1,3,2,6]
 a[2,3,2,6] a[3,3,2,6] a[0,4,2,6] a[1,4,2,6] a[2,4,2,6] a[3,4,2,6]
 a[4,4,2,6] a[0,5,2,6] a[1,5,2,6] a[2,5,2,6] a[3,5,2,6] a[4,5,2,6]
 a[5,5,2,6] a[0,6,2,6] a[1,6,2,6] a[2,6,2,6] a[3,6,2,6] a[4,6,2,6]
 a[5,6,2,6] a[6,6,2,6] a[0,7,2,6] a[1,7,2,6] a[2,7,2,6] a[3,7,2,6]
 a[4,7,2,6] a[5,7,2,6] a[6,7,2,6] a[7,7,2,6] a[0,8,2,6] a[1,8,2,6]
 a[2,8,2,6] a[3,8,2,6] a[4,8,2,6] a[5,8,2,6] a[6,8,2,6] a[7,8,2,6]
 a[8,8,2,6] a[0,9,2,6] a[1,9,2,6] a[2,9,2,6] a[3,9,2,6] a[4,9,2,6]
 a[5,9,2,6] a[6,9,2,6] a[7,9,2,6] a[8,9,2,6] a[9,9,2,6] a[0,0,2,8]
 a[0,1,2,8] a[1,1,2,8] a[0,2,2,8] a[1,2,2,8] a[2,2,2,8] a[0,3,2,8]
 a[1,3,2,8] a[2,3,2,8] a[3,3,2,8] a[0,4,2,8] a[1,4,2,8] a[2,4,2,8]
 a[3,4,2,8] a[4,4,2,8] a[0,5,2,8] a[1,5,2,8] a[2,5,2,8] a[3,5,2,8]
 a[4,5,2,8] a[5,5,2,8] a[0,6,2,8] a[1,6,2,8] a[2,6,2,8] a[3,6,2,8]
 a[4,6,2,8] a[5,6,2,8] a[6,6,2,8] a[0,7,2,8] a[1,7,2,8] a[2,7,2,8]
 a[3,7,2,8] a[4,7,2,8] a[5,7,2,8] a[6,7,2,8] a[7,7,2,8] a[0,8,2,8]
 a[1,8,2,8] a[2,8,2,8] a[3,8,2,8] a[4,8,2,8] a[5,8,2,8] a[6,8,2,8]
 a[7,8,2,8] a[8,8,2,8] a[0,9,2,8] a[1,9,2,8] a[2,9,2,8] a[3,9,2,8]
 a[4,9,2,8] a[5,9,2,8] a[6,9,2,8] a[7,9,2,8] a[8,9,2,8] a[9,9,2,8]
 a[0,0,2,9] a[0,1,2,9] a[1,1,2,9] a[0,2,2,9] a[1,2,2,9] a[2,2,2,9]
 a[0,3,2,9] a[1,3,2,9] a[2,3,2,9] a[3,3,2,9] a[0,4,2,9] a[1,4,2,9]
 a[2,4,2,9] a[3,4,2,9] a[4,4,2,9] a[0,5,2,9] a[1,5,2,9] a[2,5,2,9]
 a[3,5,2,9] a[4,5,2,9] a[5,5,2,9] a[0,6,2,9] a[1,6,2,9] a[2,6,2,9]
 a[3,6,2,9] a[4,6,2,9] a[5,6,2,9] a[6,6,2,9] a[0,7,2,9] a[1,7,2,9]
 a[2,7,2,9] a[3,7,2,9] a[4,7,2,9] a[5,7,2,9] a[6,7,2,9] a[7,7,2,9]
 a[0,8,2,9] a[1,8,2,9] a[2,8,2,9] a[3,8,2,9] a[4,8,2,9] a[5,8,2,9]
 a[6,8,2,9] a[7,8,2,9] a[8,8,2,9] a[0,9,2,9] a[1,9,2,9] a[2,9,2,9]
 a[3,9,2,9] a[4,9,2,9] a[5,9,2,9] a[6,9,2,9] a[7,9,2,9] a[8,9,2,9]
 a[9,9,2,9] a[0,0,3,3] a[0,1,3,3] a[1,1,3,3] a[0,2,3,3] a[1,2,3,3]
 a[2,2,3,3] a[0,3,3,3] a[1,3,3,3] a[2,3,3,3] a[3,3,3,3] a[0,4,3,3]
 a[1,4,3,3] a[2,4,3,3] a[3,4,3,3] a[4,4,3,3] a[0,5,3,3] a[1,5,3,3]
 a[2,5,3,3] a[3,5,3,3] a[4,5,3,3] a[5,5,3,3] a[0,6,3,3] a[1,6,3,3]
 a[2,6,3,3] a[3,6,3,3] a[4,6,3,3] a[5,6,3,3] a[6,6,3,3] a[0,7,3,3]
 a[1,7,3,3] a[2,7,3,3] a[3,7,3,3] a[4,7,3,3] a[5,7,3,3] a[6,7,3,3]
 a[7,7,3,3] a[0,8,3,3] a[1,8,3,3] a[2,8,3,3] a[3,8,3,3] a[4,8,3,3]
 a[5,8,3,3] a[6,8,3,3] a[7,8,3,3] a[8,8,3,3] a[0,9,3,3] a[1,9,3,3]
 a[2,9,3,3] a[3,9,3,3] a[4,9,3,3] a[5,9,3,3] a[6,9,3,3] a[7,9,3,3]
 a[8,9,3,3] a[9,9,3,3] a[0,0,3,7] a[0,1,3,7] a[1,1,3,7] a[0,2,3,7]
 a[1,2,3,7] a[2,2,3,7] a[0,3,3,7] a[1,3,3,7] a[2,3,3,7] a[3,3,3,7]
 a[0,4,3,7] a[1,4,3,7] a[2,4,3,7] a[3,4,3,7] a[4,4,3,7] a[0,5,3,7]
 a[1,5,3,7] a[2,5,3,7] a[3,5,3,7] a[4,5,3,7] a[5,5,3,7] a[0,6,3,7]
 a[1,6,3,7] a[2,6,3,7] a[3,6,3,7] a[4,6,3,7] a[5,6,3,7] a[6,6,3,7]
 a[0,7,3,7] a[1,7,3,7] a[2,7,3,7] a[3,7,3,7] a[4,7,3,7] a[5,7,3,7]
 a[6,7,3,7] a[7,7,3,7] a[0,8,3,7] a[1,8,3,7] a[2,8,3,7] a[3,8,3,7]
 a[4,8,3,7] a[5,8,3,7] a[6,8,3,7] a[7,8,3,7] a[8,8,3,7] a[0,9,3,7]
 a[1,9,3,7] a[2,9,3,7] a[3,9,3,7] a[4,9,3,7] a[5,9,3,7] a[6,9,3,7]
 a[7,9,3,7] a[8,9,3,7] a[9,9,3,7] a[0,0,4,4] a[0,1,4,4] a[1,1,4,4]
 a[0,2,4,4] a[1,2,4,4] a[2,2,4,4] a[0,3,4,4] a[1,3,4,4] a[2,3,4,4]
 a[3,3,4,4] a[0,4,4,4] a[1,4,4,4] a[2,4,4,4] a[3,4,4,4] a[4,4,4,4]
 a[0,5,4,4] a[1,5,4,4] a[2,5,4,4] a[3,5,4,4] a[4,5,4,4] a[5,5,4,4]
 a[0,6,4,4] a[1,6,4,4] a[2,6,4,4] a[3,6,4,4] a[4,6,4,4] a[5,6,4,4]
 a[6,6,4,4] a[0,7,4,4] a[1,7,4,4] a[2,7,4,4] a[3,7,4,4] a[4,7,4,4]
 a[5,7,4,4] a[6,7,4,4] a[7,7,4,4] a[0,8,4,4] a[1,8,4,4] a[2,8,4,4]
 a[3,8,4,4] a[4,8,4,4] a[5,8,4,4] a[6,8,4,4] a[7,8,4,4] a[8,8,4,4]
 a[0,9,4,4] a[1,9,4,4] a[2,9,4,4] a[3,9,4,4] a[4,9,4,4] a[5,9,4,4]
 a[6,9,4,4] a[7,9,4,4] a[8,9,4,4] a[9,9,4,4] a[0,0,5,5] a[0,1,5,5]
 a[1,1,5,5] a[0,2,5,5] a[1,2,5,5] a[2,2,5,5] a[0,3,5,5] a[1,3,5,5]
 a[2,3,5,5] a[3,3,5,5] a[0,4,5,5] a[1,4,5,5] a[2,4,5,5] a[3,4,5,5]
 a[4,4,5,5] a[0,5,5,5] a[1,5,5,5] a[2,5,5,5] a[3,5,5,5] a[4,5,5,5]
 a[5,5,5,5] a[0,6,5,5] a[1,6,5,5] a[2,6,5,5] a[3,6,5,5] a[4,6,5,5]
 a[5,6,5,5] a[6,6,5,5] a[0,7,5,5] a[1,7,5,5] a[2,7,5,5] a[3,7,5,5]
 a[4,7,5,5] a[5,7,5,5] a[6,7,5,5] a[7,7,5,5] a[0,8,5,5] a[1,8,5,5]
 a[2,8,5,5] a[3,8,5,5] a[4,8,5,5] a[5,8,5,5] a[6,8,5,5] a[7,8,5,5]
 a[8,8,5,5] a[0,9,5,5] a[1,9,5,5] a[2,9,5,5] a[3,9,5,5] a[4,9,5,5]
 a[5,9,5,5] a[6,9,5,5] a[7,9,5,5] a[8,9,5,5] a[9,9,5,5] a[0,0,6,6]
 a[0,1,6,6] a[1,1,6,6] a[0,2,6,6] a[1,2,6,6] a[2,2,6,6] a[0,3,6,6]
 a[1,3,6,6] a[2,3,6,6] a[3,3,6,6] a[0,4,6,6] a[1,4,6,6] a[2,4,6,6]
 a[3,4,6,6] a[4,4,6,6] a[0,5,6,6] a[1,5,6,6] a[2,5,6,6] a[3,5,6,6]
 a[4,5,6,6] a[5,5,6,6] a[0,6,6,6] a[1,6,6,6] a[2,6,6,6] a[3,6,6,6]
 a[4,6,6,6] a[5,6,6,6] a[6,6,6,6] a[0,7,6,6] a[1,7,6,6] a[2,7,6,6]
 a[3,7,6,6] a[4,7,6,6] a[5,7,6,6] a[6,7,6,6] a[7,7,6,6] a[0,8,6,6]
 a[1,8,6,6] a[2,8,6,6] a[3,8,6,6] a[4,8,6,6] a[5,8,6,6] a[6,8,6,6]
 a[7,8,6,6] a[8,8,6,6] a[0,9,6,6] a[1,9,6,6] a[2,9,6,6] a[3,9,6,6]
 a[4,9,6,6] a[5,9,6,6] a[6,9,6,6] a[7,9,6,6] a[8,9,6,6] a[9,9,6,6]
 a[0,0,6,8] a[0,1,6,8] a[1,1,6,8] a[0,2,6,8] a[1,2,6,8] a[2,2,6,8]
 a[0,3,6,8] a[1,3,6,8] a[2,3,6,8] a[3,3,6,8] a[0,4,6,8] a[1,4,6,8]
 a[2,4,6,8] a[3,4,6,8] a[4,4,6,8] a[0,5,6,8] a[1,5,6,8] a[2,5,6,8]
 a[3,5,6,8] a[4,5,6,8] a[5,5,6,8] a[0,6,6,8] a[1,6,6,8] a[2,6,6,8]
 a[3,6,6,8] a[4,6,6,8] a[5,6,6,8] a[6,6,6,8] a[0,7,6,8] a[1,7,6,8]
 a[2,7,6,8] a[3,7,6,8] a[4,7,6,8] a[5,7,6,8] a[6,7,6,8] a[7,7,6,8]
 a[0,8,6,8] a[1,8,6,8] a[2,8,6,8] a[3,8,6,8] a[4,8,6,8] a[5,8,6,8]
 a[6,8,6,8] a[7,8,6,8] a[8,8,6,8] a[0,9,6,8] a[1,9,6,8] a[2,9,6,8]
 a[3,9,6,8] a[4,9,6,8] a[5,9,6,8] a[6,9,6,8] a[7,9,6,8] a[8,9,6,8]
 a[9,9,6,8] a[0,0,6,9] a[0,1,6,9] a[1,1,6,9] a[0,2,6,9] a[1,2,6,9]
 a[2,2,6,9] a[0,3,6,9] a[1,3,6,9] a[2,3,6,9] a[3,3,6,9] a[0,4,6,9]
 a[1,4,6,9] a[2,4,6,9] a[3,4,6,9] a[4,4,6,9] a[0,5,6,9] a[1,5,6,9]
 a[2,5,6,9] a[3,5,6,9] a[4,5,6,9] a[5,5,6,9] a[0,6,6,9] a[1,6,6,9]
 a[2,6,6,9] a[3,6,6,9] a[4,6,6,9] a[5,6,6,9] a[6,6,6,9] a[0,7,6,9]
 a[1,7,6,9] a[2,7,6,9] a[3,7,6,9] a[4,7,6,9] a[5,7,6,9] a[6,7,6,9]
 a[7,7,6,9] a[0,8,6,9] a[1,8,6,9] a[2,8,6,9] a[3,8,6,9] a[4,8,6,9]
 a[5,8,6,9] a[6,8,6,9] a[7,8,6,9] a[8,8,6,9] a[0,9,6,9] a[1,9,6,9]
 a[2,9,6,9] a[3,9,6,9] a[4,9,6,9] a[5,9,6,9] a[6,9,6,9] a[7,9,6,9]
 a[8,9,6,9] a[9,9,6,9] a[0,0,7,7] a[0,1,7,7] a[1,1,7,7] a[0,2,7,7]
 a[1,2,7,7] a[2,2,7,7] a[0,3,7,7] a[1,3,7,7] a[2,3,7,7] a[3,3,7,7]
 a[0,4,7,7] a[1,4,7,7] a[2,4,7,7] a[3,4,7,7] a[4,4,7,7] a[0,5,7,7]
 a[1,5,7,7] a[2,5,7,7] a[3,5,7,7] a[4,5,7,7] a[5,5,7,7] a[0,6,7,7]
 a[1,6,7,7] a[2,6,7,7] a[3,6,7,7] a[4,6,7,7] a[5,6,7,7] a[6,6,7,7]
 a[0,7,7,7] a[1,7,7,7] a[2,7,7,7] a[3,7,7,7] a[4,7,7,7] a[5,7,7,7]
 a[6,7,7,7] a[7,7,7,7] a[0,8,7,7] a[1,8,7,7] a[2,8,7,7] a[3,8,7,7]
 a[4,8,7,7] a[5,8,7,7] a[6,8,7,7] a[7,8,7,7] a[8,8,7,7] a[0,9,7,7]
 a[1,9,7,7] a[2,9,7,7] a[3,9,7,7] a[4,9,7,7] a[5,9,7,7] a[6,9,7,7]
 a[7,9,7,7] a[8,9,7,7] a[9,9,7,7] a[0,0,8,8] a[0,1,8,8] a[1,1,8,8]
 a[0,2,8,8] a[1,2,8,8] a[2,2,8,8] a[0,3,8,8] a[1,3,8,8] a[2,3,8,8]
 a[3,3,8,8] a[0,4,8,8] a[1,4,8,8] a[2,4,8,8] a[3,4,8,8] a[4,4,8,8]
 a[0,5,8,8] a[1,5,8,8] a[2,5,8,8] a[3,5,8,8] a[4,5,8,8] a[5,5,8,8]
 a[0,6,8,8] a[1,6,8,8] a[2,6,8,8] a[3,6,8,8] a[4,6,8,8] a[5,6,8,8]
 a[6,6,8,8] a[0,7,8,8] a[1,7,8,8] a[2,7,8,8] a[3,7,8,8] a[4,7,8,8]
 a[5,7,8,8] a[6,7,8,8] a[7,7,8,8] a[0,8,8,8] a[1,8,8,8] a[2,8,8,8]
 a[3,8,8,8] a[4,8,8,8] a[5,8,8,8] a[6,8,8,8] a[7,8,8,8] a[8,8,8,8]
 a[0,9,8,8] a[1,9,8,8] a[2,9,8,8] a[3,9,8,8] a[4,9,8,8] a[5,9,8,8]
 a[6,9,8,8] a[7,9,8,8] a[8,9,8,8] a[9,9,8,8] a[0,0,8,9] a[0,1,8,9]
 a[1,1,8,9] a[0,2,8,9] a[1,2,8,9] a[2,2,8,9] a[0,3,8,9] a[1,3,8,9]
 a[2,3,8,9] a[3,3,8,9] a[0,4,8,9] a[1,4,8,9] a[2,4,8,9] a[3,4,8,9]
 a[4,4,8,9] a[0,5,8,9] a[1,5,8,9] a[2,5,8,9] a[3,5,8,9] a[4,5,8,9]
 a[5,5,8,9] a[0,6,8,9] a[1,6,8,9] a[2,6,8,9] a[3,6,8,9] a[4,6,8,9]
 a[5,6,8,9] a[6,6,8,9] a[0,7,8,9] a[1,7,8,9] a[2,7,8,9] a[3,7,8,9]
 a[4,7,8,9] a[5,7,8,9] a[6,7,8,9] a[7,7,8,9] a[0,8,8,9] a[1,8,8,9]
 a[2,8,8,9] a[3,8,8,9] a[4,8,8,9] a[5,8,8,9] a[6,8,8,9] a[7,8,8,9]
 a[8,8,8,9] a[0,9,8,9] a[1,9,8,9] a[2,9,8,9] a[3,9,8,9] a[4,9,8,9]
 a[5,9,8,9] a[6,9,8,9] a[7,9,8,9] a[8,9,8,9] a[9,9,8,9] a[0,0,9,9]
 a[0,1,9,9] a[1,1,9,9] a[0,2,9,9] a[1,2,9,9] a[2,2,9,9] a[0,3,9,9]
 a[1,3,9,9] a[2,3,9,9] a[3,3,9,9] a[0,4,9,9] a[1,4,9,9] a[2,4,9,9]
 a[3,4,9,9] a[4,4,9,9] a[0,5,9,9] a[1,5,9,9] a[2,5,9,9] a[3,5,9,9]
 a[4,5,9,9] a[5,5,9,9] a[0,6,9,9] a[1,6,9,9] a[2,6,9,9] a[3,6,9,9]
 a[4,6,9,9] a[5,6,9,9] a[6,6,9,9] a[0,7,9,9] a[1,7,9,9] a[2,7,9,9]
 a[3,7,9,9] a[4,7,9,9] a[5,7,9,9] a[6,7,9,9] a[7,7,9,9] a[0,8,9,9]
 a[1,8,9,9] a[2,8,9,9] a[3,8,9,9] a[4,8,9,9] a[5,8,9,9] a[6,8,9,9]
 a[7,8,9,9] a[8,8,9,9] a[0,9,9,9] a[1,9,9,9] a[2,9,9,9] a[3,9,9,9]
 a[4,9,9,9] a[5,9,9,9] a[6,9,9,9] a[7,9,9,9] a[8,9,9,9] a[9,9,9,9] b[0,0,0]
 b[0,0,1] b[0,0,2] b[0,1,0] b[0,1,1] b[0,1,2] b[0,2,0] b[0,2,1] b[0,2,2]
 b[0,3,0] b[0,3,1] b[0,3,2] b[0,4,0] b[0,4,1] b[0,4,2] b[0,5,0] b[0,5,1]
 b[0,5,2] b[0,6,0] b[0,6,1] b[0,6,2] b[0,7,0] b[0,7,1] b[0,7,2] b[0,8,0]
 b[0,8,1] b[0,8,2] b[0,9,0] b[0,9,1] b[0,9,2] b[2,0,0] b[2,0,1] b[2,0,2]
 b[2,1,0] b[2,1,1] b[2,1,2] b[2,2,0] b[2,2,1] b[2,2,2] b[2,3,0] b[2,3,1]
 b[2,3,2] b[2,4,0] b[2,4,1] b[2,4,2] b[2,5,0] b[2,5,1] b[2,5,2] b[2,6,0]
 b[2,6,1] b[2,6,2] b[2,7,0] b[2,7,1] b[2,7,2] b[2,8,0] b[2,8,1] b[2,8,2]
 b[2,9,0] b[2,9,1] b[2,9,2] b[6,0,0] b[6,0,1] b[6,0,2] b[6,1,0] b[6,1,1]
 b[6,1,2] b[6,2,0] b[6,2,1] b[6,2,2] b[6,3,0] b[6,3,1] b[6,3,2] b[6,4,0]
 b[6,4,1] b[6,4,2] b[6,5,0] b[6,5,1] b[6,5,2] b[6,6,0] b[6,6,1] b[6,6,2]
 b[6,7,0] b[6,7,1] b[6,7,2] b[6,8,0] b[6,8,1] b[6,8,2] b[6,9,0] b[6,9,1]
 b[6,9,2] b[8,0,0] b[8,0,1] b[8,0,2] b[8,1,0] b[8,1,1] b[8,1,2] b[8,2,0]
 b[8,2,1] b[8,2,2] b[8,3,0] b[8,3,1] b[8,3,2] b[8,4,0] b[8,4,1] b[8,4,2]
 b[8,5,0] b[8,5,1] b[8,5,2] b[8,6,0] b[8,6,1] b[8,6,2] b[8,7,0] b[8,7,1]
 b[8,7,2] b[8,8,0] b[8,8,1] b[8,8,2] b[8,9,0] b[8,9,1] b[8,9,2] b[9,0,0]
 b[9,0,1] b[9,0,2] b[9,1,0] b[9,1,1] b[9,1,2] b[9,2,0] b[9,2,1] b[9,2,2]
 b[9,3,0] b[9,3,1] b[9,3,2] b[9,4,0] b[9,4,1] b[9,4,2] b[9,5,0] b[9,5,1]
 b[9,5,2] b[9,6,0] b[9,6,1] b[9,6,2] b[9,7,0] b[9,7,1] b[9,7,2] b[9,8,0]
 b[9,8,1] b[9,8,2] b[9,9,0] b[9,9,1] b[9,9,2] b[1,0,0] b[1,0,1] b[1,0,2]
 b[1,1,0] b[1,1,1] b[1,1,2] b[1,2,0] b[1,2,1] b[1,2,2] b[1,3,0] b[1,3,1]
 b[1,3,2] b[1,4,0] b[1,4,1] b[1,4,2] b[1,5,0] b[1,5,1] b[1,5,2] b[1,6,0]
 b[1,6,1] b[1,6,2] b[1,7,0] b[1,7,1] b[1,7,2] b[1,8,0] b[1,8,1] b[1,8,2]
 b[1,9,0] b[1,9,1] b[1,9,2] b[3,0,0] b[3,0,1] b[3,0,2] b[3,1,0] b[3,1,1]
 b[3,1,2] b[3,2,0] b[3,2,1] b[3,2,2] b[3,3,0] b[3,3,1] b[3,3,2] b[3,4,0]
 b[3,4,1] b[3,4,2] b[3,5,0] b[3,5,1] b[3,5,2] b[3,6,0] b[3,6,1] b[3,6,2]
 b[3,7,0] b[3,7,1] b[3,7,2] b[3,8,0] b[3,8,1] b[3,8,2] b[3,9,0] b[3,9,1]
 b[3,9,2] b[7,0,0] b[7,0,1] b[7,0,2] b[7,1,0] b[7,1,1] b[7,1,2] b[7,2,0]
 b[7,2,1] b[7,2,2] b[7,3,0] b[7,3,1] b[7,3,2] b[7,4,0] b[7,4,1] b[7,4,2]
 b[7,5,0] b[7,5,1] b[7,5,2] b[7,6,0] b[7,6,1] b[7,6,2] b[7,7,0] b[7,7,1]
 b[7,7,2] b[7,8,0] b[7,8,1] b[7,8,2] b[7,9,0] b[7,9,1] b[7,9,2] b[4,0,0]
 b[4,0,1] b[4,0,2] b[4,1,0] b[4,1,1] b[4,1,2] b[4,2,0] b[4,2,1] b[4,2,2]
 b[4,3,0] b[4,3,1] b[4,3,2] b[4,4,0] b[4,4,1] b[4,4,2] b[4,5,0] b[4,5,1]
 b[4,5,2] b[4,6,0] b[4,6,1] b[4,6,2] b[4,7,0] b[4,7,1] b[4,7,2] b[4,8,0]
 b[4,8,1] b[4,8,2] b[4,9,0] b[4,9,1] b[4,9,2] b[5,0,0] b[5,0,1] b[5,0,2]
 b[5,1,0] b[5,1,1] b[5,1,2] b[5,2,0] b[5,2,1] b[5,2,2] b[5,3,0] b[5,3,1]
 b[5,3,2] b[5,4,0] b[5,4,1] b[5,4,2] b[5,5,0] b[5,5,1] b[5,5,2] b[5,6,0]
 b[5,6,1] b[5,6,2] b[5,7,0] b[5,7,1] b[5,7,2] b[5,8,0] b[5,8,1] b[5,8,2]
 b[5,9,0] b[5,9,1] b[5,9,2] x[0,0] x[0,1] x[0,2] x[1,0] x[1,1] x[1,2]
 x[2,0] x[2,1] x[2,2] x[3,0] x[3,1] x[3,2] x[4,0] x[4,1] x[4,2] x[5,0]
 x[5,1] x[5,2] x[6,0] x[6,1] x[6,2] x[7,0] x[7,1] x[7,2] x[8,0] x[8,1]
 x[8,2] x[9,0] x[9,1] x[9,2]
Generals
 Q[0] Q[1] Q[2] Q[3] Q[4] Q[5] Q[6] Q[7] Q[8] Q[9] z[0,0] z[0,1] z[0,2]
 z[1,0] z[1,1] z[1,2] z[2,0] z[2,1] z[2,2] z[3,0] z[3,1] z[3,2] z[4,0]
 z[4,1] z[4,2] z[5,0] z[5,1] z[5,2] z[6,0] z[6,1] z[6,2] z[7,0] z[7,1]
 z[7,2] z[8,0] z[8,1] z[8,2] z[9,0] z[9,1] z[9,2] tp[0] tp[1] tp[2] tp[3]
 tp[4] tp[5] tp[6] tp[7] tp[8] tp[9] tn[0] tn[1] tn[2] tn[3] tn[4] tn[5]
 tn[6] tn[7] tn[8] tn[9]
End
