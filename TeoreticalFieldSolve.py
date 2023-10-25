import numpy as np
import constant as const

#Count of Basis elements
N = const.CountOfBasisElements

L = const.GetPlateLong()
c = const.GetC()
hi = const.GetHi()
pi = np.pi
i = complex(0, 1)
w = const.GetAngelFrecuncy()
k = const.KLong
mu0 = const.GetMu0()
XMid = const.PlatePositionStart + const.PlateLong/2

def E0(xPosition):
    return np.exp(i * (w/c) * xPosition)

#delta function
def DF(m, n):
    if (m == n):
        return 1
    else:
        return 0

def GetEvenPolinom(polinomNumber):
    def PolinomFunction(xPosition):
        if (polinomNumber == 0):
            return np.sqrt(c * c / (L * hi(xPosition)))
        else:
            return np.sqrt(2 * c * c / (L * hi(xPosition))) * np.cos(polinomNumber * pi * xPosition / L)

    return PolinomFunction

def GetOddPolinom(polinomNumber):
    def PolinomFunction(x_position):
        return np.sqrt(2 * c * c / (L * hi(x_position))) * np.sin(polinomNumber * pi * x_position / L)

    return PolinomFunction

def GetAngelFrecuncy(frecuncyNumber):
    def FrecuncyFunction(xPosition):
        return pi * c * 2 * frecuncyNumber / (L * np.sqrt(hi(xPosition)))

    return FrecuncyFunction


# m - row in matrix, n - colum in matrix
def GetGreenMatrixComponent(m, n):
    wn = GetAngelFrecuncy(n)(XMid)
    if (m == n == 0):
        return -(i/L) * (c/w) * (1 - np.exp(i * (w / c) * L) + i * (w / c) * L)

    if (m % 2 == 0) and (n % 2 == 0):
        return i * (w ** 3) / (L * c * (wn ** 2 - w ** 2)) \
               * ((i * (w/c) * (L/2) * DF(m, n))/((w/c)**2 - (m*np.pi/L)**2)
                  - ((-1)**((n+m)/2)) * ((w/c)**2) * (np.exp(i*w*L/c) - 1)
                  * (1/((w/c)**2 - (m*np.pi/L)**2)) * (1/((w/c)**2 - (n*np.pi/L)**2)))

    if (m % 2 != 0) and (n % 2 != 0):
        return i * (w ** 3) / (L * c * (wn ** 2 - w ** 2)) \
               * ((i * (w / c) * L * DF(m, n)) / ((w / c) ** 2 - (m * np.pi / L) ** 2)
                  - ((-1) ** ((n + m) / 2)) * 2 * ((w / c) ** 2) * (np.exp(i * w * L / c) + 1)
                  * (1 / ((w / c) ** 2 - (m * np.pi / L) ** 2)) * (1 / ((w / c) ** 2 - (n * np.pi / L) ** 2)))

    if ((m % 2 != 0) & (n % 2 == 0)) or ((n % 2 != 0) & (m % 2 == 0)):
        return 0

# Fi == f_i component, where \beta_i = f_i + \sum_j G_{ij} for E_0 == \exp{i w/c * x}
def Fi(n):
    a = 1
    if n == 0: a = 1 / np.sqrt(2)

    if n % 2 == 0:
        return a * ((-1) ** (n/2)) * (w/c) * np.sqrt((2 * hi(XMid)) / (L * c**2)) * (2*np.sin((w/c) * (L/2))) / ((w / c) ** 2 - (n * np.pi / L) ** 2)
    else:
        return ((-1) ** ((n+1)/2)) * (w/c) * np.sqrt((2 * hi(XMid)) / (L * c**2)) * (2*i*np.cos((w/c) * (L/2))) / ((w / c) ** 2 - (n * np.pi / L) ** 2)

# J_m(x) where E(x) = E_0(x) + w^2/c^2 * \sum_j w^2/(w_j^2 - w^2) * b_j * J_j
def Jm(xPosition,m):
    a = 1
    if m == 0: a = 1 / np.sqrt(2)

    if m % 2 == 0:
        return np.exp(-i*w*xPosition/c) * (i*np.sqrt(2*c**2 / (L*hi(xPosition))) * ((-1)**(m/2)) * (np.sin(w*L/c/2) / ((w / c) ** 2 - (m * np.pi / L) ** 2)))
    else:
        return np.exp(-i*w*xPosition/c) * (np.sqrt(2*c**2 / (L*hi(xPosition))) * ((-1)**((m-1)/2)) * (np.cos(w*L/c/2) / ((w / c) ** 2 - (m * np.pi / L) ** 2)))

# Bi == \beta_i = f_i + \sum_j G_{ij} for E_0 == \exp{i w/c * x}
def Bi(Number):
    sumOfGElements = 0
    for count in range(0, N):
        sumOfGElements += GetGreenMatrixComponent(Number, count)

    return Fi(Number) + sumOfGElements

# Elector Field E = E_0(x) + w^2/c^2 * \sum_j w^2/(w_j^2 - w^2) * b_j * J_j
def TeorFieldE(xPosition):

    sumOfComponents = 0
    for count in range(0, N):
        sumOfComponents += (w**2) / ((GetAngelFrecuncy(count)(xPosition))**2 - w**2) * Bi(count) * Jm(xPosition, count)

    return E0(xPosition) + sumOfComponents

def TeorFieldH(xPosition):
    return (k / (mu0 * w)) * TeorFieldE(xPosition)
