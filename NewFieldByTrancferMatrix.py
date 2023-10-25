import numpy as np
import constant as const


l = const.PlateLong

epsl0 = const.GetEpsilon0()
mu0 = const.GetMu0()

i = complex(0, 1)
k = np.transpose(np.array([const.GetKLong(), 0]))
kx = k[0]
ky = k[1]

w = const.GetAngelFrecuncy()

#TE Polarisation
#magnetic field input
def GetInTE_H(x_position, isForDiffer):
    A_plus = 1
    A_minus = 1
    if (isForDiffer):
        return i * kx * A_plus * np.exp(i * kx * x_position) - i * kx * A_minus * np.exp(-i * kx * x_position)
    else:
        return A_plus * np.exp(i * kx * x_position) + A_minus * np.exp(-i * kx * x_position)

#Electrone field input , x_pos == x_position
def GetInTE_Ex(x_pos):
    return -kx / (w * epsl0 * const.GetEpsilone(x_pos)) * GetInTE_H(x_pos, False)

def GetInTE_Ey(x_pos):
    return - i / (w * epsl0 * const.GetEpsilone(x_pos)) * GetInTE_H(x_pos, True)


#получаем значение полей (Ey and H = Hz) при данной позиции x и при данной длинне l (на которую хотим посмотреть). На вход идут напряженности электрического и магнитного полей в тьочке x, а на выходе получаем значение напряженностей в точке x+l. Важно! - преобразоание в однородном слое, где epsilone == const
def GetNewFieldComponentTE(x_pos_start, transferLong, fieldH, fieldEy):

    transferMatrix = np.array([[np.cos(kx * transferLong), i * w * epsl0 * const.GetEpsilone(x_pos_start) * (1/kx) * np.sin(kx * transferLong)],
                           [i * kx * np.sin(kx * transferLong) / (w * epsl0 * const.GetEpsilone(x_pos_start)), np.cos(kx * transferLong)]])

    startFieldVector = np.transpose(np.array([fieldH, fieldEy]))

    NewEx = -kx/(w * epsl0 * const.GetEpsilone(x_pos_start)) * np.dot(transferMatrix, startFieldVector)[0]

    return {"H": np.dot(transferMatrix, startFieldVector)[0], "Ex": NewEx, "Ey": np.dot(transferMatrix, startFieldVector)[1]}




#TM Polarisation
# Electric field E == Ez input
def GetInTM_E(x_position, isForDiffer):
    A_plus = 1
    A_minus = 1
    if (isForDiffer):
        return i * kx * A_plus * np.exp(i * kx * x_position) - i * kx * A_minus * np.exp(-i * kx * x_position)
    else:
        return A_plus * np.exp(i * kx * x_position) + A_minus * np.exp(-i * kx * x_position)


# Magnetic field input
def GetInTM_Hx(x_pos):
    return ky / (w * mu0 * const.GetMu(x_pos)) * GetInTM_E(x_pos, False)


def GetInTM_Hy(x_pos):
    return i / (w * mu0 * const.GetMu(x_pos)) * GetInTE_H(x_pos, True)


# получаем значение полей при данной позиции x и при данной длинне l (на которую хотим посмотреть). На вход идут напряженности электрического и магнитного полей в тьочке x, а на выходе получаем значение напряженностей в точке x+l. Важно! - преобразоание в однородном слое, где epsilone == const
def GetNewFieldComponentTM(x_pos_start, transferLong):
    transferMatrix = np.array(
        [[np.cos(kx * transferLong), i * w * mu0 * const.GetMu(x_pos_start) * (1/kx) * np.sin(kx * transferLong)],
         [i * kx * np.sin(kx * transferLong) / (w * mu0 * const.GetMu(x_pos_start)), np.cos(kx * transferLong)]]
    )

    startFieldVector = np.transpose(
        np.array([GetInTM_E(x_pos_start, False), GetInTM_Hy(x_pos_start)]))

    NewHx = ky / (w * mu0 * const.GetMu(x_pos_start)) * np.dot(transferMatrix, startFieldVector)[0]

    return {"E": np.dot(transferMatrix, startFieldVector)[0], "Hx": NewHx,
            "Hy": np.dot(transferMatrix, startFieldVector)[1]}


