import numpy as np

PlateLong = 5
PlatePositionStart = 5
AngularFrecuncy = 1
Epsilone_0 = 1
Mu_0 = 1
EpsilonePlate = 1
MuPlate = 1
HiPlate = 1
CountOfBasisElements = 1
KLong = 1

# скорость света
def GetPlatePosition():
    positionStart = PlatePositionStart
    return {'startPosition': positionStart, 'endPosition': positionStart + GetPlateLong()}


def GetPlateLong():
    plateLong = PlateLong
    return plateLong

def GetC():
    return np.sqrt(GetMu0() * GetEpsilon0())

def GetEpsilon0():
    epsilone_0 = Epsilone_0
    return epsilone_0

def GetMu0():
    mu_0 = Mu_0
    return mu_0

def GetEpsilone(x_position):
    plateEpsilone = EpsilonePlate

    if (GetPlatePosition().get('startPosition') <= x_position
            and x_position <= GetPlatePosition().get('endPosition')):
        return plateEpsilone
    else:
        return 1 #возможно стоит заменить на ноль, но это зависит от того как используем в формулах (может вылезти ошибк аделения на ноль)

def GetMu(x_position):
    plateMu = MuPlate

    if (GetPlatePosition().get('startPosition') <= x_position
            and x_position <= GetPlatePosition().get('endPosition')):
        return plateMu
    else:
        return 1

def GetRefractiveIndex(x_position):

    return np.sqrt(GetMu(x_position) * GetEpsilone(x_position))


def GetHi():
    plateHi = HiPlate
    def HiFunction(x_position):
        if ( x_position >= GetPlatePosition().get('startPosition')
                and x_position <= GetPlatePosition().get('endPosition') ):
            return plateHi
        else:
            return 0
    return HiFunction

def GetAngelFrecuncy():
    w = AngularFrecuncy
    return w

def GetCountBasisElements():
    N = CountOfBasisElements
    return N

def GetKLong():
    k = KLong
    return k
