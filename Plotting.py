import matplotlib
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import constant as const
import TeoreticalFieldSolve as tfs
import NewFieldByTrancferMatrix as tm

#Ser Plate Parametrs
const.PlateLong = 1
const.PlatePositionStart = 20


#Set constanse value
const.AngularFrecuncy = 0.556
const.CounOfBasisElements = 20
const.KLong = 1
const.HiPlate = 1
const.MuPlate = 1
const.Mu_0 = 1
const.EpsilonePlate = 1
const.Epsilone_0 = 1


#Plotting Grafics
StartX = 1
SideX = 10
N = 1000
ResultDirectionTeorFieldReal = {"Ex": [], "Ey": [], "Hz": []}
ResultDirectionTeorFieldImag = {"Ex": [], "Ey": [], "Hz": []}

ResultDictionaryReal = {"TE_H": [], "TE_Ex": [], "TE_Ey": [], "TM_E": [], "TM_Hx": [], "TM_Hy": []}
ResultDictionaryImag = {"TE_H": [], "TE_Ex": [], "TE_Ey": [], "TM_E": [], "TM_Hx": [], "TM_Hy": []}
XPositionList = []

for count in range(0, 1000):
    currentTransferLong = (SideX - StartX) / N * count + StartX
    XPositionList.append(currentTransferLong)

    #RealComponents
    ResultDictionaryReal["TE_H"].append(GetNewFieldComponentTE(StartX, currentTransferLong)["H"].real)
    ResultDictionaryReal["TE_Ex"].append(GetNewFieldComponentTE(StartX, currentTransferLong)["Ex"].real)
    ResultDictionaryReal["TE_Ey"].append(GetNewFieldComponentTE(StartX, currentTransferLong)["Ey"].real)

    # ResultDictionaryReal["TM_E"].append(GetNewFieldComponentTM(StartX, currentTransferLong)["E"].real)
    # ResultDictionaryReal["TM_Hx"].append(GetNewFieldComponentTM(StartX, currentTransferLong)["Hx"].real)
    # ResultDictionaryReal["TM_Hy"].append(GetNewFieldComponentTM(StartX, currentTransferLong)["Hy"].real)

    #Imaginary Components
    ResultDictionaryImag["TE_H"].append(GetNewFieldComponentTE(StartX, currentTransferLong)["H"].imag)
    ResultDictionaryImag["TE_Ex"].append(GetNewFieldComponentTE(StartX, currentTransferLong)["Ex"].imag)
    ResultDictionaryImag["TE_Ey"].append(GetNewFieldComponentTE(StartX, currentTransferLong)["Ey"].imag)

    # ResultDictionaryImag["TM_E"].append(GetNewFieldComponentTM(StartX, currentTransferLong)["E"].imag)
    # ResultDictionaryImag["TM_Hx"].append(GetNewFieldComponentTM(StartX, currentTransferLong)["Hx"].imag)
    # ResultDictionaryImag["TM_Hy"].append(GetNewFieldComponentTM(StartX, currentTransferLong)["Hy"].imag)


plt.plot(XPositionList, ResultDictionaryReal["TM_Hx"], "r", label = "Re")
plt.plot(XPositionList, ResultDictionaryImag["TM_Hx"], "b", label = "Im")
plt.title("TM_Hx")
plt.legend(loc = "upper left")
plt.show()