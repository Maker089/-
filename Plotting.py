import matplotlib.pyplot as plt
import constant as const
import TeoreticalFieldSolve as tfs
import NewFieldByTrancferMatrix as tm

#Set Plate Parametrs
const.PlateLong = 5
const.PlatePositionStart = 5


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
SideX = 15
SideToPlate = const.PlatePositionStart - StartX
SideInPLate = const.PlateLong
SideOutPlate = SideX - (const.PlateLong + const.PlatePositionStart)

N = 100
# ResultDirectionTeorFieldReal = {"Ex": [], "Ey": [], "Hz": []}
# ResultDirectionTeorFieldImag = {"Ex": [], "Ey": [], "Hz": []}
ResultDictionaryTeor = {"Ex": [], "Ey": [], "Hz": []}

# ResultDictionaryReal = {"TE_H": [], "TE_Ex": [], "TE_Ey": [], "TM_E": [], "TM_Hx": [], "TM_Hy": []}
# ResultDictionaryImag = {"TE_H": [], "TE_Ex": [], "TE_Ey": [], "TM_E": [], "TM_Hx": [], "TM_Hy": []}
ResultDictionary = {"TE_H": [], "TE_Ex": [], "TE_Ey": [], "TM_E": [], "TM_Hx": [], "TM_Hy": []}
XPositionList = []
XPositionList.append(0)

#Set Start field components for not teor
ResultDictionary["TE_H"].append(tm.GetInTE_H(StartX, False))
ResultDictionary["TE_Ex"].append(tm.GetInTE_Ex(StartX))
ResultDictionary["TE_Ey"].append(tm.GetInTE_Ey(StartX))

#Set Start field components for teor
ResultDictionaryTeor["Hz"].append(tm.GetInTE_H(StartX, False))
ResultDictionaryTeor["Ex"].append(tm.GetInTE_Ex(StartX))
ResultDictionaryTeor["Ey"].append(tm.GetInTE_Ey(StartX))



#Solve fields to plate
for count in range(1, N):
    currentXPosition = (SideToPlate) / N * count + StartX
    XPositionList.append(currentXPosition)
    currentTransferLong = SideToPlate / N * count

    #NotTeorComponents
    ResultDictionary["TE_H"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong, ResultDictionary["TE_H"][0],
                                                              ResultDictionary["TE_Ey"][0])["H"])
    ResultDictionary["TE_Ex"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong, ResultDictionary["TE_H"][0],
                                                               ResultDictionary["TE_Ey"][0])["Ex"])
    ResultDictionary["TE_Ey"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong, ResultDictionary["TE_H"][0],
                                                               ResultDictionary["TE_Ey"][0])["Ey"])

    #TeorComponents
    ResultDictionaryTeor["Hz"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong, ResultDictionary["TE_H"][0],
                                                              ResultDictionary["TE_Ey"][0])["H"])
    ResultDictionaryTeor["Ex"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong, ResultDictionary["TE_H"][0],
                                                               ResultDictionary["TE_Ey"][0])["Ex"])
    ResultDictionaryTeor["Ey"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong, ResultDictionary["TE_H"][0],
                                                               ResultDictionary["TE_Ey"][0])["Ey"])



#Solve fields in plate
for count in range(0, N):
    currentXPosition1 = (SideInPLate) / N * count + StartX + SideToPlate
    currentTransferLong1 = SideInPLate / N * count
    XPositionList.append(currentXPosition1)

    # NotTeorComponents
    ResultDictionary["TE_H"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong1, ResultDictionary["TE_H"][0],
                                                              ResultDictionary["TE_Ey"][0])["H"])
    ResultDictionary["TE_Ex"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong1, ResultDictionary["TE_H"][0],
                                                               ResultDictionary["TE_Ey"][0])["Ex"])
    ResultDictionary["TE_Ey"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong1, ResultDictionary["TE_H"][0],
                                                               ResultDictionary["TE_Ey"][0])["Ey"])
    # TeorComponents
    ResultDictionaryTeor["Hz"].append(tfs.TeorFieldH(currentXPosition1))
    ResultDictionaryTeor["Ex"].append(0)
    ResultDictionaryTeor["Ey"].append(tfs.TeorFieldE(currentXPosition1))


#Parametrs field for transfers matrix
InOutTE_H = ResultDictionary["TE_H"][-1]
InOutTE_Ex = ResultDictionary["TE_Ex"][-1]
InOutTE_Ey = ResultDictionary["TE_Ey"][-1]

TeorInOut_Hz = ResultDictionaryTeor["Hz"][-1]
TeorInOut_Ex = ResultDictionaryTeor["Ex"][-1]
TeorInOut_Ey = ResultDictionaryTeor["Ey"][-1]


#Solve fields out plate
for count in range(0, N):
    currentXPosition2 = (SideOutPlate) / N * count + StartX + SideToPlate + SideInPLate
    currentTransferLong2 = SideOutPlate / N * count
    XPositionList.append(currentXPosition2)

    # NotTeorComponents
    ResultDictionary["TE_H"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong2, InOutTE_H, InOutTE_Ey)["H"])
    ResultDictionary["TE_Ex"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong2, InOutTE_H, InOutTE_Ey)["Ex"])
    ResultDictionary["TE_Ey"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong2, InOutTE_H, InOutTE_Ey)["Ey"])

    # TeorComponents
    ResultDictionaryTeor["Hz"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong2, TeorInOut_Hz, TeorInOut_Ey)["H"])
    ResultDictionaryTeor["Ex"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong2, TeorInOut_Hz, TeorInOut_Ey)["Ex"])
    ResultDictionaryTeor["Ey"].append(tm.GetNewFieldComponentTE(StartX, currentTransferLong2, TeorInOut_Hz, TeorInOut_Ey)["Ey"])


RealComponents = []
ImagComponents = []
TeorRealComponents = []
TeorImagComponents = []
for count in range(0, 3 * N):
    RealComponents.append(ResultDictionary["TE_Ey"][count].real)
    ImagComponents.append(ResultDictionary["TE_Ey"][count].imag)

    TeorRealComponents.append(ResultDictionaryTeor["Ey"][count].real)
    TeorImagComponents.append(ResultDictionaryTeor["Ey"][count].imag)


plt.plot(XPositionList, RealComponents, "r", label = "Re")
plt.plot(XPositionList, ImagComponents, "b", label = "Im")

plt.plot(XPositionList, TeorRealComponents, "y-", label = "Teor Re")
plt.plot(XPositionList, TeorImagComponents, "g-", label = "Teor Im")

plt.plot([const.GetPlatePosition()['startPosition'], const.GetPlatePosition()['startPosition']], [-2,2], 'black')
plt.plot([const.GetPlatePosition()['startPosition'] + const.GetPlateLong(), const.GetPlatePosition()['startPosition'] + const.GetPlateLong()], [-2,2], 'black')

plt.title("Ey")
plt.legend(loc = "upper left")
plt.show()