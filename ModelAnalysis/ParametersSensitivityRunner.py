# writing functions that will :
# - read a file 
# - extract the following information :
#   • cell mutation state (results.vizcelltypes) -> number of each type of cells, number of branches 
#   • center of cells ie nodes coordinates (results.viznodes) -> closest and furthest tip cell coordinates from the first fixed cell -> area and volume of the blood vessel network computed (area and volume of a half sphere calculated from the furthest tip cell coordinates which represent the radius)
#   • consecutive branches (results.branchesnumber) -> maximum and minimum of consecutive branches -> comparison branches obtained and branches expected 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib import cm
from matplotlib.ticker import LinearLocator 
from matplotlib.colors import LinearSegmentedColormap

from scipy.spatial import ConvexHull

# function reading the file 'results.vizcelltypes' and returning a list of 0, 1, 2 corresponding to the mutation states 
def MutationStates(file_cellmutation):
    # initialisation of the list
    list_cellmutation = []

    # we open the file
    f = open(file_cellmutation, 'r')

    for line in f:
        pass
    last_line = line

    list_cellmutation = [float(x) for x in last_line.split()[1:]]

    # we close the file 
    f.close()

    return list_cellmutation

# function reading the file 'results.viznodes' and returning an array containing the coordinates of all the nodes
def NodesCoordinates(file_nodescoordinates):
    # initialisation of the list
    list_nodescoordinates = []

    # we open the file
    f = open(file_nodescoordinates, 'r')

    for line in f:
        pass
    last_line = line

    list_nodescoordinates = [float(x) for x in last_line.split()[1:]]

    NumberNodes = int((len(list_nodescoordinates))/3)
    array_nodescoordinates = np.reshape(list_nodescoordinates, (NumberNodes,3))

    # we close the file 
    f.close()

    return array_nodescoordinates

# function reading the file 'results.branchesnumber' and returning a list of the consecutive branches of the cell 
def ConsecutiveBranches(file_consecutivebranches):
    # initialisation of the list
    list_consecutivebranches = []

    # we open the file
    f = open(file_consecutivebranches, 'r')

    for line in f:
        pass
    last_line = line

    list_consecutivebranches = [float(x) for x in last_line.split()[1:]]

    # we close the file 
    f.close()

    return list_consecutivebranches

def BranchesNumber(file_branchesnumber):
    # initialisation of the list
    list_branchesnumber = []

    # we open the file
    f = open(file_branchesnumber, 'r')

    for line in f:
        pass
    last_line = line

    list_branchesnumber = [float(x) for x in last_line.split()[1:]]

    # we close the file 
    f.close()

    return list_branchesnumber

def TimeFirstBranch(file_branchesnumber):
    # initialisation 
    time_firstbranch = 0
    time_increment = 0

    # we open the file
    f = open(file_branchesnumber, 'r')

    for line in f:
        list_branchesnumber = [float(x) for x in line.split()[1:]]
        if(1 in list_branchesnumber):
            time_firstbranch = time_increment 
            break
        time_increment += 1

    # we close the file 
    f.close()

    return time_firstbranch

def LongestPath(file_tortuosity):
    # initialisation of the values 
    arc = 0
    length = 0

    # we open the file
    f = open(file_tortuosity, 'r')

    for line in f:
        pass
    last_line = line

    list_longestpath = [float(x) for x in last_line.split()[1:]]

    # we close the file 
    f.close()

    list_arcs = []
    list_tipcellscoordinates = []

    # we collect the length of the arcs and the coordinates of the tip cells corresponding
    for k in range(len(list_longestpath)):
        if (k % 4 == 0):
            list_arcs.append(list_longestpath[k])
        else :
            list_tipcellscoordinates.append(list_longestpath[k])

    NbTipCells = int(len(list_tipcellscoordinates)/3)

    # we re arrange the tip cells coordinates list into an array 
    array_tipcellscoordinates = np.reshape(list_tipcellscoordinates, (NbTipCells,3))

    # we look for the longest arc 
    arc = list_arcs[0]
    LongestArcNb = 0
    for k in range(len(list_arcs)):
        if(list_arcs[k] > arc):
            arc = list_arcs[k]
            LongestArcNb = k

    # we collect the coordinates of the tip cell associated to the longest arc 
    TipCellLength_x = array_tipcellscoordinates[LongestArcNb][0]
    TipCellLength_y = array_tipcellscoordinates[LongestArcNb][1]
    TipCellLength_z = array_tipcellscoordinates[LongestArcNb][2]

    # we calculate its norm 
    length = np.sqrt(TipCellLength_x**2 + TipCellLength_y**2 + TipCellLength_z**2)

    return arc, length 


# a function to obtain the total number of cells at the end of the simulation 
def TotalNumberCells(file_cellmutation):
    list_cellmutation = MutationStates(file_cellmutation)
    return len(list_cellmutation)

# a function to obtain the total number of tip cells at the end of the simulation 
def TotalNumberTipCells(file_cellmutation):
    NumberTipCells = 0

    list_cellmutation = MutationStates(file_cellmutation)

    for elem in list_cellmutation:
        if(int(elem) == 0):
            NumberTipCells += 1

    return NumberTipCells

# a function to obtain the total number of stalk cells at the end of the simulation (including branching points)
def TotalNumberStalkCells(file_cellmutation):
    NumberStalkCells = 0

    list_cellmutation = MutationStates(file_cellmutation)

    for elem in list_cellmutation:
        if(elem != 0):
            NumberStalkCells += 1

    return NumberStalkCells

# a function to obtain the total number of branches at the end of the simulation 
def TotalNumberBranches(file_cellmutation):
    NumberBranches = 0

    list_cellmutation = MutationStates(file_cellmutation)

    for elem in list_cellmutation:
        if(elem == 2):
            NumberBranches += 1
        elif(elem == 0):
            NumberBranches += 1

    return NumberBranches

# a function to obtain the minimum consecutive branches 
def MinimumConsecutiveBranches(list_consecutivebranches):
    NumberMinimumConsecutiveBranches = 1.0

    for elem in list_consecutivebranches:
        if(elem < NumberMinimumConsecutiveBranches):
            NumberMinimumConsecutiveBranches = elem

    return NumberMinimumConsecutiveBranches 

# a function to obtain the maximum consecutive branches 
def MaximumConsecutiveBranches(list_consecutivebranches):
    NumberMaximumConsecutiveBranches = 1.0

    for elem in list_consecutivebranches:
        if(elem > NumberMaximumConsecutiveBranches):
            NumberMaximumConsecutiveBranches = elem

    return NumberMaximumConsecutiveBranches 

def BirthTimeCells(file_birthtimecells):
    # initialisation of the list
    list_birthtimecells = []

    # we open the file
    f = open(file_birthtimecells, 'r')

    for line in f:
        pass
    last_line = line

    list_birthfile_birthtimecells = [float(x) for x in last_line.split()[1:]]

    # we close the file 
    f.close()

    return list_birthfile_birthtimecells

def ComputationalConsecutiveBranches(Tcycle):
    nb_test = 1000
    t = 0
    TNcycle = int(np.floor(time/Tcycle))

    AverageNumberExpectedConsecutiveBranches = 0

    # while t < nb_test :
    #     NumberExpectedConsecutiveBranches = 0
    #     list_cellcycle = [0]
    #     list_birthtimecells = [0]
    #     for k in range(time):
    #         new_list_cellcycle = list_cellcycle
    #         for j in range(len(list_cellcycle)):
    #             proba = np.random.random()
    #             if(proba <= Psprout and list_cellcycle[j] >= 13):
    #                 new_list_cellcycle[j] = 0
    #                 list_birthtimecells.append(k)
    #                 new_list_cellcycle.append(0)

    #                 NumberExpectedConsecutiveBranches += 2
    #             elif(proba > Psprout and list_cellcycle[j] >= 13):
    #                 new_list_cellcycle[j] = 0
    #                 list_birthtimecells[j] = k
    #             else:
    #                 new_list_cellcycle[j] += 1
    #         list_cellcycle = new_list_cellcycle
    #     t += 1
    #     AverageNumberExpectedConsecutiveBranches = AverageNumberExpectedConsecutiveBranches + NumberExpectedConsecutiveBranches/1000

    while t < nb_test :
        # initialisation : at the beginning of the simulation : only 1 tip cell and hence only 1 branch 
        NbTipCells = 1 
        NbCells = 1
        NumberExpectedConsecutiveBranches = 1 

        # we reset the random seed for each simulation : 
        #np.random.seed()

        for j in range(TNcycle):
            NbCells = NbCells + NbTipCells

            for i in range(NbTipCells):
                proba = np.random.random()
                if(proba <= Psprout):
                    NbTipCells = NbTipCells + 1
                    NumberExpectedConsecutiveBranches = NumberExpectedConsecutiveBranches + 2
        
        t += 1
        AverageNumberExpectedConsecutiveBranches = AverageNumberExpectedConsecutiveBranches + NumberExpectedConsecutiveBranches/nb_test
        NumberExpectedConsecutiveBranches = 0 # re-initialisation for the next simulation

    return AverageNumberExpectedConsecutiveBranches


# a function to compare the maximum consecutive branches obtained and the one expected 
def ComparisonConsecutiveBranches(file_cellmutation, Psprout, Tcycle):
    # branches obtained 
    NumberMaximumConsecutiveBranches = TotalNumberBranches(file_cellmutation)

    # expected branches code 
    NumberExpectedConsecutiveBranches = 1

    # formula of the expectation 
    NumberExpectedConsecutiveBranches = ComputationalConsecutiveBranches(Tcycle)

    NumberFormulaExpectedConsecutiveBranches = 2*(1+Psprout)**(np.floor(time/Tcycle)) - 1
    
    return NumberMaximumConsecutiveBranches, NumberExpectedConsecutiveBranches, NumberFormulaExpectedConsecutiveBranches

# function returning the closest tip cell coordinates from the first fixed cell
def ClosestTipCell(file_nodescoordinates, file_cellmutation):
    list_cellmutation = MutationStates(file_cellmutation)
    array_nodescoordinates = NodesCoordinates(file_nodescoordinates)

    CoordinatesClosestTipCell = np.zeros((3,1))
    CoordinatesClosestTipCell[0] = array_nodescoordinates[1][0]
    CoordinatesClosestTipCell[1] = array_nodescoordinates[1][1]
    CoordinatesClosestTipCell[2] = array_nodescoordinates[1][2]
    NormClosestTipCell = np.sqrt(CoordinatesClosestTipCell[0]**2 + CoordinatesClosestTipCell[1]**2+ CoordinatesClosestTipCell[2]**2)

    NumberNodes = array_nodescoordinates.shape[0]
    for k in range(NumberNodes):
        if(list_cellmutation[k] == 0) :
            CoordinatesNextTipCell = np.zeros((3,1))
            CoordinatesNextTipCell[0] = array_nodescoordinates[k][0]
            CoordinatesNextTipCell[1] = array_nodescoordinates[k][1]
            CoordinatesNextTipCell[2] = array_nodescoordinates[k][2]
            NormNextTipCell = np.sqrt(CoordinatesNextTipCell[0]**2 + CoordinatesNextTipCell[1]**2+ CoordinatesNextTipCell[2]**2)

            if(NormNextTipCell < NormClosestTipCell):
                NormClosestTipCell = NormNextTipCell
                CoordinatesClosestTipCell = CoordinatesNextTipCell

    return CoordinatesClosestTipCell, NormClosestTipCell

# function returning the furthest tip cell coordinates from the first fixed cell
def FurthestTipCell(file_nodescoordinates, file_cellmutation):
    list_cellmutation = MutationStates(file_cellmutation)
    array_nodescoordinates = NodesCoordinates(file_nodescoordinates)

    CoordinatesFurthestTipCell = np.zeros((3,1))
    CoordinatesFurthestTipCell[0] = array_nodescoordinates[1][0]
    CoordinatesFurthestTipCell[1] = array_nodescoordinates[1][1]
    CoordinatesFurthestTipCell[2] = array_nodescoordinates[1][2]
    NormFurthestTipCell = np.sqrt(CoordinatesFurthestTipCell[0]**2 + CoordinatesFurthestTipCell[1]**2+ CoordinatesFurthestTipCell[2]**2)

    NumberNodes = array_nodescoordinates.shape[0]
    for k in range(NumberNodes):
        if(list_cellmutation[k] == 0):
            CoordinatesNextTipCell = np.zeros((3,1))
            CoordinatesNextTipCell[0] = array_nodescoordinates[k][0]
            CoordinatesNextTipCell[1] = array_nodescoordinates[k][1]
            CoordinatesNextTipCell[2] = array_nodescoordinates[k][2]
            NormNextTipCell = np.sqrt(CoordinatesNextTipCell[0]**2 + CoordinatesNextTipCell[1]**2+ CoordinatesNextTipCell[2]**2)

            if(NormNextTipCell > NormFurthestTipCell):
                NormFurthestTipCell = NormNextTipCell
                CoordinatesFurthestTipCell = CoordinatesNextTipCell

    return CoordinatesFurthestTipCell, NormFurthestTipCell

# function returning the average area of the blood vessel tree 
def AverageArea(file_nodescoordinates, file_cellmutation):
    # Method with the half sphere formula
    # furthest_tip_cell, radius = FurthestTipCell(file_nodescoordinates, file_cellmutation)
    # area = 2*np.pi*radius**2 # area of half a sphere 

    # Method with the convex hull 
    array_nodescoordinates = NodesCoordinates(file_nodescoordinates)
    hull = ConvexHull(array_nodescoordinates)

    area = hull.area

    return area 

# function returning the average volume of the blood vessel tree 
def AverageVolume(file_nodescoordinates, file_cellmutation):
    # Method with the half sphere formula 
    #furthest_tip_cell, radius = FurthestTipCell(file_nodescoordinates, file_cellmutation)
    #volume = (2/3)*np.pi*radius**3 # volume of half a sphere 

    # Method with convex hull 
    array_nodescoordinates = NodesCoordinates(file_nodescoordinates)
    hull = ConvexHull(array_nodescoordinates)

    volume = hull.volume

    return volume 

def AverageFunction(list_average, list):
    n = len(list)

    for k in range(n):
        list_average[k] = list_average[k] + list[k]/9

    return list_average

def Tortuosity(file_tortuosity):
    arc, length = LongestPath(file_tortuosity)

    tortuosity = arc/length

    return tortuosity 

def MaxFunction(list_max, list_nbcells):
    for k in range(len(list_max)):
        if(list_nbcells[k] > list_max[k]):
            list_max[k] = list_nbcells[k]

    return list_max

def MinFunction(list_min, list_nbcells):
    for k in range(len(list_min)):
        if(list_nbcells[k] < list_min[k]):
            list_min[k] = list_nbcells[k]

    return list_min