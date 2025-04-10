# This is a class that analysis quantitatively the angiogenesis model
# you can use this class to:
# - read a file 
# - extract the information 

import numpy as np
import statistics as stats

from scipy.spatial import ConvexHull

class ParametersSensitivityRunner:
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
    def NodesCoordinates(file_nodescoordinates, dim):
        # initialisation of the list
        list_nodescoordinates = []

        # we open the file
        f = open(file_nodescoordinates, 'r')

        for line in f:
            pass
        last_line = line

        list_nodescoordinates = [float(x) for x in last_line.split()[1:]]

        NumberNodes = int((len(list_nodescoordinates))/dim)
        array_nodescoordinates = np.reshape(list_nodescoordinates, (NumberNodes,dim))

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

    def LongestPath(file_tortuosity, dim):
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
        array_tipcellscoordinates = np.reshape(list_tipcellscoordinates, (NbTipCells,dim))

        # we look for the longest arc 
        arc = list_arcs[0]
        LongestArcNb = 0
        for k in range(len(list_arcs)):
            if(list_arcs[k] > arc):
                arc = list_arcs[k]
                LongestArcNb = k

        # we collect the coordinates of the tip cell associated to the longest arc 
        if(dim==3):
            TipCellLength_x = array_tipcellscoordinates[LongestArcNb][0]
            TipCellLength_y = array_tipcellscoordinates[LongestArcNb][1]
            TipCellLength_z = array_tipcellscoordinates[LongestArcNb][2]

            # we calculate its norm 
            length = np.sqrt(TipCellLength_x**2 + TipCellLength_y**2 + TipCellLength_z**2)
        elif(dim==2):
            TipCellLength_x = array_tipcellscoordinates[LongestArcNb][0]
            TipCellLength_y = array_tipcellscoordinates[LongestArcNb][1]

            # we calculate its norm 
            length = np.sqrt(TipCellLength_x**2 + TipCellLength_y**2)
        elif(dim==1):
            TipCellLength_x = array_tipcellscoordinates[LongestArcNb][0]

            # we calculate its norm 
            length = np.sqrt(TipCellLength_x**2)

        return arc, length 


    # a function to obtain the total number of cells at the end of the simulation 
    def TotalNumberCells(file_cellmutation):
        list_cellmutation = ParametersSensitivityRunner.MutationStates(file_cellmutation)
        return len(list_cellmutation)

    # a function to obtain the total number of tip cells at the end of the simulation 
    def TotalNumberTipCells(file_cellmutation):
        NumberTipCells = 0

        list_cellmutation = ParametersSensitivityRunner.MutationStates(file_cellmutation)

        for elem in list_cellmutation:
            if(int(elem) == 0):
                NumberTipCells += 1

        return NumberTipCells

    # a function to obtain the total number of stalk cells at the end of the simulation (including branching points)
    def TotalNumberStalkCells(file_cellmutation):
        NumberStalkCells = 0

        list_cellmutation = ParametersSensitivityRunner.MutationStates(file_cellmutation)

        for elem in list_cellmutation:
            if(elem != 0):
                NumberStalkCells += 1

        return NumberStalkCells

    # a function to obtain the total number of branches at the end of the simulation 
    def TotalNumberBranches(file_cellmutation):
        NumberBranches = 0

        list_cellmutation = ParametersSensitivityRunner.MutationStates(file_cellmutation)

        for elem in list_cellmutation:
            if(elem == 2):
                NumberBranches += 1
            elif(elem == 0):
                NumberBranches += 1

        return NumberBranches
    
    # a function to obtain the number of cells at each time step t
    def NbCellsAtEachTimeStep(file_cellmutation):
        # initialisation
        list_nbcellsateachtimestep = []

        # we open the file
        f = open(file_cellmutation, 'r')

        for line in f:
            list_nbcellsattimestept = [float(x) for x in line.split()[1:]]
            nbcellsattimestept = len(list_nbcellsattimestept)
            list_nbcellsateachtimestep.append(nbcellsattimestept)

        # we close the file 
        f.close()

        return list_nbcellsateachtimestep
    
    # a function to obtain the number of cells at each time step t
    def NbCellsInPlaneAtEachTimeStep(file_nodescoordinates,plane_normal_vector, dim):
        # initialisation
        list_nbcellsinplaneateachtimestep = []
        eps = 0.5

        # we open the file
        f = open(file_nodescoordinates, 'r')

        for line in f:
            array_nbcellsinplaneattimestept = []
            list_nbcellsinplaneattimestept = [float(x) for x in line.split()[1:]]
            NumberNodesTimeT = int((len(list_nbcellsinplaneattimestept))/dim)
            array_nodescoordinates = np.reshape(list_nbcellsinplaneattimestept, (NumberNodesTimeT,dim))

            for k in range(NumberNodesTimeT):
                if(array_nodescoordinates[k][0] > plane_normal_vector[0]-eps and array_nodescoordinates[k][0] < plane_normal_vector[0] + eps):
                    array_nbcellsinplaneattimestept.append(array_nodescoordinates[k])
            
            nbcellsinplaneattimestept = len(array_nbcellsinplaneattimestept)

            #nbcellsinplaneattimestept = len(array_nbcellsinplaneattimestept)/ParametersSensitivityRunner.AreaCellsPlane(file_nodescoordinates, plane_normal_vector)
            list_nbcellsinplaneateachtimestep.append(nbcellsinplaneattimestept)

        # we close the file 
        f.close()

        return list_nbcellsinplaneateachtimestep
    
    # a function to obtain all the nodes in a specific plane (take into account the radius of the cell i.e. introduce an epsilon)
    def CellsInPlane(file_nodescoordinates, plane_normal_vector, dim):
        # initialisation 
        array_nodesplanecoordinates = []
        eps = 0.5
        array_nodescoordinates = ParametersSensitivityRunner.NodesCoordinates(file_nodescoordinates, dim)
        TotalNbCells = len(array_nodescoordinates)

        # Generalised method not working for now 
        # for k in range(TotalNbCells):
        #     MinPlane = plane_normal_vector[0]*(array_nodescoordinates[k][0]-eps) + plane_normal_vector[1]*(array_nodescoordinates[k][1]-eps) + plane_normal_vector[2]*(array_nodescoordinates[k][2]-eps) + plane_normal_vector[3]
        #     MaxPlane = plane_normal_vector[0]*(array_nodescoordinates[k][0]+eps) + plane_normal_vector[1]*(array_nodescoordinates[k][1]+eps) + plane_normal_vector[2]*(array_nodescoordinates[k][2]+eps) + plane_normal_vector[3]

        #     if(MinPlane <= 1e-6 and MaxPlane >= 1e-6):
        #         array_nodesplanecoordinates.append(array_nodescoordinates[k])

        # method corresponding to the chemotactic force in the x axis (to modify later)
        for k in range(TotalNbCells):
            if(array_nodescoordinates[k][0] > plane_normal_vector[0]-eps and array_nodescoordinates[k][0] < plane_normal_vector[0] + eps):
                array_nodesplanecoordinates.append(array_nodescoordinates[k])

        return array_nodesplanecoordinates
    
    # a function to obtain the number of cells for x > x_ref (in 2D only for now!!!)
    def NbCellsAfterPlane(file_nodescoordinates, ref_point, dim):
        # Initialisation
        NbCellsAfterPlane = 0
        array_nodescoordinates = ParametersSensitivityRunner.NodesCoordinates(file_nodescoordinates, dim)
        TotalNbCells = len(array_nodescoordinates)

        # select nodes xith x coordinates inferior to ref_point
        for k in range(TotalNbCells):
            if(array_nodescoordinates[k][0] < ref_point):
                NbCellsAfterPlane += 1

        return NbCellsAfterPlane
    
    # a function to obtain the number of cells in a specific plane 
    def NbCellsPlane(file_nodescoordinates, plane_normal_vector, dim):
        # Initialisation
        NbCellsPlane = 0
        array_nodesplanecoordinates = ParametersSensitivityRunner.CellsInPlane(file_nodescoordinates, plane_normal_vector, dim)

        NbCellsPlane = len(array_nodesplanecoordinates)

        #NbCellsPlane = len(array_nodesplanecoordinates)/ParametersSensitivityRunner.AreaCellsPlane(file_nodescoordinates,plane_normal_vector)

        return NbCellsPlane
    
    # a function to obtain the number of cells at each time step t
    # ATTENTION: ONLY FOR 2D FOR THE MOMENT 
    def TimeFirstReachingPlane(file_nodescoordinates, x_ref, dim):
        # initialisation
        TimeFirstReachingPlane = -1
        eps = 0.5

        # we open the file
        f = open(file_nodescoordinates, 'r')

        for line in f:
            t = float(line.split()[0])
            list_nbcellsinplaneattimestept = [float(x) for x in line.split()[1:]]
            NumberNodesTimeT = int((len(list_nbcellsinplaneattimestept))/dim)
            array_nodescoordinates = np.reshape(list_nbcellsinplaneattimestept, (NumberNodesTimeT,dim))

            for k in range(NumberNodesTimeT):
                if(array_nodescoordinates[k][0] < x_ref-eps):
                    TimeFirstReachingPlane = t
                    # we close the file 
                    f.close()

                    return TimeFirstReachingPlane
                
        f.close()
        return TimeFirstReachingPlane
    
    # a function to obtain the closest neighbouring distance between cells in a specific plane
    def ClosestNeighbouringDistanceInPlane(file_nodescoordinates, plane_normal_vector, dim):
        # Initialisation 
        list_neighbouringdistanceinplane = []
        list_closestneighbouringdistanceinplane = []
        array_nodesplanecoordinates = ParametersSensitivityRunner.CellsInPlane(file_nodescoordinates, plane_normal_vector, dim)

        for k in range(len(array_nodesplanecoordinates)):
            for j in range(len(array_nodesplanecoordinates)):
                if(j != k):
                    neighbouringdistance = np.linalg.norm(array_nodesplanecoordinates[k]-array_nodesplanecoordinates[j])
                    list_neighbouringdistanceinplane.append(neighbouringdistance)
            closestneighbouringdistanceinplane = min(list_neighbouringdistanceinplane)
            list_closestneighbouringdistanceinplane.append(closestneighbouringdistanceinplane)
        
        AverageClosestNeighbouringDistanceInPlane = stats.mean(list_closestneighbouringdistanceinplane)*100

        return AverageClosestNeighbouringDistanceInPlane
    
    # a function to obtain the closest neighbouring distance between cells in a specific plane at each time step
    def ClosestNeighbouringDistanceInPlaneAtEachTimeStep(file_nodescoordinates, plane_normal_vector, dim):
        # Initialisation 
        list_averageclosestneighbouringdistanceinplane = []
        array_nodesplanecoordinatesateachtimestep = []
        list_neighbouringdistanceinplane = []
        list_closestneighbouringdistanceinplane = []
        eps = 0.5
        counter = 0

        # we open the file
        f = open(file_nodescoordinates, 'r')

        for line in f:
            array_nbcellsinplaneattimestept = []
            list_neighbouringdistanceinplane = []
            list_nbcellsinplaneattimestept = [float(x) for x in line.split()[1:]]
            NumberNodesTimeT = int((len(list_nbcellsinplaneattimestept))/dim)
            array_nodescoordinates = np.reshape(list_nbcellsinplaneattimestept, (NumberNodesTimeT,dim))

            for k in range(NumberNodesTimeT):
                if(abs(array_nodescoordinates[k][0]-eps) < plane_normal_vector[0]):
                    array_nodesplanecoordinatesateachtimestep.append(array_nodescoordinates[k])

            for k in range(len(array_nodesplanecoordinatesateachtimestep)):
                for j in range(len(array_nodesplanecoordinatesateachtimestep)):
                    if(j != k):
                        neighbouringdistance = np.linalg.norm(array_nodesplanecoordinatesateachtimestep[k]-array_nodesplanecoordinatesateachtimestep[j])
                        list_neighbouringdistanceinplane.append(neighbouringdistance)
                closestneighbouringdistanceinplane = min(list_neighbouringdistanceinplane)
                list_closestneighbouringdistanceinplane.append(closestneighbouringdistanceinplane)

            averageneighbouringdistanceinplaneattimestept = stats.mean(list_closestneighbouringdistanceinplane)*100
            
            list_averageclosestneighbouringdistanceinplane.append(averageneighbouringdistanceinplaneattimestept)

            counter += 1
        # we close the file 
        f.close()

        return list_averageclosestneighbouringdistanceinplane
    
        
    # a function to obtain the average neighbouring distance between cells in a specific plane
    def AverageNeighbouringDistanceInPlane(file_nodescoordinates, plane_normal_vector, dim):
        # Initialisation 
        AverageNeighbouringDistanceInPlane = 0
        list_neighbouringdistanceinplane = []
        array_nodesplanecoordinates = ParametersSensitivityRunner.CellsInPlane(file_nodescoordinates, plane_normal_vector, dim)

        for k in range(len(array_nodesplanecoordinates)):
            for j in range(len(array_nodesplanecoordinates)):
                neighbouringdistance = np.linalg.norm(array_nodesplanecoordinates[k]-array_nodesplanecoordinates[j])
                list_neighbouringdistanceinplane.append(neighbouringdistance)

        AverageNeighbouringDistanceInPlane = stats.mean(list_neighbouringdistanceinplane)

        return AverageNeighbouringDistanceInPlane
    
    # a function to obtain the average neighbouring distance between cells in a specific plane at each time step
    def AverageNeighbouringDistanceInPlaneAtEachTimeStep(file_nodescoordinates, plane_normal_vector, dim):
        # initialisation
        list_averageneighbouringdistanceinplaneateachtimestep = []
        eps = 0.5

        # we open the file
        f = open(file_nodescoordinates, 'r')

        for line in f:
            array_nbcellsinplaneattimestept = []
            list_neighbouringdistanceinplane = []
            list_nbcellsinplaneattimestept = [float(x) for x in line.split()[1:]]
            NumberNodesTimeT = int((len(list_nbcellsinplaneattimestept))/dim)
            array_nodescoordinates = np.reshape(list_nbcellsinplaneattimestept, (NumberNodesTimeT,dim))

            for k in range(NumberNodesTimeT):
                if(abs(array_nodescoordinates[k][0]-eps) < plane_normal_vector[0]):
                    array_nbcellsinplaneattimestept.append(array_nodescoordinates[k])

            for k in range(len(array_nbcellsinplaneattimestept)):
                for j in range(len(array_nbcellsinplaneattimestept)):
                    neighbouringdistance = np.linalg.norm(array_nbcellsinplaneattimestept[k]-array_nbcellsinplaneattimestept[j])
                    list_neighbouringdistanceinplane.append(neighbouringdistance)

            averageneighbouringdistanceinplaneattimestept = stats.mean(list_neighbouringdistanceinplane)
            
            list_averageneighbouringdistanceinplaneateachtimestep.append(averageneighbouringdistanceinplaneattimestept)

        # we close the file 
        f.close()

        return list_averageneighbouringdistanceinplaneateachtimestep

    # a function to obtain the max neighbouring distance between cells in a specific plane
    def MaxNeighbouringDistanceInPlane(file_nodescoordinates, plane_normal_vector, dim):
        # Initialisation 
        MaxNeighbouringDistanceInPlane = 0
        list_neighbouringdistanceinplane = []
        array_nodesplanecoordinates = ParametersSensitivityRunner.CellsInPlane(file_nodescoordinates, plane_normal_vector, dim)

        for k in range(len(array_nodesplanecoordinates)):
            for j in range(len(array_nodesplanecoordinates)):
                neighbouringdistance = np.linalg.norm(array_nodesplanecoordinates[k]-array_nodesplanecoordinates[j])
                list_neighbouringdistanceinplane.append(neighbouringdistance)

        MaxNeighbouringDistanceInPlane = max(list_neighbouringdistanceinplane)

        return MaxNeighbouringDistanceInPlane
    
    # a function to obtain the area occupied by the cells in a specific plane
    # use of convex hull
    def AreaCellsPlane(file_nodescoordinates, plane_normal_vector, dim):
        # Initialisation
        AreaCellsPlane = 0
        eps = 0.5

        array_nodesplanecoordinates = ParametersSensitivityRunner.CellsInPlane(file_nodescoordinates, plane_normal_vector, dim)

        if (len(array_nodesplanecoordinates) == 1):
            AreaCellsPlane = 1
        elif (len(array_nodesplanecoordinates) == 0):
            AreaCellsPlane = 1
        else:
            miny = min(array_nodesplanecoordinates[0][1], array_nodesplanecoordinates[1][1])-eps
            maxy = max(array_nodesplanecoordinates[0][1], array_nodesplanecoordinates[1][1])+eps
            minz = min(array_nodesplanecoordinates[0][2], array_nodesplanecoordinates[1][2])-eps
            maxz = max(array_nodesplanecoordinates[0][2], array_nodesplanecoordinates[1][2])+eps
            
            for k in range(2, len(array_nodesplanecoordinates)):
                y = array_nodesplanecoordinates[k][1]
                z = array_nodesplanecoordinates[k][2]

                if(y < miny):
                    miny = y
                if(y > maxy):
                    maxy = y
                if(z < minz):
                    minz = z
                if(z > maxz):
                    maxz = z

            AreaCellsPlane = (maxy-miny)*(maxz-minz)

        return AreaCellsPlane

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
        time = 300
        Psprout = 0.1 # WARNING
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
        time = 300 # WARNING

        # branches obtained 
        NumberMaximumConsecutiveBranches = ParametersSensitivityRunner.TotalNumberBranches(file_cellmutation)

        # expected branches code 
        NumberExpectedConsecutiveBranches = 1

        # formula of the expectation 
        NumberExpectedConsecutiveBranches = ParametersSensitivityRunner.ComputationalConsecutiveBranches(Tcycle)

        NumberFormulaExpectedConsecutiveBranches = 2*(1+Psprout)**(np.floor(time/Tcycle)) - 1
        
        return NumberMaximumConsecutiveBranches, NumberExpectedConsecutiveBranches, NumberFormulaExpectedConsecutiveBranches

    # function returning the closest tip cell coordinates from the first fixed cell
    def ClosestTipCell(file_nodescoordinates, file_cellmutation, dim):
        list_cellmutation = ParametersSensitivityRunner.MutationStates(file_cellmutation)
        array_nodescoordinates = ParametersSensitivityRunner.NodesCoordinates(file_nodescoordinates, dim)

        CoordinatesClosestTipCell = np.zeros((dim,1))
        if (dim == 3):
            CoordinatesClosestTipCell[0] = array_nodescoordinates[1][0]
            CoordinatesClosestTipCell[1] = array_nodescoordinates[1][1]
            CoordinatesClosestTipCell[2] = array_nodescoordinates[1][2]
            NormClosestTipCell = np.sqrt(CoordinatesClosestTipCell[0]**2 + CoordinatesClosestTipCell[1]**2+ CoordinatesClosestTipCell[2]**2)
        elif(dim == 2):
            CoordinatesClosestTipCell[0] = array_nodescoordinates[1][0]
            CoordinatesClosestTipCell[1] = array_nodescoordinates[1][1]
            NormClosestTipCell = np.sqrt(CoordinatesClosestTipCell[0]**2 + CoordinatesClosestTipCell[1]**2)
        elif(dim == 1):
            CoordinatesClosestTipCell[0] = array_nodescoordinates[1][0]
            NormClosestTipCell = np.sqrt(CoordinatesClosestTipCell[0]**2)

        NumberNodes = array_nodescoordinates.shape[0]
        for k in range(NumberNodes):
            if(list_cellmutation[k] == 0) :
                if(dim == 3):
                    CoordinatesNextTipCell = np.zeros((dim,1))
                    CoordinatesNextTipCell[0] = array_nodescoordinates[k][0]
                    CoordinatesNextTipCell[1] = array_nodescoordinates[k][1]
                    CoordinatesNextTipCell[2] = array_nodescoordinates[k][2]
                    NormNextTipCell = np.sqrt(CoordinatesNextTipCell[0]**2 + CoordinatesNextTipCell[1]**2+ CoordinatesNextTipCell[2]**2)
                elif(dim == 2):
                    CoordinatesNextTipCell = np.zeros((dim,1))
                    CoordinatesNextTipCell[0] = array_nodescoordinates[k][0]
                    CoordinatesNextTipCell[1] = array_nodescoordinates[k][1]
                    NormNextTipCell = np.sqrt(CoordinatesNextTipCell[0]**2 + CoordinatesNextTipCell[1]**2)
                elif(dim == 1):
                    CoordinatesNextTipCell = np.zeros((dim,1))
                    CoordinatesNextTipCell[0] = array_nodescoordinates[k][0]
                    NormNextTipCell = np.sqrt(CoordinatesNextTipCell[0]**2)

                if(NormNextTipCell < NormClosestTipCell):
                    NormClosestTipCell = NormNextTipCell
                    CoordinatesClosestTipCell = CoordinatesNextTipCell

        return CoordinatesClosestTipCell, NormClosestTipCell

    # function returning the furthest tip cell coordinates from the first fixed cell
    # ATTENTION FUNCTION FOR A FIXED FIRST CELL AT POSITION (25, 50, 50)
    def FurthestTipCell(file_nodescoordinates, file_cellmutation, dim):
        list_cellmutation = ParametersSensitivityRunner.MutationStates(file_cellmutation)
        array_nodescoordinates = ParametersSensitivityRunner.NodesCoordinates(file_nodescoordinates, dim)

        NumberNodes = array_nodescoordinates.shape[0]
        counter = 0
        while counter < NumberNodes:
            if(list_cellmutation[counter] == 0):
                if(dim==3):
                    CoordinatesFurthestTipCell = np.zeros((dim,1))
                    CoordinatesFurthestTipCell[0] = array_nodescoordinates[1][0]-25
                    CoordinatesFurthestTipCell[1] = array_nodescoordinates[1][1]-50
                    CoordinatesFurthestTipCell[2] = array_nodescoordinates[1][2]-50
                    NormFurthestTipCell = np.sqrt(CoordinatesFurthestTipCell[0]**2 + CoordinatesFurthestTipCell[1]**2+ CoordinatesFurthestTipCell[2]**2)
                elif(dim==2):
                    CoordinatesFurthestTipCell = np.zeros((dim,1))
                    CoordinatesFurthestTipCell[0] = array_nodescoordinates[1][0]-25
                    CoordinatesFurthestTipCell[1] = array_nodescoordinates[1][1]-50
                    NormFurthestTipCell = np.sqrt(CoordinatesFurthestTipCell[0]**2 + CoordinatesFurthestTipCell[1]**2)
                elif(dim==1):
                    CoordinatesFurthestTipCell = np.zeros((dim,1))
                    CoordinatesFurthestTipCell[0] = array_nodescoordinates[1][0]-25
                    NormFurthestTipCell = np.sqrt(CoordinatesFurthestTipCell[0]**2)
            counter += 1

        for k in range(NumberNodes):
            if(list_cellmutation[k] == 0):
                if(dim==3):
                    CoordinatesNextTipCell = np.zeros((dim,1))
                    CoordinatesNextTipCell[0] = array_nodescoordinates[k][0]-25
                    CoordinatesNextTipCell[1] = array_nodescoordinates[k][1]-50
                    CoordinatesNextTipCell[2] = array_nodescoordinates[k][2]-50
                    NormNextTipCell = np.sqrt(CoordinatesNextTipCell[0]**2 + CoordinatesNextTipCell[1]**2+ CoordinatesNextTipCell[2]**2)
                elif(dim==2):
                    CoordinatesNextTipCell = np.zeros((dim,1))
                    CoordinatesNextTipCell[0] = array_nodescoordinates[k][0]-25
                    CoordinatesNextTipCell[1] = array_nodescoordinates[k][1]-50
                    NormNextTipCell = np.sqrt(CoordinatesNextTipCell[0]**2 + CoordinatesNextTipCell[1]**2)
                elif(dim==1):
                    CoordinatesNextTipCell = np.zeros((dim,1))
                    CoordinatesNextTipCell[0] = array_nodescoordinates[k][0]-25
                    NormNextTipCell = np.sqrt(CoordinatesNextTipCell[0]**2)

                if(NormNextTipCell > NormFurthestTipCell):
                    NormFurthestTipCell = NormNextTipCell
                    CoordinatesFurthestTipCell = CoordinatesNextTipCell

        return CoordinatesFurthestTipCell, NormFurthestTipCell
    
    def NormFurthestTipCell(file_nodescoordinates, file_cellmutation, dim):
        list_cellmutation = ParametersSensitivityRunner.MutationStates(file_cellmutation)
        array_nodescoordinates = ParametersSensitivityRunner.NodesCoordinates(file_nodescoordinates, dim)

        NumberNodes = array_nodescoordinates.shape[0]
        counter = 0
        NormFurthestTipCell = 0
        while counter < NumberNodes:
            if(list_cellmutation[counter] == 0):
                if(dim==3):
                    CoordinatesFurthestTipCell = np.zeros(dim)
                    CoordinatesFurthestTipCell[0] = array_nodescoordinates[1][0]-25
                    CoordinatesFurthestTipCell[1] = array_nodescoordinates[1][1]-50
                    CoordinatesFurthestTipCell[2] = array_nodescoordinates[1][2]-50
                    NormFurthestTipCell = np.sqrt(CoordinatesFurthestTipCell[0]**2 + CoordinatesFurthestTipCell[1]**2+ CoordinatesFurthestTipCell[2]**2)
                elif(dim==2):
                    CoordinatesFurthestTipCell = np.zeros(dim)
                    CoordinatesFurthestTipCell[0] = array_nodescoordinates[1][0]-25
                    CoordinatesFurthestTipCell[1] = array_nodescoordinates[1][1]-50
                    NormFurthestTipCell = np.sqrt(CoordinatesFurthestTipCell[0]**2 + CoordinatesFurthestTipCell[1]**2)
                elif(dim==1):
                    CoordinatesFurthestTipCell = np.zeros(dim)
                    CoordinatesFurthestTipCell[0] = array_nodescoordinates[1][0]-25
                    NormFurthestTipCell = np.sqrt(CoordinatesFurthestTipCell[0]**2)
            counter += 1

        for k in range(NumberNodes):
            if(list_cellmutation[k] == 0):
                if(dim==3):
                    CoordinatesNextTipCell = np.zeros(dim)
                    CoordinatesNextTipCell[0] = array_nodescoordinates[k][0]-25
                    CoordinatesNextTipCell[1] = array_nodescoordinates[k][1]-50
                    CoordinatesNextTipCell[2] = array_nodescoordinates[k][2]-50
                    NormNextTipCell = np.sqrt(CoordinatesNextTipCell[0]**2 + CoordinatesNextTipCell[1]**2+ CoordinatesNextTipCell[2]**2)
                elif(dim==2):
                    CoordinatesNextTipCell = np.zeros(dim)
                    CoordinatesNextTipCell[0] = array_nodescoordinates[k][0]-25
                    CoordinatesNextTipCell[1] = array_nodescoordinates[k][1]-50
                    NormNextTipCell = np.sqrt(CoordinatesNextTipCell[0]**2 + CoordinatesNextTipCell[1]**2)
                elif(dim==1):
                    CoordinatesNextTipCell = np.zeros(dim)
                    CoordinatesNextTipCell[0] = array_nodescoordinates[k][0]-25
                    NormNextTipCell = np.sqrt(CoordinatesNextTipCell[0]**2)

                if(NormNextTipCell > NormFurthestTipCell):
                    NormFurthestTipCell = NormNextTipCell
                    CoordinatesFurthestTipCell = CoordinatesNextTipCell

        return NormFurthestTipCell

    # function returning the average area of the blood vessel tree 
    def AverageArea(file_nodescoordinates, file_cellmutation):
        # Method with the half sphere formula
        # furthest_tip_cell, radius = FurthestTipCell(file_nodescoordinates, file_cellmutation)
        # area = 2*np.pi*radius**2 # area of half a sphere 

        # Method with the convex hull 
        array_nodescoordinates = ParametersSensitivityRunner.NodesCoordinates(file_nodescoordinates, dim)
        hull = ConvexHull(array_nodescoordinates)

        area = hull.area

        return area 

    # function returning the average volume of the blood vessel tree 
    def AverageVolume(file_nodescoordinates, file_cellmutation):
        # Method with the half sphere formula 
        #furthest_tip_cell, radius = FurthestTipCell(file_nodescoordinates, file_cellmutation)
        #volume = (2/3)*np.pi*radius**3 # volume of half a sphere 

        # Method with convex hull 
        array_nodescoordinates = ParametersSensitivityRunner.NodesCoordinates(file_nodescoordinates, dim)
        hull = ConvexHull(array_nodescoordinates)

        volume = hull.volume

        return volume 

    def AverageFunction(list_average, list, TotalTestNb):
        n = len(list)

        for k in range(n):
            list_average[k] = list_average[k] + list[k]/TotalTestNb

        return list_average

    def Tortuosity(file_tortuosity):
        arc, length = ParametersSensitivityRunner.LongestPath(file_tortuosity, dim)

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