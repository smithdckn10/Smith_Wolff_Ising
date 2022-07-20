import matplotlib.pyplot as plt
from numpy import sqrt
updates = 0
dataE = []
dataM = []
dataAvgE = []
dataAvgE2 = []
dataAvgM = []
dataAvgM2 = []
count = 0
##expE = 0 
##expM = 0 
thermCut = 200  ##This tells the program what sweep to start calculating average values and stuff
  ##200 should be plenty of updates for the system to thermalize

with open('/Users/Smith/Desktop/UnderGradProject/MINE/Smith_Wolff_Energies.txt') as in_file:
    lines = in_file.readlines()
    line = lines[0].split(",", -1)
    print("Size: ", line[0], " Updates: ", line[1], " Temps: ") ##This just takes the first line from the text file and pulls out the size, number of updates, and the temperatures the updates were performed at
    updates = int(line[1])
    count = (len(line) - 2) ##This gets the number of different temperatures
    name = ["" for k in range(count)]
    for x in range(count):
        name[x] = "J = " + str(line[x+2])
    if(updates <= thermCut):
        print("The number of updates is less than the thermCut! No analysis can be performed.")

    #print("count: ",count)
    dataEA = [[0 for x in range(updates+1)] for y in range(count+1)] ##Array that stores the energies after each sweep, along with the "temperature" the sweep was performed at
    dataMA = [[0 for x in range(updates+1)] for y in range(count+1)] ##Array that stores the magnetizations of each sweep, along with the "temp"
    expE = [0 for y in range(count)]  ##This array stores the expectation value of the energy 
    expM = [0 for y in range(count)]  ##This array stores the expectation value of the magnetzation
    expE2 = [0 for y in range(count)]  ##This array stores the expectation value of the energy^2 
    expM2 = [0 for y in range(count)]  ##This array stores the expectation value of the magnetzation^2 
    errE = [0 for y in range(count)] ##This is the error for the estimate of <E>
    stdE = [0 for y in range(count)] ##This is the standard deviation for the estimate of <E>
    errM = [0 for y in range(count)] ##This is the error for the estimate of <M>
    stdM = [0 for y in range(count)] ##This is the standard deviation for the estimate of <M>
    sweepArr = [x for x in range(1,updates+1)] ##Just an array to serve as the x component of the data
    tempArr = [0 for x in range(count)] ##This array stores the temperatures the sims were performed at, and serves as x component of data

    del lines[0] ##Delete the first line in the text file, leaving just the data
    i = 0
    j = 1
    k = 1
    temp = ""
    check = False
    while i < count:
        temp = lines[updates*i].split(",",3)
        while j <= updates:
            if i == count - 1 and j == updates:
                ##print("temp: ", temp[0])
                pass

            line = lines[updates * i + (j-1)].split(",",7)
            ##print(line)
            dataEA[i][j] = (float(line[2]))
            dataMA[i][j] = (float(line[3]))

            ##print("dataEA[",i,"][",j,"]: ",dataEA[i][j],"  dataEA[500][18]: ",dataEA[18][500])
            
            ##print("lol")
            j+=1
        ##print("i: ",i)
        dataEA[i][0] = float(temp[0])
        dataMA[i][0] = float(temp[0])
        

        tempArr[i] = float(temp[0])



        j = 1
        i+=1

    i = 0
    while i < count: ##This loop is for calculating the expected values and their errors
        j = thermCut
        while j <= updates:
            expE[i] = expE[i] + dataEA[i][j] # <E>
            expE2[i] = expE2[i] + dataEA[i][j]*dataEA[i][j] # <E^2>
            expM[i] = expM[i] + dataMA[i][j] # <M> 
            expM2[i] = expM[i] + dataMA[i][j]*dataMA[i][j] # <M^2>
            j+=1

        expE[i] = expE[i] / (updates - thermCut - 1)   ##This makes expE[i][1] = <E> at a given temp (expE[i][0] stores the temp). This is an unbiased estimate, since we are ignoring the first ~200 updates, then adding the energies of the updates, then dividing by the number of updates - 200 - 1. If we didn't include that minus one, it would be a biased estimate I think
        expM[i] = expM[i] / (updates - thermCut - 1)   ##This is <M>, but everything else is the same as expE
        expE2[i] = expE2[i] / (updates - thermCut - 1) ##This is the expectation value of the energy^2 or <E^2>
        expM2[i] = expM2[i] / (updates - thermCut - 1) ##This is <M^2>
        print('expE2[',i,']: ', expE2[i], ', expE[',i,']*expE[',i,']: ', expE[i]*expE[i])
        stdE[i] = sqrt(abs(expE2[i] - (expE[i]*expE[i]))) ##I think theres supposed to be an abs() in here, but I'm not sure
        stdM[i] = sqrt(abs(expM2[i] - (expM[i]*expM[i])))
        errE[i] = stdE[i] / sqrt((updates - thermCut - 1))  ##I'm not sure if were supposed to subtract one like we did for the expectation value, but I'll leave it for now
        errM[i] = stdM[i] / sqrt((updates - thermCut - 1))
        i+=1
    i = 0
    
    ##print("dataEA[1]: ",dataEA[1])
    ##print("dataEA[18][500]: ",dataEA[18][500])
    fig = plt.figure(figsize = (10,8))
    eGraf = fig.add_subplot(221)
    plt.ylabel('Energy Density')
    plt.xlabel('Update')
    eAx = fig.add_subplot(223)
    plt.xlabel('Temp Modeled ( J )')
    plt.ylabel('<E>')
    mGraf = fig.add_subplot(222)
    plt.ylabel('Magnetization')
    plt.xlabel('Update')
    mAx = fig.add_subplot(224)
    plt.xlabel('Temp Modeled ( J )')
    plt.ylabel('<M>')
    
    for i in range(count):
        del dataEA[i][0]
        del dataMA[i][0] 
        eGraf.plot(sweepArr, dataEA[i], label = name[i])
        mGraf.plot(sweepArr, dataMA[i], label = name[i])
        #eGraf.plot(sweepArr[thermCut:], dataEA[i][thermCut:], label = name[i])
        #mGraf.plot(sweepArr[thermCut:], dataMA[i][thermCut:], label = name[i])
    i = 0
    



    eGraf.legend(ncol = 3, loc = "upper center", fontsize = 'small', borderaxespad = -4.0)
    mGraf.legend(ncol = 3, loc = "upper center", fontsize = 'small', borderaxespad = -4.0)
    fig.suptitle('Energy Density & Magnetization for Varying J (Using Wolff Cluster Algorithm)') ##Fix This
    eAx.scatter(tempArr, expE)
    eAx.errorbar(tempArr, expE, yerr= errE, fmt='o--')
    #eAx.ylabel('<E>')
    mAx.scatter(tempArr, expM)
    mAx.errorbar(tempArr, expM, yerr= errM, fmt='o--')
    #mAx.xlabel('Temp J')
    #mAx.ylabel('<M>')

    for x,y in zip(tempArr,expE):

        labelE = f"σ = {round(errE[i],5)}"
        eAx.annotate(labelE, (x,y), textcoords="offset points", xytext=(0,12), ha='center')
        i+=1
    i=0
    for x,y in zip(tempArr,expM):

        labelM = f"σ = {round(errM[i],5)}"
        mAx.annotate(labelM, (x,y), textcoords="offset points", xytext=(0,12), ha='center')
        i+=1
    i=0
    plt.figtext(.5,0.02, 'Note: The expectation values were calculated ignoring the first ' + str(thermCut) + ' updates', transform = fig.transFigure, horizontalalignment = 'center')
    plt.show()

    """
    fig = plt.figure(figsize = (10,8))
    eGraf = fig.add_subplot(221)
    eAx = fig.add_subplot(223)
    mGraf = fig.add_subplot(224)
    mAx = fig.add_subplot(122)
    fig.suptitle('Horizontally stacked subplots')
    eAx.scatter(tempArr, expE)
    eAx.errorbar(tempArr, expE, yerr= errE, fmt='o--')
    mAx.scatter(tempArr, expM)
    mAx.errorbar(tempArr, expM, yerr= errM, fmt='o--')
    for x,y in zip(tempArr,expE):

        labelE = f"σ = {round(errE[i],5)}"
        eAx.annotate(labelE, (x,y), textcoords="offset points", xytext=(0,12), ha='center')
        i+=1
    i=0
    for x,y in zip(tempArr,expM):

        labelM = f"σ = {round(errM[i],5)}"
        mAx.annotate(labelM, (x,y), textcoords="offset points", xytext=(0,12), ha='center')
        i+=1
    i=0
    plt.show()
    """
    """
    plt.figure(0)
    plt.figure(figsize = (8,5))
    plt.scatter(tempArr, expE)
    plt.errorbar(tempArr, expE, yerr= errE, fmt='o--')
    # zip joins x and y coordinates in pairs
    for x,y in zip(tempArr,expE):

        label = f"σ = {round(errE[i],5)}"

        plt.annotate(label, # this is the text
                    (x,y), # these are the coordinates to position the label
                    textcoords="offset points", # how to position the text
                    xytext=(0,12), # distance from text to points (x,y)
                    ha='center')
        i+=1
    plt.title('Expectation value of E for varying J')
    plt.xlabel('J value the sim was performed at')
    plt.ylabel('Expectation Value of E')
    #leg = plt.legend()
    i = 0
    plt.figure(1)
    plt.figure(figsize = (8,5))
    plt.scatter(tempArr, expM)
    plt.errorbar(tempArr, expM, yerr= errM, fmt='o--')
    # zip joins x and y coordinates in pairs
    for x,y in zip(tempArr,expM):

        label = f"σ = {round(errM[i],5)}"

        plt.annotate(label, # this is the text
                    (x,y), # these are the coordinates to position the label
                    textcoords="offset points", # how to position the text
                    xytext=(0,12), # distance from text to points (x,y)
                    ha='center')
        i+=1
    plt.title('Expectation value of M for varying J')
    plt.xlabel('J value the sim was performed at')
    plt.ylabel('Expectation Value of M')

    plt.show()

    """
    """"
    plt.figure(figsize=(10,8))
    for i in range(count):
        del dataEA[i][0]
        del dataMA[i][0] 
        plt.plot(sweepArr, dataMA[i], label = name[i]) ##Just change the dataEA to dataMA to plot the magnetizations (and change the title of the plot)
    ##plt.plot(sweepArr, dataEA[18])
    ##plt.plot(dataM)
    plt.title('Comparing Magnetization For J in [.01 , .5]')
    plt.xlabel('Sweep')
    plt.ylabel('Magnetization After Each Sweep[idk units]')
    leg = plt.legend(ncol = 4) ##I did some digging on how to make the legend not so awful, but the only thing I was actually able to implement was changing the number of columns in the legend :(

    plt.show()
    """