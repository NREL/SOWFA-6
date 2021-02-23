# Function for processing SOWFA source term data.
#
#
# Matt Churchfield 
# National Renewable Energy Laboratory
# 22 February 2021










# Figure out how many time directories there are within a directory, and put
# them in numerical order.
def getOutputTimes(dir):
  # Import necessary modules
  import os
  
  
  data = os.listdir(dir)
  outputTimesI = []
  outputTimes = []
  nTimes = len(data)
  ii = 0
  for i in range(nTimes):
     if data[i][0].isnumeric():
        outputTimesI.append(data[i])
        ii = ii + 1
      
  nTimes = len(outputTimesI)
  
  outputTimesIndex = 0
  outputTimesSort = []
  for i in range(nTimes):
     outputTimesSort.append([i,float(outputTimesI[i])])
     
  outputTimesSort = sorted(outputTimesSort, key=lambda index: index[1])
  
  for i in range(nTimes):
     outputTimes.append(outputTimesI[outputTimesSort[i][0]])
  
  
  return nTimes, outputTimes

  
  
  
  
  
  
  
  
  
# Assemble a complete source history.
def assembleSourceHistory(inputDir):
    # Import necessary modules
    import numpy as np
  
  
    # Get the number of time directories and their names.
    [nTimes, outputTimes] = getOutputTimes(inputDir)
    
    
    # Initialize the big arrays.
    timeMomentum = []
    timeTemperature = []
    sourceMomentum = []
    sourceTemperature = []
  
  
    # Loop through the time directories and get the source information.
    for n in range(nTimes):
        sourceName = ['SourceMomentumHistory','SourceTemperatureHistory']
          
        inputFile = inputDir + '/' + outputTimes[n] + '/' + sourceName[0]
  
        heightMomentum, timeMomentumI, sourceMomentumI = readSourceHistoryFile(inputFile)
           
        if (n == 0):
            timeMomentum = timeMomentumI
            sourceMomentum = sourceMomentumI
        else:
            startTime = timeMomentumI[0]
  
            l = len(timeMomentum)
  
            if (timeMomentum[l-1] > startTime):
                indEnd = (np.where(timeMomentum >= startTime))[0][0] - 1
  
            else:
                indEnd = l-1
                 
            timeMomentum = np.append(timeMomentum[0:indEnd],timeMomentumI)
            sourceMomentum = np.append(sourceMomentum[0:indEnd][:],sourceMomentumI,axis=0)
               
            timeMomentumI = []
            sourceMomentumI = []
               
  
        inputFile = inputDir + '/' + outputTimes[n] + '/' + sourceName[1]

        [heightTemperature, timeTemperatureI, sourceTemperatureI] = readSourceHistoryFile(inputFile)  
  
        if (n == 0):
            timeTemperature = timeTemperatureI
            sourceTemperature = sourceTemperatureI
        else:
            startTime = timeTemperatureI[0]
  
            l = len(timeTemperature)
  
            if (timeTemperature[l-1] > startTime):
                indEnd = (np.where(timeTemperature >= startTime))[0][0] - 1
  
            else:
                indEnd = l-1
  
            timeTemperature = np.append(timeTemperature[0:indEnd],timeTemperatureI)
            sourceTemperature = np.append(sourceTemperature[0:indEnd][:],sourceTemperatureI,axis=0)
  
            timeTemperatureI = []
            sourceTemperatureI = []
    
  
    return heightMomentum,heightTemperature,timeMomentum,timeTemperature,sourceMomentum,sourceTemperature
  
  
  
  
  
  

  
  
  
# Read a single source history file.
def readSourceHistoryFile(inputFile):
    import numpy as np
  
    # Open the file.
    f = open(inputFile,'r')

    line = f.readline().split()

    if line[0].startswith('Time'):
        heights = np.zeros(1)
    elif line[0].startswith('Heights'):
        heights = [ float(val) for val in line[2:] ]
        f.readline()
    else:
        print('Error: Expected first line to start with "Time" or "Heights", but instead read',line[0])
        return


    data = []
    for line in f:
        line = [ float(val) for val in line.replace('(','').replace(')','').split() ]
        data.append(line)

    data = np.array(data)
    
    
    time = data[:,0]
    source = data[:,2:]

    return heights, time, source
  
  
  
  
  
  
  
  
  

# Write out a source file that will become SOWFA input.
def writeSourceForInput(fileName,height,time,source,sourceName,verbose):
    
    nTimes,nComponents = source.shape 
    nHeights = len(height)
    component = ['X','Y','Z']

    if (verbose):
        print('Writing '+fileName+':')
        print('   Source name = ',sourceName,' with ',nComponents,' components')
        if (nComponents > 1):
            for j in range(nComponents):
                print('   Source - ',component[j],' range = ',source[0,j],'..',source[-1,j])
        else:
            print('   Source range = ',source[0,0],'..',source[-1,0])
        print('   Number of heights = ',nHeights)
        print('   Height range = ',height[0],'..',height[-1],'m')
        print('   Number of times = ',nTimes)
        print('   Time range = ',time[0],'..',time[-1],'s')
        print('\n')
    
    # Open the file.
    fid = open(fileName,'w')
    
    
    # Write the momentum source height list.
    fid.write('sourceHeights'+sourceName+'\n')
    fid.write('(\n')
    for i in range(len(height)):
       fid.write('    ' + str(height[i]) + '\n')
       
    fid.write(');\n\n')
  
  
    # Write the source tables
    for j in range(nComponents):
        if (nComponents > 1):
            fid.write('sourceTable'+sourceName+component[j]+'\n')
        else:
            fid.write('sourceTable'+sourceName+'\n')
        fid.write('(\n')
        for n in range(len(time)):
            textStr = '    (' + str(time[n])
            for i in range(len(height)):
                textStr = textStr + ' ' + str(source[n][3*i+j])
            
            textStr = textStr + ')\n'
            fid.write(textStr)
                
        fid.write(');\n')
  
  
    # Close the file.
    fid.close()
    
