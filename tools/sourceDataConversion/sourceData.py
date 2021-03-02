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
def assembleSourceHistory(inputDir,sourceName):
    # Import necessary modules
    import numpy as np
  
  
    # Get the number of time directories and their names.
    [nTimes, outputTimes] = getOutputTimes(inputDir)
    
    
    # Initialize the big arrays.
    time = []
    source = []
  
  
    # Loop through the time directories and get the source information.
    for n in range(nTimes):
        print('\n' + 'Reading ' + sourceName)
        print('   -reading time: ' + str(outputTimes[n]) + ' s  [' + str(n+1) + '/' + str(nTimes) + ']...')
          
        inputFile = inputDir + '/' + outputTimes[n] + '/' + sourceName
  
        height, timeI, sourceI = readSourceHistoryFile(inputFile)
           
        if (n == 0):
            time = timeI
            source = sourceI
        else:
            startTime = timeI[0]
  
            l = len(time)
  
            if (time[l-1] >= startTime):
                indEnd = (np.where(time >= startTime))[0][0] - 1
  
            else:
                indEnd = l
                 
            time = np.append(time[0:indEnd],timeI)
            source = np.append(source[0:indEnd][:],sourceI,axis=0)
               
            timeI = []
            sourceI = []
               
    
    return height,time,source
  
  
  
  
  
  

  
  
  
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
def writeSourceForInput(fileName,height,time,source,sourceName,verbose=True,timeRange=[]):
    import numpy as np
    
    nHeights = len(height)
    nTimes,nComponents = source.shape 
    nComponents = int(nComponents/nHeights)
    component = ['X','Y','Z']

    indMin = 0
    indMax = nTimes-1
    if (bool(timeRange)):
        indMin = np.argmax(time >= timeRange[0])
        indMax = np.argmax(time >= timeRange[1])
        if (indMax < indMin):
            indMax = nTimes-1
        nTimes = indMax - indMin + 1

    if (verbose):
        print('\n' + 'Writing '+fileName+':')
        if (nComponents > 1):
            print('   Source name = ',sourceName,' with ',nComponents,' components')
        else:
            print('   Source name = ',sourceName,' with ',nComponents,' component')

        if (nComponents > 1):
            for j in range(nComponents):
                print('   Source - ' + str(component[j]) + ' range = ' + str(source[indMin,j]) + '..' + str(source[indMax,j]))
        else:
            print('   Source range = ' + str(source[indMin,0]) + '..' + str(source[indMax,0]))

        print('   Number of heights = ',nHeights)
        print('   Height range = ' + str(height[0]) + '..' + str(height[-1]) + ' m')
        print('   Number of times = ',nTimes)
        print('   Time range = ' + str(time[indMin]) + '..' + str(time[indMax]) + ' s')
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
        for n in range(indMin,indMax+1):
            textStr = '    (' + str(time[n])
            for i in range(nHeights):
                textStr = textStr + ' ' + str(source[n][nComponents*i+j])
            
            textStr = textStr + ')\n'
            fid.write(textStr)
                
        fid.write(');\n')
  
  
    # Close the file.
    fid.close()
    
