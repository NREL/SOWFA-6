import sourceData as sd
import matplotlib.pyplot as plt

  
  
  
  
  
  
# Specify the directory where the source history files reside.
inputDir = './sourceHistory'
outputFileMomentum = './givenSourceU'
outputFileTemperature = './givenSourceT'
timeRange = [20000.0,110000.00]
plotOn = False


# Assemble the source information.
[heightMomentum,timeMomentum,sourceMomentum] = sd.assembleSourceHistory(inputDir,'SourceMomentumHistory')
[heightTemperature,timeTemperature,sourceTemperature] = sd.assembleSourceHistory(inputDir,'SourceTemperatureHistory')


# Write the source file for input to the solver.
#sd.writeSourceForInput(outputFileMomentum,heightMomentum,timeMomentum,sourceMomentum,'Momentum',True,timeRange)
#sd.writeSourceForInput(outputFileTemperature,heightTemperature,timeTemperature,sourceTemperature,'Temperature',True,timeRange)
sd.writeSourceForInput(outputFileMomentum,heightMomentum,timeMomentum,sourceMomentum,'Momentum',True)
sd.writeSourceForInput(outputFileTemperature,heightTemperature,timeTemperature,sourceTemperature,'Temperature',True)


# Plot if desired.
if (plotOn):
    plt.figure(1)
    plt.plot(timeMomentum,sourceMomentum)
    plt.legend(['x','y','z'])
    plt.xlabel('time (s)')
    plt.ylabel('momentum source (m/s^2)')
    plt.show()

    plt.figure(2)
    plt.plot(timeTemperature,sourceTemperature,'k-')
    plt.xlabel('time (s)')
    plt.ylabel('temperature source (K/s)')
    plt.show()
