import sourceData as sd
import matplotlib.pyplot as plt

  
  
  
  
  
  
# Specify the directory where the source history files reside.
inputDir = './sourceHistory'
outputFileMomentum = './givenSourceU'
outputFileTemperature = './givenSourceT'
plotOn = True


# Assemble the source information.
[heightMomentum,heightTemperature,timeMomentum,timeTemperature,sourceMomentum,sourceTemperature] = sd.assembleSourceHistory(inputDir)


# Write the source file for input to the solver.
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
