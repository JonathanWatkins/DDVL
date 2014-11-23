; geometry 0=channel  1=tube  2=custom
;
[Header]
outputType=1
geometry=1	
simulationTime=1000
temp=0.0001
lorentzForce=0.0

[Geometry]
sourceBfield=0.25
sinkBfield=0.24
bathLength=10
bathWidth=8
channelLength=60
channelWidth=8

[InputData]
altPosFile=false
altPosFileName=posdata.txt
altPinsFile=false
altPinsFileName=gs.txt

[GeneralParameters]
a0 = 1
binSize = 5
cellSize=6.66
pi = 3.14159265358979
forceRange=6.66
eta=1.0
kB=1.0
Ap=1
dt=0.01
tau=1
triangulationInterval=5
framedataInterval=100
drawCoordinateGrid=false
showParticleTracker=false
thermostat=Anderson
disorderDensity=0
disorderStrength=1e-11
disorderRange=0.2e-7

[BathParameters]
applyBathVelocities=false
applyStiffBath=false

[Interactions]
vvForce=1
Phi=0.2067
lambda=1.11

[Job]
jobtag=tube
