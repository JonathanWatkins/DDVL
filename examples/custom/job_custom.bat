; geometry 0=channel  1=tube  2=custom
;
[Header]
outputType=1
geometry=2	
simulationTime=10000
temp=0.01
lorentzForce=0.0

[Geometry]
periodicity=x

[InputData]
altPosFile=true
altPosFileName=posdata.txt

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

[Interactions]
vvForce=1
Phi=0.2067
lambda=1.11

[Job]
jobtag=custom
