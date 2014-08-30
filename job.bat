[Header]
runtype=2
outputType=1
geometry=0
sourceBfield=0.25
sinkBfield=0.14
bathLength=20
bathWidth=8
channelLength=60
channelWidth=8
simulationTime=40000
relaxtime=0
temp=0.000001
lorentzForce=0.0

[InputData]
altPosFile=false
altPosFileName=gs.txt
altPinsFile=false
altPinsFileName=gs.txt

[GeneralParameters]
a0 = 1
binSize = 5
cellSize=6.66
vfieldBinSize = 1
pi = 3.14159265358979
forceRange=6.66
eta=1.0
kB=1.0
Ap=1
dt=0.01
tau=1
drawInterval=5
triangulationInterval=5
framedataInterval=100
drawCoordinateGrid=false
calcTrajectories=true
showParticleTracker=false
thermostat=Anderson
disorderDensity=0
disorderStrength=1e-11
disorderRange=0.2e-7

[Annealing]
annealing=false
annealingtime=2000
annealingfactor=0.5
annealingendtime=5000
annealingT=0.03

[BathParameters]
applyBathVelocities=false
applyStiffBath=false
flatChannelEnds=false
reflectedChannelEnds=false

[WallParameters]
applyBounceBack=false

[Interactions]
vvForce=1
Phi=0.2067
lambda=1.11

[Bubble]
makeBubble=false
centerx=45
centery=7
radius=4

[BurgersVector]
calculateBurgersVector=false
centerx=45
centery=7
radius=7
tracked=true

[Job]
jobtag=jobname
