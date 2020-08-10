from Gaudi.Configuration import *
from Configurables import ApplicationMgr

from Configurables import HelloWorld
alg = HelloWorld()

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=[ 'gem_tracker_disc.xml'], OutputLevel = DEBUG)

#from Configurables import TestCellCounting
#cells = TestCellCounting("cells", readoutName="ECalHits",
#                         fieldNames=["system"],
#                         fieldValues=[0],
#                         volumeMatchName="BoxECal",
#                         OutputLevel = DEBUG)
# ApplicationMgr
ApplicationMgr(EvtSel='NONE',
               EvtMax=1,
               TopAlg=[alg],
               ExtSvc=[geoservice],
               OutputLevel=DEBUG)
