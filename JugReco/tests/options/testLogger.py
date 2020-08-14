from Gaudi.Configuration import *
from Configurables import ApplicationMgr



from Configurables import TestACTSLogger
logtest = TestACTSLogger()


#  auto trackFindingCfg = Options::readTrackFindingConfig(vm);
#trackFindingCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
#trackFindingCfg.inputInitialTrackParameters =
#    particleSmearingCfg.outputTrackParameters;
#trackFindingCfg.outputTrajectories = "trajectories";
#trackFindingCfg.findTracks = TrackFindingAlgorithm::makeTrackFinderFunction(
#    trackingGeometry, magneticField, logLevel);
#sequencer.addAlgorithm(
#    std::make_shared<TrackFindingAlgorithm>(trackFindingCfg, logLevel));




ApplicationMgr(EvtSel='NONE',
               EvtMax=1,
               TopAlg=[logtest],
               ExtSvc=[],
               OutputLevel=DEBUG)
