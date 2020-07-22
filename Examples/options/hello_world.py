from Gaudi.Configuration import *
from Configurables import HelloWorld

alg = HelloWorld()

ApplicationMgr(
    EvtMax = 10,
    EvtSel = 'NONE',
    HistogramPersistency = 'NONE',
    TopAlg = [alg],
)
