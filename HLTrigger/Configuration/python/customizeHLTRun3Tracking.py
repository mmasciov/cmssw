import copy
import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.CUDACore.SwitchProducerCUDA import SwitchProducerCUDA
from HLTrigger.Configuration.common import *
from HLTrigger.Configuration.customizeHLTforPatatrack import *
from Configuration.Eras.Modifier_run3_common_cff import run3_common

def customizeHLTRun3Tracking(process):

    process = customizeHLTforPatatrack(process)    
    if hasattr(process,'hltPixelTracksCUDA'):
        process.hltPixelTracksCUDA.includeJumpingForwardDoublets = cms.bool(True)
        process.hltPixelTracksCUDA.minHitsPerNtuplet             = cms.uint32(3)
        process.hltPixelTracksCUDA.idealConditions = cms.bool(True)
    if hasattr(process,'hltPixelTracksSoA'):
        process.hltPixelTracksSoA.cpu.includeJumpingForwardDoublets = cms.bool(True)
        process.hltPixelTracksSoA.cpu.minHitsPerNtuplet             = cms.uint32(3)
        process.hltPixelTracksSoA.cpu.idealConditions = cms.bool(True)

    if hasattr(process,'HLTIter0PSetTrajectoryFilterIT'):
        process.HLTIter0PSetTrajectoryFilterIT.minHitsMinPt        = cms.int32(3)
        process.HLTIter0PSetTrajectoryFilterIT.minimumNumberOfHits = cms.int32(3)

    if hasattr(process,'hltIter0PFLowPixelSeedsFromPixelTracks'):
        process.hltIter0PFLowPixelSeedsFromPixelTracks.includeFourthHit = cms.bool(True)

    if hasattr(process,'hltSiStripRawToClustersFacility'):
        process.hltSiStripRawToClustersFacility.onDemand = cms.bool( False )

    if hasattr(process,'hltIter0PFlowTrackCutClassifier'):
        delattr(process,'hltIter0PFlowTrackCutClassifier')
        process.hltIter0PFlowTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
            src = cms.InputTag("hltIter0PFlowCtfWithMaterialTracks"),
            beamspot = cms.InputTag("hltOnlineBeamSpot"),
            vertices = cms.InputTag("hltTrimmedPixelVertices"),
            qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
            mva = cms.PSet(
                minPixelHits = cms.vint32(0, 0, 0),
                maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
                dr_par = cms.PSet(
                    d0err = cms.vdouble(0.003, 0.003, 0.003),
                    dr_par2 = cms.vdouble(3.40282346639e+38, 0.6, 0.6),
                    dr_par1 = cms.vdouble(3.40282346639e+38, 0.8, 0.8),
                    dr_exp = cms.vint32(4, 4, 4),
                    d0err_par = cms.vdouble(0.001, 0.001, 0.001)
                ),
                maxLostLayers = cms.vint32(1, 1, 1),
                min3DLayers = cms.vint32(0, 0, 0),
                dz_par = cms.PSet(
                    dz_par1 = cms.vdouble(3.40282346639e+38, 0.75, 0.75),
                    dz_par2 = cms.vdouble(3.40282346639e+38, 0.5, 0.5),
                    dz_exp = cms.vint32(4, 4, 4)
                ),
                minNVtxTrk = cms.int32(3),
                maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
                minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
                maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
                maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
                maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
                minLayers = cms.vint32(3, 3, 3)
            ),
            ignoreVertices = cms.bool(False)
        )
    
    if hasattr(process,'hltMergedTracks'):
        delattr(process,'hltMergedTracks')
        process.hltMergedTracks = process.hltIter0PFlowTrackSelectionHighPurity.clone()

    if hasattr(process,'HLTIterativeTrackingIteration0'):
        delattr(process,'HLTIterativeTrackingIteration0')
        process.HLTIterativeTrackingIteration0 = cms.Sequence( process.hltIter0PFLowPixelSeedsFromPixelTracks + process.hltIter0PFlowCkfTrackCandidates + process.hltIter0PFlowCtfWithMaterialTracks + process.hltIter0PFlowTrackCutClassifier + process.hltMergedTracks )
    
    if hasattr(process,'HLTIterativeTrackingIter02'):
        delattr(process,'HLTIterativeTrackingIter02')
        process.HLTIterativeTrackingIter02 = cms.Sequence( process.HLTIterativeTrackingIteration0 )
    
    if hasattr(process,'MC_ReducedIterativeTracking_v12'):
        delattr(process,'MC_ReducedIterativeTracking_v12')
        process.MC_ReducedIterativeTracking_v12 = cms.Path( process.HLTBeginSequence + process.hltPreMCReducedIterativeTracking + process.HLTDoLocalPixelSequence + process.HLTRecopixelvertexingSequence + process.HLTDoLocalStripSequence + process.HLTIterativeTrackingIter02 + process.HLTEndSequence )

    return process



def customizeHLTRun3TrackingAllPixelVertices(process):

    process.hltPixelTracksClean = cms.EDProducer(
        "TrackWithVertexSelector",
        # the track collection
        src = cms.InputTag('hltPixelTracks'),
        # kinematic cuts  (pT in GeV)
        etaMin = cms.double(0.0),
        etaMax = cms.double(5.0),
        ptMin = cms.double(0.0),
        ptMax = cms.double(999999.),
        # impact parameter cut (in cm)
        d0Max = cms.double(999.),
        dzMax = cms.double(999.),
        # quality cuts (valid hits, normalized chi2)
        normalizedChi2 = cms.double(999999.),
        numberOfValidHits = cms.uint32(0),
        numberOfValidPixelHits = cms.uint32(0),
        numberOfValidHitsForGood = cms.uint32(4),
        numberOfValidPixelHitsForGood = cms.uint32(4),
        numberOfLostHits = cms.uint32(999), ## at most 999 lost hits
        ptErrorCut = cms.double(999999.), ## [pTError/pT]*max(1,normChi2) <= ptErrorCut
        quality = cms.string("any"), # quality cut as defined in reco::TrackBase
        # compatibility with a vertex ?
        useVtx = cms.bool(True),
        #vertexTag = cms.InputTag('hltTrimmedPixelVertices'),
        vertexTag = cms.InputTag('hltPixelVertices'),
        timesTag = cms.InputTag(''),
        timeResosTag = cms.InputTag(''),
        nVertices = cms.uint32(100), ## how many vertices to look at before dropping the track
        vtxFallback = cms.bool(True), ## falback to beam spot if there are no vertices
        # uses vtx=(0,0,0) with deltaZeta=15.9, deltaRho = 0.2
        zetaVtx = cms.double(0.3),
        rhoVtx = cms.double(0.1), ## tags used by b-tagging folks
        zetaVtxScale = cms.double(5.0),
        rhoVtxScale = cms.double(2.5), 
        zetaVtxSig = cms.double(2.0),
        rhoVtxSig = cms.double(2.0),
        nSigmaDtVertex = cms.double(0),
        # should _not_ be used for the TrackWithVertexRefSelector
        copyExtras = cms.untracked.bool(True), ## copies also extras and rechits on RECO
        copyTrajectories = cms.untracked.bool(False), # don't set this to true on AOD!
    )
    
    if hasattr(process,'hltIter0PFLowPixelSeedsFromPixelTracks'):
        process.hltIter0PFLowPixelSeedsFromPixelTracks.InputVertexCollection = cms.InputTag( "hltPixelVertices" )
        process.hltIter0PFLowPixelSeedsFromPixelTracks.InputCollection = cms.InputTag("hltPixelTracksClean"),

    if hasattr(process,'hltIter0PFlowTrackCutClassifier'):
        process.hltIter0PFlowTrackCutClassifier.vertices = cms.InputTag("hltPixelVertices")

    if hasattr(process,'HLTIterativeTrackingIteration0'):
        delattr(process,'HLTIterativeTrackingIteration0')
        process.HLTIterativeTrackingIteration0 = cms.Sequence( process.hltPixelTracksClean + process.hltIter0PFLowPixelSeedsFromPixelTracks + process.hltIter0PFlowCkfTrackCandidates + process.hltIter0PFlowCtfWithMaterialTracks + process.hltIter0PFlowTrackCutClassifier + process.hltMergedTracks )
            
    return process
