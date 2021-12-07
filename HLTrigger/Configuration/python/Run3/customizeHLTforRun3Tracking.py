import copy
import FWCore.ParameterSet.Config as cms
from HLTrigger.Configuration.common import *
from HLTrigger.Configuration.customizeHLTforPatatrack import *
from Configuration.ProcessModifiers.pixelNtupletFit_cff import pixelNtupletFit
#from Configuration.ProcessModifiers.gpu_cff import gpu

def customizeHLTforRun3Tracking(process):
    
    process.extend(pixelNtupletFit)
    #process.extend(gpu)

    process = customizeHLTforPatatrackTriplets(process)    
    if hasattr(process,'hltPixelTracksCUDA'):
        process.hltPixelTracksCUDA.includeJumpingForwardDoublets = cms.bool(True)
        process.hltPixelTracksCUDA.idealConditions               = cms.bool(False)
        process.hltPixelTracksCUDA.fillStatistics                = cms.bool(True)
        process.hltPixelTracksCUDA.useSimpleTripletCleaner       = cms.bool(False)
    if hasattr(process,'hltPixelTracksSoA'):
        process.hltPixelTracksSoA.cpu.includeJumpingForwardDoublets = cms.bool(True)
        process.hltPixelTracksSoA.cpu.idealConditions               = cms.bool(False)
        process.hltPixelTracksSoA.cpu.fillStatistics                = cms.bool(True)
        process.hltPixelTracksSoA.cpu.useSimpleTripletCleaner       = cms.bool(False)

    if hasattr(process,'hltPixelTracks'):
        process.hltPixelTracks.minNumberOfHits = cms.int32(0)
        process.hltPixelTracks.minQuality = cms.string('loose')

    if hasattr(process,'HLTIter0PSetTrajectoryFilterIT'):
        process.HLTIter0PSetTrajectoryFilterIT.minHitsMinPt        = cms.int32(3)
        process.HLTIter0PSetTrajectoryFilterIT.minimumNumberOfHits = cms.int32(3)

    if hasattr(process,'hltSiStripRawToClustersFacility'):
        process.hltSiStripRawToClustersFacility.onDemand = cms.bool( True )

    if hasattr(process,'hltIter0PFLowPixelSeedsFromPixelTracks'):
        process.hltIter0PFLowPixelSeedsFromPixelTracks.includeFourthHit = cms.bool(True)

    if hasattr(process,'hltIter0PFlowTrackCutClassifier'):
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
        process.hltMergedTracks = process.hltIter0PFlowTrackSelectionHighPurity.clone()

    process.HLTIterativeTrackingIteration0Task = cms.Sequence(
        process.hltIter0PFLowPixelSeedsFromPixelTracks +
        process.hltIter0PFlowCkfTrackCandidates +
        process.hltIter0PFlowCtfWithMaterialTracks +
        process.hltIter0PFlowTrackCutClassifier +
        process.hltMergedTracks
    )
    if hasattr(process,'HLTIterativeTrackingIteration0'):
        process.HLTIterativeTrackingIteration0 = cms.Sequence( process.HLTIterativeTrackingIteration0Task )
    
    if hasattr(process,'HLTIterativeTrackingIter02'):
        process.HLTIterativeTrackingIter02 = cms.Sequence( process.HLTIterativeTrackingIteration0 )
    
    if hasattr(process,'MC_ReducedIterativeTracking_v12'):
        process.MC_ReducedIterativeTracking_v12 = cms.Path( 
            process.HLTBeginSequence +
            process.hltPreMCReducedIterativeTracking +
            process.HLTDoLocalPixelSequence +
            process.HLTRecopixelvertexingSequence +
            process.HLTDoLocalStripSequence +
            process.HLTIterativeTrackingIter02 +
            process.HLTEndSequence
        )

    ### Cleanup
    if hasattr(process, 'HLTIterativeTrackingIteration1'):
        delattr(process, 'HLTIterativeTrackingIteration1')
    if hasattr(process, 'hltIter1ClustersRefRemoval'):
        delattr(process, 'hltIter1ClustersRefRemoval')
    if hasattr(process, 'hltIter1MaskedMeasurementTrackerEvent'):
        delattr(process, 'hltIter1MaskedMeasurementTrackerEvent')
    if hasattr(process, 'hltIter1PixelLayerQuadruplets'):
        delattr(process, 'hltIter1PixelLayerQuadruplets')
    if hasattr(process, 'hltIter1PFlowPixelTrackingRegions'):
        delattr(process, 'hltIter1PFlowPixelTrackingRegions')
    if hasattr(process, 'hltIter1PFlowPixelClusterCheck'):
        delattr(process, 'hltIter1PFlowPixelClusterCheck')
    if hasattr(process, 'hltIter1PFlowPixelHitDoublets'):
        delattr(process, 'hltIter1PFlowPixelHitDoublets')
    if hasattr(process, 'hltIter1PFlowPixelHitQuadruplets'):
        delattr(process, 'hltIter1PFlowPixelHitQuadruplets')
    if hasattr(process, 'hltIter1PixelTracks'):
        delattr(process, 'hltIter1PixelTracks')
    if hasattr(process, 'hltIter1PFLowPixelSeedsFromPixelTracks'):
        delattr(process, 'hltIter1PFLowPixelSeedsFromPixelTracks')
    if hasattr(process, 'hltIter1PFlowCkfTrackCandidates'):
        delattr(process, 'hltIter1PFlowCkfTrackCandidates')
    if hasattr(process, 'hltIter1PFlowCtfWithMaterialTracks'):
        delattr(process, 'hltIter1PFlowCtfWithMaterialTracks')
    if hasattr(process, 'hltIter1PFlowTrackCutClassifierPrompt'):
        delattr(process, 'hltIter1PFlowTrackCutClassifierPrompt')
    if hasattr(process, 'hltIter1PFlowTrackCutClassifierDetached'):
        delattr(process, 'hltIter1PFlowTrackCutClassifierDetached')
    if hasattr(process, 'hltIter1PFlowTrackCutClassifierMerged'):
        delattr(process, 'hltIter1PFlowTrackCutClassifierMerged')
    if hasattr(process, 'hltIter1PFlowTrackSelectionHighPurity'):
        delattr(process, 'hltIter1PFlowTrackSelectionHighPurity')

    if hasattr(process, 'HLTIter1TrackAndTauJets4Iter2Sequence'):
        delattr(process, 'HLTIter1TrackAndTauJets4Iter2Sequence')
    if hasattr(process, 'hltIter1TrackRefsForJets4Iter2'):
        delattr(process, 'hltIter1TrackRefsForJets4Iter2')
    if hasattr(process, 'hltAK4Iter1TrackJets4Iter2'):
        delattr(process, 'hltAK4Iter1TrackJets4Iter2')
    if hasattr(process, 'hltIter1TrackAndTauJets4Iter2'):
        delattr(process, 'hltIter1TrackAndTauJets4Iter2')

    if hasattr(process, 'HLTIterativeTrackingIteration2'):
        delattr(process, 'HLTIterativeTrackingIteration2')
    if hasattr(process, 'hltIter2ClustersRefRemoval'):
        delattr(process, 'hltIter2ClustersRefRemoval')
    if hasattr(process, 'hltIter2MaskedMeasurementTrackerEvent'):
        delattr(process, 'hltIter2MaskedMeasurementTrackerEvent')
    if hasattr(process, 'hltIter2PixelLayerTriplets'):
        delattr(process, 'hltIter2PixelLayerTriplets')
    if hasattr(process, 'hltIter2PFlowPixelTrackingRegions'):
        delattr(process, 'hltIter2PFlowPixelTrackingRegions')
    if hasattr(process, 'hltIter2PFlowPixelClusterCheck'):
        delattr(process, 'hltIter2PFlowPixelClusterCheck')
    if hasattr(process, 'hltIter2PFlowPixelHitDoublets'):
        delattr(process, 'hltIter2PFlowPixelHitDoublets')
    if hasattr(process, 'hltIter2PFlowPixelHitTriplets'):
        delattr(process, 'hltIter2PFlowPixelHitTriplets')
    if hasattr(process, 'hltIter2PFlowPixelSeeds'):
        delattr(process, 'hltIter2PFlowPixelSeeds')
    if hasattr(process, 'hltIter2PFlowCkfTrackCandidates'):
        delattr(process, 'hltIter2PFlowCkfTrackCandidates')
    if hasattr(process, 'hltIter2PFlowCtfWithMaterialTracks'):
        delattr(process, 'hltIter2PFlowCtfWithMaterialTracks')
    if hasattr(process, 'hltIter2PFlowTrackCutClassifier'):
        delattr(process, 'hltIter2PFlowTrackCutClassifier')
    if hasattr(process, 'hltIter2PFlowTrackSelectionHighPurity'):
        delattr(process, 'hltIter2PFlowTrackSelectionHighPurity')

    if hasattr(process, 'hltDoubletRecoveryClustersRefRemoval'):
        process.hltDoubletRecoveryClustersRefRemoval.trajectories = cms.InputTag( "hltMergedTracks" )
        process.hltDoubletRecoveryClustersRefRemoval.oldClusterRemovalInfo = cms.InputTag( "" )

    return process



def enableDoubletRecovery(process):

    if hasattr(process, 'hltMergedTracks'):
        process.hltIter0PFlowTrackSelectionHighPurity = process.hltMergedTracks.clone()

    if hasattr(process, 'hltDoubletRecoveryClustersRefRemoval'):
        process.hltDoubletRecoveryClustersRefRemoval.trajectories = cms.InputTag( "hltIter0PFlowTrackSelectionHighPurity" )
        process.hltDoubletRecoveryClustersRefRemoval.oldClusterRemovalInfo = cms.InputTag( "" )

    process.hltMergedTracks = cms.EDProducer(
        "TrackListMerger",
        ShareFrac = cms.double( 0.19 ),
        FoundHitBonus = cms.double( 5.0 ),
        LostHitPenalty = cms.double( 20.0 ),
        MinPT = cms.double( 0.05 ),
        Epsilon = cms.double( -0.001 ),
        MaxNormalizedChisq = cms.double( 1000.0 ),
        MinFound = cms.int32( 3 ),
        TrackProducers = cms.VInputTag( 'hltIter0PFlowTrackSelectionHighPurity','hltDoubletRecoveryPFlowTrackSelectionHighPurity' ),
        hasSelector = cms.vint32( 0, 0 ),
        indivShareFrac = cms.vdouble( 1.0, 1.0 ),
        selectedTrackQuals = cms.VInputTag( 'hltIter0PFlowTrackSelectionHighPurity','hltDoubletRecoveryPFlowTrackSelectionHighPurity' ),
        setsToMerge = cms.VPSet( 
            cms.PSet(  pQual = cms.bool( False ),
                       tLists = cms.vint32( 0, 1 )
                   )
        ),
        trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
        allowFirstHitShare = cms.bool( True ),
        newQuality = cms.string( "confirmed" ),
        copyExtras = cms.untracked.bool( True ),
        writeOnlyTrkQuals = cms.bool( False ),
        copyMVA = cms.bool( False )
    )

    process.HLTIterativeTrackingIteration0Task = cms.Sequence(
        process.hltIter0PFLowPixelSeedsFromPixelTracks +
        process.hltIter0PFlowCkfTrackCandidates +
        process.hltIter0PFlowCtfWithMaterialTracks +
        process.hltIter0PFlowTrackCutClassifier +
        process.hltIter0PFlowTrackSelectionHighPurity
    )
    process.HLTIterativeTrackingIteration0 = cms.Sequence( process.HLTIterativeTrackingIteration0Task )
    
    process.HLTIterativeTrackingIter02 = cms.Sequence( process.HLTIterativeTrackingIteration0 + process.HLTIterativeTrackingDoubletRecovery + process.hltMergedTracks )

    return process
