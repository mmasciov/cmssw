import FWCore.ParameterSet.Config as cms

def customizeHLTPixelTrackingCAHardCurvCut(process):

    # if any of the following objects does not exist, do not apply any customisation
    for objLabel in [
        'hltPixelTracksSoA',
    ]:
        if not hasattr(process, objLabel):
            print(f'# WARNING: customizeHLTPixelTrackingCAHardCurvCut failed (object with label "{objLabel}" not found) - no customisation applied !')
            return process

    process.hltPixelTracksSoA.hardCurvCut = cms.double( 0.0328407225 )

    if hasattr(process, 'hltPixelTracksSoASerialSync'):
        process.hltPixelTracksSoASerialSync.hardCurvCut = cms.double( 0.0328407225 )

    return process
