#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/RecoAlgos/interface/TrackFullCloneSelectorBase.h"
#include "CommonTools/RecoAlgos/interface/TrackWithVertexSelector.h"
#include "RecoPixelVertexing/PixelVertexFinding/interface/PVClusterComparer.h"

namespace reco {
  namespace modules {

    typedef TrackFullCloneSelectorBase< ::TrackWithVertexSelector> TrackWithVertexSelector;

    DEFINE_FWK_MODULE(TrackWithVertexSelector);

  }  // namespace modules
}  // namespace reco
