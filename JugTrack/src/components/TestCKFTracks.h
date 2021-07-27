#ifndef TESTRECONSTRUCTION_TESTACTSLOGGER_H
#define TESTRECONSTRUCTION_TESTACTSLOGGER_H

// GAUDI
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Property.h"

// FCCSW
#include "JugBase/DataHandle.h"


/** CKF track finding/fitting.
 *
 * \ingroup track
 */
class TestCKFTracks : public GaudiAlgorithm {
public:

  explicit TestCKFTracks(const std::string&, ISvcLocator*);
  virtual ~TestCKFTracks();

  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;

  /**  Execute.
   *   @return status code
   */
  virtual StatusCode execute() final;

  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

};
#endif /* TESTRECONSTRUCTION_TESTACTSLOGGER_H */
