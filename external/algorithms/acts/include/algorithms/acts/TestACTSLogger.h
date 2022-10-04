// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#ifndef TESTRECONSTRUCTION_TESTACTSLOGGER_H
#define TESTRECONSTRUCTION_TESTACTSLOGGER_H

// GAUDI
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Property.h"

// FCCSW
#include "JugBase/DataHandle.h"

/** Logging for ACTS.
 *
 * \ingroup tracking
 */
class TestACTSLogger : public GaudiAlgorithm {
public:
  explicit TestACTSLogger(const std::string&, ISvcLocator*);
  virtual ~TestACTSLogger();
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
