// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#ifndef TESTRECONSTRUCTION_TESTACTSLOGGER_H
#define TESTRECONSTRUCTION_TESTACTSLOGGER_H

// GAUDI
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

// FCCSW
#include <k4FWCore/DataHandle.h>

/** Logging for ACTS.
 *
 * \ingroup tracking
 */
class TestACTSLogger : public Gaudi::Algorithm {
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
  virtual StatusCode execute(const EventContext&) const final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;
};
#endif /* TESTRECONSTRUCTION_TESTACTSLOGGER_H */
