#ifndef IPileUpTool_h
#define IPileUpTool_h

#include "GaudiKernel/IAlgTool.h"

/** @class IPileUpTool IPileUpTool.h "Generation/IPileUpTool.h"
 *
 *  Abstract interface to pile up tools. Generates the number of pile-up
 *  interactions to generate for each event.
 *
 *  @author Patrick Robbe
 *  @date   2005-08-17
 */

namespace Jug::Base {

static const InterfaceID IID_IPileUpTool("IPileUpTool", 3, 0);

class IPileUpTool : virtual public IAlgTool {
public:
  /** Computes the number of pile-up interactions in the event.
   *  @param[out] currentLuminosity  Luminosity of the current event.
   *  @return Number of pile-up interactions to generate.
   */
  virtual unsigned int numberOfPileUp() = 0;

  virtual double getMeanPileUp() = 0;

  /// Print various counters at the end of the job
  virtual void printPileUpCounters() = 0;

};

}

#endif  // IPileUpTool_h
