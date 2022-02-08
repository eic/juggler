#include <algorithm>

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/UniqueID.h"

// Tensorflow headers
#include "tensorflow/lite/interpreter.h"
#include "tensorflow/lite/kernels/register.h"
#include "tensorflow/lite/model.h"
#include "tensorflow/lite/optional_debug_tools.h"

// Event Model related classes
#include "eicd/ClusterCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  /** Simple cluster identification algorithm using ML.
   *
   * \ingroup reco
   */
  class ClusterIdentification : public GaudiAlgorithm, AlgorithmIDMixin<> {
  public:
    DataHandle<eic::ClusterCollection> m_inputClusterCollection{"inputClusterCollection", Gaudi::DataHandle::Reader, this};

    Gaudi::Property<std::string> m_modelTFLiteFile{this, "modelTFLiteFile", ""};

    // interpreter
    std::unique_ptr<tflite::Interpreter> m_interpreter;

    ClusterIdentification(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc)
      , AlgorithmIDMixin<>(name, info()) {
      declareProperty("inputClusterCollection", m_inputClusterCollection, "");
      declareProperty("modelTFLiteFile", m_modelTFLiteFile, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }

      // load model from file
      std::unique_ptr<tflite::FlatBufferModel> model =
        tflite::FlatBufferModel::BuildFromFile(m_modelTFLiteFile.value().data());

      // build interpreter from model
      tflite::ops::builtin::BuiltinOpResolver resolver;
      tflite::InterpreterBuilder builder(*model, resolver);
      builder(&m_interpreter);

      // allocate tensors for interpreter
      m_interpreter->AllocateTensors();

      // debug
      printf("=== Pre-invoke Interpreter State ===\n");
      tflite::PrintInterpreterState(m_interpreter.get());

      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto& clusters = *m_inputClusterCollection.get();

      // fill input tensors
      double* input0 = m_interpreter->typed_input_tensor<double>(0);

      // run inference
      if (m_interpreter->Invoke() != kTfLiteOk) return StatusCode::FAILURE;

      // debug
      printf("\n\n=== Post-invoke Interpreter State ===\n");
      tflite::PrintInterpreterState(m_interpreter.get());

      // read output tensors
      double* output = m_interpreter->typed_output_tensor<double>(0);

      return StatusCode::SUCCESS;
    }
  };

  DECLARE_COMPONENT(ClusterIdentification)

} // namespace Jug::Reco
