#pragma once

#include "core.hpp"
#include "coefficients.hpp"

namespace surge::sd3d {

struct NearFieldOptions {
  int maxPrintedM = 11;
  const CoefficientTable* highOrderTable = nullptr;
};

DenseBlocks assemblePairwiseResistance3D(const System& system, const NearFieldOptions& options = {});

}  // namespace surge::sd3d
