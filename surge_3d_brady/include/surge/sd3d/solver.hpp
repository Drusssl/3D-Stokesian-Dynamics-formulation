#pragma once

#include "core.hpp"
#include "ewald.hpp"
#include "near_field.hpp"

namespace surge::sd3d {

enum class FarFieldModel {
  None,
  OpenSpaceRPY,
  PeriodicStokesletEwald
};

enum class HydrodynamicAssembly {
  PairwiseResistanceOnly,
  FarFieldOnly,
  HybridFarPlusPairwise
};

struct SolverOptions {
  HydrodynamicAssembly assembly = HydrodynamicAssembly::PairwiseResistanceOnly;
  FarFieldModel farField = FarFieldModel::OpenSpaceRPY;
  NearFieldOptions near{};
  EwaldOptions ewald{};
};

DenseBlocks assembleHydrodynamicBlocks(const System& system, const SolverOptions& options);

Vector hydrodynamicForce(const DenseBlocks& blocks, const Vector& Ustar, const Mat3& E);

SolverResult solveForceBalance(
    const System& system,
    const Mat3& E,
    const Vector& nonHydroForce,
    const SolverOptions& options);

}  // namespace surge::sd3d
