#include "surge/sd3d/solver.hpp"

#include <stdexcept>

namespace surge::sd3d {
namespace {

DenseBlocks zeroBlocks(int n) {
  DenseBlocks out;
  out.n = n;
  out.R = Matrix::Zero(6 * n, 6 * n);
  out.Phi = Matrix::Zero(6 * n, 9);
  return out;
}

DenseBlocks assembleFar(const System& system, const SolverOptions& options) {
  switch (options.farField) {
    case FarFieldModel::None:
      return zeroBlocks(static_cast<int>(system.particles.size()));
    case FarFieldModel::OpenSpaceRPY:
      return farFieldResistanceFromMobilityTT(system, assembleOpenSpaceRPYMobilityTT(system));
    case FarFieldModel::PeriodicStokesletEwald:
      return farFieldResistanceFromMobilityTT(system, assemblePeriodicStokesletEwaldMobilityTT(system, options.ewald));
  }
  throw std::invalid_argument("unknown far-field model");
}

}  // namespace

DenseBlocks assembleHydrodynamicBlocks(const System& system, const SolverOptions& options) {
  const int n = static_cast<int>(system.particles.size());
  if (n == 0) return zeroBlocks(0);

  switch (options.assembly) {
    case HydrodynamicAssembly::PairwiseResistanceOnly:
      return assemblePairwiseResistance3D(system, options.near);

    case HydrodynamicAssembly::FarFieldOnly:
      return assembleFar(system, options);

    case HydrodynamicAssembly::HybridFarPlusPairwise: {
      DenseBlocks out = assembleFar(system, options);
      DenseBlocks near = assemblePairwiseResistance3D(system, options.near);

      // This is the practical modular form. A fully rigorous SD correction is:
      // R = R_far + R_2B_exact - R_2B_far.
      // The current implementation adds non-translational near-field and Phi
      // information immediately; the exact two-body far subtraction is the next
      // validation task once high-order coefficient tables are installed.
      for (int row = 0; row < 6 * n; ++row) {
        for (int col = 0; col < 6 * n; ++col) {
          const bool translational = (row % 6 < 3) && (col % 6 < 3);
          if (!translational) out.R(row, col) += near.R(row, col);
        }
      }
      out.Phi += near.Phi;
      return out;
    }
  }

  throw std::invalid_argument("unknown hydrodynamic assembly");
}

Vector hydrodynamicForce(const DenseBlocks& blocks, const Vector& Ustar, const Mat3& E) {
  const int n6 = 6 * blocks.n;
  if (Ustar.size() != n6) throw std::invalid_argument("Ustar size mismatch");
  return blocks.R * Ustar + blocks.Phi * flattenStrain(E);
}

SolverResult solveForceBalance(
    const System& system,
    const Mat3& E,
    const Vector& nonHydroForce,
    const SolverOptions& options) {
  const int n6 = 6 * static_cast<int>(system.particles.size());
  if (nonHydroForce.size() != n6) throw std::invalid_argument("nonHydroForce size must be 6N");

  SolverResult result;
  result.blocks = assembleHydrodynamicBlocks(system, options);
  const Vector rhs = -(result.blocks.Phi * flattenStrain(E) + nonHydroForce);

  Eigen::LDLT<Matrix> ldlt(result.blocks.R);
  if (ldlt.info() == Eigen::Success) {
    result.Ustar = ldlt.solve(rhs);
  } else {
    result.Ustar = result.blocks.R.colPivHouseholderQr().solve(rhs);
  }
  result.hydroForce = hydrodynamicForce(result.blocks, result.Ustar, E);
  return result;
}

}  // namespace surge::sd3d
