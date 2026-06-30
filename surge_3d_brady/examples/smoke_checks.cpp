#include <cmath>
#include <iostream>
#include <stdexcept>

#include "surge/sd3d/projection.hpp"
#include "surge/sd3d/time_integrator.hpp"

namespace {

void require(bool ok, const char* message) {
  if (!ok) throw std::runtime_error(message);
}

bool finiteMatrix(const surge::sd3d::Matrix& m) {
  return m.array().isFinite().all();
}

bool finiteVector(const surge::sd3d::Vector& v) {
  return v.array().isFinite().all();
}

surge::sd3d::System twoSphereOpenSystem() {
  surge::sd3d::System system;
  system.viscosity = 1.0;
  system.particles = {
      {{-1.25, 0.25, 0.0}, 1.0},
      {{1.25, -0.25, 0.0}, 1.0},
  };
  return system;
}

}  // namespace

int main() {
  using namespace surge::sd3d;

  const System openSystem = twoSphereOpenSystem();
  SolverOptions pairwise;
  pairwise.assembly = HydrodynamicAssembly::PairwiseResistanceOnly;
  pairwise.farField = FarFieldModel::None;
  pairwise.near.maxPrintedM = 11;

  const DenseBlocks near = assembleHydrodynamicBlocks(openSystem, pairwise);
  require(near.n == 2, "near-field block particle count mismatch");
  require(near.R.rows() == 12 && near.R.cols() == 12, "near-field R size mismatch");
  require(near.Phi.rows() == 12 && near.Phi.cols() == 9, "near-field Phi size mismatch");
  require(finiteMatrix(near.R), "near-field R contains non-finite values");
  require(finiteMatrix(near.Phi), "near-field Phi contains non-finite values");
  require((near.R - near.R.transpose()).norm() < 1e-8, "near-field R is not symmetric");

  SolverOptions far;
  far.assembly = HydrodynamicAssembly::FarFieldOnly;
  far.farField = FarFieldModel::OpenSpaceRPY;
  const DenseBlocks rpy = assembleHydrodynamicBlocks(openSystem, far);
  require(finiteMatrix(rpy.R), "open-space far-field R contains non-finite values");
  require((rpy.R - rpy.R.transpose()).norm() < 1e-8, "open-space far-field R is not symmetric");

  System periodic = openSystem;
  periodic.box.length = {12.0, 12.0, 12.0};
  periodic.box.periodicX = true;
  periodic.box.periodicY = true;
  periodic.box.periodicZ = true;
  SolverOptions ewald;
  ewald.assembly = HydrodynamicAssembly::FarFieldOnly;
  ewald.farField = FarFieldModel::PeriodicStokesletEwald;
  ewald.ewald.xi = 0.35;
  ewald.ewald.realCutImages = 1;
  ewald.ewald.kMax = 3;
  const DenseBlocks periodicFar = assembleHydrodynamicBlocks(periodic, ewald);
  require(finiteMatrix(periodicFar.R), "periodic Ewald R contains non-finite values");

  const Matrix P = monolayerProjectionMatrix(2);
  const DenseBlocks projected = projectBlocksToMonolayer(near);
  require(P.rows() == 6 && P.cols() == 12, "monolayer projection size mismatch");
  require(projected.R.rows() == 6 && projected.R.cols() == 6, "projected R size mismatch");
  require(finiteMatrix(projected.R), "projected R contains non-finite values");

  Mat3 E = Mat3::Zero();
  E(0, 1) = 0.5;
  E(1, 0) = 0.5;
  Vector nonHydro = Vector::Zero(12);
  const SolverResult solved = solveForceBalance(openSystem, E, nonHydro, pairwise);
  require(finiteVector(solved.Ustar), "force-balance velocity contains non-finite values");
  require(finiteVector(solved.hydroForce), "hydrodynamic force contains non-finite values");

  TimeIntegratorOptions time;
  time.dt = 1e-3;
  time.steps = 4;
  const auto frames = integrateRK4(openSystem, E, {}, pairwise, time);
  require(frames.size() == 5, "RK4 frame count mismatch");
  for (const auto& frame : frames) {
    for (const Vec3& x : frame.positions) require(x.array().isFinite().all(), "trajectory contains non-finite position");
  }

  std::cout << "smoke checks passed\n";
  return 0;
}
