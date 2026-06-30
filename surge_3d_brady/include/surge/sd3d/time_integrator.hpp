#pragma once

#include <functional>

#include "solver.hpp"

namespace surge::sd3d {

using ForceCallback = std::function<Vector(const System& system, double time)>;

struct TimeIntegratorOptions {
  double dt = 1e-3;
  int steps = 1;
};

struct TrajectoryFrame {
  double time = 0.0;
  std::vector<Vec3> positions;
};

std::vector<TrajectoryFrame> integrateEuler(
    System initial,
    const Mat3& E,
    ForceCallback nonHydroForce,
    const SolverOptions& solverOptions,
    const TimeIntegratorOptions& timeOptions);

std::vector<TrajectoryFrame> integrateRK4(
    System initial,
    const Mat3& E,
    ForceCallback nonHydroForce,
    const SolverOptions& solverOptions,
    const TimeIntegratorOptions& timeOptions);

}  // namespace surge::sd3d
