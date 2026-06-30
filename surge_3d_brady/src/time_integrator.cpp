#include "surge/sd3d/time_integrator.hpp"

#include <stdexcept>

namespace surge::sd3d {
namespace {

TrajectoryFrame frameFromSystem(const System& system, double time) {
  TrajectoryFrame f;
  f.time = time;
  f.positions.reserve(system.particles.size());
  for (const Particle& p : system.particles) f.positions.push_back(p.x);
  return f;
}

std::vector<Vec3> translationalVelocity(const Vector& Ustar, const Mat3& E, const System& system) {
  std::vector<Vec3> v(system.particles.size(), Vec3::Zero());
  for (int i = 0; i < static_cast<int>(system.particles.size()); ++i) {
    v[static_cast<std::size_t>(i)] =
        Ustar.segment<3>(6 * i) + E * system.particles[static_cast<std::size_t>(i)].x;
  }
  return v;
}

void applyPositionIncrement(System& system, const std::vector<Vec3>& dx, double scale) {
  for (int i = 0; i < static_cast<int>(system.particles.size()); ++i) {
    system.particles[static_cast<std::size_t>(i)].x += scale * dx[static_cast<std::size_t>(i)];
  }
}

std::vector<Vec3> velocityRhs(
    const System& system,
    double time,
    const Mat3& E,
    const ForceCallback& nonHydroForce,
    const SolverOptions& solverOptions) {
  const Vector fp = nonHydroForce ? nonHydroForce(system, time)
                                  : Vector::Zero(6 * static_cast<int>(system.particles.size()));
  SolverResult solved = solveForceBalance(system, E, fp, solverOptions);
  return translationalVelocity(solved.Ustar, E, system);
}

}  // namespace

std::vector<TrajectoryFrame> integrateEuler(
    System initial,
    const Mat3& E,
    ForceCallback nonHydroForce,
    const SolverOptions& solverOptions,
    const TimeIntegratorOptions& timeOptions) {
  if (timeOptions.dt <= 0.0 || timeOptions.steps < 0) throw std::invalid_argument("bad time options");
  std::vector<TrajectoryFrame> frames;
  frames.reserve(static_cast<std::size_t>(timeOptions.steps + 1));
  System system = initial;
  double time = 0.0;
  frames.push_back(frameFromSystem(system, time));

  for (int step = 0; step < timeOptions.steps; ++step) {
    const std::vector<Vec3> v = velocityRhs(system, time, E, nonHydroForce, solverOptions);
    applyPositionIncrement(system, v, timeOptions.dt);
    time += timeOptions.dt;
    frames.push_back(frameFromSystem(system, time));
  }
  return frames;
}

std::vector<TrajectoryFrame> integrateRK4(
    System initial,
    const Mat3& E,
    ForceCallback nonHydroForce,
    const SolverOptions& solverOptions,
    const TimeIntegratorOptions& timeOptions) {
  if (timeOptions.dt <= 0.0 || timeOptions.steps < 0) throw std::invalid_argument("bad time options");
  std::vector<TrajectoryFrame> frames;
  frames.reserve(static_cast<std::size_t>(timeOptions.steps + 1));
  System system = initial;
  double time = 0.0;
  frames.push_back(frameFromSystem(system, time));

  for (int step = 0; step < timeOptions.steps; ++step) {
    const double h = timeOptions.dt;
    const std::vector<Vec3> k1 = velocityRhs(system, time, E, nonHydroForce, solverOptions);

    System s2 = system;
    applyPositionIncrement(s2, k1, 0.5 * h);
    const std::vector<Vec3> k2 = velocityRhs(s2, time + 0.5 * h, E, nonHydroForce, solverOptions);

    System s3 = system;
    applyPositionIncrement(s3, k2, 0.5 * h);
    const std::vector<Vec3> k3 = velocityRhs(s3, time + 0.5 * h, E, nonHydroForce, solverOptions);

    System s4 = system;
    applyPositionIncrement(s4, k3, h);
    const std::vector<Vec3> k4 = velocityRhs(s4, time + h, E, nonHydroForce, solverOptions);

    for (int i = 0; i < static_cast<int>(system.particles.size()); ++i) {
      system.particles[static_cast<std::size_t>(i)].x +=
          (h / 6.0) * (k1[static_cast<std::size_t>(i)]
              + 2.0 * k2[static_cast<std::size_t>(i)]
              + 2.0 * k3[static_cast<std::size_t>(i)]
              + k4[static_cast<std::size_t>(i)]);
    }
    time += h;
    frames.push_back(frameFromSystem(system, time));
  }
  return frames;
}

}  // namespace surge::sd3d
