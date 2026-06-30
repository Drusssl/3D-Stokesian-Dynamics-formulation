#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "surge/sd3d/time_integrator.hpp"

int main(int argc, char** argv) {
  std::string outPath = "two_sphere_trajectory.csv";
  if (argc > 1) outPath = argv[1];

  surge::sd3d::System system;
  system.viscosity = 1.0;
  system.particles = {
      {{-1.25, 0.25, 0.0}, 1.0},
      {{ 1.25, -0.25, 0.0}, 1.0},
  };

  // Simple shear strain-rate tensor. The vorticity/background rotation is not
  // included here yet; this first driver validates force/strain coupling.
  surge::sd3d::Mat3 E = surge::sd3d::Mat3::Zero();
  E(0, 1) = 0.5;
  E(1, 0) = 0.5;

  surge::sd3d::SolverOptions solver;
  solver.assembly = surge::sd3d::HydrodynamicAssembly::PairwiseResistanceOnly;
  solver.farField = surge::sd3d::FarFieldModel::None;
  solver.near.maxPrintedM = 11;

  surge::sd3d::TimeIntegratorOptions time;
  time.dt = 1e-3;
  time.steps = 200;

  auto zeroForce = [](const surge::sd3d::System& s, double) {
    return surge::sd3d::Vector::Zero(6 * static_cast<int>(s.particles.size()));
  };

  const auto frames = surge::sd3d::integrateRK4(system, E, zeroForce, solver, time);

  std::ofstream out(outPath);
  out << std::setprecision(16);
  out << "time,x1,y1,z1,x2,y2,z2\n";
  for (const auto& frame : frames) {
    out << frame.time;
    for (const auto& x : frame.positions) out << ',' << x.x() << ',' << x.y() << ',' << x.z();
    out << '\n';
  }

  std::cout << "wrote " << outPath << "\n";
}
