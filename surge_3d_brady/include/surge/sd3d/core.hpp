#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace surge::sd3d {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

struct Particle {
  Vec3 x = Vec3::Zero();
  double radius = 1.0;
};

struct PeriodicBox {
  Vec3 length = Vec3::Zero();
  bool periodicX = false;
  bool periodicY = false;
  bool periodicZ = false;
};

struct System {
  double viscosity = 1.0;
  std::vector<Particle> particles;
  PeriodicBox box{};
};

struct DenseBlocks {
  int n = 0;
  Matrix R;     // 6N x 6N
  Matrix Phi;   // 6N x 9, row contracts with row-major vec(E)
};

struct SolverResult {
  Vector Ustar;       // 6N
  Vector hydroForce;  // 6N
  DenseBlocks blocks;
};

inline int dof(int particle, int component, bool angular) {
  return 6 * particle + (angular ? 3 : 0) + component;
}

inline int strainIndex(int i, int j) {
  return 3 * i + j;
}

inline Eigen::Vector<double, 9> flattenStrain(const Mat3& E) {
  Eigen::Vector<double, 9> out;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) out[strainIndex(i, j)] = E(i, j);
  }
  return out;
}

inline Vec3 minimumImage(Vec3 dr, const PeriodicBox& box) {
  auto wrap = [](double x, double L, bool periodic) {
    if (!periodic || L <= 0.0) return x;
    return x - L * std::round(x / L);
  };
  return {
      wrap(dr.x(), box.length.x(), box.periodicX),
      wrap(dr.y(), box.length.y(), box.periodicY),
      wrap(dr.z(), box.length.z(), box.periodicZ),
  };
}

inline double boxVolume(const PeriodicBox& box) {
  return box.length.x() * box.length.y() * box.length.z();
}

}  // namespace surge::sd3d
