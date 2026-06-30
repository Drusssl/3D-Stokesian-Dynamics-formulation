#include "surge/sd3d/projection.hpp"

#include <stdexcept>

namespace surge::sd3d {

Matrix monolayerProjectionMatrix(int n) {
  Matrix P = Matrix::Zero(3 * n, 6 * n);
  for (int i = 0; i < n; ++i) {
    P(3 * i + 0, dof(i, 0, false)) = 1.0;  // Ux
    P(3 * i + 1, dof(i, 1, false)) = 1.0;  // Uy
    P(3 * i + 2, dof(i, 2, true)) = 1.0;   // Omega_z
  }
  return P;
}

DenseBlocks projectBlocksToMonolayer(const DenseBlocks& blocks) {
  const Matrix P = monolayerProjectionMatrix(blocks.n);
  DenseBlocks out;
  out.n = blocks.n;
  out.R = P * blocks.R * P.transpose();
  out.Phi = P * blocks.Phi;
  return out;
}

Vector liftMonolayerVelocity(const Vector& qdot3n) {
  if (qdot3n.size() % 3 != 0) throw std::invalid_argument("monolayer velocity size must be 3N");
  const int n = static_cast<int>(qdot3n.size() / 3);
  Vector full = Vector::Zero(6 * n);
  for (int i = 0; i < n; ++i) {
    full[dof(i, 0, false)] = qdot3n[3 * i + 0];
    full[dof(i, 1, false)] = qdot3n[3 * i + 1];
    full[dof(i, 2, true)] = qdot3n[3 * i + 2];
  }
  return full;
}

}  // namespace surge::sd3d
