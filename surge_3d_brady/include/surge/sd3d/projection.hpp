#pragma once

#include "core.hpp"

namespace surge::sd3d {

// Project 6N generalized vector onto Brady monolayer DOFs:
// per particle [Ux, Uy, Omega_z].
Matrix monolayerProjectionMatrix(int n);

DenseBlocks projectBlocksToMonolayer(const DenseBlocks& blocks);

Vector liftMonolayerVelocity(const Vector& qdot3n);

}  // namespace surge::sd3d
