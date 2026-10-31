#pragma once

#include "modal_solver.hpp"

namespace BellOptimizer {

// Dependency-free axisymmetric Q4 modal solver operating directly on the in-memory mesh.
class InternalModalSolver final : public ModalSolver {
public:
    explicit InternalModalSolver(bool quiet = true, int requested_modes = 10, bool capture_mode_shapes = false,
                                 int sparse_solver_threshold = 128, bool use_rcm_ordering = true)
        : quiet_(quiet), requested_modes_(requested_modes), capture_mode_shapes_(capture_mode_shapes),
          sparse_solver_threshold_(sparse_solver_threshold < 0 ? 0 : sparse_solver_threshold),
          use_rcm_ordering_(use_rcm_ordering) {}

    bool solve(const ModalProblem& problem, ModalData& result) override;

    const char* name() const noexcept override { return "InternalAxisymmetricFEM"; }

private:
    bool quiet_ = true;
    int requested_modes_ = 10;
    bool capture_mode_shapes_ = false;
    int sparse_solver_threshold_ = 128;
    bool use_rcm_ordering_ = true;
};

} // namespace BellOptimizer
