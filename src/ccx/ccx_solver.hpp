#pragma once

#include "modal_solver.hpp"
#include "ccx_parser.hpp"
#include <string>

namespace BellOptimizer {

class CalculiXSolver final : public ModalSolver {
public:
    explicit CalculiXSolver(std::string executable =
#ifdef BELL_CCX_EXECUTABLE
        BELL_CCX_EXECUTABLE
#else
        "ccx"
#endif
        , bool quiet = true);

    bool solve(const ModalProblem& problem, ModalData& result) override;
    const char* name() const noexcept override { return "CalculiX"; }

private:
    std::string executable_;
    bool quiet_ = true;
};

} // namespace BellOptimizer
