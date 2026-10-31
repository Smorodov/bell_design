#pragma once

#include "bell_types.h"
#include "modal_types.hpp"
#include "modal_influence.hpp"
#include <filesystem>
#include <string>

namespace BellOptimizer {

struct ModalProblem {
    const BellMeshData& mesh;
    const BellMaterial& material;
    std::filesystem::path work_directory = ".";
    std::string job_name = "bell";
    const ModalInfluenceProvider* influence_provider = nullptr;
};

class ModalSolver {
public:
    virtual ~ModalSolver() = default;
    virtual bool solve(const ModalProblem& problem, ModalData& result) = 0;
    virtual const char* name() const noexcept = 0;
};

} // namespace BellOptimizer
