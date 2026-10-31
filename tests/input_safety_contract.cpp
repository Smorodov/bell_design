#include "core/simple_json.hpp"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"
#include "ccx/ccx_solver.hpp"

#include <limits>
#include <string>

int main() {
    using namespace BellOptimizer;

    // A malformed/non-contiguous node map must fail safely before matrix assembly.
    BellMeshData mesh;
    mesh.nodes = {{1, 10.0, 0.0}, {2, 10.0, 1.0}, {4, 11.0, 0.0}};
    mesh.elements = {{1, 2, 4, 1}};
    BellMaterial material;
    ModalData result;
    InternalModalSolver internal(true, 3);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    ModalProblem problem{mesh, material, ".", "invalid_mesh", &modal_influence_1};
    if (internal.solve(problem, result)) return 1;

    // Out-of-range numeric-to-int conversions must return the fallback.
    const auto large = bell::json::Parser("2147483648").parse();
    if (large.integer(17) != 17) return 2;
    const auto fractional = bell::json::Parser("12.75").parse();
    if (fractional.integer(17) != 12) return 3;

    // Numeric overflow and excessive nesting are rejected, not accepted as infinities
    // or allowed to recurse indefinitely.
    bool rejected_overflow = false;
    try { (void)bell::json::Parser("1e9999").parse(); }
    catch (const std::runtime_error&) { rejected_overflow = true; }
    if (!rejected_overflow) return 4;

    std::string deeply_nested(300, '[');
    deeply_nested += '0';
    deeply_nested.append(300, ']');
    bool rejected_depth = false;
    try { (void)bell::json::Parser(deeply_nested).parse(); }
    catch (const std::runtime_error&) { rejected_depth = true; }
    if (!rejected_depth) return 5;

    // Invalid job names must be rejected before writing files or launching a shell.
    BellMeshData empty_mesh;
    CalculiXSolver ccx("definitely-not-launched", true);
    ModalProblem bad_job{empty_mesh, material, ".", "bad;touch_injected"};
    if (ccx.solve(bad_job, result)) return 6;

    return 0;
}
