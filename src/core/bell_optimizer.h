// bell_optimizer.h
#ifndef BELL_OPTIMIZER_H
#define BELL_OPTIMIZER_H

#include "bell_types.h"
#include "bell_generator.h"
#include "modal_solver.hpp"
#include "modal_tracking.hpp"
#include "models/bell_model_profile.hpp"
#include <memory>
#include <vector>
#include <random>
#include <string>

using namespace BellOptimizer;

class EvolutionaryOptimizer {
public:
    EvolutionaryOptimizer(const BellMaterial& material,
                          const TargetData& targets,
                          const EvolutionParams& evo_params,
                          std::unique_ptr<ModalSolver> solver,
                          std::unique_ptr<Models::GeometrySeedProvider> seed_provider = nullptr);
    
    void run();
    
    const BellGeometry& getBestGeometry() const { return best_geometry_; }
    double getBestError() const { return best_error_; }
    
    void setInitialGuess(const BellGeometry& initial);
    bool loadBestFromFile(const std::string& filename);

    // ============================================================
    // СТРУКТУРА КАНДИДАТА (ПУБЛИЧНАЯ)
    // ============================================================
    struct Candidate {
        BellGeometry geometry;
        double error = BellConstants::ERROR_INVALID;
        ModalData modal_data;
        bool valid = false;
    };

private:
    BellMaterial material_;
    TargetData targets_;
    EvolutionParams evo_params_;
    std::unique_ptr<Models::GeometrySeedProvider> seed_provider_;
    BellGeometry initial_guess_;
    bool has_initial_guess_ = false;
    
    BellGeometry best_geometry_;
    double best_error_ = BellConstants::ERROR_INVALID;
    ModalData best_modal_data_;
    ModalData modal_reference_;
    bool has_modal_reference_ = false;
    
    std::random_device rd_;
    std::mt19937 rng_;
    std::unique_ptr<ModalSolver> solver_;
    
    static BellGeometry getDefaultGeometry();
    std::vector<BellGeometry> getReasonableGeometries();
    
    void saveBestResult();
    void copyFile(const std::string& src, const std::string& dst);
    
    BellGeometry generateRandomGeometry();
    Candidate createRandomCandidate();
    Candidate createNearCandidate(const BellGeometry& base, double error);
    
    double evaluateCandidate(BellGeometry& geometry, ModalData& out_data);
    double computeError(const ModalData& modal_data, const BellGeometry& geometry);
    
    std::vector<Candidate> initializePopulation();
    void evaluatePopulation(std::vector<Candidate>& population);
    std::vector<Candidate> selectElite(const std::vector<Candidate>& population);
    Candidate crossover(const Candidate& p1, const Candidate& p2);
    void mutate(BellGeometry& geom);
    static BellGeometry canonicalizeShape(const BellGeometry& geometry, double reference_height = 1000.0);
    
    std::vector<double> generateMonotonicRadii(double Rmouth, int num_internal = 3);
    
    double getMinRmouth(double H) const { return BellConstants::MIN_RMOUTH_RATIO * H; }
    double getMaxRmouth(double H) const { return BellConstants::MAX_RMOUTH_RATIO * H; }
    
    void printHeader();
    void printGeneration(int gen, const Candidate& best);
    void printFinal(const Candidate& best);
};

#endif // BELL_OPTIMIZER_H