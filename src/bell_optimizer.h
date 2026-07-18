// bell_optimizer.h
#ifndef BELL_OPTIMIZER_H
#define BELL_OPTIMIZER_H

#include "bell_types.h"
#include "bell_generator.h"
#include "ccx_parser.hpp"
#include <vector>
#include <random>
#include <string>

using namespace BellOptimizer;

class EvolutionaryOptimizer {
public:
    EvolutionaryOptimizer(const BellMaterial& material, 
                          const TargetData& targets,
                          const EvolutionParams& evo_params);
    
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
        CCXData ccx_data;
        bool valid = false;
    };

private:
    BellMaterial material_;
    TargetData targets_;
    EvolutionParams evo_params_;
    BellGeometry initial_guess_;
    bool has_initial_guess_ = false;
    
    BellGeometry best_geometry_;
    double best_error_ = BellConstants::ERROR_INVALID;
    CCXData best_ccx_data_;
    
    std::random_device rd_;
    std::mt19937 rng_;
    
    static BellGeometry getDefaultGeometry();
    static std::vector<BellGeometry> getReasonableGeometries();
    
    void saveBestResult();
    void copyFile(const std::string& src, const std::string& dst);
    
    BellGeometry generateRandomGeometry();
    Candidate createRandomCandidate();
    Candidate createNearCandidate(const BellGeometry& base, double error);
    
    double evaluateCandidate(BellGeometry& geometry, CCXData& out_data);
    double computeError(const CCXData& ccx_data, const BellGeometry& geometry);
    
    std::vector<Candidate> initializePopulation();
    void evaluatePopulation(std::vector<Candidate>& population);
    std::vector<Candidate> selectElite(const std::vector<Candidate>& population);
    Candidate crossover(const Candidate& p1, const Candidate& p2);
    void mutate(BellGeometry& geom);
    
    std::vector<double> generateMonotonicRadii(double Rmouth, int num_internal = 3);
    
    double getMinRmouth(double H) const { return BellConstants::MIN_RMOUTH_RATIO * H; }
    double getMaxRmouth(double H) const { return BellConstants::MAX_RMOUTH_RATIO * H; }
    
    void printHeader();
    void printGeneration(int gen, const Candidate& best);
    void printFinal(const Candidate& best);
};

#endif // BELL_OPTIMIZER_H