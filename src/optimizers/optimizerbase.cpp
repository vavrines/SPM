#include "optimizers/optimizerbase.hpp"
#include "common/config.hpp"
#include "entropies/entropybase.hpp"
#include "optimizers/newtonoptimizer.hpp"

OptimizerBase::OptimizerBase( Config* settings ) {
    _entropy  = EntropyBase::Create( settings );
    _settings = settings;
}

OptimizerBase::~OptimizerBase() { delete _entropy; }

OptimizerBase* OptimizerBase::Create( Config* settings ) {
    switch( settings->GetOptimizerName() ) {
        case NEWTON: return new NewtonOptimizer( settings );
        default: return new NewtonOptimizer( settings );
    }
}
