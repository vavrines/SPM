#include "quadratures/quadraturebase.hpp"
#include "common/config.hpp"
#include "quadratures/qgausschebyshev1D.hpp"
#include "quadratures/qgausslegendre1D.hpp"
#include "quadratures/qgausslegendretensorized.hpp"
#include "quadratures/qldfesa.hpp"
#include "quadratures/qlebedev.hpp"
#include "quadratures/qlevelsymmetric.hpp"
#include "quadratures/qmontecarlo.hpp"
#include "quadratures/qproduct.hpp"
#include "toolboxes/errormessages.hpp"

QuadratureBase::QuadratureBase( Config* settings ) {
    _settings = settings;
    _order    = settings->GetQuadOrder();
}

QuadratureBase::QuadratureBase( unsigned order ) : _order( order ) { _settings = nullptr; }

QuadratureBase* QuadratureBase::Create( Config* settings ) {
    QUAD_NAME name = settings->GetQuadName();

    switch( name ) {
        case QUAD_MonteCarlo: return new QMonteCarlo( settings );
        case QUAD_GaussLegendreTensorized: return new QGaussLegendreTensorized( settings );
        case QUAD_GaussLegendre1D: return new QGaussLegendre1D( settings );
        case QUAD_GaussLegendreTensorized2D: return new QGaussLegendreTensorized2D( settings );

        case QUAD_GaussChebyshev1D: return new QGaussChebyshev1D( settings );
        case QUAD_LevelSymmetric: return new QLevelSymmetric( settings );
        case QUAD_LDFESA: return new QLDFESA( settings );
        case QUAD_Lebedev: return new QLebedev( settings );
        case QUAD_Product: return new QProduct( settings );
        default: ErrorMessages::Error( "Creator for the chosen quadrature does not yet exist. This is the fault of the coder!", CURRENT_FUNCTION );
    }
    return nullptr;
}

QuadratureBase* QuadratureBase::Create( QUAD_NAME name, unsigned quadOrder ) {

    switch( name ) {
        case QUAD_MonteCarlo: return new QMonteCarlo( quadOrder );
        case QUAD_GaussLegendreTensorized:
            ErrorMessages::Error( "This quadrature must be initialized with a settings constructor!", CURRENT_FUNCTION );
            break;
        case QUAD_GaussLegendre1D: return new QGaussLegendre1D( quadOrder );
        case QUAD_GaussChebyshev1D: return new QGaussChebyshev1D( quadOrder );
        case QUAD_LevelSymmetric: return new QLevelSymmetric( quadOrder );
        case QUAD_LDFESA: return new QLDFESA( quadOrder );
        case QUAD_Lebedev: return new QLebedev( quadOrder );
        case QUAD_Product: return new QProduct( quadOrder );
        default: ErrorMessages::Error( "Creator for the chosen quadrature does not yet exist. This is the fault of the coder!", CURRENT_FUNCTION );
    }
    return nullptr;
}

double QuadratureBase::Integrate( double ( *f )( double, double, double ) ) {
    double result = 0.0;
    for( unsigned i = 0; i < _nq; i++ ) {
        result += _weights[i] * f( _points[i][0], _points[i][1], _points[i][2] );
    }
    return result;
}

double QuadratureBase::IntegrateSpherical( double ( *f )( double, double ) ) {
    double result = 0.0;
    for( unsigned i = 0; i < _nq; i++ ) {
        result += _weights[i] * f( _pointsSphere[i][0], _pointsSphere[i][1] );
    }
    return result;
}

std::vector<double> QuadratureBase::Integrate( std::vector<double> ( *f )( double, double, double ), unsigned len ) {
    std::vector<double> result( len, 0.0 );
    std::vector<double> funcEval( len, 0.0 );

    for( unsigned i = 0; i < _nq; i++ ) {
        funcEval = f( _points[i][0], _points[i][1], _points[i][2] );
        for( unsigned idx_len = 0; idx_len < len; idx_len++ ) {
            result[idx_len] += _weights[i] * funcEval[idx_len];
        }
    }
    return result;
}

double QuadratureBase::SumUpWeights() { return sum( _weights ); }

void QuadratureBase::PrintWeights() {
    auto log = spdlog::get( "event" );
    for( unsigned i = 0; i < _nq; i++ ) {
        log->info( _weights[i] );
    }
}

void QuadratureBase::PrintPoints() {
    auto log = spdlog::get( "event" );
    for( unsigned i = 0; i < _nq; i++ ) {
        log->info( "{0}, {1}, {2}", _points[i][0], _points[i][1], _points[i][2] );
    }
}
void QuadratureBase::PrintPointsAndWeights() {
    auto log = spdlog::get( "event" );
    for( unsigned i = 0; i < _nq; i++ ) {
        log->info( "{0}, {1}, {2}, {3}", _points[i][0], _points[i][1], _points[i][2], _weights[i] );
    }
}

void QuadratureBase::ScaleWeights( double leftBound, double rightBound ) {
    if( !_settings ) {
        ErrorMessages::Error( "This function is only available with an active settings file.", CURRENT_FUNCTION );
    }
    if( _settings->GetDim() != 1 ) {
        ErrorMessages::Error( "This function is only available in 1 spatial dimension.", CURRENT_FUNCTION );
    }

    _weights = _weights * ( ( rightBound - leftBound ) / 2.0 );
}
