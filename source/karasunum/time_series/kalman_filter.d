/**
Kalman filter module.
*/
module karasunum.time_series.kalman_filter;

import std.traits : isDynamicArray, isAssociativeArray;
import std.typecons : Rebindable;

import karasunum.differential :
    Parameter,
    Differentiable,
    constant,
    square,
    log,
    exp;

@safe:

/**
Kalman filter.

Params:
    F = functional type.
    P = parameter type.
*/
struct KalmanFilter(F, P)
{
    /**
    Kalman filter parameters.
    */
    struct Parameters
    {
        P drift;
        P tension;
        P cons;
        P logMeasureVariance;
        P logStateVariance;
        size_t likelihoodSkipCount = 0;
    }

    @disable this();

    /**
    Construct from parameters.

    Params:
        parameters = Kalman filter parameters.
        initState = initial state value.
        initVariance = initial kalman filter variance.
    */
    this(Parameters parameters, F initState, F initVariance) nothrow pure scope
    {
        this.params_ = parameters;
        this.tensionSquare_ = parameters.tension.square;
        this.stateVariance_ = exp(parameters.logStateVariance);
        this.measureVariance_ = exp(parameters.logMeasureVariance);
        this.state_ = initState;
        this.variance_ = initVariance;
    }

    /**
    Estimate next value.

    Params:
        x = x value.
    Returns:
        estimated value.
    */
    F estimate(F x) nothrow pure return scope
    {
        estimateState_ = params_.drift + params_.tension * state_;
        estimateMeasure_ = params_.cons + estimateState_ * x;
        estimateVariance_ = tensionSquare_ * variance_ + stateVariance_;
        return estimateMeasure_;
    }

    /**
    Filtering kalman filter.

    Params:
        x = x value.
        y = y value.
    Returns:
        filtered state value.
    */
    F filtering(F x, F y) nothrow pure scope
    {
        auto error = y - estimateMeasure_;
        auto errorVariance = estimateVariance_ * x.square + measureVariance_;
        auto currentLikelihood = log(errorVariance) + error.square / errorVariance;

        if (time_ == params_.likelihoodSkipCount)
        {
            likelihood_ = currentLikelihood;
        }
        else if (time_ > params_.likelihoodSkipCount)
        {
            likelihood_ = likelihood_ + currentLikelihood;
        }

        ++time_;

        auto xev = x * estimateVariance_;
        auto k = xev / (x * xev + measureVariance_);
        state_ = estimateState_ + k * error;
        variance_ = estimateVariance_ - xev * k;
        return state_;
    }

    /**
    Returns:
        likelihood value.
    */
    @property F likelihood() const @nogc nothrow pure scope
    {
        return likelihood_;
    }

private:
    Parameters params_;
    F tensionSquare_;
    F stateVariance_;
    F measureVariance_;
    RebindableType!F state_;
    RebindableType!F variance_;

    RebindableType!F estimateState_;
    RebindableType!F estimateVariance_;
    RebindableType!F estimateMeasure_;
    RebindableType!F likelihood_;
    size_t time_;
}

///
@nogc nothrow pure @safe unittest
{
    import std.math : isClose;

    alias KF = KalmanFilter!(double, double);

    immutable KF.Parameters parameters = {
        drift: 1.0,
        tension: 2.0,
        cons: 3.0,
        logMeasureVariance: log(10.0),
        logStateVariance: log(10.0),
        likelihoodSkipCount: 0,
    };

    auto kalmanFilter = KF(parameters, 1.0, 2.0);
    assert(kalmanFilter.estimate(3.0).isClose(12.0));

    immutable filtered = kalmanFilter.filtering(3.0, 12.5);
    assert(filtered.isClose(3.156977, 1e-5, 1e-5));
    assert(kalmanFilter.likelihood.isClose(5.14895, 1e-5, 1e-5));

    assert(kalmanFilter.estimate(1.0).isClose(10.314, 1e-5, 1e-5));
    immutable filtered2 = kalmanFilter.filtering(3.0, 10.0);
    assert(filtered2.isClose(7.2169, 1e-5, 1e-5));
    assert(kalmanFilter.likelihood.isClose(10.0746, 1e-5, 1e-5));
}

private:

enum isReferenceType(T)
    = is(T == class)
    || is(T == interface)
    || isDynamicArray!T
    || isAssociativeArray!T;

template RebindableType(T)
{
    static if (isReferenceType!T)
    {
        alias RebindableType = Rebindable!T;
    }
    else
    {
        alias RebindableType = T;
    }
}

