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
    log;

@safe:

private template RebindableType(T)
{
    static if (is(T == class) || is(T == interface) || isDynamicArray!T || isAssociativeArray!T)
    {
        alias RebindableType = Rebindable!T;
    }
    else
    {
        alias RebindableType = T;
    }
}

struct KalmanFilter(F, P)
{
    struct Parameters
    {
        P drift;
        P tension;
        P cons;
        P measureVariance;
        P stateVariance;
        size_t likelihoodSkipCount = 0;
    }

    @disable this();

    this(Parameters parameters, F initState, F initVariance, F one) nothrow pure scope
    {
        this.params_ = parameters;
        this.state_ = initState;
        this.variance_ = initVariance;
        this.one_ = one;
    }

    F estimate(F x) nothrow pure return scope
    {
        estimateState_ = params_.drift + params_.tension * state_;
        estimateMeasure_ = params_.cons + estimateState_ * x;
        estimateVariance_ = params_.tension.square * variance_ + params_.stateVariance;
        return estimateMeasure_;
    }

    F filtering(F x, F y) nothrow pure scope
    {
        auto error = y - estimateMeasure_;
        auto errorVariance = estimateVariance_ * x.square + params_.measureVariance;
        auto currentLikelihood = (log(errorVariance) + error.square / errorVariance);

        if (time_ >= params_.likelihoodSkipCount)
        {
            likelihood_ = (likelihood_) ? (likelihood_ + currentLikelihood) : currentLikelihood;
        }
        ++time_;

        auto k = (x * estimateVariance_)
            / (x.square * estimateVariance_ + params_.measureVariance);
        state_ = estimateState_ + k * error;
        variance_ = (one_ - x * k) * estimateVariance_;
        return state_;
    }

    @property F likelihood() const @nogc nothrow pure scope
    {
        return likelihood_;
    }

private:
    Parameters params_;
    F one_;
    RebindableType!F state_;
    RebindableType!F variance_;

    RebindableType!F estimateState_;
    RebindableType!F estimateVariance_;
    RebindableType!F estimateMeasure_;
    RebindableType!F likelihood_;
    size_t time_;
}

