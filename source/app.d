import std.stdio;

import std.algorithm : map;
import std.array : array;
import std.math : isClose;

import karasunum.time_series.kalman_filter : KalmanFilter;
import karasunum.differential :
    Parameter,
    param,
    constant,
    Differentiable,
    one,
    zero,
    diffContext,
    evalContext;
import karasunum.optimize : NewtonMethod;

immutable PRICES = [
    106.138, 106.140, 106.138, 106.142, 106.142,
    106.146, 106.150, 106.140, 106.140, 106.140,
    106.139, 106.137, 106.144, 106.145, 106.148,
    106.148, 106.141, 106.146, 106.146, 106.144,
    106.146, 106.148, 106.146, 106.148, 106.150,
    106.148, 106.149, 106.150, 106.149, 106.152,
    106.152, 106.150, 106.150, 106.136, 106.130,
    106.118, 106.122, 106.124, 106.120, 106.116,
    106.124, 106.116, 106.114, 106.115, 106.126,
    106.134, 106.130, 106.129, 106.120, 106.126,
    106.131, 106.136, 106.133, 106.132, 106.130,
    106.142, 106.146, 106.144, 106.152, 106.157,
    106.152, 106.144, 106.148, 106.143, 106.144,
    106.142, 106.140, 106.149, 106.146, 106.144,
    106.144, 106.146, 106.146, 106.144, 106.144,
    106.151, 106.158, 106.157, 106.156, 106.157,
    106.166, 106.162, 106.158, 106.152, 106.152,
    106.145, 106.144, 106.160, 106.150, 106.160,
    106.155, 106.150, 106.141, 106.156, 106.149,
    106.144, 106.132, 106.127, 106.132, 106.138,
    106.135, 106.128, 106.126, 106.136, 106.125,
    106.127, 106.130, 106.140, 106.140, 106.141,
    106.140, 106.138, 106.140, 106.140, 106.144,
    106.148, 106.140, 106.138, 106.133, 106.135,
];

void main()
{
    alias KF = KalmanFilter!(const(Differentiable!double), const(Differentiable!double));
    auto drift = param!double(0.0);
    auto tension = param!double(1.0);
    auto measureVariance = param!double(0.01);
    auto stateVariance = param!double(0.01);
    auto initState = constant!double(106.138);
    auto initVariance = constant!double(1.0);
    immutable oneValue = one!double();
    auto variables = [drift, tension, measureVariance, stateVariance];

    KF.Parameters parameters = {
        drift: drift,
        tension: tension,
        cons: zero!double(),
        measureVariance: measureVariance,
        stateVariance: stateVariance,
        likelihoodSkipCount: 1,
    };
    scope kalmanFilter = KF(parameters, initState, initVariance, oneValue);

    foreach (price; PRICES[0 .. 60])
    {
        kalmanFilter.estimate(oneValue);
        kalmanFilter.filtering(oneValue, constant(price));
    }

    auto lf = kalmanFilter.likelihood;
    auto dVariables = variables.map!((v) => lf.differentiate(v.diffContext)).array;
    scope newton = new NewtonMethod!(double, 4)(dVariables, variables);
    double likelihood = double.nan;
    foreach (i; 0 .. 20)
    {
        newton.step();

        scope ec = evalContext!double();
        immutable newLikelihood = ec.evaluate(lf);
        if (likelihood.isClose(newLikelihood))
        {
            break;
        }

        debug
        {
            import std.stdio : writefln;
            writefln("%s", likelihood);
        }

        likelihood = newLikelihood;
   }
}

