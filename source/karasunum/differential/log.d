/**
Differentiable log function module.
*/
module karasunum.differential.log;

import std.math : mathLog = log;
import std.traits : isFloatingPoint, isNumeric;

import karasunum.differential.differentiable :
    Differentiable,
    DiffContext,
    EvalContext;
import karasunum.differential.mul : mul;
import karasunum.differential.div : div;

@safe:

/**
Differentiable log class.

Params:
    R = result type.
*/
private final class Log(R) : Differentiable!R
    if (isFloatingPoint!R)
{
    this(const(Differentiable!R) x) const @nogc nothrow pure scope
        in (x)
    {
        this.x_ = x;
    }

    override R opCall() const @nogc nothrow pure return scope
    {
        return mathLog(x_());
    }

    override R evaluate(scope EvalContext!R context) const nothrow pure
    {
        return mathLog(context.evaluate(x_));
    }

    const(Differentiable!R) differentiate(scope DiffContext!R context) const nothrow pure return scope
        in (false)
    {
        auto xDiff = context.diff(x_);
        auto dlog = context.div(context.one, x_);
        return context.mul(dlog, xDiff);
    }

private:
    const(Differentiable!R) x_;
}

const(Log!R) log(R)(const(Differentiable!R) x) nothrow pure
    if (isFloatingPoint!R)
{
    return new const(Log!R)(x);
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : diffContext;
    import karasunum.differential.parameter : param;
    import karasunum.differential.pow : square;

    auto p = param(3.0);
    auto plog = p.log();
    assert(plog().isClose(mathLog(3.0)));

    auto context = p.diffContext;
    auto plogd = plog.differentiate(context);
    assert(plogd().isClose(1.0/3.0));

    auto plog2x = log(p.square);
    auto plog2xd = plog2x.differentiate(context);
    assert(plog2xd().isClose(2.0/3.0));
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : evalContext;
    import karasunum.differential.parameter : param;
    import karasunum.differential.pow : square;

    auto p = param(3.0);
    auto plog = p.log();
    auto context = evalContext!double();
    assert(context.evaluate(plog).isClose(mathLog(3.0)));
    assert(context.callCount == 2);
    assert(context.evaluateCount == 2);
    assert(context.cacheHitCount == 0);

    assert(context.evaluate(plog).isClose(mathLog(3.0)));
    assert(context.callCount == 3);
    assert(context.evaluateCount == 2);
    assert(context.cacheHitCount == 1);
}

/**
log for numeric.

Params:
    x = value
Returns:
    log value.
*/
T log(T)(T x) @nogc nothrow pure
    if (isNumeric!T)
{
    return mathLog(x);
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;
    assert(log(4.0).isClose(mathLog(4.0)));
}

