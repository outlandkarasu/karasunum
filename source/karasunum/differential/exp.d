/**
Differentiable exp function module.
*/
module karasunum.differential.exp;

import std.math : mathExp = exp;
import std.traits : isFloatingPoint, isNumeric;

import karasunum.differential.differentiable :
    Differentiable,
    DiffContext,
    EvalContext;
import karasunum.differential.mul : mul;

@safe:

/**
Differentiable exp class.

Params:
    R = result type.
*/
private final class Exp(R) : Differentiable!R
    if (isFloatingPoint!R)
{
    this(const(Differentiable!R) x) const @nogc nothrow pure scope
        in (x)
    {
        this.x_ = x;
    }

    override R opCall() const @nogc nothrow pure return scope
    {
        return mathExp(x_());
    }

    override R evaluate(scope EvalContext!R context) const nothrow pure
    {
        return mathExp(context.evaluate(x_));
    }

    const(Differentiable!R) differentiate(scope DiffContext!R context) const nothrow pure return scope
        in (false)
    {
        auto xDiff = context.diff(x_);
        return context.mul(this, xDiff);
    }

private:
    const(Differentiable!R) x_;
}

const(Exp!R) exp(R)(const(Differentiable!R) x) nothrow pure
    if (isFloatingPoint!R)
{
    return new const(Exp!R)(x);
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : diffContext;
    import karasunum.differential.parameter : param;
    import karasunum.differential.pow : square;

    auto p = param(3.0);
    auto pexp = p.exp();
    assert(pexp().isClose(mathExp(3.0)));

    auto context = p.diffContext;
    auto pexpd = pexp.differentiate(context);
    assert(pexpd().isClose(mathExp(3.0)));

    auto pexpx2 = exp(p.square);
    auto pexpx2d = pexpx2.differentiate(context);
    assert(pexpx2d().isClose(mathExp(9.0) * 6.0));
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : evalContext;
    import karasunum.differential.parameter : param;
    import karasunum.differential.pow : square;

    auto p = param(3.0);
    auto pexp = p.exp();
    auto context = evalContext!double();
    assert(context.evaluate(pexp).isClose(mathExp(3.0)));
    assert(context.callCount == 2);
    assert(context.evaluateCount == 2);
    assert(context.cacheHitCount == 0);

    assert(context.evaluate(pexp).isClose(mathExp(3.0)));
    assert(context.callCount == 3);
    assert(context.evaluateCount == 2);
    assert(context.cacheHitCount == 1);
}

/**
exp for numeric.

Params:
    x = value
Returns:
    exp value.
*/
T exp(T)(T x) @nogc nothrow pure
    if (isNumeric!T)
{
    return mathExp(x);
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;
    assert(exp(4.0).isClose(mathExp(4.0)));
}

