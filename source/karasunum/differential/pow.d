/**
Differentiable power function module.
*/
module karasunum.differential.pow;

import std.math : mathLog = log, mathPow = pow;
import std.traits : isNumeric;

import karasunum.differential.differentiable :
    Differentiable,
    DiffContext,
    EvalContext;
import karasunum.differential.mul : mul;
import karasunum.differential.div : div;
import karasunum.differential.add_sub : add;
import karasunum.differential.log : log;

@safe:

/**
Differentiable square class.

Params:
    R = result type.
*/
private final class Square(R) : Differentiable!R
{
    this(const(Differentiable!R) x) const @nogc nothrow pure scope
        in (x)
    {
        this.x_ = x;
    }

    override R opCall() const @nogc nothrow pure return scope
    {
        auto x = x_();
        return x * x;
    }

    override R evaluate(scope EvalContext!R context) const nothrow pure
    {
        auto x = context.evaluate(x_);
        return x * x;
    }

    const(Differentiable!R) differentiate(scope DiffContext!R context) const nothrow pure return scope
    {
        auto xDiff = context.diff(x_);
        return context.mul(context.mul(context.two, x_), xDiff);
    }

private:
    const(Differentiable!R) x_;
}

const(Square!R) square(R)(const(Differentiable!R) x) nothrow pure
{
    return new const(Square!R)(x);
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : diffContext;
    import karasunum.differential.parameter : param;

    auto p = param(3.0);
    auto p2 = p.square();
    assert(p2().isClose(9.0));

    auto p2d = p2.differentiate(p.diffContext);
    assert(p2d().isClose(6.0));

    auto p2dd = p2d.differentiate(p.diffContext);
    assert(p2dd().isClose(2.0));
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : evalContext;
    import karasunum.differential.parameter : param;

    auto p = param(3.0);
    auto p2 = p.square();
    auto context = evalContext!double();
    assert(context.evaluate(p2).isClose(9.0));
    assert(context.callCount == 2);
    assert(context.evaluateCount == 2);
    assert(context.cacheHitCount == 0);

    assert(context.evaluate(p2).isClose(9.0));
    assert(context.callCount == 3);
    assert(context.evaluateCount == 2);
    assert(context.cacheHitCount == 1);
}

/**
Square for numeric.

Params:
    x = value
Returns:
    squared value.
*/
T square(T)(T x) @nogc nothrow pure if(isNumeric!T)
{
    return x ^^ T(2);
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;
    assert(square(4.0).isClose(16.0));
}

/**
Differentiable power class.

Params:
    R = result type.
*/
final class Power(R) : Differentiable!R
{
    this(const(Differentiable!R) lhs, const(Differentiable!R) rhs) const @nogc nothrow pure scope
        in (lhs && rhs)
    {
        this.lhs_ = lhs;
        this.rhs_ = rhs;
    }

    override R opCall() const @nogc nothrow pure return scope
    {
        return mathPow(lhs_(), rhs_());
    }

    override R evaluate(scope EvalContext!R context) const nothrow pure
    {
        return mathPow(context.evaluate(lhs_), context.evaluate(rhs_));
    }

    const(Differentiable!R) differentiate(scope DiffContext!R context) const nothrow pure return scope
        in (false)
    {
        auto lhsDiff = context.diff(lhs_);
        auto rhsDiff = context.diff(rhs_);
        auto ld = context.mul(lhsDiff, context.div(rhs_, lhs_));
        auto rd = context.mul(rhsDiff, log(lhs_));
        return context.mul(this, context.add(ld, rd));
    }

private:
    const(Differentiable!R) lhs_;
    const(Differentiable!R) rhs_;
}

const(Power!R) pow(R)(const(Differentiable!R) lhs, const(Differentiable!R) rhs) nothrow pure
{
    return new const(Power!R)(lhs, rhs);
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : diffContext;
    import karasunum.differential.parameter : param;

    auto p1 = param(2.0);
    auto p2 = param(3.0);
    auto m = p1.pow(p2);
    assert(m().isClose(8.0));

    auto p1d = m.differentiate(p1.diffContext);
    assert(p1d().isClose(12.0));

    auto p2d = m.differentiate(p2.diffContext);
    assert(p2d().isClose(8.0 * mathLog(2.0)));
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : evalContext;
    import karasunum.differential.parameter : param;

    auto p1 = param(2.0);
    auto p2 = param(3.0);
    auto m = p1.pow(p2);
    auto context = evalContext!double();
    assert(context.evaluate(m).isClose(8.0));
    assert(context.callCount == 3);
    assert(context.evaluateCount == 3);
    assert(context.cacheHitCount == 0);

    assert(context.evaluate(m).isClose(8.0));
    assert(context.callCount == 4);
    assert(context.evaluateCount == 3);
    assert(context.cacheHitCount == 1);
}

/**
Power for numeric.

Params:
    lhs = lhs value
    rhs = rhs value
Returns:
    powered value.
*/
T pow(T)(T lhs, T rhs) @nogc nothrow pure if(isNumeric!T)
{
    return mathPow(lhs, rhs);
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;
    assert(pow(4.0, 2.0).isClose(16.0));
}

