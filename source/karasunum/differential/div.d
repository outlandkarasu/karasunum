/**
Differentiable division function module.
*/
module karasunum.differential.div;

import karasunum.differential.differentiable :
    Differentiable,
    DiffContext,
    EvalContext;
import karasunum.differential.add_sub : sub;
import karasunum.differential.mul : mul;
import karasunum.differential.pow : square;

@safe:

/**
Differentiable division class.

Params:
    R = result type.
*/
final class Division(R) : Differentiable!R
{
    this(const(Differentiable!R) lhs, const(Differentiable!R) rhs) const @nogc nothrow pure scope
        in (lhs && rhs)
    {
        this.lhs_ = lhs;
        this.rhs_ = rhs;
    }

    override R opCall() const @nogc nothrow pure return scope
    {
        return lhs_() / rhs_();
    }

    override R evaluate(scope EvalContext!R context) const nothrow pure
    {
        return context.evaluate(lhs_) / context.evaluate(rhs_);
    }

    const(Differentiable!R) differentiate(scope DiffContext!R context) const nothrow pure return scope
        in (false)
    {
        auto lhsDiff = context.diff(lhs_);
        auto rhsDiff = context.diff(rhs_);
        auto ldy = context.mul(lhsDiff, rhs_);
        auto rdy = context.mul(lhs_, rhsDiff);
        auto numerator = context.sub(ldy, rdy);
        return context.div(numerator, rhs_.square);
    }

private:
    const(Differentiable!R) lhs_;
    const(Differentiable!R) rhs_;
}

const(Division!R) div(R)(const(Differentiable!R) lhs, const(Differentiable!R) rhs) nothrow pure
{
    return new const(Division!R)(lhs, rhs);
}

pure nothrow unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : diffContext;
    import karasunum.differential.parameter : param;

    auto p1 = param(2.0);
    auto p2 = param(3.0);
    auto m = p1.div(p2);
    assert(m().isClose(2.0/3.0));

    auto p1d = m.differentiate(p1.diffContext);
    assert(p1d().isClose(3.0/9.0));

    auto p2d = m.differentiate(p2.diffContext);
    assert(p2d().isClose(-2.0/9.0));
}

pure nothrow unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : evalContext;
    import karasunum.differential.parameter : param;

    auto p1 = param(2.0);
    auto p2 = param(3.0);
    auto m = p1.div(p2);
    auto context = evalContext!double();
    assert(context.evaluate(m).isClose(2.0/3.0));
    assert(context.callCount == 3);
    assert(context.evaluateCount == 3);
    assert(context.cacheHitCount == 0);

    assert(context.evaluate(m).isClose(2.0/3.0));
    assert(context.callCount == 4);
    assert(context.evaluateCount == 3);
    assert(context.cacheHitCount == 1);
}

const(Differentiable!R) div(R)(scope DiffContext!R context, const(Differentiable!R) lhs, const(Differentiable!R) rhs) nothrow pure
{
    if (context.isOne(rhs))
    {
        return lhs;
    }

    return new const(Division!R)(lhs, rhs);
}

pure nothrow unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : diffContext;
    import karasunum.differential.parameter : param;

    auto p1 = param(2.0);
    auto p2 = param(3.0);
    auto context = diffContext(p1);
    assert(context.div(p1, p2)().isClose(2.0/3.0));
    assert(context.div(p1, context.one) is p1);
}

