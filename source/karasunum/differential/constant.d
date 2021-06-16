/**
Differentiable constant type.
*/
module karasunum.differential.constant;

import karasunum.differential.differentiable :
    Differentiable,
    DiffContext,
    diffContext,
    EvalContext;

@safe:

/**
Differentiable static constant class.

Params:
    R = result type.
    value = constant value.
*/
final class StaticConstant(R, R value) : Differentiable!R
{
    private this() immutable @nogc nothrow pure return scope {}

    override R opCall() const nothrow pure return scope
    {
        return R(value);
    }

    override R evaluate(scope EvalContext!R context) const nothrow pure
        in (false)
    {
        return R(value);
    }

    override const(Differentiable!R) differentiate(scope DiffContext!R context) const nothrow pure return scope
    {
        return context.zero;
    }
}

/**
Create zero constant.

Returns:
    Zero constant.
*/
auto zero(R)() nothrow pure
    out (r; r)
{
    return new immutable(StaticConstant!(R, 0))();
}

///
nothrow pure unittest
{
    import std.math : isClose;

    // evaluated is zero.
    auto z = zero!real();
    assert(z().isClose(0.0));

    // differentiated.
    auto context = diffContext(z);
    auto d = z.differentiate(context);
    assert(d is context.zero);
}

/**
Create one constant.

Returns:
    One constant.
*/
auto one(R)() nothrow pure
    out (r; r)
{
    return new immutable(StaticConstant!(R, 1))();
}

///
nothrow pure unittest
{
    import std.math : isClose;

    // evaluated is zero.
    auto o = one!real();
    assert(o().isClose(1.0));

    // differentiate
    auto context = diffContext(o);
    auto d = o.differentiate(context);
    assert(d is context.zero);
}

/**
Create two constant.

Returns:
    Two constant.
*/
auto two(R)() nothrow pure
    out (r; r)
{
    return new immutable(StaticConstant!(R, 2))();
}

///
nothrow pure unittest
{
    import std.math : isClose;

    // evaluated is two.
    auto o = two!real();
    assert(o().isClose(2.0));

    // differentiate
    auto context = diffContext(o);
    auto d = o.differentiate(context);
    assert(d is context.zero);
}

/**
Differentiable constant class.

Params:
    R = result type.
*/
final class Constant(R) : Differentiable!R
{
    override R opCall() const nothrow pure return scope
    {
        return value_;
    }

    override R evaluate(scope EvalContext!R context) const nothrow pure
        in (false)
    {
        return value_;
    }

    override const(Differentiable!R) differentiate(scope DiffContext!R context) const nothrow pure return scope
        in (false)
    {
        return context.zero;
    }

private:
    R value_;

    this()(auto return scope ref const(R) value) immutable nothrow pure return scope
    {
        this.value_ = value;
    }
}

/**
Create a constant.

Returns:
    A constant.
*/
immutable(Constant!R) constant(R)(auto return scope ref const(R) value) nothrow pure
    out (r; r)
{
    return new immutable(Constant!R)(value);
}

///
nothrow pure unittest
{
    import std.math : isClose;

    // evaluated is constant value.
    auto c = constant(1.234);
    assert(c().isClose(1.234));

    // differentiate
    auto context = diffContext(c);
    auto d = c.differentiate(context);
    assert(d is context.zero);
}

