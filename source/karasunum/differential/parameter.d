/**
Differentiable parameter type.
*/
module karasunum.differential.parameter;

import karasunum.differential.differentiable :
    Differentiable,
    DiffContext,
    EvalContext;

@safe:

/**
Differentiable plus minus class.

Params:
    R = result type.
*/
final class Parameter(R) : Differentiable!R
{
    this()(auto return scope ref const(R) value) nothrow pure return scope
    {
        bind(value);
    }

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
        return (context.target is this) ? context.one : context.zero;
    }

    /**
    Get current value.

    Returns:
        Current value.
    */
    @property R value() const @nogc nothrow pure scope
    {
        return this.value_;
    }

    /**
    Bind new value.

    Params:
        value = new value.
    */
    void bind()(auto return scope ref const(R) value) nothrow pure return scope
    {
        this.value_ = value;
    }

private:
    R value_;
}

Parameter!R param(R)(auto return scope ref const(R) value) nothrow pure
{
    return new Parameter!R(value);
}

///
nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.differentiable : diffContext;

    auto p = param(1.0);
    assert(p().isClose(1.0));
    assert(p.value.isClose(1.0));

    auto context = diffContext(p);
    auto d = p.differentiate(context);
    assert(d is context.one);

    auto p2 = param(1.1);
    auto d2 = p2.differentiate(context);
    assert(d2 is context.zero);

    p.bind(100.0);
    assert(p.value.isClose(100.0));
}

