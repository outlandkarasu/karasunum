/**
Differentiable type.
*/
module karasunum.differential.differentiable;

import std.typecons : Rebindable, rebindable;
import std.exception : assumeWontThrow;

import karasunum.differential.add_sub : Addition, add, Subtraction, sub;
import karasunum.differential.constant : zero, one, two, constant;
import karasunum.differential.div : Division, div;
import karasunum.differential.mul : Multiply, mul;
import karasunum.differential.pow : Power, pow;

@safe:

/**
Differentiable interface.

Params:
    R = result type.
*/
interface Differentiable(R)
{
    /**
    Evaluate function.

    Returns;
        function result.
    */
    R opCall() const @nogc nothrow pure return scope;

    /**
    Evaluate function.

    Params:
        context = evaluate context.
    Returns;
        function result.
    */
    R evaluate(scope EvalContext!R context) const nothrow pure;

    /**
    Differentiate function.

    Params:
        context = differentiate context.
    Returns:
        Differentiate result.
    */
    const(Differentiable!R) differentiate(scope DiffContext!R context) const nothrow pure return scope
        in (context)
        out (r; r);

    /**
    Add operator.

    Params:
        rhs = right hand side.
    Returns:
        add expression.
    */
    final const(Addition!R) opBinary(string op)(const(Differentiable!R) rhs) const nothrow pure if (op == "+")
        in (rhs)
        out (r; r)
    {
        return this.add(rhs);
    }

    /**
    Subtract operator.

    Params:
        rhs = right hand side.
    Returns:
        subtract expression.
    */
    final const(Subtraction!R) opBinary(string op)(const(Differentiable!R) rhs) const nothrow pure if (op == "-")
        in (rhs)
        out (r; r)
    {
        return this.sub(rhs);
    }

    /**
    Multiply operator.

    Params:
        rhs = right hand side.
    Returns:
        multiply expression.
    */
    final const(Multiply!R) opBinary(string op)(const(Differentiable!R) rhs) const nothrow pure if (op == "*")
        in (rhs)
        out (r; r)
    {
        return this.mul(rhs);
    }

    /**
    Divide operator.

    Params:
        rhs = right hand side.
    Returns:
        divite expression.
    */
    final const(Division!R) opBinary(string op)(const(Differentiable!R) rhs) const nothrow pure if (op == "/")
        in (rhs)
        out (r; r)
    {
        return this.div(rhs);
    }

    /**
    Power operator.

    Params:
        rhs = right hand side.
    Returns:
        power expression.
    */
    final const(Power!R) opBinary(string op)(const(Differentiable!R) rhs) const nothrow pure if (op == "^^")
        in (rhs)
        out (r; r)
    {
        return this.pow(rhs);
    }

    /**
    Binary operator.

    Params:
        rhs = right hand side.
    Returns:
        binary operator expression.
    */
    final auto opBinary(string op)(auto scope ref const(R) rhs) const nothrow pure
        if (op == "+" || op == "-" || op == "*" || op == "/" || op == "^^")
    {
        return opBinary!op(rhs.constant);
    }

    /**
    Binary operator.

    Params:
        lhs = left hand side.
    Returns:
        binary operator expression.
    */
    final auto opBinaryRight(string op)(auto scope ref const(R) lhs) const nothrow pure
        if (op == "+" || op == "-" || op == "*" || op == "/" || op == "^^")
    {
        return lhs.constant.opBinary!op(this);
    }
}

/**
Differentiate context.

Params:
    R = result type.
*/
final class DiffContext(R)
{
    this(const(Differentiable!R) target) nothrow pure scope return
        in (target)
    {
        this.target_ = target;
        this.zero_ = .zero!R();
        this.one_ = .one!R();
        this.two_ = .two!R();
    }

    @property const @nogc nothrow pure @safe scope
    {
        const(Differentiable!R) target() { return target_; }
        const(Differentiable!R) zero() { return zero_; }
        const(Differentiable!R) one() { return one_; }
        const(Differentiable!R) two() { return two_; }
    }

    const(Differentiable!R) diff(const(Differentiable!R) f) nothrow pure @safe scope
        in (f)
    {
        return assumeWontThrow(memo_.require(f, f.differentiate(this).rebindable));
    }

    bool isZero(scope const(Differentiable!R) value) @nogc nothrow pure @safe scope
    {
        return value is zero_;
    }

    bool isOne(scope const(Differentiable!R) value) @nogc nothrow pure @safe scope
    {
        return value is one_;
    }

    bool isTwo(scope const(Differentiable!R) value) @nogc nothrow pure @safe scope
    {
        return value is two_;
    }

private:
    const(Differentiable!R) target_;
    const(Differentiable!R) zero_;
    const(Differentiable!R) one_;
    const(Differentiable!R) two_;
    Rebindable!(const(Differentiable!R))[const(Differentiable!R)] memo_;
}

/**
Create differentiate context.

Params:
    R = result type.
    target = differentiate target.
Returns:
    Differentiate context.
*/
DiffContext!R diffContext(R)(const(Differentiable!R) target) nothrow pure
    in (target)
{
    return new DiffContext!R(target);
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.parameter : param;

    auto p = param(1.0);
    auto context = diffContext(p);
    assert(context.target is p);
    assert(context.zero()().isClose(0.0));
    assert(context.one()().isClose(1.0));
    assert(context.two()().isClose(2.0));

    assert(context.isZero(context.zero));
    assert(!context.isZero(context.one));
    assert(context.isOne(context.one));
    assert(!context.isOne(context.zero));
    assert(context.isTwo(context.two));
    assert(!context.isTwo(context.one));
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.parameter : param;

    auto p = param(1.0);
    auto context = diffContext(p);
    auto f = p * p;
    auto df = context.diff(f);
    assert(df().isClose(2.0));

    // cached for same function.
    assert(context.diff(f) is df);
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.parameter : param;

    auto p1 = param(1.0);
    auto p2 = param(2.0);
    assert((p1 + p2)().isClose(3.0));
    assert((p1 - p2)().isClose(-1.0));
    assert((p1 * p2)().isClose(2.0));
    assert((p1 / p2)().isClose(0.5));
    assert((p1 ^^ p2)().isClose(1.0));
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.parameter : param;

    auto p1 = param(1.0);
    assert((p1 + 2.0)().isClose(3.0));
    assert((p1 - 2.0)().isClose(-1.0));
    assert((p1 * 2.0)().isClose(2.0));
    assert((p1 / 2.0)().isClose(0.5));
    assert((p1 ^^ 2.0)().isClose(1.0));
}

nothrow pure unittest
{
    import std.math : isClose;
    import karasunum.differential.parameter : param;

    auto p1 = param(2.0);
    assert((1.0 + p1)().isClose(3.0));
    assert((1.0 - p1)().isClose(-1.0));
    assert((1.0 * p1)().isClose(2.0));
    assert((1.0 / p1)().isClose(0.5));
    assert((1.0 ^^ p1)().isClose(1.0));
}

/**
Evaluate context.

Params:
    R = result type.
*/
final class EvalContext(R)
{
    R evaluate(const(Differentiable!R) f) nothrow pure @safe scope
        in (f)
    {
        ++callCount_;
        return assumeWontThrow(memo_.require(
            f, { ++evaluateCount_; return f.evaluate(this); }()));
    }

    @property @nogc nothrow pure @safe const scope
    {
        size_t callCount() { return callCount_; }
        size_t evaluateCount() { return evaluateCount_; }
        size_t cacheHitCount() { return callCount_ - evaluateCount_; }
    }

private:
    R[const(Differentiable!R)] memo_;
    size_t callCount_;
    size_t evaluateCount_;
}

EvalContext!R evalContext(R)() nothrow pure @safe
    out(r; r)
{
    return new EvalContext!R();
}

unittest
{
    import std.conv : to;
    import std.math : isClose;
    import karasunum.differential.parameter : param;

    auto p1 = param(1.0);
    auto p2 = param(2.0);
    auto add = p1 + p2;
    auto context = evalContext!double();
    assert(context.evaluate(add).isClose(3.0));
    assert(context.callCount == 3);
    assert(context.evaluateCount == 3);
    assert(context.cacheHitCount == 0);

    assert(context.evaluate(add).isClose(3.0));
    assert(context.callCount == 4);
    assert(context.evaluateCount == 3);
    assert(context.cacheHitCount == 1);

    assert(context.evaluate(add).isClose(3.0));
    assert(context.callCount == 5);
    assert(context.evaluateCount == 3);
    assert(context.cacheHitCount == 2);
}

