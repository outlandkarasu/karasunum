/**
Newton's method implementation module.
*/
module karasunum.optimize.newton;

import karasunum.differential.differentiable :
    Differentiable,
    DiffContext,
    EvalContext;
import karasunum.differential.parameter : Parameter;

import karasunum.linear_algebra : Matrix, luDecomposition, Vector;

@safe:

/**
Newton's method implementation.

Params:
    R = element type.
    dimensions = parameter and result count.
*/
final class NewtonMethod(R, size_t dimensions)
{
    /**
    Initialize by functions and parameters.

    Params:
        functions = target functions.
        parameters = target parameters.
    */
    this(const(Differentiable!R)[] functions, Parameter!(R)[] parameters) nothrow pure scope
        in (functions.length == dimensions)
        in (parameters.length == dimensions)
    {
        this.functions_[] = functions[];
        this.parameters_[] = parameters[];

        foreach (i, f; functions)
        {
            foreach (j, p; parameters)
            {
                scope context = new DiffContext!R(p);
                this.jacobian_ ~= context.diff(f);
            }
        }
    }

    /**
    Update newton step.
    */
    void step() nothrow pure scope
    {
        scope context = new EvalContext!R();

        Vector!(dimensions, R) currentResult;
        foreach (i, f; functions_)
        {
            currentResult[i] = context.evaluate(f);
        }

        Matrix!(dimensions, dimensions, R) jacobian;
        foreach (i, f; jacobian_)
        {
            jacobian[i / dimensions, i % dimensions] = context.evaluate(f);
        }

        Vector!(dimensions, R) offset;
        luDecomposition(jacobian).solve(currentResult, offset);

        foreach (i; 0 .. dimensions)
        {
            auto p = parameters_[i];
            p.bind(p.value - offset[i]);
        }
    }

private:
    const(Differentiable!R)[dimensions] functions_;
    Parameter!R[dimensions] parameters_;
    const(Differentiable!R)[] jacobian_;
}

///
nothrow pure @safe unittest
{
    import std.math : isClose;
    import karasunum.differential.parameter : param;

    auto p1 = param(1.0);
    auto p2 = param(1.0);
    auto f1 = p1 ^^ 2.0;
    auto f2 = p2 ^^ 2.0;

    scope newton = new NewtonMethod!(double, 2)([f1, f2], [p1, p2]);
    foreach (i; 0 .. 20)
    {
        newton.step();
    }
    assert(p1.value.isClose(0.0, 1e-6, 1e-6));
    assert(p2.value.isClose(0.0, 1e-6, 1e-6));
}

