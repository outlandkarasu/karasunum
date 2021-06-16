/**
Matrix LU decomposition module.
*/
module karasunum.linear_algebra.lu;

import std.algorithm : swap;
import std.traits : isNumeric;

import karasunum.linear_algebra.matrix : Matrix;
import karasunum.linear_algebra.vector : Vector;

@safe:

/**
LU decomposition matrices.
*/
struct LUDecomposition(size_t N, E) if (isNumeric!E)
{
    /**
    Matrix type.
    */
    alias Mat = Matrix!(N, N, E);

    /**
    Vector type.
    */
    alias Vec = Vector!(N, E);

    @disable this();

    /**
    Initialize by matrix.

    Params:
        m = target matrix.
    */
    this()(auto scope ref const(Mat) m) @nogc nothrow pure scope
    {
        m.createPivots(this.pivot_);

        Mat swapped;
        foreach (i, pos; pivot_)
        {
            foreach (j; 0 .. N)
            {
                swapped[pos, j] = m[i, j];
            }
        }

        swapped.createLUDecomposition(this.l_, this.u_);
    }

    /**
    Create inverse matrix.

    Params:
        dest = destination matrix.
    Returns:
        destination matrix reference.
    */
    ref Mat createInverse(return scope ref Mat dest) const @nogc nothrow pure scope
    {
        Mat inverseL;
        l_.inverseLMatrix(inverseL);

        Mat inverseU;
        u_.inverseUMatrix(inverseU);

        Mat inverseLUP;
        inverseLUP.mul(inverseU, inverseL);

        foreach (i; 0 .. N)
        {
            foreach (j; 0 .. N)
            {
                dest[i, j] = inverseLUP[i, pivot_[j]];
            }
        }

        return dest;
    }

    /**
    Create inverse matrix.

    Returns:
        destination matrix reference.
    */
    Mat createInverse() const @nogc nothrow pure scope
    {
        Mat dest;
        return createInverse(dest);
    }

    /**
    Solve equation.

    Params:
        y = equation Y.
        x = solved X.
    Returns:
        reference to X.
    */
    ref Vec solve()(auto scope ref const(Vec) y, return scope ref Vec x) const @nogc nothrow pure scope
    {
        Vec swapped;
        foreach (i, pos; pivot_)
        {
            swapped[i] = y[pos];
        }

        Vec a;
        l_.solveByLMatrix(swapped, a);
        u_.solveByUMatrix(a, x);
        return x;
    }

    /**
    Solve equation.

    Params:
        y = equation Y.
    Returns:
        solved X.
    */
    Vec solve()(auto scope ref const(Vec) y) const @nogc nothrow pure scope
    {
        Vec x;
        return solve(y, x);
    }

private:
    Mat l_;
    Mat u_;
    size_t[N] pivot_;
}

/**
Construct LU decomposition matrices.

Params:
    N = dimensions.
    E = element type.
    m = target matrix.
Returns:
    LU decomposition matrices.
*/
auto luDecomposition(size_t N, E)(auto scope ref const(Matrix!(N, N, E)) m) @nogc nothrow pure
    if (isNumeric!E)
{
    return LUDecomposition!(N, E)(m);
}

///
@nogc nothrow pure @safe unittest
{
    import karasunum.linear_algebra.matrix : isClose;
    enum N = 4;

    immutable m = Matrix!(N, N).fromRows([
        [5.0, 6.0, 7.0, 8.0],
        [10.0, 21.0, 24.0, 27.0],
        [15.0, 54.0, 73.0, 81.0],
        [25.0, 84.0, 179.0, 211.0],
    ]);
    immutable lu = m.luDecomposition();

    immutable inverse = lu.createInverse();
    auto result = typeof(lu).Mat();
    result.mul(inverse, m);
    assert(result.isUnitMatrix);
}

///
@nogc nothrow pure @safe unittest
{
    import karasunum.linear_algebra.vector : isClose;
    enum N = 4;

    immutable m = Matrix!(N, N).fromRows([
        [5.0, 6.0, 7.0, 8.0],
        [10.0, 21.0, 24.0, 27.0],
        [15.0, 54.0, 73.0, 81.0],
        [25.0, 84.0, 179.0, 211.0],
    ]);
    immutable lu = m.luDecomposition();

    immutable y = typeof(lu).Vec([3.0, 5.0, 4.0, -5.0]);
    immutable x = lu.solve(y);

    typeof(lu).Vec result;
    result.mul(m, x);
    assert(result.isClose(y));
}

private:

version(unittest)
{
    bool isUnitMatrix(E, size_t N)(auto scope ref const(Matrix!(N, N, E)) m) @nogc nothrow pure @safe
    {
        import std.math : isClose;

        foreach (i; 0 .. N)
        {
            foreach (j; 0 .. N)
            {
                if (!m[i, j].isClose((i == j) ? 1.0 : 0.0, 1e-4, 1e-4))
                {
                    return false;
                }
            }
        }

        return true;
    }
}

/**
Create matrix pivots.

Params:
    m = target matrix.
    pivots = pivots indicies.
*/
void createPivots(size_t N, E)(
    auto scope ref const(Matrix!(N, N, E)) m,
    scope ref size_t[N] pivots)
    if (isNumeric!E)
{
    foreach (i, ref e; pivots)
    {
        e = i;
    }

    foreach (i; 0 .. N)
    {
        size_t maxPos = i;
        E maxValue = m[maxPos, i];
        foreach (j; (i + 1) .. N)
        {
            immutable value = m[pivots[j], i];
            if (maxValue < value)
            {
                maxPos = j;
                maxValue = value;
            }
        }

        swap(pivots[i], pivots[maxPos]);
    }
}

///
@nogc nothrow pure @safe unittest
{
    enum N = 2;
    alias Mat = Matrix!(N, N);

    size_t[N] pivots;
    auto m = Mat.fromRows([
        [1.0f, 2.0f],
        [2.0f, 1.0f]]);

    m.createPivots(pivots);
    assert(pivots == [1, 0]);

    m = Mat.fromRows([
        [1.0f, 1.0f],
        [2.0f, 3.0f]]);

    m.createPivots(pivots);
    assert(pivots == [1, 0]);
}

///
@nogc nothrow pure @safe unittest
{
    enum N = 3;
    alias Mat = Matrix!(N, N);

    size_t[N] pivots;
    auto m = Mat.fromRows([
        [1.0f, 2.0f, 3.0f],
        [2.0f, 1.0f, 4.0f],
        [3.0f, 4.0f, 0.0f],
    ]);

    m.createPivots(pivots);
    assert(pivots == [2, 0, 1]);
}

/**
LU decomposition.

Params:
    m = target matrix.
    l = L destination matrix.
    u = U destination matrix.
*/
void createLUDecomposition(size_t N, E)(
    auto scope ref const(Matrix!(N, N, E)) m,
    scope ref Matrix!(N, N, E) l,
    scope ref Matrix!(N, N, E) u)
    if (isNumeric!E)
{
    // initialize L first row.
    l[0, 0] = E(1);
    foreach (column; 1 .. N)
    {
        l[0, column] = E(0);
    }

    // initialize U first row.
    foreach (column; 0 .. N)
    {
        u[0, column] = m[0, column];
    }

    auto m00 = m[0, 0];
    foreach (row; 1 .. N)
    {
        // clear L fixed cells.
        l[row, row] = E(1);
        foreach (column; row + 1 .. N)
        {
            l[row, column] = E(0);
        }

        // clear R fiexd cells.
        foreach (column; 0 .. row)
        {
            u[row, column] = E(0);
        }

        // setting up L first column.
        l[row, 0] = m[row, 0] / m00;

        // calculate L columns.
        foreach (column; 1 .. row)
        {
            E sum = E(0);
            foreach (lc; 0 .. column)
            {
                sum += l[row, lc] * u[lc, column];
            }
            l[row, column] = (m[row, column] - sum) / u[column, column];
        }

        // calculate U columns.
        foreach (column; row .. N)
        {
            E sum = E(0);
            foreach (lc; 0 .. row)
            {
                sum += l[row, lc] * u[lc, column];
            }
            u[row, column] = m[row, column] - sum;
        }
    }
}

///
@nogc nothrow pure @safe unittest
{
    import karasunum.linear_algebra.matrix : isClose;

    immutable m = Matrix!(2, 2).fromRows([
        [4.0f, 3.0f],
        [6.0f, 3.0f]]);
    auto l = Matrix!(2, 2)();
    auto u = Matrix!(2, 2)();
    m.createLUDecomposition(l, u);

    assert(l.isClose(Matrix!(2, 2).fromRows([
        [1.0, 0.0],
        [1.5, 1.0]
    ])));
    assert(u.isClose(Matrix!(2, 2).fromRows([
        [4.0, 3.0],
        [0.0, -1.5]
    ])));

    auto mul = Matrix!(2, 2)();
    mul.mul(l, u);
    assert(mul.isClose(m));
}

@nogc nothrow pure @safe unittest
{
    import karasunum.linear_algebra.matrix : isClose;

    immutable m = Matrix!(3, 3).fromRows([
        [5.0, 6.0, 7.0],
        [10.0, 20.0, 23.0],
        [15.0, 50.0, 67.0],
    ]);
    auto l = Matrix!(3, 3)();
    auto u = Matrix!(3, 3)();
    m.createLUDecomposition(l, u);
    assert(l.isClose(Matrix!(3, 3).fromRows([
        [1.0,  0.0,  0.0],
        [2.0,  1.0,  0.0],
        [3.0,  4.0,  1.0],
    ])));
    assert(u.isClose(Matrix!(3, 3).fromRows([
        [5.0, 6.0,  7.0],
        [0.0, 8.0,  9.0],
        [0.0, 0.0, 10.0],
    ])));

    auto mul = Matrix!(3, 3)();
    mul.mul(l, u);
    assert(mul.isClose(m));
}

@nogc nothrow pure @safe unittest
{
    import karasunum.linear_algebra.matrix : isClose;

    immutable m = Matrix!(4, 4).fromRows([
        [5.0, 6.0, 7.0, 8.0],
        [10.0, 21.0, 24.0, 27.0],
        [15.0, 54.0, 73.0, 81.0],
        [25.0, 84.0, 179.0, 211.0],
    ]);
    auto l = Matrix!(4, 4)();
    auto u = Matrix!(4, 4)();
    m.createLUDecomposition(l, u);
    assert(l.isClose(Matrix!(4, 4).fromRows([
        [1.0,  0.0,  0.0, 0.0],
        [2.0,  1.0,  0.0, 0.0],
        [3.0,  4.0,  1.0, 0.0],
        [5.0,  6.0,  7.0, 1.0],
    ])));
    assert(u.isClose(Matrix!(4, 4).fromRows([
        [5.0, 6.0,  7.0, 8.0],
        [0.0, 9.0, 10.0, 11.0],
        [0.0, 0.0, 12.0, 13.0],
        [0.0, 0.0,  0.0, 14.0],
    ])));

    auto mul = Matrix!(4, 4)();
    mul.mul(l, u);
    assert(mul.isClose(m));
}

/**
Inverse lower triangle matrix.

Params:
    N = matrix dimensions.
    E = matrix element.
    l = lower triangle matrix.
    inverse = inverse matrix.
*/
void inverseLMatrix(size_t N, E)(
    scope ref const(Matrix!(N, N, E)) l,
    scope ref Matrix!(N, N, E) inverse)
    if (isNumeric!E)
{
    // for each inverse element.
    foreach (i; 0 .. N)
    {
        foreach (j; 0 .. i)
        {
            auto sum = E(0);
            foreach (k; j .. i)
            {
                sum += l[i, k] * inverse[k, j];
            }

            inverse[i, j] = -(sum / l[i, i]);
        }

        // diagonal element.
        inverse[i, i] = E(1) / l[i, i];

        // fill uppder triangle elements to 0.
        foreach (j; (i + 1) .. N)
        {
            inverse[i, j] = E(0);
        }
    }
}

@nogc nothrow pure @safe unittest
{
    enum N = 1;
    alias Mat = Matrix!(N, N);

    immutable m = Mat.fromRows([[5.0]]);
    auto inverse = Mat();
    m.inverseLMatrix(inverse);

    auto result = Mat();
    result.mul(inverse, m);
    assert(result.isUnitMatrix);
}

@nogc nothrow pure @safe unittest
{
    enum N = 2;
    alias Mat = Matrix!(N, N);

    immutable m = Mat.fromRows([
        [5.0, 0.0],
        [6.0, 7.0],
    ]);
    auto inverse = Mat();
    m.inverseLMatrix(inverse);

    auto result = Mat();
    result.mul(inverse, m);
    assert(result.isUnitMatrix);
}

@nogc nothrow pure @safe unittest
{
    enum N = 3;
    alias Mat = Matrix!(N, N);

    immutable m = Mat.fromRows([
        [5.0, 0.0,  0.0],
        [6.0, 7.0,  0.0],
        [8.0, 9.0, 10.0],
    ]);
    auto inverse = Mat();
    m.inverseLMatrix(inverse);

    auto result = Mat();
    result.mul(inverse, m);
    assert(result.isUnitMatrix);
}

@nogc nothrow pure @safe unittest
{
    enum N = 4;
    alias Mat = Matrix!(N, N);

    immutable m = Mat.fromRows([
        [5.0,  0, 0, 0],
        [10.0, 21.0, 0, 0],
        [15.0, 54.0, 73.0, 0],
        [25.0, 84.0, 179.0, 211.0],
    ]);
    auto inverse = Mat();
    m.inverseLMatrix(inverse);

    auto result = Mat();
    result.mul(inverse, m);
    assert(result.isUnitMatrix);
}

/**
Inverse upper triangle matrix.

Params:
    N = matrix dimensions.
    E = matrix element.
    l = upper triangle matrix.
    inverse = inverse matrix.
*/
void inverseUMatrix(size_t N, E)(
    scope ref const(Matrix!(N, N, E)) u,
    scope ref Matrix!(N, N, E) inverse)
    if (isNumeric!E)
{
    // for each inverse element in reverse order.
    foreach_reverse (i; 0 .. N)
    {
        // fill lower triangle elements to 0.
        foreach (j; 0 .. i)
        {
            inverse[i, j] = E(0);
        }

        // diagonal element.
        inverse[i, i] = E(1) / u[i, i];

        foreach (j; (i + 1) .. N)
        {
            auto sum = E(0);
            foreach_reverse (k; (i + 1) .. (j + 1))
            {
                sum += u[i, k] * inverse[k, j];
            }

            inverse[i, j] = -(sum / u[i, i]);
        }
    }
}

@nogc nothrow pure @safe unittest
{
    enum N = 1;
    alias Mat = Matrix!(N, N);

    immutable m = Mat.fromRows([[5.0]]);
    auto inverse = Mat();
    m.inverseUMatrix(inverse);

    auto result = Mat();
    result.mul(inverse, m);
    assert(result.isUnitMatrix);
}

@nogc nothrow pure @safe unittest
{
    enum N = 2;
    alias Mat = Matrix!(N, N);

    immutable m = Mat.fromRows([
        [5.0, 6.0],
        [0.0, 7.0],
    ]);
    auto inverse = Mat();
    m.inverseUMatrix(inverse);

    auto result = Mat();
    result.mul(inverse, m);
    assert(result.isUnitMatrix);
}

@nogc nothrow pure @safe unittest
{
    enum N = 3;
    alias Mat = Matrix!(N, N);

    immutable m = Mat.fromRows([
        [5.0, 4.0,  3.0],
        [0.0, 7.0,  8.0],
        [0.0, 0.0, 10.0],
    ]);
    auto inverse = Mat();
    m.inverseUMatrix(inverse);

    auto result = Mat();
    result.mul(inverse, m);

    assert(result.isUnitMatrix);
}

@nogc nothrow pure @safe unittest
{
    enum N = 4;
    alias Mat = Matrix!(N, N);

    immutable m = Mat.fromRows([
        [2.0, 5.0, 4.0,  3.0],
        [0.0, 7.0, 8.0,  1.0],
        [0.0, 0.0, 10.0, 5.0],
        [0.0, 0.0, 0.0,  5.0],
    ]);
    auto inverse = Mat();
    m.inverseUMatrix(inverse);

    auto result = Mat();
    result.mul(inverse, m);
    assert(result.isUnitMatrix);
}

/**
Solve equation by L matrix.
*/
void solveByLMatrix(size_t N, E)(
    auto scope ref const(Matrix!(N, N, E)) l,
    auto scope ref const(Vector!(N, E)) y,
    ref Vector!(N, E) x)
{
    foreach (i; 0 .. N)
    {
        E sum = E(0);
        foreach (j; 0 .. i)
        {
            sum += l[i, j] * x[j];
        }

        x[i] = (y[i] - sum) / l[i, i];
    }
}

@nogc nothrow pure @safe unittest
{
    import karasunum.linear_algebra.vector : isClose;

    enum N = 1;
    alias Mat = Matrix!(N, N, double);
    alias Vec = Vector!(N, double);

    immutable m = Mat.fromRows([
        [3.0],
    ]);
    immutable y = Vec([12.0]);
    Vec x;
    m.solveByLMatrix(y, x);

    assert(x.isClose(Vec([4.0])));
}

@nogc nothrow pure @safe unittest
{
    import karasunum.linear_algebra.vector : isClose;

    enum N = 2;
    alias Mat = Matrix!(N, N, double);
    alias Vec = Vector!(N, double);

    immutable m = Mat.fromRows([
        [3.0, 0.0],
        [4.0, 5.0],
    ]);
    immutable y = Vec([9, 3]);
    Vec x;
    m.solveByLMatrix(y, x);

    assert(x.isClose(Vec([3.0, -9.0 / 5.0])));
}

/**
Solve equation by U matrix.
*/
void solveByUMatrix(size_t N, E)(
    auto scope ref const(Matrix!(N, N, E)) l,
    auto scope ref const(Vector!(N, E)) y,
    ref Vector!(N, E) x)
{
    foreach_reverse (i; 0 .. N)
    {
        E sum = E(0);
        foreach_reverse (j; (i + 1) .. N)
        {
            sum += l[i, j] * x[j];
        }

        x[i] = (y[i] - sum) / l[i, i];
    }
}

@nogc nothrow pure @safe unittest
{
    import karasunum.linear_algebra.vector : isClose;

    enum N = 2;
    alias Mat = Matrix!(N, N, double);
    alias Vec = Vector!(N, double);

    immutable m = Mat.fromRows([
        [3.0, 4.0],
        [0.0, 5.0],
    ]);
    immutable y = Vec([9, 3]);
    Vec x;
    m.solveByUMatrix(y, x);

    assert(x.isClose(Vec([2.2, 0.6])));
}

