/**
Matrix module.
*/
module karasunum.linear_algebra.matrix;

import std.algorithm : swap;
import std.math : cos, sin;
import std.traits : isNumeric;

import karasunum.linear_algebra.vector : Vector;

@safe:

/**
Matrix structure.

Params:
    ROWS = matrix rows.
    COLS = matrix columns.
    E = element type.
*/
struct Matrix(size_t ROWS, size_t COLS, E = float)
{
    static assert(ROWS > 0);
    static assert(COLS > 0);
    static assert(isNumeric!E);

    /**
    Initialize by row major elements.

    Params:
        elements = matrix row major elements.
    Returns:
        initialized matrix.
    */
    static typeof(this) fromRows(scope const(E)[COLS][ROWS] elements)
    {
        auto m = typeof(this)();
        foreach (j; 0 .. COLS)
        {
            foreach (i; 0 .. ROWS)
            {
                m.elements_[j][i] = elements[i][j];
            }
        }
        return m;
    }

    static if(COLS == ROWS)
    {
        /**
        Initialize unit matrix.

        Returns:
            unit matrix;
        */
        static typeof(this) unit()
        {
            auto m = typeof(this)();
            foreach (j; 0 .. COLS)
            {
                foreach (i; 0 .. ROWS)
                {
                    m.elements_[j][i] = cast(E)((i == j) ? 1 : 0);
                }
            }
            return m;
        }

        /**
        Create scale matrix.

        Params:
            factors = scale factors.
        Returns:
            scale matrix.
        */
        static typeof(this) scale(Factors...)(Factors factors)
        {
            static assert(factors.length == COLS - 1);
            auto m = typeof(this)();
            m.fill(cast(E) 0);
            foreach (i, f; factors)
            {
                m[i, i] = cast(E) f;
            }
            m[ROWS - 1, COLS - 1] = cast(E) 1.0;
            return m;
        }

        /**
        Create translate matrix.

        Params:
            factors = translate factors.
        Returns:
            translate matrix.
        */
        static typeof(this) translate(Factors...)(Factors factors)
        {
            static assert(factors.length == COLS - 1);

            auto m = typeof(this).unit();
            foreach (i, f; factors)
            {
                m[i, COLS - 1] = cast(E) f;
            }
            return m;
        }

        static if (ROWS == 4 && COLS == 4)
        {
            /**
            Create rotation X matrix.

            Params:
                theta = rotate theta.
            Returns:
                rotation matrix.
            */
            static typeof(this) rotateX(E theta)
            {
                auto m = typeof(this).unit();
                m[1, 1] = cos(theta);
                m[1, 2] = -sin(theta);
                m[2, 1] = sin(theta);
                m[2, 2] = cos(theta);
                return m;
            }

            /**
            Create rotation Y matrix.

            Params:
                theta = rotate theta.
            Returns:
                rotation matrix.
            */
            static typeof(this) rotateY(E theta)
            {
                auto m = typeof(this).unit();
                m[0, 0] = cos(theta);
                m[0, 2] = sin(theta);
                m[2, 0] = -sin(theta);
                m[2, 2] = cos(theta);
                return m;
            }

            /**
            Create rotation Z matrix.

            Params:
                theta = rotate theta.
            Returns:
                rotation matrix.
            */
            static typeof(this) rotateZ(E theta)
            {
                auto m = typeof(this).unit();
                m[0, 0] = cos(theta);
                m[0, 1] = -sin(theta);
                m[1, 0] = sin(theta);
                m[1, 1] = cos(theta);
                return m;
            }
        }
    }

    @property const scope
    {
        size_t rows() { return ROWS; }
        size_t columns() { return COLS; }
    }

    /**
    Get an element.

    Params:
        i = row index.
        j = column index.
    Returns:
        element value.
    */
    ref const(E) opIndex(size_t i, size_t j) const return scope
    in (i < ROWS)
    in (j < COLS)
    {
        return elements_[j][i];
    }

    /**
    Set an element.

    Params:
        value = element value.
        i = row index.
        j = column index.
    Returns:
        assigned element value.
    */
    ref const(E) opIndexAssign()(auto ref const(E) value, size_t i, size_t j) return scope
    in (i < ROWS)
    in (j < COLS)
    {
        return elements_[j][i] = value;
    }

    /**
    operation and assign an element.

    Params:
        op = operator.
        value = element value.
        i = row index.
        j = column index.
    Returns:
        assigned element value.
    */
    ref const(E) opIndexOpAssign(string op)(auto ref const(E) value, size_t i, size_t j) return scope
    in (i < ROWS)
    in (j < COLS)
    {
        return mixin("elements_[j][i] " ~ op ~ "= value");
    }

    /**
    Operation and assign other vector.

    Params:
        value = other vetor value.
    Returns:
        this vector.
    */
    ref typeof(this) opOpAssign(string op)(auto ref const(typeof(this)) value) return scope
    {
        foreach (j, ref column; elements_)
        {
            foreach (i, ref v; column)
            {
                mixin("v " ~ op ~ "= value[i, j];");
            }
        }
        return this;
    }

    /**
    Matrix multiplication.

    Params:
        lhs = left hand side matrix.
        rhs = right hand side matrix.
    Returns:
        calculated this matrix.
    */
    ref typeof(this) mul(size_t N, E1, E2)(
            auto ref const(Matrix!(ROWS, N, E1)) lhs,
            auto ref const(Matrix!(N, COLS, E2)) rhs) return scope
    {
        foreach (j, ref column; elements_)
        {
            foreach (i, ref v; column)
            {
                v = cast(E) 0;
                foreach (k; 0 .. N)
                {
                    v += lhs[i, k] * rhs[k, j];
                }
            }
        }
        return this;
    }

    /**
    Fill elements.

    Params:
        value = filler value.
    */
    ref typeof(this) fill()(auto ref const(E) value) return scope
    {
        foreach (ref column; elements_)
        {
            column[] = value;
        }
        return this;
    }

    /**
    Matrix pointer.

    Returns:
        matrix pointer.
    */
    @property const(E)* ptr() const return scope
    out (r; r != null)
    {
        return &elements_[0][0];
    }

    /**
    Swap 2 rows.

    Params:
        row1 = swap target row 1.
        row2 = swap target row 2.
    */
    void swapRows(size_t row1, size_t row2) scope
        in (row1 < ROWS)
        in (row2 < ROWS)
    {
        if (row1 == row2)
        {
            return;
        }

        foreach (column; 0 .. COLS)
        {
            swap(elements_[column][row1], elements_[column][row2]);
        }
    }

private:
    E[ROWS][COLS] elements_;
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;

    immutable m = Matrix!(2, 3).fromRows([
        [1, 2, 3],
        [4, 5, 6],
    ]);
    assert(m.rows == 2);
    assert(m.columns == 3);

    assert(m[0, 0].isClose(1));
    assert(m[0, 1].isClose(2));
    assert(m[0, 2].isClose(3));
    assert(m[1, 0].isClose(4));
    assert(m[1, 1].isClose(5));
    assert(m[1, 2].isClose(6));
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;

    immutable m = Matrix!(3, 3).unit;
    assert(m.rows == 3);
    assert(m.columns == 3);

    assert(m[0, 0].isClose(1));
    assert(m[0, 1].isClose(0));
    assert(m[0, 2].isClose(0));
    assert(m[1, 0].isClose(0));
    assert(m[1, 1].isClose(1));
    assert(m[1, 2].isClose(0));
    assert(m[2, 0].isClose(0));
    assert(m[2, 1].isClose(0));
    assert(m[2, 2].isClose(1));
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;

    auto m = Matrix!(2, 2).fromRows([
        [1, 2],
        [3, 4]
    ]);
    m[0, 0] = 3.0f;
    m[0, 1] = 4.0f;
    m[1, 0] = 5.0f;
    m[1, 1] = 6.0f;

    assert(m[0, 0].isClose(3));
    assert(m[0, 1].isClose(4));
    assert(m[1, 0].isClose(5));
    assert(m[1, 1].isClose(6));
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;

    auto m = Matrix!(2, 2).fromRows([
        [1, 2],
        [3, 4]
    ]);
    m[0, 0] += 1.0f;
    m[0, 1] += 1.0f;
    m[1, 0] += 1.0f;
    m[1, 1] += 1.0f;

    assert(m[0, 0].isClose(2));
    assert(m[0, 1].isClose(3));
    assert(m[1, 0].isClose(4));
    assert(m[1, 1].isClose(5));
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;

    auto m = Matrix!(2, 2).fromRows([
        [1, 2],
        [3, 4]
    ]);
    immutable t = Matrix!(2, 2).fromRows([
        [3, 4],
        [5, 6]
    ]);

    m += t;

    assert(m[0, 0].isClose(4));
    assert(m[0, 1].isClose(6));
    assert(m[1, 0].isClose(8));
    assert(m[1, 1].isClose(10));
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;

    auto result = Matrix!(2, 2)();
    immutable lhs = Matrix!(2, 3).fromRows([
        [3, 4, 5],
        [6, 7, 8],
    ]);
    immutable rhs = Matrix!(3, 2).fromRows([
        [3, 4],
        [6, 7],
        [8, 9],
    ]);

    result.mul(lhs, rhs);

    assert(result[0, 0].isClose(3 * 3 + 4 * 6 + 5 * 8));
    assert(result[0, 1].isClose(3 * 4 + 4 * 7 + 5 * 9));
    assert(result[1, 0].isClose(6 * 3 + 7 * 6 + 8 * 8));
    assert(result[1, 1].isClose(6 * 4 + 7 * 7 + 8 * 9));
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;

    auto m = Matrix!(2, 2)();
    m.fill(1.0);

    assert(m[0, 0].isClose(1.0));
    assert(m[0, 1].isClose(1.0));
    assert(m[1, 0].isClose(1.0));
    assert(m[1, 1].isClose(1.0));
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;

    immutable m = Matrix!(2, 2)([[1, 2], [3, 4]]);
    assert(isClose(*(m.ptr), 1.0));
}

///
@nogc nothrow pure unittest
{
    import std.math : isClose;

    immutable m = Matrix!(4, 4).scale(2.0, 3.0, 4.0);
    assert(m[0, 0].isClose(2.0));
    assert(m[1, 1].isClose(3.0));
    assert(m[2, 2].isClose(4.0));
    assert(m[3, 3].isClose(1.0));

    foreach (i; 0 .. 4)
    {
        foreach (j; 0 .. 4)
        {
            if (i != j)
            {
                assert(m[i, j].isClose(0.0));
            }
        }
    }
}

///
@nogc nothrow pure unittest
{
    import karasunum.linear_algebra.vector : isClose;

    immutable m = Matrix!(4, 4).translate(2.0, 3.0, 4.0);
    immutable v = Vector!4([1.0, 2.0, 3.0, 1.0]);
    auto result = Vector!4();
    result.mul(m, v);

    assert(result.isClose(Vector!4([3, 5, 7, 1])));
}

///
@nogc nothrow pure unittest
{
    import karasunum.linear_algebra.vector : isClose;

    immutable m = Matrix!(4, 4).rotateX(0.5);
    immutable x = Vector!4([1.0, 0.0, 0.0, 1.0]);
    immutable y = Vector!4([0.0, 1.0, 0.0, 1.0]);
    immutable z = Vector!4([0.0, 0.0, 1.0, 1.0]);

    auto result = Vector!4();
    result.mul(m, x);
    assert(result.isClose(x));
    result.mul(m, y);
    assert(result.isClose(Vector!4([0.0, cos(0.5), sin(0.5), 1.0])));
    result.mul(m, z);
    assert(result.isClose(Vector!4([0.0, -sin(0.5), cos(0.5), 1.0])));
}

///
@nogc nothrow pure unittest
{
    import karasunum.linear_algebra.vector : isClose;

    immutable m = Matrix!(4, 4).rotateY(0.5);
    immutable x = Vector!4([1.0, 0.0, 0.0, 1.0]);
    immutable y = Vector!4([0.0, 1.0, 0.0, 1.0]);
    immutable z = Vector!4([0.0, 0.0, 1.0, 1.0]);

    auto result = Vector!4();
    result.mul(m, x);
    assert(result.isClose(Vector!4([cos(0.5), 0.0, -sin(0.5), 1.0])));
    result.mul(m, y);
    assert(result.isClose(y));
    result.mul(m, z);
    assert(result.isClose(Vector!4([sin(0.5), 0.0, cos(0.5), 1.0])));
}

///
@nogc nothrow pure unittest
{
    import karasunum.linear_algebra.vector : isClose;

    immutable m = Matrix!(4, 4).rotateZ(0.5);
    immutable x = Vector!4([1.0, 0.0, 0.0, 1.0]);
    immutable y = Vector!4([0.0, 1.0, 0.0, 1.0]);
    immutable z = Vector!4([0.0, 0.0, 1.0, 1.0]);

    auto result = Vector!4();
    result.mul(m, x);
    assert(result.isClose(Vector!4([cos(0.5), sin(0.5), 0.0, 1.0])));
    result.mul(m, y);
    assert(result.isClose(Vector!4([-sin(0.5), cos(0.5), 0.0, 1.0])));
    result.mul(m, z);
    assert(result.isClose(z));
}

/**
isClose for matrix.

Params:
    R = rows.
    C = columns.
    a = matrix.
    b = other matrix.
Returns:
    true if both matrix are close.
*/
bool isClose(size_t R, size_t C, E)(
    auto scope ref const(Matrix!(R, C, E)) a,
    auto scope ref const(Matrix!(R, C, E)) b) @nogc nothrow pure 
{
    import std.math : mathIsClose = isClose;

    foreach (c; 0 .. C)
    {
        foreach (r; 0 .. R)
        {
            if (!mathIsClose(a[r, c], b[r, c]))
            {
                return false;
            }
        }
    }

    return true;
}

///
@nogc nothrow pure @safe unittest
{
    immutable a = Matrix!(4, 4).fromRows([
        [1, 2, 3, 4],
        [5, 6, 7, 8],
        [9, 10, 11, 12],
        [13, 14, 15, 16],
    ]);

    assert(a.isClose(a));

    auto b = Matrix!(4, 4)();
    b = a;
    assert(a.isClose(b));
    assert(b.isClose(a));

    b[0, 0] = 100.0;
    assert(!a.isClose(b));
    assert(!b.isClose(a));
}

///
@nogc nothrow pure @safe unittest
{
    auto a = Matrix!(4, 4).fromRows([
        [1, 2, 3, 4],
        [5, 6, 7, 8],
        [9, 10, 11, 12],
        [13, 14, 15, 16],
    ]);

    immutable expected = Matrix!(4, 4).fromRows([
        [9, 10, 11, 12],
        [5, 6, 7, 8],
        [1, 2, 3, 4],
        [13, 14, 15, 16],
    ]);

    a.swapRows(0, 2);
    assert(a.isClose(expected));

    a.swapRows(3, 3);
    assert(a.isClose(expected));
}

