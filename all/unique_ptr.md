# `unique_ptr`
C++11 introduced the [`unique_ptr`][1] template in the `<memory>` header.  It
provides a convenient way of managing memory (and more generally, lifetimes of
arbitrary objects).  It's a bit odd to use, so this article provides a very
basic introduction to its semantics.

The article will not attempt explain its features in detail.  You can find
[more info in the API documentation][1].

(The astute reader will notice that I'm using the same ownership terminology
as [Rust][2].)

## Scalar vs array

There are two kinds of `unique_ptr`, one for scalars (i.e. non-arrays) and one
for arrays:

  - `unique_ptr<double>` can hold a scalar of type `double`;
  - `unique_ptr<double[]>` can hold an array of `double` values with an
    unknown number of elements.

Of course, one may substitute `double` with any arbitrary type.

## Construction

A `std::unique_ptr` variable can constructed using [`make_unique`][3]:

~~~c++
auto p1 = std::make_unique<double>(3.14);
auto p2 = std::make_unique<double[]>(n);
~~~

  - The first line creates a `unique_ptr` to a `double` value `3.14` and saves
    it in the `p1` variable.
  - The second line creates a `unique_ptr` to a `double` array of `n` elements
    and saves it to in the `p2` variable.

Note that `make_unique` requires C++14.  If you want to stick to plain C++11
you can substitute those with:

~~~c++
std::unique_ptr<double>   p1(new double(3.14));
std::unique_ptr<double[]> p2(new double[n]());
~~~

## Ownership and uniqueness

A `unique_ptr` variable is said to be the **owner** of the object it points
to: if the owner dies without passing its ownership to another variable, the
object is deleted (and hence the memory is deallocated), thereby invalidating
all pointers and references to it.

There can only be one owner at any one time: one does not simply **copy** a
`unique_ptr`.  For example:

~~~c++
// compile with: c++ -std=c++14
#include <memory>

void taker(std::unique_ptr<double[]>);

void uniqueness1()
{
    // p is created and owns a section of memory containing 42 numbers
    std::unique_ptr<double[]> p = std::make_unique<double[]>(42);

    // compile error! (attempted to copy p to q)
    std::unique_ptr<double[]> q = p;

    // compile error! (attempted to copy p to the argument of taker)
    taker(p);
}
~~~

Typical errors when this happens:

    (g++)
    error: use of deleted function ‘std::unique_ptr<…>::unique_ptr(const std::unique_ptr<…>&)
    error: use of deleted function ‘… std::unique_ptr<…>::operator=(const std::unique_ptr<…>&)

    (clang++)
    error: call to deleted constructor of 'std::unique_ptr<…>'
    error: overload resolution selected deleted operator '='

If an owner wants to **lend** the contents of `unique_ptr` to another function
without relinquishing its ownership, it can give a pointer/reference to it:

~~~c++
// compile with: c++ -std=c++14
#include <memory>

void borrower1(double*);

void borrower2(std::unique_ptr<double[]>*);

void borrower3(std::unique_ptr<double[]>&);

void uniqueness2()
{
    // p is created and owns a section of memory containing 42 numbers
    std::unique_ptr<double[]> p = std::make_unique<double[]>(42);

    // these are all okay but the first method is the most general
    // as it accepts any kind of pointer, not just unique_ptr
    borrower1(p.get());
    borrower2(&p);
    borrower3(p);

    // p dies, so the array is deleted automatically
}
~~~

If an owner wants to yield its ownership of a `unique_ptr` to another function
or variable, it can perform a [**`move`**][6]:

~~~c++
// compile with: c++ -std=c++14
#include <memory>
#include <utility> // required for std::move

void taker(std::unique_ptr<double[]>);

void uniqueness3()
{
    // p is created and owns a section of memory containing 42 numbers
    std::unique_ptr<double[]> p = std::make_unique<double[]>(42);

    // move p to q
    std::unique_ptr<double[]> q = std::move(p);

    // moved q into the argument of taker
    taker(std::move(q));

    // p and q both die, but they no longer own anything so nothing happens here
}
~~~

## Resetting

Although normally a `unique_ptr` will automatically delete its object once it
dies (or when it gains ownership of something else), one can hasten the
process by manually calling [`.reset()`][4]:

~~~c++
p.reset();
~~~

This will immediately delete the object pointed to by `p` and cause the `p` to
become **empty**.

## Emptiness

A `unique_ptr` can be empty, in which case it owns nothing – the analog of a
`nullptr`.  This is the default state if is not explicitly initialized with
something.  It is also the state after it loses its ownership of something, or
is forcifully emptied via `.reset()`.

## Releasing

Normally, a `unique_ptr` automatically loses its ownership when it is *moved*
to another `unique_ptr`.  Alternatively, one can force it give up its
ownership in the form of a raw pointer via [`.release()`][5]:

~~~c++
double* raw_ptr = p.release();
~~~

After releasing, the pointer will not be freed automatically by the
`unique_ptr` as it is now empty.  The programmer is then responsible for
deleting the pointer manually.

## Example

Here's a more complete example demonstrating the use of `unique_ptr`:

~~~c++
// compile with: c++ -std=c++14 -lopenblas
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory> // required for std::unique_ptr, std::make_unique, std::size_t
#include <utility> // required for std::move
#include <cblas.h> // required for cblas_daxpy

// create an owned vector initialized to zero
std::unique_ptr<double[]> new_vector(std::size_t n)
{
    return std::make_unique<double[]>(n);
}

// borrow one vector (v1), seize another (v2), and return an owned vector
// equal to their sum
std::unique_ptr<double[]> destructively_add_vectors(
    std::size_t n,
    const double* v1,
    std::unique_ptr<double[]> v2)
{
    // sum <- v2
    std::unique_ptr<double[]> sum = std::move(v2);

    // sum <- v1 + sum
    cblas_daxpy(n, 1., v1, 1, sum.get(), 1);

    // for obscure reasons (C++11 §12.8/32), using an explicit std::move here
    // is actually optional, but we use it anyway for consistency and clarity;
    // see also: https://stackoverflow.com/a/14856553
    return std::move(sum);
}

// borrow two vectors and return an owned vector equal to their sum
std::unique_ptr<double[]> add_vectors(
    std::size_t n,
    const double* v1,
    const double* v2)
{
    // v2_copy <- 0
    std::unique_ptr<double[]> v2_copy = new_vector(n);

    // v2_copy <- v2
    cblas_dcopy(n, v2, 1, v2_copy.get(), 1);

    return destructively_add_vectors(n, v1, std::move(v2_copy));
}

double fibby(double i)
{
    using std::pow;
    using std::sqrt;
    const double a = (1. + sqrt(5.)) / 2.;
    const double b = (1. - sqrt(5.)) / 2.;
    return (pow(a, i) - pow(b, i)) / sqrt(5);
}

// create an owned vector initialized to something fun
std::unique_ptr<double[]> example_vector(std::size_t n)
{
    std::unique_ptr<double[]> v = new_vector(n);
    for (std::size_t i = 0; i != n; ++i) {
        v[i] = fibby(i);
    }
    return v;
}

// borrow a vector and check that the result is correct
void check_result(std::size_t n, const double* v)
{
    for (std::size_t i = 0; i != n; ++i) {
        // always use !(difference < epsilon)
        // rather than (difference >= epsilon)
        // so NaN can't sneak past the check
        if (!(std::abs(v[i] - fibby(i) * 2.) < 1e-8)) {
            std::cerr << "what a cruel, cruel world" << std::endl;
            std::abort();
        }
    }
}

int main()
{
    const std::size_t n = 1024;

    std::unique_ptr<double[]> v1 = example_vector(n);
    std::unique_ptr<double[]> v2 = example_vector(n);

    // lend v1 and v2 to add_vectors
    std::unique_ptr<double[]> sum = add_vectors(n, v1.get(), v2.get());

    // lend sum to check_result
    check_result(n, sum.get());

    // reset sum (actually optional, since the following line will
    // automatically cause this to happen anyway)
    sum.reset();

    // lend v1 and relinquish v2 to destructively_add_vectors
    sum = destructively_add_vectors(n, v1.get(), std::move(v2));

    // lend sum to check_result
    check_result(n, sum.get());

    // v1 and sum are deleted automatically;
    // yay no memory leaks \o/
}
~~~

[1]: http://en.cppreference.com/w/cpp/memory/unique_ptr
[2]: https://www.rust-lang.org
[3]: http://en.cppreference.com/w/cpp/memory/unique_ptr/make_unique
[4]: http://en.cppreference.com/w/cpp/memory/unique_ptr/reset
[5]: http://en.cppreference.com/w/cpp/memory/unique_ptr/release
[6]: http://en.cppreference.com/w/cpp/utility/move
