> I'm a little confused about the roles of `f_ctx`, and `f`. `f` seems to be the void function that takes `y`, and advances forward at time `t` to `yp`. What is the point of `f_ctx`? The header file just says it is an arbitrary pointer passed as the first argument of `f`. But why?

It's a trick in low-level languages like C where functions are represented using a raw function pointer without any enclosing environments.  It is not needed in higher-level languages with [first-class functions](https://en.wikipedia.org/wiki/First-class_function) like Python.  In these languages, functions can be nested and thus store not just the instructions to execute the function but also their enclosing environments (the combination is known as a [closure](https://en.wikipedia.org/wiki/Closure_(computer_programming))).

The need for `f_ctx` arises when you want to write a function that depends on external parameters not known at compile time.  The `f_ctx` parameter allows you to smuggle these external parameters into `f` however you like.

It might be best to illustrate this with an example.  Consider a 1-dimensional numerical integrator like this:

~~~c
double integrate_1(
    double (*f)(double x), /* function to be integrated */
    double x1,
    double x2
);
~~~

This works fine if you know the complete form of the function `f` ahead of time.  But what if this is not the case -- what if the function requires parameters?  Say we want to calculate the [gamma function](https://en.wikipedia.org/wiki/Gamma_function) using an integral:

~~~c
double integrand(double x)
{
    double t = /* where do we get "t" from?? */;
    return pow(x, t - 1.0) * exp(-x);
}

double gamma_function(double t)
{
    /* how do we send the value of "t" into "integrand"? */
    return integrate_1(&integrand, 0.0, INFINITY) / M_PI;
}
~~~

Using `integrand_1` there are only three ways to do this:

 1. Store `t` into a global variable, sacrificing [thread safety](https://en.wikipedia.org/wiki/Thread_safety).  It would be bad to simultaneously call `gamma_function` from different threads as they will both attempt to use the same global variable.

 2. Use a thread-local variable, a feature not available until [C11](https://en.wikipedia.org/wiki/C11_(C_standard_revision)).  At least it is thread-safe now, but it is still not [reentrant](https://en.wikipedia.org/wiki/Reentrancy_(computing)).

 3. Write raw machine code to create an integrand on the fly.  This can be implemented in a thread-safe and reentrant manner, but it is both inefficient, [unportable](https://en.wikipedia.org/wiki/Software_portability), and inhibits compiler optimizations.

However, if the numerical integrator were to be re-designed like this:

~~~c
double integrate_2(
    double (*f)(void *f_ctx, double x),
    void *f_ctx, /* passed into every invocation of "f" */
    double x1,
    double x2
);
~~~

Then there is a much simpler solution that avoids all of these problems:

~~~c
double integrand(void *ctx, double x)
{
    double t = *(double *)ctx;
    return pow(x, t - 1.0) * exp(-x);
}

double gamma_function(double t)
{
    return integrate_2(&integrand, &t, 0.0, INFINITY) / M_PI;
}
~~~

This is thread-safe, reentrant, efficient, and portable.

Note that this problem does not exist in languages like Python where functions can be nested inside other functions:

~~~py
def gamma_function(t):
    def integrand(x):
        return x ** (t - 1.0) * math.exp(-x)
    return integrate(integrand, 0.0, math.inf) / math.pi
~~~
