---
layout: post
title: "Compile-time primes"
description: "Anyone who has seen C++ template meta-programming must have come across the compile-time prime number generator by Erwin Unruh. Since then the language has evolved and provides a very readible alternative. In this post I'll touch on both the classical and modern approach."
tags: [cpp]
---

The building blocks of C++ template meta-programming are usually `struct`s without any data members encapsulating a value stored in an `enum`. As a simple example consider

{% highlight c++ %}
template <int x> struct MyInteger
{
  enum { value = x };
};
{% endhighlight %}

This can be used to store an integer as a plain value. It has no address at run-time, since it is not a data member. Storing and retrieving the value is done as follows:

{% highlight c++ %}
int main()
{
  std::cout << MyInteger<100>::value << '\n';
  // prints 100.
}
{% endhighlight %}

Storing and retrieving values is not very exciting, but actually computing something using types *is*. The simplest way to do so is to exploit recursive definition of types. An instructive example is the Fibonacci sequence:

{% highlight c++ %}
template <size_t x> struct Fibonacci
{
  enum { 
    value = Fibonacci<x - 1>::value + Fibonacci<x - 2>::value
  };
};

int main()
{
  std::cout << Fibonacci<10>::value << '\n';
}
{% endhighlight %}

Compiling the above will throw errors at us, because the compiler gets stuck in an infinite recursion. This is clearly reflected in the error message, where it shows the compiler quits at depth 256:

{% highlight none %}
fib.cc:6:13: fatal error: recursive template instantiation exceeded maximum depth of 256
    value = Fibonacci<x - 1>::value + Fibonacci<x - 2>::value
            ^
fib.cc:6:13: note: in instantiation of template class 'Fibonacci<18446744073709551370>'
      requested here
fib.cc:6:13: note: in instantiation of template class 'Fibonacci<18446744073709551371>'
      requested here
[...]
{% endhighlight %}

Note the overflow in the template parameter, which happens because `size_t` is unsigned.

To stop the recursion, we must implement a base case for x = 0 and x = 1. This is done via specialization:

{% highlight c++ %}
template<> struct Fibonacci<0>
{
  enum { value = 0 };
};

template<> struct Fibonacci<1>
{
  enum { value = 1 };
};
{% endhighlight %}

Now the compiler does its job and evaluates `Fibonacci<10>::value` to 55, which happens to be the tenth Fibonacci number.

## Primality test using template meta-programming
A trivial primality test for a number p, is to do a linear scan over the numbers 2 to p - 1 and check whether p is divisible by none of them. Of course there are all kinds of improvements to make, but for the sake of simplicity we will not pursue these. The primality test is implemented using the following recursion

{% highlight c++ %}
template <size_t p, size_t i> struct PrimalityTest
{
  enum { value = p % i != 0 && PrimalityTest<p, i - 1>::value }; 
}; 

template <size_t p> struct PrimalityTest<p, 1>
{
  enum { value = 1 };
};
{% endhighlight %}

The specialization `PrimalityTest<p, 1>` guarantees the compiler will not get stuck in infinite recursion limbo. Furthermore, any prime number is of course divisible by 1, so `PrimalityTest<p, 1>::value == 1`. For any other `i` the value of `PrimalityTest<p, i>::value` will be 1 if `p` is neither divisible by `i` nor by any number less than `i`. The latter part is handled recursively.

To call the primality test we can create another struct that does exactly this:

{% highlight c++ %}
template <size_t p> struct IsPrime
{
  enum { value = PrimalityTest<p, p - 1>::value };
};
{% endhighlight %}

It simply states that a number `p` is prime if it is not divisible by any number smaller than `p`. A quick sanity check shows that this indeed works:

{% highlight c++ %}
int main()
{
  std::cout << IsPrime<13>::value << ' ' 
            << IsPrime<24>::value << '\n';
  // outputs 1 0 correctly.
}
{% endhighlight %}

As a concluding remark for this section, a very handy tool to do compile-time testing is `static_assert(bool_constexpr, message)` which arrived in C++11. Not only can it be used for testing, but as well to convince yourself it is really the compiler which evaluates an expression. In our example we can use it as follows:

{% highlight c++ %}
int main()
{
  static_assert(IsPrime<13>::value == 1, "13 is prime");
  static_assert(IsPrime<24>::value == 0, "24 is not prime");
}
{% endhighlight %}

## A modern alternative in C++14
Template meta-programming is very powerful, yet as a programming language it generates an enormous cognitive load. Although it sometimes has the looks of a quite mature functional programming language, it also gives the feeling of abuse of what probably is a C++ by-product. Indeed, historically Turing-completeness of the templating system was not a goal, but a welcome surprise.

To make life easier, C++11 introduced the notion of `constexpr` functions and variables. This identifier tells the compiler that the value of the function or variable can be evaluated at compile-time. As a simple example, our recursive Fibonacci function can be implemented as a one-liner `constexpr` function:

{% highlight c++ %}
constexpr size_t fibonacci(size_t n)
{
  return n < 2 ? n : fibonacci(n - 1) + fibonacci(n - 2);
}

int main()
{
  static_assert(fibonacci(0) == 0, "f(0) == 0");
  static_assert(fibonacci(1) == 1, "f(1) == 1");
  static_assert(fibonacci(10) == 55, "f(10) == 55");
}
{% endhighlight %}

This piece of code is almost indistinguishable from the code you would normally write in C++, except for its overdone compactness and recursive inefficiency.

A more efficient implementation would be non-recursive:

{% highlight c++ %}
constexpr size_t fibonacci(size_t n)
{
  if (n < 2)
    return n;

  size_t previous = 0;
  size_t current = 1;

  for (size_t idx = 2; idx <= n; ++idx)
  {
    size_t temp = current;
    current = previous + current;
    previous = temp;
  }

  return current;
}

int main()
{
  static_assert(fibonacci(0) == 0, "f(0) == 0");
  static_assert(fibonacci(1) == 1, "f(1) == 1");
  static_assert(fibonacci(10) == 55, "f(10) == 55");
}
{% endhighlight %}

However, this piece of code only compiles in C++14. Why? In C++11 `constexpr` functions cannot contain variable declarations and should in fact be as simple as a return statement. C++14 has greatly extended the allowed contents of a `constexpr` function body, which makes it possible to write code identical to what you would write for run-time functions.

## A compile-time prime sieve for counting prime numbers
As a final example to appreciate the beauty and extend of what `constexpr` can do, let's implement a compile-time prime counting function using the sieve of Eratosthenes. The function `count_primes<N>()` returns the number of primes strictly smaller than `N`. Using some well-known prime sieving trickery, a relatively efficient implementation is:

{% highlight c++ %}
template<size_t N> constexpr size_t count_primes()
{
  size_t total = 0;
  bool is_prime[N]{};

  // Assume all numbers are prime.
  for (size_t n = 0; n < N; ++n)
    is_prime[n] = true;

  // Sieve primes starting at 2.
  for (size_t n = 2; n * n < N; ++n)
    if (is_prime[n])
      for (size_t multiple = n * n; multiple < N; multiple += n)
        is_prime[multiple] = false;

  // Count primes.
  for (size_t n = 2; n < N; ++n)
    if (is_prime[n])
      ++total;
  
  return total;
}
{% endhighlight %}

Let's throw some compile-time tests at it for `N` as large as 10000:

{% highlight c++ %}
int main()
{
  static_assert(count_primes<0>() == 0, "");
  static_assert(count_primes<1>() == 0, "");
  static_assert(count_primes<2>() == 0, "");
  static_assert(count_primes<3>() == 1, "");
  static_assert(count_primes<100>() == 25, "");
  static_assert(count_primes<10000>() == 1229, "");
}
{% endhighlight %}

The above code compiles just fine, which can only mean that *the compiler can count prime numbers efficiently in just 22 human-readable lines of code*.

Unfortunately we cannot pass `N` as a function parameter rather than a template parameter, as `bool is_prime[N]` would then be interpreted as a variable-sized object.

### Can we improve further? What are the limitations?
One idea is to use STL algorithms to initialize the array with `true` as values. The problem however is that `std::fill` is not (and probably will never be) a `constexpr` function. This is because it will resort to `memset` to initialize the memory whenever it can.

We could also try to use `std::array` as an abstraction, but this has a `constexpr` implementation for `operator[](size_t) const` only from C++17 onwards. At the time of writing this seems not implemented in Clang 8.0.

A very obvious limitation is of course the possibility of a stack overflow for large values of N, but compile-time computations of that magnitude make probably little sense.

## Concluding remarks
C++ has come a long way from the discovery of template meta-programming to the very versatile `constexpr` identifier. Even if you're sceptical about its practical usefulness (which you shouldn't), the power of the language and in particular its compiler can leave you stunned. 


