---
layout: post
title: "Orthogonalization performance"
description: "Methods like GMRES and Jacobi-Davidson construct an orthogonal basis for their search subspace. It is well-known that classical Gram-Schmidt (CGS) is vulnerable for loss of orthogonality due to rounding errors. Modified Gram-Schmidt (MGS) is usually the fix for this, yet it is not free of rounding errors and is memory-bound when it comes to performance. In this post we'll look into iterative or refined orthogonalization methods and their performance."
tags: [gsoc, julia]
---

In GMRES it is necessary to build an orthonormal basis for the Krylov subspace, because the basis vectors {% katex %}v, Av, A^2v, \cdots{% endkatex %} become quickly linearly dependent in finite precision. In Jacobi-Davidson we require orthonormality of the search subspace and the converged Schur vectors as well, so that converged Schur vector do not re-enter the search subspace. Also, throughout both algorithms we simplify expressions via equalities like {% katex %}V^*V = I{% endkatex %}. Due to rounding errors this identity might however not hold exactly. Therefore it is important to have some guarantee of orthogonality.

The benchmarks in this post are performed on an Intel i5-4460 @ 3.20GHz using 16GB of DDR3 memory at 1600MHz.

## Classical and Modified Gram-Schmidt

Given a matrix {% katex %}V \in \mathbb{C}^{n \times m}{% endkatex %} where {% katex %}n \gg m{% endkatex %} and a vector {% katex %}w \in \mathbb{C}^n{% endkatex %}, our goal is to remove the components of {% katex %}w{% endkatex %} spanned by the columns of {% katex %}V{% endkatex %}.

This comes down to updating {% katex %}w := (I - VV^*)w{% endkatex %} and subsequently normalizing {% katex %}w := w / \|w\|{% endkatex %}. In GMRES we want to store the projection {% katex %}V^*w{% endkatex %} as well.

Classical Gram-Schmidt in Julia would look like this:

{% highlight julia %}
import Base.LinAlg.BLAS: gemv!, gemv

function classical_gram_schmidt!{T}(V::StridedMatrix{T}, w::StridedVector{T})
    
    # Orthogonalize
    h = gemv('T', one(T), V, w)
    gemv!('N', -one(T), V, h, one(T), w)

    # Normalize
    nrm = norm(w)
    scale!(w, one(T) / nrm)

    h, nrm
end
{% endhighlight %}

It is well-known that due to rounding errors the above procedure is not stable. Modified Gram-Schmidt tries to "solve" this problem by performing the projection sequentially:

{% highlight julia %}
import Base.LinAlg.BLAS: axpy!

function modified_gram_schmidt!{T}(V::StridedMatrix{T}, w::StridedVector{T})
    k = size(V, 2)
    h = zeros(T, k)

    # Orthogonalize
    for i = 1 : k
        column = @view(V[:, i])
        h[i] = dot(column, w)
        axpy!(-h[i], column, w)
    end

    # Normalize
    nrm = norm(w)
    scale!(w, one(T) / nrm)

    h, nrm
end
{% endhighlight %}

Performance-wise this for-loop is a major bottleneck. The problem is that we have gone from BLAS-2 (matrix-vector product) to BLAS-1 (inner products).

If we bench the above with `m = 20` and `n = 10_000_000` we'll quickly see the problem:

{% highlight julia %}
using BenchmarkTools

function bench(; m = 20)
    srand(1)
    
    suite = BenchmarkGroup()
    
    for threads = [1, 2, 4]
        BLAS.set_num_threads(threads)
        thread_key = string(threads)
        suite[thread_key] = BenchmarkGroup()
        for n = [1_000_000, 10_000_000]
            V, _ = qr(rand(n, m))
            new_vec = rand(n)
            key = string(n)
            suite[thread_key][key] = BenchmarkGroup()
            suite[thread_key][key]["cgs"] = @benchmark classical_gram_schmidt!($V, a) setup=(a = copy($new_vec))
            suite[thread_key][key]["mgs"] = @benchmark modified_gram_schmidt!($V, a) setup=(a = copy($new_vec))
        end
    end

    suite
end
{% endhighlight %}

This outputs:

{% highlight text %}
  "1" => 2-element BenchmarkTools.BenchmarkGroup:
          "1000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "mgs" => Trial(39.492 ms)
                  "cgs" => Trial(19.027 ms)
          "10000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "mgs" => Trial(414.534 ms)
                  "cgs" => Trial(199.866 ms)
  "2" => 2-element BenchmarkTools.BenchmarkGroup:
          "1000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "mgs" => Trial(35.880 ms)
                  "cgs" => Trial(16.462 ms)
          "10000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "mgs" => Trial(393.617 ms)
                  "cgs" => Trial(168.942 ms)
  "4" => 2-element BenchmarkTools.BenchmarkGroup:
          "1000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "mgs" => Trial(35.084 ms)
                  "cgs" => Trial(15.886 ms)
          "10000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "mgs" => Trial(392.145 ms)
                  "cgs" => Trial(159.393 ms)
{% endhighlight %}

which shows that MGS is more than twice as slow. The reason for this is that BLAS-1 is memory bound: all data is touched only once, and there is a lot of data in memory. For this reason it is almost useless to have threading enabled.

## Iterative / refined orthogonalization
Even though MGS is an easy trick to ensure better orthogonality, it does not necessarily generate vectors that are orthogonal up to machine precision. For instance, if we do

{% highlight text %}
srand(1)
V, _ = qr(rand(1_000_000, 10))
w = rand(1_000_000)
modified_gram_schmidt!(V, w)
@show norm(V' * w)
modified_gram_schmidt!(V, w)
@show norm(V' * w)
{% endhighlight %}

it outputs

{% highlight text %}
norm(V' * w) = 1.0289082309190917e-15
norm(V' * w) = 4.513770868406318e-17
{% endhighlight %}

which basically means that the projection is not idempotent. After the first pass of MGS there are still components in the direction of {% katex %}V{% endkatex %} and only in the second pass these are removed almost up to machine precision.

If we really need to maintain orthogonality, one idea is to repeat MGS twice or even more. In Krylov subspace there is the rule of thumb that "twice is enough" due to Kahan and Parlett, but we can just as well decide whether to repeat orthogonalization by using the Daniel-Gragg-Kaufman-Stewart criterion. This criterion states that we should repeat orthogonalization if the tangent between {% katex %}w{% endkatex %} and {% katex %}span(V){% endkatex %} is smaller than {% katex %}\eta{% endkatex %}. In ARPACK it seems that this parameter is chosen {% katex %}\eta = 1 / \sqrt{2}{% endkatex %}.

I'm going to assume reorthogonalization is only necessary once, but the `if` can be replaced by a `while`.

{% highlight julia %}
function repeated_classical_gram_schmidt!{T}(V::StridedMatrix{T}, w::StridedVector{T}; η = one(T) / √2)
    
    # Orthogonalize
    h = gemv('T', one(T), V, w)
    gemv!('N', -one(T), V, h, one(T), w)

    # Repeat
    nrm = norm(w)

    if nrm < η * norm(h)
        increment = gemv('T', one(T), V, w)
        gemv!('N', -one(T), V, increment, one(T), w)
        axpy!(one(T), increment, h)
        nrm = norm(w)
    end

    # Normalize
    scale!(w, one(T) / nrm)

    h, nrm
end
{% endhighlight %}

Similarly for MGS:

{% highlight julia %}
function repeated_modified_gram_schmidt!{T}(V::StridedMatrix{T}, w::StridedVector{T}; η = one(T) / √2)
    k = size(V, 2)
    h = zeros(T, k)

    nrm = norm(w)

    # Orthogonalize
    for i = 1 : k
        column = @view(V[:, i])
        h[i] = dot(column, w)
        axpy!(-h[i], column, w)
    end

    if norm(w) ≤ η * nrm
        for i = 1 : k
            column = @view(V[:, i])
            correction = dot(column, w)
            h[i] += correction
            axpy!(-correction, column, w)
        end

        nrm = norm(w)
    end

    # Normalize
    scale!(w, one(T) / nrm)

    h, nrm
end
{% endhighlight %}

We benchmark these methods as well:

{% highlight julia %}
function bench(; m = 20)
    srand(1)
    
    suite = BenchmarkGroup()
    
    for threads = [1, 2, 4]
        BLAS.set_num_threads(threads)
        thread_key = string(threads)
        suite[thread_key] = BenchmarkGroup()
        for n = [1_000_000, 10_000_000]
            V, _ = qr(rand(n, m))
            new_vec = rand(n)
            key = string(n)
            suite[thread_key][key] = BenchmarkGroup()
            suite[thread_key][key]["rcgs"] = @benchmark repeated_classical_gram_schmidt!($V, a) setup=(a = copy($new_vec))
            suite[thread_key][key]["rmgs"] = @benchmark repeated_modified_gram_schmidt!($V, a) setup=(a = copy($new_vec))
        end
    end

    suite
end
{% endhighlight %}

and their results are:

{% highlight text %}
  "1" => 2-element BenchmarkTools.BenchmarkGroup:
          "1000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "rcgs" => Trial(37.474 ms)
                  "rmgs" => Trial(79.268 ms)
          "10000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "rcgs" => Trial(393.802 ms)
                  "rmgs" => Trial(829.219 ms)
  "2" => 2-element BenchmarkTools.BenchmarkGroup:
          "1000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "rcgs" => Trial(32.237 ms)
                  "rmgs" => Trial(72.151 ms)
          "10000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "rcgs" => Trial(331.603 ms)
                  "rmgs" => Trial(787.374 ms)
  "4" => 2-element BenchmarkTools.BenchmarkGroup:
          "1000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "rcgs" => Trial(31.030 ms)
                  "rmgs" => Trial(69.855 ms)
          "10000000" => 2-element BenchmarkTools.BenchmarkGroup:
                  "rcgs" => Trial(323.154 ms)
                  "rmgs" => Trial(785.251 ms)
{% endhighlight %}

The good news is that repeated classical Gram-Schmidt is still faster than modified Gram-Schmidt without repetition.

## Conclusion
If orthogonality is so important that repeated orthogonalization is necessary, then repeated classical Gram-Schmidt seems a good candidate. Even if repeated orthogonalization is not a necessity, repeated classical Gram-Schmidt can even be faster than modified Gram-Schmidt without repetition.

Note that this implies that the orthonormal basis {% katex %}V{% endkatex %} is indeed best represented as a `Matrix` rather than a `Vector` of `Vector`s.

### Remark

In the benchmarks above I normalize a vector by first computing `t = 1 / norm(w)` and only then computing `w *= t`. The good thing is that multiplication is a lot faster than division, but the downside is that we round numbers twice. I don't know if the latter has significant consequences for orthogonality.

