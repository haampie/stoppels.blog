---
layout: post
title: "Jacobi-Davidson example"
description: "Let's take a look how the Jacobi-Davidson algorithm targets interior eigenvalues of a Poisson matrix."
tags: [gsoc, julia]
---

I've put Jacobi-Davidson to the test on the simple Poisson matrix from the previous post: it has {% katex %}\begin{bmatrix}-1 & 2 & -1\end{bmatrix}{% endkatex %} on the diagonals and is {% katex %}100 \times 100{% endkatex %}. Its eigenvalues are all real and located in the interval {% katex %}[0, 4]{% endkatex %}.

## Using the exact solver
When solving the correction equation with the "exact" solver, we should be able to see quadratic convergence. I picked {% katex %}m_{min} = 10{% endkatex %} and {% katex %}m_{max} = 15{% endkatex %} as minimum and maximum dimension of the search space. The target is {% katex %}\tau = 3{% endkatex %}, and we're after 10 eigenvalues closest to it. The API for this is currently just

{% highlight julia %}
Q, R, ritz_hist, conv_hist, residuals = harmonic_ritz_test(
  A,
  exact_solver(),
  pairs = 10,
  min_dimension = 10,
  max_dimension = 15,
  max_iter = 200,
  ɛ = 1e-8,
  τ = 3.0 + 0im
)
{% endhighlight %}

A sanity check is to see whether the partial Schur decomposition {% katex %}AQ = QR{% endkatex %} is correct:

{% highlight julia %}
> norm(A * Q - Q * R)
7.323237382364171e-9
> norm(Q' * Q - eye(10))
3.844239740211255e-9
{% endhighlight %}

That seems reasonable for our prescribed tolerance. Next we look at the convergence history.

<figure class="full_width_fig">
    {% img jd/convergence.svg %}

    <p>Green: harmonic Ritz values. Yellow: converged Ritz values. Black: actual eigenvalues. Horizontal blue line: target.</p>
</figure>

A couple things to note from the convergence plot: after 15 iterations the search space becomes too big, and clearly only the best approximate Ritz values are retained at restart. Another thing is that as soon as a first Ritz value converges, the others converge very quickly as well. This is the main reason why we want a search space in the first place.

## Using a couple steps of GMRES
Using the exact solver for the correction equation is basically cheating; Jacobi-Davidson should be seen as the method that targets interior eigenvalues without solving systems exactly at each iteration.

What we can do instead is use just 7 steps of GMRES *without* preconditioning, and see what happens then.

{% highlight julia %}
Q, R, ritz_hist, conv_hist, residuals = harmonic_ritz_test(
  A,
  gmres_solver(iterations = 7),
  pairs = 10,
  min_dimension = 10,
  max_dimension = 15,
  max_iter = 500,
  ɛ = 1e-8,
  τ = 3.0 + 0im
)
{% endhighlight %}

Again we do the sanity check:

{% highlight julia %}
> norm(Q' * Q - eye(10))
5.995057809828708e-8
> norm(A * Q - Q * R)
1.3111792852021023e-8
{% endhighlight %}

And we inspect the convergence plot

<figure class="full_width_fig">
    {% img jd/convergence_gmres.svg %}

    <p>Green: harmonic Ritz values. Yellow: converged Ritz values. Black: actual eigenvalues. Horizontal blue line: target.</p>
</figure>

Now we need a whole lot of iterations more to let a Ritz value converge, but the convergence seems quite monotonically and is *not* hindered by restarts, which is great. Also, it's clear that once the first Ritz value has converged, the others follow quickly as well.

The convergence can probably be improved if we would add a preconditioner, but this is something I have yet to implement.

## Loss of orthogonality
It seems that during the iterations a simple modified Gram-Schmidt is not always enough to guarantee orthogonality of matrices {% katex %}V{% endkatex %} and {% katex %}W{% endkatex %}. A solution is to [do Gram-Schmidt *twice*](http://www.netlib.org/utk/people/JackDongarra/etemplates/node138.html#alg:rgs), or use [iterative refinement](http://www.netlib.org/utk/people/JackDongarra/etemplates/node221.html). This seems rather costly, but maybe it's necessary.