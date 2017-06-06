---
layout: post
title: "Interior eigenvalues in Jacobi-Davidson"
description: "From the Arnoldi method we know that Ritz values tend to converge to exterior eigenvalues; similar behaviour is observed in Jacobi-Davidson. This is an issue when restarting the method: a Ritz value close to a specified interior target is not necessarily a good approximation of an eigenvalue in that neighbourhood; the Ritz value might as well be on its way to an exterior eigenvalue, and in that case we would not want the Ritz vector to be retained at restart. An alternative is to use harmonic Ritz pairs."
tags: [gsoc, julia]
---

Ritz pairs turn out to be good approximate eigenpairs when one is interested in exterior eigenpairs. If we are interested in *interior* eigenpairs (somewhere inside the convex hull of the spectrum), we must ensure that we do not restart Jacobi-Davidson with a set of Ritz vectors that correspond to Ritz values converging to *exterior* eigenvalues. The simple criterion of Ritz values being close to our target, might give us a low-quality search space at restart. Maybe we could make a better judgement about which Ritz vectors to keep if we would inspect their residuals, but that strategy would be far too expensive.

Shift-and-invert is a well-known trick to solve the above problem: rather than solving the eigenvalue problem {% katex %}Ax = \lambda x{% endkatex %} for {% katex %}\lambda{% endkatex %} close to {% katex %}\tau{% endkatex %}, we reformulate the problem to solving {% katex %}(A - \tau I)^{-1}x = \theta x{% endkatex %} with {% katex %}\theta = (\lambda - \tau)^{-1}{% endkatex %}. The eigenvectors remain unchanged, but the eigenvalues close to our target are transformed to exterior eigenvalues. The message is that we do not change our methods, rather we change the problem itself.

However, changing the problem like this can be very costly: iterating with a matrix {% katex %}(A - \tau I)^{-1}{% endkatex %} requires solving a system each iteration. In fact, the Arnoldi method requires accurate solves because otherwise the underlying Krylov subspace would be inexact.

## Harmonic Ritz values
The main question is whether we can get higher quality approximations to eigenvalues, so that we can restart with a higher quality search space (without partly converged exterior eigenvectors), without doing too much work. And in fact we can do this more or less.

If we look at the Ritz values of {% katex %}(A - \tau I)^{-1}{% endkatex %} with respect to the search and test space {% katex %}(A - \tau I)V{% endkatex %}, then we see that {% katex %}(A - \tau I)^{-1}{% endkatex %} drops out:

{% katex display %}
(A - \tau I)^{-1}x - \theta x \perp (A - \tau I)V \text{ for } x = (A - \tau I)Vy
{% endkatex %}

is equivalent to the requirement

{% katex display %}
(A - \tau I)z - \theta^{-1} z \perp (A - \tau I)V \text{ for } z = Vy.
{% endkatex %}

The latter is a Petrov-Galerkin projection: the search space is {% katex %}V{% endkatex %} and the test space is {% katex %}(A - \tau I)V{% endkatex %}.

My intuition was that we were going to solve the last problem for {% katex %}(\theta^{-1}, z){% endkatex %} and then work with {% katex %}(\theta, x){% endkatex %} as approximate eigenpairs for {% katex %}(A - \tau I)^{-1}{% endkatex %}. Turns out we stick to {% katex %}z{% endkatex %}, since working with {% katex %}x{% endkatex %} requires another matrix-vector product.

Why bother at all then? Well, the point is that this {% katex %}\theta{% endkatex %} is the same in both settings. And since it is a Ritz value of {% katex %}(A - \tau I)^{-1}{% endkatex %}, it will probably be quite good.

Let's write {% katex %}\tilde{W} := (A - \tau I)V{% endkatex %} so that the last equation can be rewritten as the generalized eigenvalue problem

{% katex display %}
\tilde{W}^*Vy = \theta \tilde{W}^*\tilde{W}y.
{% endkatex %}

For stability reasons we like to work with orthonormal matrices. The Gram-Schmidt procedure allows us to perform the QR-decomposition of {% katex %}\tilde{W}{% endkatex %} into {% katex %}\tilde{W} = WM_A{% endkatex %} with {% katex %}W{% endkatex %} orthonormal and {% katex %}M_A{% endkatex %} upper triangular. If we also define {% katex %}M := W^*V{% endkatex %} similarly to what we had before, then the generalized eigenvalue problem becomes

{% katex display %}
My = \theta M_Ay
{% endkatex %}

assuming that {% katex %}M_A{% endkatex %} is invertible (which is the case a long as the columns of {% katex %}W{% endkatex %} are independent).

So all in all we can "just" obtain harmonic Ritz pairs by computing eigenpairs {% katex %}(\theta, z){% endkatex %} from the small generalized eigenvalue problem, and transform them back to approximate eigenpairs by setting {% katex %}\lambda = \tau + \theta^{-1}{% endkatex %} and {% katex %}x = (A - \tau I)Vz{% endkatex %}.

In the implementation we will explicitly build these two small matrices {% katex %}M{% endkatex %} and {% katex %}M_A{% endkatex %}.

## Schur vectors vs eigenvectors
We have used a Schur decomposition in favour of an eigendecomposition of the low-dimensional problem to make it easier to shrink the search subspace. Let's see how this works with harmonic Ritz values.

We compute the generalized Schur (QZ) decomposition {% katex %}M_AS^R = S^LT^A{% endkatex %} and {% katex %}MS^R = S^LT{% endkatex %} for unitary matrices {% katex %}S^R, S^L{% endkatex %} and upper triangular matrices {% katex %}T^A, T{% endkatex %}. The eigenvalues {% katex %}\theta{% endkatex %} are now {% katex %}T_{ii}^A / T_{ii}{% endkatex %} and can be reordered. It is not hard to verify that the first column {% katex %}s_1^L{% endkatex %} of {% katex %}S^L{% endkatex %} is an eigenvector with {% katex %}T_{11}^A / T_{11}{% endkatex %} as eigenvalue.

