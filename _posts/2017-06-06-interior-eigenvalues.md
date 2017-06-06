---
layout: post
title: "Interior eigenvalues in Jacobi-Davidson"
description: "From the Arnoldi method we know that Ritz values tend to converge to exterior eigenvalues; similar behaviour is observed in Jacobi-Davidson. This is an issue when restarting the method: a Ritz value close to a specified interior target is not necessarily a good approximation of an eigenvalue in that neighbourhood; the Ritz value might as well be on its way to an exterior eigenvalue, and in that case we would not want the Ritz vector to be retained at restart. An alternative is to use harmonic Ritz pairs."
tags: [gsoc, julia]
---

Ritz pairs turn out to be good approximate eigenpairs when one is interested in exterior eigenpairs. If we are interested in *interior* eigenpairs (somewhere inside the convex hull of the spectrum), we must ensure that we do not restart Jacobi-Davidson with a set of Ritz vectors that correspond to Ritz values converging to *exterior* eigenvalues. The simple criterion of Ritz values being close to our target, might give us a low-quality search space at restart. Maybe we could make a better judgement about which Ritz vectors to keep if we would inspect their residuals, but that strategy would be far too expensive.

Shift-and-invert is a well-known trick to solve the above problem: rather than solving the eigenvalue problem {% katex %}Ax = \lambda x{% endkatex %} for {% katex %}\lambda{% endkatex %} close to {% katex %}\tau{% endkatex %}, we reformulate the problem to solving {% katex %}(A - \tau I)^{-1}x = \theta^{-1} x{% endkatex %} with {% katex %}\theta = \lambda - \tau{% endkatex %}. The eigenvectors remain unchanged, but the eigenvalues close to our target are transformed to exterior eigenvalues. The message is that we do not change our methods, rather we change the problem itself.

However, changing the problem like this can be very costly: iterating with a matrix {% katex %}(A - \tau I)^{-1}{% endkatex %} requires solving a system each iteration. In fact, the Arnoldi method requires accurate solves because otherwise the underlying Krylov subspace would be inexact.

## Harmonic Ritz values
The main question is whether we can get higher quality approximations to eigenvalues (without doing much work), so that we can restart with a search space that does not contain partly converged exterior eigenvectors.

One way to go is to look at the Ritz values of {% katex %}(A - \tau I)^{-1}{% endkatex %} with respect to a specially chosen subspace {% katex %}(A - \tau I)V{% endkatex %}, so that the {% katex %}(A - \tau I)^{-1}{% endkatex %} term drops out:

{% katex display %}
(A - \tau I)^{-1}x - \theta^{-1} x \perp (A - \tau I)V \text{ for } x = (A - \tau I)Vy
{% endkatex %}

is equivalent to the requirement

{% katex display %}
(A - \tau I)z - \theta z \perp (A - \tau I)V \text{ for } z = Vy.
{% endkatex %}

The latter is a Petrov-Galerkin projection: the search space is {% katex %}V{% endkatex %} and the test space is {% katex %}(A - \tau I)V{% endkatex %}. We refer to a pair {% katex %}(\theta, Vy){% endkatex %} as a **harmonic Ritz pair**.

The important thing is that {% katex %}\theta{% endkatex %} is identical in both expressions. The second shows a practical means to compute it, while the first shows its interpretation as being a Ritz value of {% katex %}(A - \tau I)^{-1}{% endkatex %}. 

Now the line of thought is as follows: if the Ritz values {% katex %}\theta^{-1}{% endkatex %} converge from within the convex hull of the spectrum of {% katex %}(A - \tau I)^{-1}{% endkatex %} towards the exterior eigenvalues, then conversely its reciprocal {% katex %}\theta{% endkatex %} becomes small only when it approximates an eigenvalue of {% katex %}(A - \tau I){% endkatex %} well. This in turn means that {% katex %}\theta + \tau{% endkatex %} can only be close to an eigenvalue of {% katex %}A{% endkatex %} near {% katex %}\tau{% endkatex %} if it is nearly converged. And that is what we wanted to ensure at restart.

All in all we work with the practical Petrov-Galerkin projection. Let's write {% katex %}\tilde{W} := (A - \tau I)V{% endkatex %} so that the last equation can be rewritten as the generalized eigenvalue problem

{% katex display %}
\tilde{W}^*\tilde{W}y = \theta \tilde{W}^*Vy.
{% endkatex %}

For stability reasons we like to work with orthonormal matrices. The Gram-Schmidt procedure allows us to perform the QR-decomposition of {% katex %}\tilde{W}{% endkatex %} into {% katex %}\tilde{W} = WM_A{% endkatex %} with {% katex %}W{% endkatex %} orthonormal and {% katex %}M_A{% endkatex %} upper triangular. If we also define {% katex %}M := W^*V{% endkatex %} similarly to what we had before, then the generalized eigenvalue problem becomes

{% katex display %}
M_Ay = \theta My
{% endkatex %}

assuming that {% katex %}M_A{% endkatex %} is invertible (which is the case a long as the columns of {% katex %}W{% endkatex %} are independent).

So we can just obtain harmonic Ritz pairs by computing eigenpairs {% katex %}(\theta, y){% endkatex %} from the small generalized eigenvalue problem, and transform them back to approximate eigenpairs by setting {% katex %}\lambda = \tau + \theta{% endkatex %} and {% katex %}z = Vz{% endkatex %}.

In the implementation we will explicitly build these two small matrices {% katex %}M{% endkatex %} and {% katex %}M_A{% endkatex %}.

## Schur vectors vs eigenvectors
We have used a Schur decomposition in favour of an eigendecomposition of the low-dimensional problem to make it easier to shrink the search subspace. Let's see how this works with harmonic Ritz values.

We compute the generalized Schur (QZ) decomposition {% katex %}M_AS^R = S^LT^A{% endkatex %} and {% katex %}MS^R = S^LT{% endkatex %} for unitary matrices {% katex %}S^R, S^L{% endkatex %} and upper triangular matrices {% katex %}T^A, T{% endkatex %}. The eigenvalues {% katex %}\theta{% endkatex %} are now {% katex %}T_{ii}^A / T_{ii}{% endkatex %} and can be reordered. In Julia it is the familiar function `ordschur!`. We must order them from small to large in absolute magnitude. It is not hard to verify that the first column {% katex %}s_1^R{% endkatex %} of {% katex %}S^R{% endkatex %} is an eigenvector with {% katex %}T_{11}^A / T_{11}{% endkatex %} as eigenvalue.

Now if the {% katex %}m_{min}{% endkatex %} harmonic Ritz values we want to keep are moved to the front, we can easily restart by updating both our search space as {% katex %}V \leftarrow VS^R[:, 1 : m]{% endkatex %} and our test space as {% katex %}W \leftarrow WS^L[:, 1 : m]{% endkatex %}.

## Rayleigh quotient vs harmonic Ritz value
Harmonic Ritz values are an excellent choice to determine which Ritz vectors to use at restart. They are however suboptimal for expansion. A few posts ago we wanted the residual {% katex %}r{% endkatex %} to be perpendicular to the approximate eigenvector {% katex %}u{% endkatex %}, so that when solving the correction equation {% katex %}(I - uu^*)(A - \phi I)t = -r{% endkatex %} with a Krylov method, {% katex %}t{% endkatex %} would automatically be perpendicular to {% katex %}u{% endkatex %}.

With harmonic Ritz values this might not be the case. Therefore we can choose to use harmonic Ritz values {% katex %}\theta{% endkatex %} when it comes to shrinking the search space, yet use the Rayleigh quotients {% katex %}\phi{% endkatex %} when it comes to expanding it.

We obtain it by setting {% katex %}(u, r) = 0 {% endkatex %}; in our case {% katex %}u = Vs_1^R{% endkatex %} and {% katex %}r = (A - \tau I)u - \phi u{% endkatex %}. This gives after some manipulations that {% katex %}\phi = \bar{T}_{11} T_{11}^A{% endkatex %}, where the bar denotes complex conjugation.   