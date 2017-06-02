---
layout: post
title: "Jacobi-Davidson"
description: "As part of Google Summer of Code I'm implementing several iterative methods for linear systems of equations and eigenvalue problems in Julia. Up to now I have studied Jacobi-Davidson and have a poor man's version working. The first post is about my notes on the method."
tags: [julia,gsoc]
---

Jacobi-Davidson is an iterative method for approximately solving the generalized eigenvalue problem 

{% katex display %}
Ax = \lambda Bx
{% endkatex %}

for a handful of eigenpairs. The matrices are huge and can be general. For ease of discussion we take {% katex %}B = I{% endkatex %}.

Its major selling point is that it can compute interior eigenvalues without iterating with {% katex %}A^{-1}{% endkatex %} exactly; approximate solves are enough. For very large matrices it is typically infeasible to do an exact solve per iteration. 

Other advantages include the possibility to target any point in the complex plane; freedom in the selection procedure (ordinary / harmonic / refined Ritz pairs); and the possibility to specify a good initial search subspace when one knows approximately the directions of the eigenvectors.

Its main disadvantage is expensive iterations in comparison to Arnoldi and the like.

## Subspace methods & Ritz pairs

Jacobi-Davidson is a subspace method, which means that it selects an approximate eigenvector from a low-dimensional **search subspace** {% katex %}\mathcal{V} \subset \mathbb{C}^n{% endkatex %}. Clearly one cannot expect to find the exact eigenvector in a small subspace; therefore one requires the residual to be perpendicular to a **test space** {% katex %}\mathcal{W} \subset \mathbb{C}^n{% endkatex %}. For now we will assume that these two coincide: {% katex %}\mathcal{V} = \mathcal{W}{% endkatex %}. This is called a Galerkin condition.

If the columns of an orthonormal matrix {% katex %}V \in \mathbb{C}^{n \times m}{% endkatex %} form a basis for {% katex %}\mathcal{V}{% endkatex %}, then the above can be summarized as

{% katex display %}
Au - \theta u \perp \mathcal{V} \text{ s.t. } u \in \mathcal{V} \Longleftrightarrow V^*AVy = \theta y
{% endkatex %}

for {% katex %}y \in \mathbb{C}^m{% endkatex %}. Finding {% katex %}(y, \theta){% endkatex %} is cheap, because the matrix {% katex %}V^*AV{% endkatex %} is small. An approximate eigenpair is then obtained as {% katex %}(\theta, Vy){% endkatex %} and is called a Ritz pair. We pick the approximate eigenpair that is of most interest; for instance with a Ritz value closest to a specified target. 

Note that iterative methods build these search and test spaces incrementally: each iteration the spaces are expanded by a new basis vector.

## Jacobi-Davidson's search space

Now what really defines Jacobi-Davidson is the way in which it expands the search subspace: it extends it by a vector that would roughly correct the error in the best approximate eigenvector so far. Suppose we find a Ritz pair {% katex %}(\theta, u){% endkatex %} in our search subspace, then we can improve on this as follows: find a vector {% katex %}t \perp u{% endkatex %} and a scalar {% katex %}\varepsilon{% endkatex %} so that the actual eigenvector {% katex %}x = u + t {% endkatex %} and the actual eigenvalue {% katex %}\lambda = \theta + \varepsilon{% endkatex %}. Let {% katex %}r := Au - \theta u{% endkatex %} be the residual, then it is not hard to show (noting {% katex %} u \perp r{% endkatex %}) that we have the identity

{% katex display %}
r + \varepsilon u + (A - \theta I)t - \varepsilon t = 0.
{% endkatex %}

The term {% katex %}\varepsilon t{% endkatex %} should be really small as it is the product of two errors, so we disregard that. Premultiplying this with the projection matrix {% katex %}(I - uu^*){% endkatex %} which removes components in the direction of {% katex %}u{% endkatex %} we find an equation only in terms of {% katex %}t{% endkatex %}:

{% katex display %}
(I - uu^*)(A - \theta I)t = -r \text{ with } t \perp u.
{% endkatex %}

This is the **correction equation**. Having solved it for {% katex %}t{% endkatex %}, we could in principle decide to just update our approximate eigenvector {% katex %}u{% endkatex%}. However, the point of Jacobi-Davidson is rather to extend the search subspace with {% katex %}t{% endkatex%} and repeat the whole process untill the residual is small enough. 

Lastly, note that if our matrix {% katex %}A{% endkatex %} is Hermetian, we should try and preserve that. One way is to apply the projection matrix to the other side as well: let {% katex display %}\tilde{A} := (I - uu^*)(A - \theta I)(I - uu^*){% endkatex %} and rewrite the correction equation as

{% katex display %}
\tilde{A}t = -r \text{ with } t \perp u.
{% endkatex %}

Since {% katex %}t \perp u{% endkatex %} anyway, this will not influence the solution.

## Solving the correction equation (approximately)
Remember the selling point of Jacobi-Davidson? We aren't even going to solve the correction equation exactly, but only *approximately*, using a handful of iterations of an iterative method.

It turns out that Krylov methods are an excellent choice for this. The condition {% katex %}t \perp u{% endkatex %} is automatically satisfied, since the Krylov subspace 

{% katex display %}
\mathcal{K}_s(\tilde{A}, r) = \text{span} \{ r, \tilde{A}r, \cdots, \tilde{A}^{s-1}r \}
{% endkatex %}

as a whole is automatically orthogonal to {% katex %}u{% endkatex %}, simply because {% katex %}r \perp u {% endkatex %}. This means you should take {% katex %}0{% endkatex %} as the initial guess of the solution.

## Solving the correction equation for testing purposes
Julia's multiple dispatch allows some inversion of control when solving the correction equation. In the Jacobi-Davidson routine itself we don't bother about *how* the correction {% katex %}t{% endkatex %} is retrieved, but we simply specify that is should do that. So we pass the correction equation solver as an argument to the Jacobi-Davidson routine, and simply call it when needed.

It is useful to have an exact correction equation solver, so you can isolate the components of the algorithm when writing tests. The way to go is to use a direct solver for the augmented linear system

{% katex display %}
\begin{bmatrix}
A - \theta I & u \\
u^* & 0
\end{bmatrix}
\begin{bmatrix}
t \\
z
\end{bmatrix}
=
\begin{bmatrix}
-r \\
0
\end{bmatrix}
{% endkatex %}

to find {% katex %}t{% endkatex %}.

## What's next?
So far we've seen that, given a search and test subspace, we can cheaply compute a Ritz pair as an approximate eigenpair. Jacobi-Davidson is characterized by roughly computing a correction to this approximation and use it to extend the search space.

What we haven't seen so far is:

- How to find the next eigenpair once a previous one is converged. *(This is really easy: simply remove it from the search space and subsequently update the correction equation by requiring orthogonality to the converged vectors as well.)*
- Implementation tricks for solving the small-dimensional eigenvalue problem. *(Actually you want to use a Schur decomposition rather than an eigenvalue decomposition.)*
- How to restart the method when the search subspace becomes to big. *(If you use Schur vectors, then this is identical to IRAM's method.)*
- How to incorporate preconditioners when solving the correction equation *(You need to apply deflation to the preconditioner as well, but surprisingly it is not required to perform a lot of matrix-vector multiplications.)*
- How to extract different types of Ritz pairs from the search subspace. *(It turns out harmonic Ritz pairs are better for interior eigenvalue problems. The problem with ordinary Ritz values is that they are typically on their way converging to an exterior eigenvalue, and therefore we could get useless approximate eigenvectors.)*
- Actual Julia code. *(At this point my code is somewhat Matlab-ish, just coding the ideas to verify I get everything.)*

I hope to write more about these things very soon and actually show code.
