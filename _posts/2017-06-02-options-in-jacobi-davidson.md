---
layout: post
title: "Options in Jacobi-Davidson"
description: "In my previous post I showed the defining property for Jacobi-Davidson, but it turns out that there is a lot of freedom in the remaining parts of the algorithm."
tags: [gsoc, julia]
---

Let me point out some of the choices we can make when implementing Jacobi-Davidson.

## Initial search space
We could initialize the search space {% katex %}\mathcal{V}{% endkatex %} with a single random vector {% katex %}v_0{% endkatex %}; extract an approximate eigenvector from {% katex %}\mathcal{V}{% endkatex %}; solve the correction equation approximately for {% katex %}t{% endkatex %}; expand the search space by {% katex %}t{% endkatex %}; and repeat the process.

Clearly this makes little sense at the start. The correction equation was motivated by the assumption that the search space already contains a reasonable approximation to an eigenvector, and this is not the case at the very first iteration. In fact the first approximate eigenvector is the initial random vector {% katex %}v_0{% endkatex %}.

Now we have a couple options:

1. As long as the residual of a Ritz pair {% katex %}(\theta, u){% endkatex %} is large, we could replace {% katex %}\theta{% endkatex %} in the correction equation with our target {% katex %}\tau{% endkatex %}.
2. In the first few iterations we could simply add the residual {% katex %}r{% endkatex %} to the search subspace.
3. Option 2 is actually equivalent to "solving" the correction equation with a single iteration of a Krylov solver, so a generalized approach would be to incrementally increase the accuracy level of the Krylov solver at each Jacobi-Davidson iteration.
4. We could also initialize the search space as a Krylov subspace; so Jacobi-Davidson would start out as the Arnoldi method.
5. In some applications we already have a good guess for the eigenvectors and we could initialize the search subspace with them. For instance, in bifurcation analysis people typically compute eigenvalues of a Jacobian matrix multiple times near a bifurcation point. This Jacobian matrix only varies slightly, so one would expect the eigenpairs to vary slightly as well. Hence, you could reuse the eigenvectors from the first solve, to obtain a very good initial search subspace for the second solve.

It seems sensible to use option 3 as a default when there is no prior information about the eigenpairs. Option 4 and 5 can be implemented by letting the user provide an initial search subspace.

## Hermetian problems and shrinking the search space
Hermetian matrices have orthogonal eigenvectors for distinct eigenvalues. This is useful, because we would like to work with an orthonormal basis for our search space. In fact this orthogonality property is most pronounced when we remove converged eigenvectors from our search space, or when we do a restart because the search space becomes to big.

Suppose {% katex %}A{% endkatex %} is Hermetian and the columns of the orthonormal matrix {% katex %}V{% endkatex %} span the search subspace. Then the Galerkin projection of {% katex %}A{% endkatex %} onto the search subspace gives rise to a Hermetian matrix {% katex %}M = V^*AV{% endkatex %}. So if its eigenvalues are distinct, its has an eigendecomposition {% katex %}MS = S\Lambda{% endkatex %} where {% katex %}S{% endkatex %} is unitary and {% katex %}\Lambda = \textrm{diag}(\theta_1, \cdots, \theta_m){% endkatex %}.

Suppose that the Ritz pair {% katex %}(\theta_1, Vs_1){% endkatex %} is converged, then we can remove it from our search space by choosing a new basis {% katex %}V{% endkatex %} in terms of all Ritz vectors excluding {% katex %}Vs_1{% endkatex %}. This is equivalent as updating {% katex %}V{% endkatex %} as

{% katex display%}
V \leftarrow VS[:, 2 : m].
{% endkatex %}

Similarly for restarts. If the search subspace becomes to large so that {% katex %}m = m_{max}{% endkatex %}, we decide to keep only the most promising {% katex %}m_{min}{% endkatex %} Ritz vectors. That would for instance mean updating {% katex %}V{% endkatex %} as

{% katex display %}
V \leftarrow VS[:, 1 : m_{min}].
{% endkatex %}

## Shrinking the search space without temporaries
There are multiple reasons to restart Jacobi-Davidson. One being that orthogonalization becomes too computationally expensive, the other being that the search space must be kept small due to memory constraints. In the latter case we cannot afford to store both the old basis {% katex %}V{% endkatex %} and the new basis {% katex %}VS{% endkatex %}. So somehow we must update {% katex %}V{% endkatex %} in-place.

One way to go is (assuming {% katex %}S{% endkatex %} is square) to compute the LU decomposition {% katex %}S = LU{% endkatex %}, and compute the product

{% katex display %}
V \leftarrow (V * L) * U
{% endkatex %}

column-wise from left to right with {% katex %}L{% endkatex %} and from right to left with {% katex %}U{% katex %}. If we must compute the product with only a few columns of {% katex %}S{% endkatex %}, then we could perform only a partial LU decomposition. Lastly, there is of course the trouble with pivotting for stable LU decompositions.

## Non-Hermetian problems
Non-Hermatian matrices do not necessarily have an orthonormal eigenvectors, and since we really want to work with orthonormal vectors, an alternative is to use the Schur decomposition. In fact our goal simply changes somewhat: rather than returning the eigenvectors, we compute a partial QR decomposition.

