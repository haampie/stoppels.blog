---
layout: post
title: "Shrinking the search subspace in Jacobi-Davidson"
description: "Once an approximate eigenvector has converged, we must remove it from the search space, so that other eigenpairs can be found. Also, when a search subspace becomes too large, a restart is necessary. This also requires removing vectors from the search space. Fortunately, shrinking the search subspace is easier in comparison to IRAM."
tags: [gsoc, julia]
---

Once the residual norm of a Ritz pair is small enough, we say a Ritz pair has converged. Until now I have not mentioned how we can find the next eigenpair, but it is clear that the converged eigenvector should not be part of the search space anymore, so that we can find eigenvector corresponding to eigenvalues farther away from our target.

Another reason to shrink the search subspace is when it simply becomes too big. A large search space not only requires much memory, but is computationally expensive as well. We want to keep the computational costs of orthogonalization and extraction fixed, and therefore we simply impose a maximum dimension {% katex %} m_{max} {% endkatex %}.

## Hermetian problems and shrinking the search space
Hermetian matrices have orthogonal eigenvectors for distinct eigenvalues. This is useful, because we work with an orthonormal basis for our search space as well. In fact this orthogonality property is most pronounced when we shrink the search space.

Suppose {% katex %}A{% endkatex %} is Hermetian and the columns of the orthonormal matrix {% katex %}V{% endkatex %} span the search subspace. Then the Galerkin projection of {% katex %}A{% endkatex %} onto the search subspace gives rise to a Hermetian matrix {% katex %}M = V^*AV{% endkatex %}. So if its eigenvalues are distinct, its has an eigendecomposition {% katex %}MS = S\Lambda{% endkatex %} where {% katex %}S{% endkatex %} is unitary and {% katex %}\Lambda = diag(\theta_1, \cdots, \theta_m){% endkatex %}.

Suppose that the Ritz pair {% katex %}(\theta_1, Vs_1){% endkatex %} is converged, then we can remove it from our search space by choosing a new basis {% katex %}V{% endkatex %} in terms of all Ritz vectors excluding {% katex %}Vs_1{% endkatex %}. This is equivalent as updating {% katex %}V{% endkatex %} as

{% katex display%}
V \leftarrow VS[:, 2 : m].
{% endkatex %}

Similarly for restarts. If the search subspace becomes to large so that {% katex %}m = m_{max}{% endkatex %}, we decide to keep only the most promising {% katex %}m_{min}{% endkatex %} Ritz vectors. That would for instance mean updating {% katex %}V{% endkatex %} as

{% katex display %}
V \leftarrow VS[:, 1 : m_{min}].
{% endkatex %}

Each time we update our basis of the search space {% katex %}V{% endkatex %}, we must also update our Galerkin approximation of {% katex %}A{% endkatex %}. Since our new basis consists of Ritz vectors, our Galerkin approximation will be diagonal:

{% katex display %}
(VS)^*A(VS) = S^*MS =  S^*S\Lambda = \Lambda.
{% endkatex %}

Note that we could also restart with our current best approximate eigenvector (that is {% katex %}m_{min} = 1{% endkatex %}), but usually it's a shame to throw away so many good directions. IRAM is another eigensolver that restarts with a small basis of Ritz vectors, but this restart is more subtle, as one must also have that the shrunken search space is a Krylov subspace.

## Change of basis without temporaries
The following is probably not a real issue these days, but if a restart is necessary because of memory limitations, you really would not want to store temporarily both the old basis {% katex %}V{% endkatex %} together with the new basis {% katex %}VS{% endkatex %}. It is however possible to compute this matrix-matrix product in-place.

To do so (assuming {% katex %}S{% endkatex %} is square), compute the LU decomposition {% katex %}S = LU{% endkatex %}, and compute the product

{% katex display %}
V \leftarrow (V * L) * U
{% endkatex %}

column-wise from left to right with {% katex %}L{% endkatex %} and from right to left with {% katex %}U{% endkatex %}. If we must compute the product with only a few columns of {% katex %}S{% endkatex %}, then we could perform only a partial LU decomposition. Lastly, for a stable LU decomposition we need pivoting, but I will not go into that detail.

## Non-Hermetian problems
Non-Hermatian matrices do not necessarily have orthonormal eigenvectors, and since we really want to work with orthonormal vectors, an alternative is to use the Schur decomposition. More about this in a next post.