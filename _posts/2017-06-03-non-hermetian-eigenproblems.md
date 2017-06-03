---
layout: post
title: "Non-Hermetian eigenproblems in Jacobi-Davidson"
description: "Non-Hermetian matrices do not necessarily have orthogonal eigenvectors for distinct eigenvalues. So far we have relied upon this property when shrinking the search subspace. Fortunately we can work with the Schur decomposition to tackle this problem."
tags: [gsoc, julia]
---

In the previous post we replaced our orthonormal basis {% katex %}V{% endkatex %} by a basis of Ritz vectors, so that it was easy to remove a Ritz vector from the search subspace, or even discard multiple of the least interesting Ritz vectors as a way to restart Jacobi-Davidson. The new basis of Ritz vectors is orthonormal because {% katex %}M = V^*AV{% endkatex %} was Hermetian. If our matrix {% katex %}A{% endkatex %} is non-Hermetian, we lose the orthogonality property, so we must come up with something better.

Fortunately we can work with the Schur decomposition of {% katex %}M{% endkatex %}

{% katex display %}
MS = SU
{% endkatex %}

where {% katex %}S{% endkatex %} is unitary as desired and {% katex %}U{% endkatex %} is upper triangular with the eigenvalues of {% katex %}M{% endkatex %} (the Ritz values of {% katex %}A{% endkatex %}) on the diagonal. Clearly the eigendecomposition of a Hermetian matrix is a special case of the Schur decomposition.

The downside is that we do not have direct access to all Ritz vectors. Only the "first" Ritz vector {% katex %}Vs_1{% endkatex %} corresponding to the Ritz value {% katex %}U_{11}{% endkatex %} is directly accessible. If the Ritz value closest to our target is not in {% katex %}U_{11}{% endkatex %}, we must do some work to fetch its Ritz vector. This can be done quite elegantly, but I won't bother about the details for now. The point is that the diagonal elements of {% katex %}U{% endkatex %} can be (partially) reorded with computationally inexpensive Given's rotations, so that we can retrieve any Ritz vector we like. In Julia this can be done via `ordschur`.

So basically not much changes when we're after a single eigenpair of a general matrix in Jacobi-Davidson:

- Since {% katex %}M = V^*AV{% endkatex %} is now non-Hermetian, we need to compute both its new column and new row when the subspace is expanded. In the Hermetian case we only have to compute a new colum.
- In the Hermetian case we can immediately select any Ritz pair {% katex %}(\theta_i, Vs_i){% endkatex %} from the eigendecomposition of {% katex %}M{% endkatex %}. In the non-Hermetian case we must reorder the Schur form to bring the Ritz value of interest to the upper left corner of {% katex %}U{% endkatex %}.
- Restarts in the non-Hermetian case require us to bring {% katex %}m_{min}{% endkatex %} of the Ritz values we want to keep to the first entries on the diagonal of {% katex %}U{% endkatex %}.

In summary we now know how to prune the search space in the case of non-Hermetian matrices. Yet we have not seen how to find the next eigenpair after a first has converged. This will be my next post.