---
layout: post
title: "Towards the next eigenpair in Jacobi-Davidson"
description: "After removing a converged Ritz vector from the search space, we must ensure that new expansions do not bring it back. In practice this means updating the correction equation so that these directions are deflated."
tags: [gsoc, julia]
---

In the Hermetian case it is very easy to find the next approximate eigenpair once another is converged. We simply require our correction {% katex %}t{% endkatex %} to be orthogonal not only to {% katex %}u{% endkatex %}, but as well to all previously converged Ritz vectors. If {% katex %}X{% endkatex %} contains the converged Ritz vectors, define the orthonormal matrix {% katex %}Q = [X, u]{% endkatex %}, so that {% katex %}(I - QQ^*){% endkatex %} is a projection. Our new correction equation will then be

{% katex display %}
(I - QQ^*)A(I - QQ^*)t = -r \text{ and } t \perp u.
{% endkatex %}

In the non-Hermetian case this strategy will not work as {% katex %}(I - QQ^*){% endkatex %} will fail to be a projection if the Ritz vectors are not orthogonal. So what can we do? 

## The non-Hermetian case

To fix this problem, we simply change our goal somewhat. Rather than computing and storing eigenvectors, we compute Schur vectors, so that we can guarantee orthonormality of {% katex %}Q{% endkatex %} and therefore use the new correction equation.

However, if we choose to work with Schur vectors rather than eigenvectors, we must also come up with a new convergence criterion. Suppose the converged Schur vectors are stored in {% katex %}Q{% endkatex %}. Previously we used {% katex %}\|Au - \theta u\|_2 < \varepsilon{% endkatex %} for a Ritz pair {% katex %}(\theta, u){% endkatex %}. But if {% katex %}u{% endkatex %} is a Schur vector it lacks directions of the eigenvector spanned by the columns of {% katex %}Q{% endkatex %}, so the residual will always be large. To remedy this, we inspect the size of the residual with the directions of {% katex %}Q{% endkatex %} removed:

{% katex display %}
\|(I - QQ^*)(Au - \theta u)\|_2 = \|r - Q(Q^*r)\|_2 < \varepsilon
{% endkatex %}

When a few Schur vectors are converged we get a decomposition {% katex %}AQ = QR{% endkatex %} with {% katex %}R{% endkatex %} upper triangular and {% katex %}Q{% endkatex %} unitary. The next converged Schur vector {% katex %}u{% endkatex %} is appended to {% katex %}Q{% endkatex %} as {% katex %}Q \leftarrow \begin{bmatrix}Q & u\end{bmatrix}{% endkatex %}. Updating {% katex %}R{% endkatex %} is easy as well, as it must be of the form

{% katex display %}
R \leftarrow \begin{bmatrix}R & z \\ 0 & \theta\end{bmatrix}.
{% endkatex %}

Equating {% katex %}AQ = QR{% endkatex %} implies {% katex %}Qz = r{% endkatex %}, so therefore {% katex %}z := Q^*r{% endkatex %}. Note that this quantity is already computed when we validate the convergence criterion.

