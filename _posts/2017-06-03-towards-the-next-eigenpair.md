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

To fix this problem, we simply change our goal somewhat. Rather than computing and storing eigenvectors, we compute Schur vectors, so that we can guarantee orthonormality of {% katex %}Q{% endkatex %} and therefore use the new correction equation.