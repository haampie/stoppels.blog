---
layout: post
title: "Options in Jacobi-Davidson"
description: "In my previous post I showed the defining property for Jacobi-Davidson, but it turns out that there is a lot of freedom in the remaining parts of the algorithm. For instance in how to initialize it; this post will be about exactly that."
tags: [gsoc, julia]
---

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


