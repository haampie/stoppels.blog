---
layout: post
title: "Jacobi-Davidson example"
description: "Let's take a look how the Jacobi-Davidson algorithm targets interior eigenvalues of a Poisson matrix."
tags: [gsoc, julia]
---

I've put Jacobi-Davidson to the test on the simple Poisson matrix from the previous post: it has {% katex %}\begin{bmatrix}-1 & 2 & -1\end{bmatrix}{% endkatex %} on the diagonals and is {% katex %}100 \times 100{% endkatex %}. Its eigenvalues are all real and located in the interval {% katex %}[0, 4]{% endkatex %}.

## Using the exact solver
When solving the correction equation with the "exact" solver, we should be able to see quadratic convergence.

<figure class="full_width_fig">
    {% img jd/convergence.svg %}
</figure>

