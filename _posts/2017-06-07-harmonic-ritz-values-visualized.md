---
layout: post
title: "Harmonic Ritz values visualized"
description: "A hopefully insightful visualizion of harmonic Ritz values versus ordinary Ritz values."
tags: [gsoc, julia]
---

I thought it would be nice to clarify the previous post with some plots. Take a simple Poisson matrix {% katex %}A{% endkatex %} of size {% katex %}n \times n{% endkatex %} with entries {% katex %}\begin{bmatrix}-1 & 2 & -1\end{bmatrix}{% endkatex %} on the diagonals. We know its eigenvalues are all concentrated in the interval {% katex %}[0, 4]{% endkatex %} on the real line. 

As a preliminary note: I'm not exploiting the symmetry of the matrix in my code. I want to test it on general matrices as well.

At any rate, in Julia we generate this Poisson matrix as

{% highlight julia %}
my_matrix(n) = SymTridiagonal(fill(2.0, n), fill(-1.0, n - 1))
{% endhighlight %}

Sometimes Jacobi-Davidson's search space is initialized with the search space of Arnoldi's method, namely a Krylov subspace with a random initial vector {% katex %}v_0{% endkatex %}:

{% katex display %}
\mathcal{V} = \mathcal{K}_k(A, v_0) = span\{v_0, Av_0, \cdots, A^{k-1}v_0\}.
{% endkatex %}

The well-known Arnoldi decomposition {% katex %}AV_k = V_{k+1}H_{k+1,k}{%endkatex%} gives us an orthonormal matrix {% katex %}V_k{% endkatex %} whose {% katex %}k{% endkatex %} columns span the search space, and a Hessenberg matrix of size {% katex %}(k + 1) \times k {% endkatex %}. In Julia we can generate these using

{% highlight julia %}
function init(n, dimension)
    V = zeros(n, dimension + 1)
    H = zeros(dimension + 1, dimension)
    v0 = rand(n)
    V[:, 1] = v0 / norm(v0)

    V, H
end

function expand!(A, V, k)
    V[:, k + 1] = A * @view(V[:, k])
end

function orthogonalize!(V, H, k)
    for i = 1 : k
        H[i, k] = dot(@view(V[:, i]), @view(V[:, k + 1]))
        V[:, k + 1] -= H[i, k] * @view(V[:, i])
    end
    H[k + 1, k] = norm(@view(V[:, k + 1]))
    V[:, k + 1] /= H[k + 1, k]
end

function arnoldi(A, m)
    n = size(A, 1)
    V, H = init(n, m)

    for k = 1 : m
        expand!(A, V, k)
        orthogonalize!(V, H, k)
    end

    V, H
end
{% endhighlight %}

Now the ordinary Ritz values {% katex %}\theta{% endkatex %} are obtained by solving the small eigenvalue problem

{% katex display %}
V_k^TAV_ky = H_{k,k}y = \theta y.
{% endkatex %}

Here the square matrix {% katex %}H_{k,k}{% endkatex %} is {% katex %}H_{k+1,k}{% endkatex %} with its last row removed. We can compute and plot the Ritz values for each iteration as follows:

{% highlight julia %}
using Plots

ritz(H) = eigvals(H)

function ritz_example(; n = 60, m = 30)
    A = my_matrix(n)
    V, H = arnoldi(A, m)
    p = plot(;legend = :none)

    for k = 1 : m
        ritz_values = ritz(@view(H[1 : k, 1 : k]))
        scatter!(real(ritz_values), k * ones(ritz_values), marker = (:x, :green))
    end
    p
end
{% endhighlight %}

What we get is something like

<figure class="full_width_fig">
    {% img ritz/ritz.svg %}
</figure>

On the vertical axis we have the iteration (or size of our search space), and horizontally we have the real line. The black stars are the actual eigenvalues, while the green crosses the Ritz values. Now if we were to target interior eigenvalues near {% katex %}\tau = 3{% endkatex %}, we can clearly see that many Ritz values start out near this target, but actually converge away from it. This is undesirable for restarts.

Enter harmonic Ritz values. As shown in the previous post harmonic Ritz values satisfy the Petrov-Galerkin condition

{% katex display %}
(A - \tau I)z - \theta z \perp (A - \tau I)V_k \text{ for } z = V_k y.
{% endkatex %}

Now it's not hard to verify we can write {% katex display %}(A - \tau I)V_k = V_{k+1}H_{k+1,k}^\tau{% endkatex %} for {% katex %}H_{k+1,k}^\tau := H_{k+1,k} - \tau I{% endkatex %}. Note that the Hessenberg matrix is not square, but by abuse of notation I mean that {% katex %}\tau{% endkatex %} must be substracted from the diagonal of the Hessenberg matrix. Some more work shows that the Petrov-Galerkin condition is equivalent to solving the generalized eigenvalue problem

{% katex display %}
(H^\tau_{k+1,k})^* H^\tau_{k+1,k} y = \theta (H^\tau_{k,k})^* y.
{% endkatex %}

This can be implemented using:

{% highlight julia %}
function harmonic_ritz(H, f, τ)
    Hτ = H - UniformScaling(τ)
    A = Hτ' * Hτ
    A[end, end] += conj(f) * f
    eigvals(A, Hτ') + τ
end

function harmonic_ritz_example(; n = 60, m = 30)
    A = my_matrix(n)
    V, H = arnoldi(A, m)
    p = plot(;legend = :none)

    for k = 1 : m
        Hk = @view(H[1 : k, 1 : k]
        harmonic = harmonic_ritz(Hk, H[k + 1, k], τ))
        scatter!(real(harmonic), k * ones(harmonic), marker = (:x, :green))
    end
    p
end
{% endhighlight %}

And produces harmonic Ritz values like this

<figure class="full_width_fig">
    {% img ritz/harmonic.svg %}
</figure>

The point is that harmonic Ritz values come close to our target {% katex %}\tau{% endkatex %} only when they converge to an eigenvalue near {% katex %}\tau{% endkatex %}.

For restarts this means that retaining harmonic Ritz vectors with harmonic Ritz values close to {% katex %}\tau{% endkatex %} will be a good strategy. In sharp contrast with ordinary Ritz values, where lots of Ritz values near our target are useless to retain, since they are on their way converging to exterior eigenvalues.