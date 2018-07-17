# Mike's Notes: July 17 2018
- Finally getting back to paper.
- Paper goal: present a simple coupled ODE model of ribosome loading and a novel iterative solution for its equilibrium state.
- This new goal would make the PDE version seem unnecessary, but that may be due to not thinking about this much.
- Would be useful to illustrate model behavior graphically (extreme cases as well as reasonable ones).
- The above point suggests it would be worthwhile non-dimensionalizing the main equations by $\tau_0$.
- May want to use subscripts for \kappa and $\tau$ instead of $(i)$

## Questions for Nate
     - How general is the approach?
     Specifically, if I used generic functions for $\mu(i)$ and $\tau(i)$ would our approach still work?
    - If the answer above is, "No," would any linear function work?
    - Seems like the $s_i$ function is equal to $\kappa(i-1)$, if so we should just use $k$ instead of introducing another function.
    - Similarly, $p_i$ seems equivalent to $\tau(i)$
