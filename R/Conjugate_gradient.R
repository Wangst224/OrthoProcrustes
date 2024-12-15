conjugate_gradient = function(theta, A, B, V0, d, tol, max_iter) {

    optim_res = optim(theta, Obj_func, gr = Obj_func_gradient,
                      method = "BFGS", A = A, B = B, V0 = V0, d = d,
                      control = list(maxit = max_iter, abstol = tol, reltol = tol))

    G = generate_G(optim_res$par, sin(optim_res$par), cos(optim_res$par))

    return(list(
        V = V0 %*% G,
        U = V0 %*% G %*% diag(d),
        value = optim_res$value,
        convergence = optim_res$convergence
    ))
}
