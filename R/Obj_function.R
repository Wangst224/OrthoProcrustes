generate_G = function(theta, sn, cs) {
    n = length(theta)
    p = (1 + sqrt(1 + 8*n)) / 2
    j = 1

    G = diag(p)

    for (k in 2:p) {
        for (i in 1:(k-1)) {
            G[, c(i, k)] = G[, c(i, k)] %*% matrix(c(cs[j], -sn[j], sn[j], cs[j]), nrow=2)
            j = j + 1
        }
    }

    return(G)
}

Obj_func = function(theta, A, B, V0, d) {
    # Evaluate the objective function ||A-BU||^2 as a function of theta
    # U = DV
    # V = V0 * G(theta)

    sn = sin(theta)
    cs = cos(theta)
    G = generate_G(theta, sn, cs)

    V = V0 %*% G
    U = V %*% diag(d)
    diff = A - B %*% U

    return(sum(diff^2))
}
