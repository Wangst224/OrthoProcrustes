generate_G = function(theta, sn, cs) {
    n = length(theta)
    p = (1 + sqrt(1 + 8*n)) / 2 # n = p(p-1)/2

    j = 1
    G = diag(p)

    for (k in 2:p) {
        for (i in 1:(k-1)) {
            G[, c(i, k)] = G[, c(i, k)] %*% matrix(c(cs[j], sn[j], -sn[j], cs[j]), nrow = 2)
            j = j + 1
        }
    }

    return(G)
}

Obj_func = function(theta, A, B, V0, d) {
    # Evaluate the objective function ||A-BU||^2 as a function of theta

    sn = sin(theta)
    cs = cos(theta)
    G = generate_G(theta, sn, cs)

    V = V0 %*% G
    U = V %*% diag(d)
    diff = A - B %*% U

    return(sum(diff^2))
}

Obj_func_gradient = function(theta, A, B, V0, d) {
    # Evaluate the gradient of the objective function ||A-BU||^2 as a function of theta
    n = length(theta)
    p = (1 + sqrt(1 + 8*n)) / 2 # n = p(p-1)/2

    sn = sin(theta)
    cs = cos(theta)
    G = generate_G(theta, sn, cs)

    Hd = t(A) %*% B %*% V0
    Tl = G %*% diag(d)

    X = B %*% V0
    Y = G %*% diag(d)
    Fr = X %*% Y
    N = nrow(Fr)
    M = ncol(Fr)
    dFr = matrix(NA, N, M)

    j = 1
    grad = rep(NA, n)
    for (k in 2:p) {
        for (i in 1:(k-1)) {

            # Contribution from -2Tr{A'BU}

            ## Remove G_j from Tl
            Tl[c(i,k),] = matrix(c(cs[j], -sn[j], sn[j], cs[j]), nrow = 2) %*% Tl[c(i,k),]

            grad_pt1 = -2 * sum(-Hd[,i] * Tl[i,] * sn[j] +
                                 Hd[,k] * Tl[i,] * cs[j] -
                                 Hd[,i] * Tl[k,] * cs[j] -
                                 Hd[,k] * Tl[k,] * sn[j])

            ## Add G_j to Hd
            Hd[,c(i,k)] = Hd[,c(i,k)] %*% matrix(c(cs[j], sn[j], -sn[j], cs[j]), nrow = 2)

            # Contribution from Tr{D'V'B'BVD}

            ## Remove G_j from Y
            Y[c(i,k),] = matrix(c(cs[j], -sn[j], sn[j], cs[j]), nrow = 2) %*% Y[c(i,k),]

            ## Compute dFr
            for (s in 1:N) {
                for (t in 1:M) {
                    dFr[s, t] = -X[s, i] * Y[i, t] * sn[j] +
                                X[s, k] * Y[i, t] * cs[j] -
                                X[s, i] * Y[k, t] * cs[j] -
                                X[s, k] * Y[k, t] * sn[j]
                }
            }

            grad_pt2 = 2*sum(dFr * Fr)

            ## Add G_j to X
            X[,c(i,k)] = X[,c(i,k)] %*% matrix(c(cs[j], sn[j], -sn[j], cs[j]), nrow = 2)

            grad[j] = grad_pt1 + grad_pt2
            j = j + 1
        }
    }

    return(grad)
}
