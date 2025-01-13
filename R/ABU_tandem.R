#' Orthogonal Procrustes Problem
#'
#' Orthogonal but not necessarily orthonormal Procrustes problem. Numerically solved using tandem algorithm through conjugate gradient.
#'
#' @details This functions solves the orthogonal Procrustes problem.
#' Given matrices A and B, the orthogonal procrustes seeks a matrix U, which has orthogonal columns (U'U = D^2), such that ||A - B U||^2 is minimized.
#' The solution is found by a tandem algorithm using conjugate gradient method.
#' The current implementation assumes that the diagonal matrix D is known.
#'
#' @param A An m by p matrix. The target matrix.
#' @param B An m by p matrix. The matrix to be aligned.
#' @param d A vector of length p. The diagonal matrix D.
#'
#' @return TBD
#'
#' @export
#'
#' @examples
#' ABU_tandem(A, B, d)
ABU_tandem = function(A, B, d, tol = 1e-6, max_iter = 1000) {
    # Given matrices A and B, find the V with orthonormal columns and d that minimize
    # rho = ||A - B V diag(d) ||^2

    # Minimization stops when the relative improvement in rho from one
    # iteration to the next is relatively less than tol or itmax iterations are
    # exceeded.
    #
    # Approximations of V and d are found in tandem.  When d is fixed, V is
    # found by a conjugate gradient method; when V is fixed d is found by
    # differentiating rho wrt to d and setting the result equal to zero.
    #
    # For current implementation, d is assumed to be known.

    # Check inputs
    if (any(dim(A) != dim(B))) {
        stop("A and B must have the same dimensions")
    }

    m = dim(A)[1]
    p = dim(A)[2]

    if (length(d) != p) {
        stop("d must have length equal to the number of columns of A and B")
    }

    # Initialization
    svd_init = svd(t(B) %*% A)
    V_init = svd_init$u %*% t(svd_init$v)
    theta_init = rep(0, p*(p-1)/2)

    # For fixed V, find d through differentiation. Not implemented yet.
    # Therefore there is no iteration in this function now.

    # For fixed d, find V through conjugate gradient method.
    cg_res = conjugate_gradient(theta_init, A, B, V_init, d, tol, max_iter)

    if (cg$res != 0) warning("Conjugate gradient did not converge")

    return(cg_res)
}
