# import LsqFit: estimate_covar, LsqFitResult

# # NOTE: Fixes LsqFit's estimate_covar, which crashes with a LAPACK exception
# # when R is not invertible
# function estimate_covar(fit::LsqFitResult)
#     # computes covariance matrix of fit parameters
#     J = fit.jacobian

#     if isempty(fit.wt)
#         r = fit.resid

#         # compute the covariance matrix from the QR decomposition
#         Q, R = qr(J)
        
#         if det(R) == 0.0
#             covar = similar(R)
#             covar .= +Inf
#         else
#             Rinv = inv(R)
#             covar = Rinv*Rinv'*mse(fit)
#         end

#     else
#         covar = inv(J'*J)
#     end

#     return covar
# end