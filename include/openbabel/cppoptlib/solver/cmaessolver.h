// CppNumericalSolver
#ifndef CMAES_H_
#define CMAES_H_

#include <random>
#include <Eigen/Dense>
#include "isolver.h"

namespace cppoptlib {

/**
 * @brief Covariance Matrix Adaptation
 */
template<typename TProblem>
class CMAesSolver : public ISolver<TProblem, 1> {
public:
    using Super = ISolver<TProblem, 1>;
    using typename Super::Scalar;
    using typename Super::TVector;
    using typename Super::THessian;
    using typename Super::TCriteria;
    using TMatrix = Eigen::Matrix<Scalar, TProblem::Dim, Eigen::Dynamic>;
    using TVarVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

protected:
    std::mt19937 gen;
    Scalar m_stepSize;
    /*
     * @brief Create a vector sampled from a normal distribution
     *
     */
    TVector normDist(const int n) {
        TVector sample = TVector::Zero(n);
        std::normal_distribution<Scalar> d(0, 1);
        for (int i = 0; i < n; ++i) {
          sample[i] = d(gen);
        }
        return sample;
    }

    std::vector<size_t> index_partial_sort(const TVarVector &x, Eigen::ArrayXd::Index N)
    {
        eigen_assert(x.size() >= N);
        std::vector<size_t> allIndices(x.size()), indices(N);
        for(size_t i = 0; i < allIndices.size(); i++) {
            allIndices[i] = i;
        }
        partial_sort(allIndices.begin(), allIndices.begin() + N, allIndices.end(), [&x](size_t i1, size_t i2) { return x[i1] < x[i2]; });
        //std::cout << "SORTED: ";
        for (Eigen::ArrayXd::Index i = 0; i < N; i++) {
            indices[i] = allIndices[i];
            //std::cout << indices[i] << ",";
        }
        //std::cout << std::endl;
        return indices;
    }

public:
    CMAesSolver() : gen((std::random_device())()) {
        m_stepSize = 0.5;
        // Set some sensible defaults for the stop criteria
        Super::m_stop.iterations = 1e5;
        Super::m_stop.gradNorm = 0; // Switch this off
        Super::m_stop.condition = 1e14;
        Super::m_stop.xDelta = 1e-7;
        Super::m_stop.fDelta = 1e-9;
    }

    void minimize(TProblem &objFunc, TVector &x0) {
        TVector var0 = TVector::Ones(x0.rows());
        this->minimize(objFunc, x0, var0);
    }

    /**
    * @brief minimize
    * @details [long description]
    *
    * @param objFunc [description]
    */
    void minimize(TProblem &objFunc, TVector &x0, const TVector &var0) {
        const int n = x0.rows();
        eigen_assert(x0.rows() == var0.rows());
        int la = ceil(4 + round(3 * log(n)));
        const int mu = floor(la / 2);
        TVarVector w = TVarVector::Zero(mu);
        for (int i = 0; i < mu; ++i) {
          w[i] = log(mu+1/2)-log(i+1);
        }
        w /= w.sum();
        const Scalar mu_eff = (w.sum()*w.sum()) / w.dot(w);
        la = std::max<int>(16, la); // Increase to 16 samples for very small populations, but AFTER calcaulting mu_eff
        const Scalar cc     = (4. + mu_eff / n) / (n + 4. + 2.*mu_eff/n);
        const Scalar cs     = (mu_eff + 2.) / (n + mu_eff + 5.);
        const Scalar c1     = 2. / (pow(n + 1.3, 2.) + mu_eff);
        const Scalar cmu    = std::min(1. - c1, 2.*(mu_eff - 2. + 1./mu_eff) / (pow(n+2, 2.) + mu_eff));
        const Scalar ds     = 1. + cs + 2.*std::max(0., sqrt((mu_eff - 1.) / (n + 1.)) - 1.);
        const Scalar chi    = sqrt(n) * (1. - 1./(4.*n) + 1./(21.*n*n));
        const Scalar hsig_thr = (1.4 + 2 / (n + 1.)) * chi;

        TVector pc = TVector::Zero(n);
        TVector ps = TVector::Zero(n);
        THessian B = THessian::Identity(n, n);
        THessian D = THessian::Identity(n, n);
        THessian C = B*D*(B*D).transpose();

        Scalar sigma = m_stepSize;

        TVector xmean = x0;
        TVector zmean = TVector::Zero(n);
        TMatrix arz(n, la);
        TMatrix arx(n, la);
        TVarVector costs(la);
        Scalar prevCost = objFunc.value(x0);
        // CMA-ES Main Loop
        size_t eigen_last_eval = 0;
        size_t eigen_next_eval = std::max<Scalar>(1, 1/(10*n*(c1+cmu)));
        this->m_current.reset();
        if (Super::m_debug >= DebugLevel::Low) {
            std::cout << "CMA-ES Initial Config" << std::endl;
            std::cout << "n " << n << " la " << la << " mu " << mu << " mu_eff " << mu_eff << " sigma " << sigma << std::endl;
            std::cout << "cs " << cs << " ds " << ds << " chi " << chi << " cc " << cc << " c1 " << c1 << " cmu " << cmu
                      << " hsig_thr " << hsig_thr << std::endl;
            std::cout << "C" << std::endl << C << std::endl;
            std::cout << "Hessian will be updated every " << eigen_next_eval << " iterations." << std::endl;
            std::cout << "Iteration: " << this->m_current.iterations << " best cost " << prevCost << " sigma " << sigma
                      << " cond " << this->m_current.condition << " xmean " << x0.transpose() << std::endl;
        }
        do {
            for (int k = 0; k < la; ++k) {
              arz.col(k) = normDist(n);
              arx.col(k) = xmean + sigma * B*D*arz.col(k);
              costs[k] = objFunc(arx.col(k));
            }
            if (Super::m_debug >= DebugLevel::High) {
                std::cout << "arz" << std::endl << arz << std::endl;
                std::cout << "arx" << std::endl << arx << std::endl;
                std::cout << "costs " << costs.transpose() << std::endl;
            }
            std::vector<size_t> indices = index_partial_sort(costs, mu);
            xmean = TVector::Zero(n);
            zmean = TVector::Zero(n);
            for (int k = 0; k < mu; k++) {
              zmean += w[k]*arz.col(indices[k]);
              xmean += w[k]*arx.col(indices[k]);
            }
            // Update evolution paths
            ps = (1. - cs)*ps + sqrt(cs*(2. - cs)*mu_eff) * B*zmean;
            Scalar hsig = (ps.norm()/sqrt(pow(1 - (1. - cs), (2 * (this->m_current.iterations + 1)))) < hsig_thr) ? 1.0 : 0.0;
            pc = (1. - cc)*pc + hsig*sqrt(cc*(2. - cc)*mu_eff)*(B*D*zmean);
            // Adapt covariance matrix
            C = (1 - c1 - cmu)*C
                + c1*(pc*pc.transpose() + (1 - hsig)*cc*(2 - cc)*C);
            for (int k = 0; k < mu; ++k) {
                TVector temp = B*D*arz.col(k);
                C += cmu*w(k)*temp*temp.transpose();
            }
            sigma = sigma * exp((cs/ds) * ((ps).norm()/chi - 1.));

            ++Super::m_current.iterations;
            if ((Super::m_current.iterations - eigen_last_eval) == eigen_next_eval) {
                // Update B and D
                eigen_last_eval = Super::m_current.iterations;
                Eigen::SelfAdjointEigenSolver<THessian> eigenSolver(C);
                B = eigenSolver.eigenvectors();
                D.diagonal() = eigenSolver.eigenvalues().array().sqrt();
                if (Super::m_debug >= DebugLevel::High) {
                    std::cout << "Updated hessian." << std::endl;
                    std::cout << "C" << std::endl << C << std::endl;
                    std::cout << "B" << std::endl << B << std::endl;
                    std::cout << "D" << std::endl << D << std::endl;
                }
            }
            Super::m_current.condition = D.diagonal().maxCoeff() / D.diagonal().minCoeff();
            Super::m_current.xDelta = (xmean - x0).norm();
            x0 = xmean;
            Super::m_current.fDelta = fabs(costs[indices[0]] - prevCost);
            prevCost = costs[indices[0]];
            if (Super::m_debug >= DebugLevel::Low) {
                std::cout << "Iteration: " << this->m_current.iterations << " best cost " << costs[indices[0]] << " sigma " << sigma
                          << " cond " << this->m_current.condition << " xmean " << xmean.transpose() << std::endl;
            }
            if (fabs(costs[indices[0]] - costs[indices[mu-1]]) < 1e-6) {
                sigma = sigma * exp(0.2+cs/ds);
                if (Super::m_debug >= DebugLevel::Low) {
                    std::cout << "Flat fitness " << costs[indices[0]] << " " << costs[indices[mu-1]] << std::endl;
                }
            }
            Super::m_status = checkConvergence(this->m_stop, this->m_current);
        } while (objFunc.callback(this->m_current, x0) && (this->m_status == Status::Continue));
        // Return the best evaluated solution
        x0 = xmean;
        m_stepSize = sigma;
        if (Super::m_debug >= DebugLevel::Low) {
            std::cout << "Stop" << std::endl;
            this->m_stop.print(std::cout);
            std::cout << "Current" << std::endl;
            this->m_current.print(std::cout);
            std::cout << "Reason: " << Super::m_status << std::endl;
        }
    }
};

} /* namespace cppoptlib */

#endif /* CMAES_H_ */
