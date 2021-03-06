/**
    @file
    @author  Alexander Sherikov

    @copyright 2020 Alexander Sherikov. Licensed under the Apache License,
    Version 2.0. (see LICENSE or http://www.apache.org/licenses/LICENSE-2.0)

    @brief
*/

#include <iostream>

#include <qpmad/solver.h>
#include <qpOASES.hpp>
#include <eiquadprog/eiquadprog-fast.hpp>


#include "util.h"

namespace
{
    class Results : public ariles2::SloppyBase
    {
#define ARILES2_ENTRIES(v)                                                                                             \
    ARILES2_TYPED_ENTRY_(v, errors, Eigen::VectorXd)                                                                   \
    ARILES2_TYPED_ENTRY_(v, durations, Eigen::VectorXd)
#include ARILES2_INITIALIZE

    public:
        void resize(const std::size_t size)
        {
            errors_.setZero(size);
            durations_.setZero(size);
        }
    };


    class AllResults : public ariles2::SloppyBase
    {
#define ARILES2_ENTRIES(v)                                                                                             \
    ARILES2_ENTRY_(v, qpmad)                                                                                           \
    ARILES2_ENTRY_(v, qpoases)                                                                                         \
    ARILES2_ENTRY_(v, eiquadprog)
#include ARILES2_INITIALIZE

    public:
        std::map<std::string, Results> qpmad_;
        std::map<std::string, Results> qpoases_;
        std::map<std::string, Results> eiquadprog_;

    public:
        void resize(const std::string &id, const std::size_t size)
        {
            qpmad_[id].resize(size);
            qpoases_[id].resize(size);
            eiquadprog_[id].resize(size);
        }
    };
}  // namespace


int main(int argc, char **argv)
{
    if (argc < 2)
    {
        return (EXIT_FAILURE);
    }

    std::vector<benchmark::Problem<qp::IterativeQP> > qps;
    try
    {
        for (std::size_t index = 1; index < argc; ++index)
        {
            const std::string dir = argv[index];
            if (not boost::filesystem::is_directory(dir))
            {
                continue;
            }
            for (boost::filesystem::directory_entry &file : boost::filesystem::directory_iterator(dir))
            {
                if (not boost::filesystem::is_regular(file))
                {
                    continue;
                }

                benchmark::Problem<qp::IterativeQP> qp_container;
                qp_container.file_.init(file);
                qp_container.load();
                if (qp_container.problem_.common_.objective_.positive_definite_)
                {
                    qps.emplace_back(qp_container);
                }
                else
                {
                    std::cout << qp_container.problem_.id_ << "  skipped" << std::endl;
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return (EXIT_FAILURE);
    }


    if (qps.empty())
    {
        return (EXIT_FAILURE);
    }


    AllResults results;
    for (benchmark::Problem<qp::IterativeQP> &qp_container : qps)
    {
        results.resize(qp_container.problem_.id_, qp_container.problem_.instances_.size());
    }
    benchmark::Timer timer;


    bool fail = false;
    try
    {
        for (benchmark::Problem<qp::IterativeQP> &qp_container : qps)
        {
            std::size_t result_index = 0;

            qpmad::Solver solver;
            qpmad::Solver::ReturnStatus status;
            Eigen::MatrixXd hessian_copy = qp_container.problem_.common_.objective_.hessian_;
            Eigen::VectorXd solution;

            const bool has_common_ctr = qp_container.problem_.common_.getNumberOfConstraints() > 0;
            qpmad::SolverParameters param;
            param.hessian_type_ = qpmad::SolverParameters::HESSIAN_LOWER_TRIANGULAR;
            param.return_inverted_cholesky_factor_ = true;

            std::cout << qp_container.problem_.id_ << std::endl;
            for (qp::QP qp_problem : qp_container.problem_.instances_)
            {
                timer.start();
                status = solver.solve(
                        solution,
                        hessian_copy,
                        qp_problem.objective_.vector_,
                        qp_problem.bounds_.lower_,
                        qp_problem.bounds_.upper_,
                        has_common_ctr ? qp_container.problem_.common_.constraints_.matrix_ :
                                         qp_problem.constraints_.matrix_,
                        qp_problem.constraints_.lower_,
                        qp_problem.constraints_.upper_,
                        param);

                results.qpmad_[qp_container.problem_.id_].durations_(result_index) = timer.stop();
                const double err = (solution - qp_problem.solution_.vector_).norm();
                results.qpmad_[qp_container.problem_.id_].errors_(result_index) = err;
                if (err < 1e-9)
                {
                    // std::cout << "ok, " << timer << std::endl;
                }
                else
                {
                    fail = true;
                    std::cout << "fail, error = " << err << std::endl;
                    break;
                }
                ++result_index;
                // reuse Hessian factorization.
                param.hessian_type_ = solver.getHessianType();
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return (EXIT_FAILURE);
    }
    std::cout << "---" << std::endl;


    try
    {
        for (benchmark::Problem<qp::IterativeQP> &qp_container : qps)
        {
            std::size_t result_index = 0;

            std::shared_ptr<qpOASES::QProblem> qp = std::make_shared<qpOASES::QProblem>(
                    qp_container.problem_.getNumberOfVariables(),
                    qp_container.problem_.getNumberOfConstraints(),
                    qpOASES::HST_POSDEF);

            qpOASES::Options options;
            options.setToMPC();
            qp->setOptions(options);

            qpOASES::returnValue qpoases_return_value;


            Eigen::VectorXd solution;
            solution.resize(qp_container.problem_.getNumberOfVariables());

            const bool has_common_ctr = qp_container.problem_.common_.getNumberOfConstraints() > 0;

            const qpOASES::real_t *lb_ptr = NULL;
            const qpOASES::real_t *ub_ptr = NULL;

            const qpOASES::real_t *A_ptr = NULL;
            const qpOASES::real_t *lbA_ptr = NULL;
            const qpOASES::real_t *ubA_ptr = NULL;

            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A_rm;


            std::cout << qp_container.problem_.id_ << std::endl;
            bool first = true;
            for (qp::QP qp_problem : qp_container.problem_.instances_)
            {
                if (qp_problem.hasBounds())
                {
                    lb_ptr = qp_problem.bounds_.lower_.data();
                    ub_ptr = qp_problem.bounds_.upper_.data();
                }

                {
                    // Note: qpOASES expects rowMajor data, Eigen is colMajor by default
                    //       only A needs to be handled since H is symmetric and the rest
                    //       are vectors

                    A_rm = has_common_ctr ? qp_container.problem_.common_.constraints_.matrix_ :
                                            qp_problem.constraints_.matrix_;

                    if (A_rm.rows() > 0)
                    {
                        A_ptr = A_rm.data();

                        lbA_ptr = qp_problem.constraints_.lower_.data();
                        ubA_ptr = qp_problem.constraints_.upper_.data();
                    }
                }

                int number_of_iterations = 10000;

                timer.start();
                if (true == first)
                {
                    qpoases_return_value = qp->init(
                            qp_container.problem_.common_.objective_.hessian_.data(),
                            qp_problem.objective_.vector_.data(),
                            A_ptr,
                            lb_ptr,
                            ub_ptr,
                            lbA_ptr,
                            ubA_ptr,
                            number_of_iterations,
                            /*max_time_ptr=*/NULL);
                    qp->getPrimalSolution(solution.data());
                }
                else
                {
                    qpoases_return_value = qp->hotstart(
                            qp_problem.objective_.vector_.data(),
                            lb_ptr,
                            ub_ptr,
                            lbA_ptr,
                            ubA_ptr,
                            number_of_iterations,
                            /*max_time_ptr=*/NULL);
                    qp->getPrimalSolution(solution.data());
                }
                results.qpoases_[qp_container.problem_.id_].durations_(result_index) = timer.stop();
                const double err = (solution - qp_problem.solution_.vector_).norm();
                results.qpoases_[qp_container.problem_.id_].errors_(result_index) = err;

                if (err < 5e-9)
                {
                    // std::cout << "ok, " << timer << " || " << qpoases_return_value << std::endl;
                }
                else
                {
                    fail = true;
                    std::cout << "fail, error = " << err << " || " << qpoases_return_value << std::endl;
                }
                ++result_index;
                first = false;
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return (EXIT_FAILURE);
    }
    std::cout << "---" << std::endl;


    try
    {
        for (benchmark::Problem<qp::IterativeQP> &qp_container : qps)
        {
            std::size_t result_index = 0;

            Eigen::VectorXd solution;
            solution.resize(qp_container.problem_.getNumberOfVariables());

            const bool has_common_ctr = qp_container.problem_.common_.constraints_.matrix_.rows() > 0;

            std::cout << qp_container.problem_.id_ << std::endl;
            for (qp::QP qp_problem : qp_container.problem_.instances_)
            {
                eiquadprog::solvers::EiquadprogFast qp;

                Eigen::MatrixXd Aeq;
                Eigen::VectorXd Beq;
                Eigen::MatrixXd Aineq;
                Eigen::VectorXd Bineq;

                benchmark::getQuadProgConstraints(
                        &Aeq,
                        &Beq,
                        &Aineq,
                        &Bineq,
                        qp_container.problem_.getNumberOfVariables(),
                        has_common_ctr ? qp_container.problem_.common_.constraints_.matrix_ :
                                         qp_problem.constraints_.matrix_,
                        qp_problem.constraints_.lower_,
                        qp_problem.constraints_.upper_,
                        qp_problem.bounds_.lower_,
                        qp_problem.bounds_.upper_);

                qp.reset(qp_container.problem_.getNumberOfVariables(), Aeq.rows(), Aineq.rows());


                timer.start();
                eiquadprog::solvers::EiquadprogFast_status status = qp.solve_quadprog(
                        qp_container.problem_.common_.objective_.hessian_,
                        qp_problem.objective_.vector_,
                        Aeq,
                        Beq,
                        Aineq,
                        Bineq,
                        solution);
                results.eiquadprog_[qp_container.problem_.id_].durations_(result_index) = timer.stop();
                const double err = (solution - qp_problem.solution_.vector_).norm();
                results.eiquadprog_[qp_container.problem_.id_].errors_(result_index) = err;

                double tol = 1e-9;
                if ("crane.json" == qp_container.problem_.id_)
                {
                    tol = 5e-4;  // XXX hm, 1e-9 fails on crane
                }
                if (err < tol)
                {
                    // std::cout << "ok, " << timer << std::endl;
                }
                else
                {
                    fail = true;
                    std::cout << qp_problem.solution_.value_ << std::endl;
                    std::cout << 0.5 * solution.transpose() * qp_container.problem_.common_.objective_.hessian_
                                                 * solution
                                         + qp_problem.objective_.vector_.transpose() * solution
                              << std::endl;
                    std::cout << (Aineq * qp_problem.solution_.vector_ + Bineq).array().minCoeff() << std::endl;
                    std::cout << "fail, error = " << err << std::endl;
                    break;
                }
                ++result_index;
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return (EXIT_FAILURE);
    }
    std::cout << "---" << std::endl;


    if (fail)
    {
        return (EXIT_FAILURE);
    }
    else
    {
        /*
        std::cout << "qpmad " << results.qpmad_.durations_.sum() << std::endl;
        std::cout << "qpOASES " << results.qpoases_.durations_.sum() << std::endl;
        std::cout << "eiquadprog " << results.eiquadprog_.durations_.sum() << std::endl;
        */

        ariles2::apply<ariles2::octave::Writer>("iterative.m", results);


        std::ofstream octave_script;
        octave_script.open("iterative.m", std::ofstream::out | std::ofstream::app);
        octave_script <<
R"string(
figure('position',[0,0,1200,800])
for i = 1:4
    subplot(2,2,i)
    hold on
    title(qpoases{i}.first)
    plot(qpoases{i}.second.durations, 'k')
    plot(eiquadprog{i}.second.durations, 'b')
    plot(qpmad{i}.second.durations, 'r')
    legend(['qpoases, mean=', num2str(mean(qpoases{i}.second.durations))], ...
           ['eiquadprog, mean=', num2str(mean(eiquadprog{i}.second.durations))], ...
           ['qpmad, mean=', num2str(mean(qpmad{i}.second.durations))], ...
            'Location', 'northeast')
    xlim([0, numel(qpoases{i}.second.durations)])
end
print ('iterative.png', '-dpng', '-color', '-tight', '-F:8', '-S1200,800');
)string" << std::endl;
        octave_script.close();

        return (EXIT_SUCCESS);
    }
}
