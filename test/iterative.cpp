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
    class Results : public ariles2::DefaultBase
    {
#define ARILES2_ENTRIES(v)                                                                                             \
    ARILES2_TYPED_ENTRY_(v, errors, Eigen::VectorXd)                                                                   \
    ARILES2_TYPED_ENTRY_(v, durations, Eigen::VectorXd)
#include ARILES2_INITIALIZE

    public:
        Results(const std::size_t size)
        {
            resize(size);
        }

        void resize(const std::size_t size)
        {
            errors_.setZero(size);
            durations_.setZero(size);
        }
    };


    class AllResults : public ariles2::DefaultBase
    {
#define ARILES2_ENTRIES(v)                                                                                             \
    ARILES2_TYPED_ENTRY_(v, qpmad, Results)                                                                            \
    ARILES2_TYPED_ENTRY_(v, qpoases, Results)                                                                          \
    ARILES2_TYPED_ENTRY_(v, eiquadprog, Results)
#include ARILES2_INITIALIZE

    public:
        AllResults(const std::size_t size) : qpmad_(size), qpoases_(size), eiquadprog_(size)
        {
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


    std::size_t results_size = 0;
    for (benchmark::Problem<qp::IterativeQP> &qp_container : qps)
    {
        results_size += qp_container.problem_.instances_.size();
    }
    AllResults results(results_size);
    bool fail = false;

    try
    {
        benchmark::Timer timer;
        std::size_t result_index = 0;

        for (benchmark::Problem<qp::IterativeQP> &qp_container : qps)
        {
            qpmad::Solver solver;
            qpmad::Solver::ReturnStatus status;
            Eigen::MatrixXd hessian_copy = qp_container.problem_.common_.objective_.hessian_;
            Eigen::VectorXd solution;

            const bool has_common_ctr = qp_container.problem_.common_.constraints_.matrix_.rows() > 0;
            qpmad::SolverParameters param;
            param.hessian_type_ = qpmad::SolverParameters::HESSIAN_LOWER_TRIANGULAR;

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
                results.qpmad_.durations_(result_index) = timer.stop();

                results.qpmad_.errors_(result_index) = (solution - qp_problem.solution_.vector_).norm();
                if (results.qpmad_.errors_(result_index) < 1e-9)
                {
                    std::cout << "ok, " << timer << std::endl;
                }
                else
                {
                    fail = true;
                    std::cout << "fail, error = " << results.qpmad_.errors_(result_index) << std::endl;
                    break;
                }
                ++result_index;
                // reuse Hessian factorization.
                param.hessian_type_ = qpmad::SolverParameters::HESSIAN_CHOLESKY_FACTOR;
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return (EXIT_FAILURE);
    }
    std::cout << "---" << std::endl;


    /*
    try
    {
        benchmark::Timer timer;
        std::size_t result_index = 0;

        for (benchmark::Problem<qp::IterativeQP> &qp_container : qps)
        {
            std::shared_ptr<qpOASES::QProblem> qp = std::make_shared<qpOASES::QProblem>(
                    qp_container.problem_.getNumberOfVariables(),
                    qp_container.problem_.getNumberOfConstraints(),
                    qpOASES::HST_POSDEF);


            qpOASES::real_t *max_time_ptr = NULL;

            const qpOASES::real_t *lb_ptr = NULL;
            const qpOASES::real_t *ub_ptr = NULL;
            if (qp_container.problem_.hasBounds())
            {
                lb_ptr = qp_container.problem_.bounds_.lower_.data();
                ub_ptr = qp_container.problem_.bounds_.upper_.data();
            }

            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A_rm;
            const qpOASES::real_t *A_ptr = NULL;
            const qpOASES::real_t *lbA_ptr = NULL;
            const qpOASES::real_t *ubA_ptr = NULL;
            if (qp_container.problem_.getNumberOfConstraints() > 0)
            {
                // Note: qpOASES expects rowMajor data, Eigen is colMajor by default
                //       only A needs to be handled since H is symmetric and the rest
                //       are vectors

                A_rm = qp_container.problem_.constraints_.matrix_;
                A_ptr = A_rm.data();

                lbA_ptr = qp_container.problem_.constraints_.lower_.data();
                ubA_ptr = qp_container.problem_.constraints_.upper_.data();
            }

            const qpOASES::real_t *solution_guess_ptr = NULL;

            qpOASES::Bounds *active_set_bounds_ptr = NULL;
            qpOASES::Constraints *active_set_constraints_ptr = NULL;

            qpOASES::returnValue qpoases_return_value;

            Eigen::VectorXd solution;
            solution.resize(qp_container.problem_.getNumberOfVariables());

            int number_of_iterations = 10000;

            std::cout << qp_container.problem_.id_ << "  ";
            timer.start();
            qpoases_return_value = qp->init(
                    qp_container.problem_.objective_.hessian_.data(),
                    qp_container.problem_.objective_.vector_.data(),
                    A_ptr,
                    lb_ptr,
                    ub_ptr,
                    lbA_ptr,
                    ubA_ptr,
                    number_of_iterations,
                    max_time_ptr,
                    solution_guess_ptr,
                    NULL,  // Optimal dual solution vector
                    active_set_bounds_ptr,
                    active_set_constraints_ptr);
            qp->getPrimalSolution(solution.data());
            results.qpoases_.durations_(result_index) = timer.stop();

            results.qpoases_.errors_(result_index) = (solution - qp_container.problem_.solution_.vector_).norm();
            if (results.qpoases_.errors_(result_index) < 1e-9)
            {
                std::cout << "ok, " << timer << std::endl;
            }
            else
            {
                fail = true;
                std::cout << "fail, error = " << results.qpoases_.errors_(result_index) << std::endl;
            }
            ++result_index;
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
        benchmark::Timer timer;
        std::size_t result_index = 0;

        for (benchmark::Problem<qp::IterativeQP> &qp_container : qps)
        {
            eiquadprog::solvers::EiquadprogFast qp;

            std::size_t total_ctr_num = 0;
            Eigen::MatrixXd constraints;
            Eigen::VectorXd constraints_lb;
            Eigen::VectorXd constraints_ub;
            if (qp_container.problem_.hasBounds())
            {
                total_ctr_num =
                        qp_container.problem_.getNumberOfVariables() + qp_container.problem_.getNumberOfConstraints();

                constraints.resize(total_ctr_num, qp_container.problem_.getNumberOfVariables());
                constraints_lb.resize(total_ctr_num);
                constraints_ub.resize(total_ctr_num);


                constraints.topRows(qp_container.problem_.getNumberOfConstraints()) =
                        qp_container.problem_.constraints_.matrix_;
                constraints_lb.head(qp_container.problem_.getNumberOfConstraints()) =
                        qp_container.problem_.constraints_.lower_;
                constraints_ub.head(qp_container.problem_.getNumberOfConstraints()) =
                        qp_container.problem_.constraints_.upper_;

                constraints.bottomRows(qp_container.problem_.getNumberOfVariables()).setIdentity();
                constraints_lb.tail(qp_container.problem_.getNumberOfVariables()) =
                        qp_container.problem_.bounds_.lower_;
                constraints_ub.tail(qp_container.problem_.getNumberOfVariables()) =
                        qp_container.problem_.bounds_.upper_;
            }
            else
            {
                total_ctr_num = qp_container.problem_.getNumberOfConstraints();

                constraints = qp_container.problem_.constraints_.matrix_;
                constraints_lb = qp_container.problem_.constraints_.lower_;
                constraints_ub = qp_container.problem_.constraints_.upper_;
            }


            std::vector<std::size_t> equalities;
            std::vector<std::size_t> inequalities;
            for (std::size_t i = 0; i < total_ctr_num; ++i)
            {
                if (std::abs(constraints_lb(i) - constraints_ub(i)) < 1e-9)
                {
                    equalities.push_back(i);
                }
                else
                {
                    inequalities.push_back(i);
                }
            }

            const std::size_t num_eq = equalities.size();
            const std::size_t num_ineq = total_ctr_num - equalities.size();


            Eigen::MatrixXd Aeq;
            Eigen::VectorXd Beq;
            Eigen::MatrixXd Aineq;
            Eigen::VectorXd Bineq;

            Aeq.resize(num_eq, qp_container.problem_.getNumberOfVariables());
            Beq.resize(num_eq);
            for (std::size_t i = 0; i < num_eq; ++i)
            {
                Aeq.row(i) = constraints.row(equalities[i]);
                Beq(i) = -constraints_lb(equalities[i]);
            }

            Aineq.resize(num_ineq * 2, qp_container.problem_.getNumberOfVariables());
            Bineq.resize(num_ineq * 2);
            for (std::size_t i = 0; i < num_ineq; ++i)
            {
                Aineq.row(i) = constraints.row(inequalities[i]);
                Bineq(i) = -constraints_lb(inequalities[i]);
                Aineq.row(num_ineq + i) = -constraints.row(inequalities[i]);
                Bineq(num_ineq + i) = constraints_ub(inequalities[i]);
            }

            qp.reset(qp_container.problem_.getNumberOfVariables(), num_eq, num_ineq * 2);


            Eigen::VectorXd solution;
            solution.resize(qp_container.problem_.getNumberOfVariables());

            std::cout << qp_container.problem_.id_ << "  ";
            timer.start();
            eiquadprog::solvers::EiquadprogFast_status status = qp.solve_quadprog(
                    qp_container.problem_.objective_.hessian_,
                    qp_container.problem_.objective_.vector_,
                    Aeq,
                    Beq,
                    Aineq,
                    Bineq,
                    solution);
            results.eiquadprog_.durations_(result_index) = timer.stop();

            results.eiquadprog_.errors_(result_index) = (solution - qp_container.problem_.solution_.vector_).norm();
            if (results.eiquadprog_.errors_(result_index) < 1e-9)
            {
                std::cout << "ok, " << timer << std::endl;
            }
            else
            {
                fail = true;
                std::cout << "fail, error = " << results.eiquadprog_.errors_(result_index) << std::endl;
            }
            ++result_index;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return (EXIT_FAILURE);
    }
    std::cout << "---" << std::endl;
    */


    if (fail)
    {
        return (EXIT_FAILURE);
    }
    else
    {
        std::cout << "qpmad " << results.qpmad_.durations_.mean() << std::endl;
        std::cout << "qpOASES " << results.qpoases_.durations_.mean() << std::endl;
        std::cout << "eiquadprog " << results.eiquadprog_.durations_.mean() << std::endl;

        // ariles2::apply<ariles2::octave::Writer>(std::cout, results);
        ariles2::apply<ariles2::octave::Writer>("oneshot.m", results);
        return (EXIT_SUCCESS);
    }
}
