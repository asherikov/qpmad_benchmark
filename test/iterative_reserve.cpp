/**
    @file
    @author  Alexander Sherikov

    @copyright 2020 Alexander Sherikov. Licensed under the Apache License,
    Version 2.0. (see LICENSE or http://www.apache.org/licenses/LICENSE-2.0)

    @brief
*/

#define EIGEN_RUNTIME_NO_MALLOC

#include <iostream>

#include <qpmad/solver.h>

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
    ARILES2_ENTRY_(v, qpmad_reserve)
#include ARILES2_INITIALIZE

    public:
        std::map<std::string, Results> qpmad_;
        std::map<std::string, Results> qpmad_reserve_;

    public:
        void resize(const std::string &id, const std::size_t size)
        {
            qpmad_[id].resize(size);
            qpmad_reserve_[id].resize(size);
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
                solver.reserve(
                        hessian_copy.rows(),
                        qp_problem.bounds_.lower_.rows(),
                        has_common_ctr ? qp_container.problem_.common_.constraints_.matrix_.rows():
                                         qp_problem.constraints_.matrix_.rows());

                timer.start();
                Eigen::internal::set_is_malloc_allowed(false);
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
                Eigen::internal::set_is_malloc_allowed(true);
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

                results.qpmad_reserve_[qp_container.problem_.id_].durations_(result_index) = timer.stop();
                const double err = (solution - qp_problem.solution_.vector_).norm();
                results.qpmad_reserve_[qp_container.problem_.id_].errors_(result_index) = err;
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




    if (fail)
    {
        return (EXIT_FAILURE);
    }
    else
    {
        ariles2::apply<ariles2::octave::Writer>("iterative_reserve.m", results);


        std::ofstream octave_script;
        octave_script.open("iterative_reserve.m", std::ofstream::out | std::ofstream::app);
        octave_script <<
R"string(
figure('position',[0,0,1200,800])
hold on
i=1
title(qpmad{i}.first)
plot(qpmad{i}.second.durations, 'k')
plot(qpmad_reserve{i}.second.durations, 'r')
legend(['qpmad, mean=', num2str(mean(qpmad{i}.second.durations))], ...
       ['qpmad_reserve, mean=', num2str(mean(qpmad_reserve{i}.second.durations))], ...
        'Location', 'northeast')
xlim([0, numel(qpmad{i}.second.durations)])
print ('iterative_reserve.png', '-dpng', '-color', '-tight', '-F:8', '-S1200,800');
)string" << std::endl;
        octave_script.close();

        return (EXIT_SUCCESS);
    }
}
