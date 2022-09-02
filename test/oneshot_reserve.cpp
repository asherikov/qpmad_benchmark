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
    ARILES2_TYPED_ENTRY_(v, qpmad, Results)
#include ARILES2_INITIALIZE

    public:
        AllResults(const std::size_t size) : qpmad_(size)
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

    std::vector<benchmark::Problem<qp::QP> > qps;
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

                benchmark::Problem<qp::QP> qp_container;
                qp_container.file_.init(file);
                qp_container.load();
                if (qp_container.problem_.objective_.positive_definite_)
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


    AllResults results(qps.size());

    bool fail = false;
    try
    {
        benchmark::Timer timer;
        std::size_t result_index = 0;

        for (benchmark::Problem<qp::QP> &qp_container : qps)
        {
            qpmad::Solver solver;
            qpmad::Solver::ReturnStatus status;
            Eigen::MatrixXd hessian_copy = qp_container.problem_.objective_.hessian_;
            Eigen::VectorXd solution;

            solver.reserve(
                    hessian_copy.rows(),
                    qp_container.problem_.bounds_.lower_.rows(),
                    qp_container.problem_.constraints_.matrix_.rows());

            std::cout << qp_container.problem_.id_ << "  ";
            timer.start();
            Eigen::internal::set_is_malloc_allowed(false);
            status = solver.solve(
                    solution,
                    hessian_copy,
                    qp_container.problem_.objective_.vector_,
                    qp_container.problem_.bounds_.lower_,
                    qp_container.problem_.bounds_.upper_,
                    qp_container.problem_.constraints_.matrix_,
                    qp_container.problem_.constraints_.lower_,
                    qp_container.problem_.constraints_.upper_);
                Eigen::internal::set_is_malloc_allowed(true);
            results.qpmad_.durations_(result_index) = timer.stop();

            results.qpmad_.errors_(result_index) = (solution - qp_container.problem_.solution_.vector_).norm();
            if (results.qpmad_.errors_(result_index) < 1e-9)
            {
                std::cout << "ok, " << timer << std::endl;
            }
            else
            {
                fail = true;
                std::cout << "fail, error = " << results.qpmad_.errors_(result_index) << std::endl;
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


    if (fail)
    {
        return (EXIT_FAILURE);
    }
    else
    {
        std::cout << "qpmad " << results.qpmad_.durations_.sum() << std::endl;

        // ariles2::apply<ariles2::octave::Writer>(std::cout, results);
        ariles2::apply<ariles2::octave::Writer>("oneshot.m", results);

        return (EXIT_SUCCESS);
    }
}
