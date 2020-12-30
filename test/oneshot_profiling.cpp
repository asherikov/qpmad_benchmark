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
                if (qp_container.problem_.objective_.positive_definite_ and "HS118.json" == qp_container.problem_.id_)
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


    benchmark::Timer timer;

    try
    {
        for (benchmark::Problem<qp::QP> &qp_container : qps)
        {
            qpmad::Solver solver;
            qpmad::Solver::ReturnStatus status;
            Eigen::MatrixXd hessian_copy = qp_container.problem_.objective_.hessian_;
            Eigen::VectorXd solution;

            timer.start();
            status = solver.solve(
                    solution,
                    hessian_copy,
                    qp_container.problem_.objective_.vector_,
                    qp_container.problem_.bounds_.lower_,
                    qp_container.problem_.bounds_.upper_,
                    qp_container.problem_.constraints_.matrix_,
                    qp_container.problem_.constraints_.lower_,
                    qp_container.problem_.constraints_.upper_);
            const double duration = timer.stop();

            const double err = (solution - qp_container.problem_.solution_.vector_).norm();
            std::cout << "error = " << err << "  " << timer << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return (EXIT_FAILURE);
    }


    try
    {
        for (benchmark::Problem<qp::QP> &qp_container : qps)
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
                    qp_container.problem_.constraints_.matrix_,
                    qp_container.problem_.constraints_.lower_,
                    qp_container.problem_.constraints_.upper_,
                    qp_container.problem_.bounds_.lower_,
                    qp_container.problem_.bounds_.upper_);


            qp.reset(qp_container.problem_.getNumberOfVariables(), Aeq.rows(), Aineq.rows());


            Eigen::VectorXd solution;
            solution.resize(qp_container.problem_.getNumberOfVariables());

            timer.start();
            eiquadprog::solvers::EiquadprogFast_status status = qp.solve_quadprog(
                    qp_container.problem_.objective_.hessian_,
                    qp_container.problem_.objective_.vector_,
                    Aeq,
                    Beq,
                    Aineq,
                    Bineq,
                    solution);
            const double duration = timer.stop();

            const double err = (solution - qp_container.problem_.solution_.vector_).norm();
            std::cout << "error = " << err << "  " << timer << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return (EXIT_FAILURE);
    }


    return (EXIT_SUCCESS);
}
