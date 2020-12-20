/**
    @file
    @author  Alexander Sherikov

    @copyright 2020 Alexander Sherikov. Licensed under the Apache License,
    Version 2.0. (see LICENSE or http://www.apache.org/licenses/LICENSE-2.0)

    @brief
*/

#include <iostream>

#include <boost/filesystem/operations.hpp>

#include <qpmad/solver.h>
#include <qpOASES.hpp>

#include "qp.h"
#include "timer.h"

namespace
{
    struct QPFile
    {
        std::string path_;
        std::string name_;
        std::string extension_;
    };
}  // namespace


int main(int argc, char **argv)
{
    if (argc < 2)
    {
        return (EXIT_FAILURE);
    }

    std::vector<QPFile> files;
    for (std::size_t index = 1; index < argc; ++index)
    {
        const std::string dir = argv[index];
        if (boost::filesystem::is_directory(dir))
        {
            for (boost::filesystem::directory_entry &file : boost::filesystem::directory_iterator(dir))
            {
                if (boost::filesystem::is_regular(file))
                {
                    QPFile qp_file;

                    qp_file.path_ = file.path().string();
                    qp_file.name_ = file.path().filename().string();
                    qp_file.extension_ = file.path().extension().string();

                    files.emplace_back(qp_file);
                }
            }
        }
    }

    if (files.empty())
    {
        return (EXIT_FAILURE);
    }


    std::vector<qp::QP> qp_problems;
    try
    {
        for (const QPFile &qp_file : files)
        {
            std::cout << "Reading " << qp_file.path_ << std::endl;
            qp::QP qp_problem;
            ariles2::apply<ariles2::rapidjson::Reader>(qp_file.path_, qp_problem);
            // ariles2::apply<ariles2::rapidjson::Writer>(std::cout, qp);
            // std::cout << qp_problem.objective_.hessian_.rows() << std::endl;
            if (qp_problem.id_.empty())
            {
                qp_problem.id_ = qp_file.name_;
            }
            qp_problems.push_back(qp_problem);
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return (EXIT_FAILURE);
    }


    bool fail = false;
    try
    {
        benchmark::Timer timer;

        for (const qp::QP &qp_problem : qp_problems)
        {
            qpmad::Solver solver;
            qpmad::Solver::ReturnStatus status;
            Eigen::MatrixXd hessian_copy = qp_problem.objective_.hessian_;
            Eigen::VectorXd solution;

            std::cout << qp_problem.id_ << "  ";
            if (qp_problem.objective_.positive_definite_)
            {
                timer.start();
                status = solver.solve(
                        solution,
                        hessian_copy,
                        qp_problem.objective_.vector_,
                        qp_problem.bounds_.lower_,
                        qp_problem.bounds_.upper_,
                        qp_problem.constraints_.matrix_,
                        qp_problem.constraints_.lower_,
                        qp_problem.constraints_.upper_);
                timer.stop();

                const double err = (solution - qp_problem.solution_.vector_).norm();
                if (err < 1e-9)
                {
                    std::cout << "ok, " << timer << std::endl;
                }
                else
                {
                    fail = true;
                    std::cout << "fail, error = " << err << std::endl;
                }
            }
            else
            {
                std::cout << "skipped" << std::endl;
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return (EXIT_FAILURE);
    }

    if (fail)
    {
        return (EXIT_FAILURE);
    }
    else
    {
        return (EXIT_SUCCESS);
    }
}
