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


namespace
{
    struct QPFile
    {
        std::string path_;
        std::string name_;
        std::string extension_;
    };
}


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
        for (const QPFile & qp_file : files)
        {
            std::cout << "Reading " << qp_file.path_ << std::endl;
            qp::QP qp;
            ariles2::apply<ariles2::rapidjson::Reader>(qp_file.path_, qp);
            //ariles2::apply<ariles2::rapidjson::Writer>(std::cout, qp);
            //std::cout << qp.objective_.hessian_.rows() << std::endl;
            qp_problems.push_back(qp);
        }
    }
    catch(const std::exception & e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return (EXIT_FAILURE);
    }

    return (EXIT_SUCCESS);
}
