/**
    @file
    @author  Alexander Sherikov

    @copyright 2020 Alexander Sherikov, Licensed under the Apache License, Version 2.0.
    (see @ref LICENSE or http://www.apache.org/licenses/LICENSE-2.0)

    @brief
*/

#include <boost/filesystem/operations.hpp>
#include "qp.h"
#include "timer.h"


namespace benchmark
{
    class InputFileDescription
    {
    public:
        std::string path_;
        std::string name_;
        std::string extension_;

    public:
        void init(const boost::filesystem::directory_entry &file)
        {
            path_ = file.path().string();
            name_ = file.path().filename().string();
            extension_ = file.path().extension().string();
        }
    };


    template <class t_Problem>
    class Problem
    {
    public:
        benchmark::InputFileDescription file_;
        t_Problem problem_;

    public:
        void load()
        {
            std::cout << "Reading " << file_.path_ << std::endl;
            ariles2::apply<ariles2::rapidjson::Reader>(file_.path_, problem_);
            if (problem_.id_.empty())
            {
                problem_.id_ = file_.name_;
            }
        }
    };
}  // namespace qp
