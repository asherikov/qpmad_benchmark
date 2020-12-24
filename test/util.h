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


    inline void getQuadProgConstraints(
            Eigen::MatrixXd *Aeq,
            Eigen::VectorXd *Beq,
            Eigen::MatrixXd *Aineq,
            Eigen::VectorXd *Bineq,
            const std::size_t num_var,
            const Eigen::MatrixXd &ctr,
            const Eigen::VectorXd &ctr_lb,
            const Eigen::VectorXd &ctr_ub,
            const Eigen::VectorXd &lb,
            const Eigen::VectorXd &ub)
    {
        std::size_t total_ctr_num = 0;
        Eigen::MatrixXd constraints;
        Eigen::VectorXd constraints_lb;
        Eigen::VectorXd constraints_ub;
        if (lb.rows() > 0)
        {
            total_ctr_num = lb.rows() + ctr_lb.rows();

            constraints.resize(total_ctr_num, lb.rows());
            constraints_lb.resize(total_ctr_num);
            constraints_ub.resize(total_ctr_num);


            if (ctr.rows() > 0)
            {
                constraints.topRows(ctr.rows()) = ctr;
                constraints_lb.head(ctr.rows()) = ctr_lb;
                constraints_ub.head(ctr.rows()) = ctr_ub;
            }

            constraints.bottomRows(lb.rows()).setIdentity();
            constraints_lb.tail(lb.rows()) = lb;
            constraints_ub.tail(lb.rows()) = ub;
        }
        else
        {
            total_ctr_num = ctr.rows();

            constraints = ctr;
            constraints_lb = ctr_lb;
            constraints_ub = ctr_ub;
        }


        Eigen::MatrixXd Aeq_tmp(total_ctr_num, num_var);
        Eigen::VectorXd Beq_tmp(total_ctr_num);
        Eigen::MatrixXd Aineq_tmp(total_ctr_num * 2, num_var);
        Eigen::VectorXd Bineq_tmp(total_ctr_num * 2);


        std::size_t num_eq = 0;
        std::size_t num_ineq = 0;
        for (std::size_t i = 0; i < total_ctr_num; ++i)
        {
            if (std::abs(constraints_lb(i) - constraints_ub(i)) < 1e-9)  // equality
            {
                Aeq_tmp.row(num_eq) = constraints.row(i);
                Beq_tmp(num_eq) = -constraints_lb(i);
                ++num_eq;
            }
            else  // inequality
            {
                if (constraints_lb(i) > -1e10)  // finite lower bound
                {
                    Aineq_tmp.row(num_ineq) = constraints.row(i);
                    Bineq_tmp(num_ineq) = -constraints_lb(i);
                    ++num_ineq;
                }

                if (constraints_ub(i) < 1e10)  // finite upper bound
                {
                    Aineq_tmp.row(num_ineq) = -constraints.row(i);
                    Bineq_tmp(num_ineq) = constraints_ub(i);
                    ++num_ineq;
                }
            }
        }

        *Aeq = Aeq_tmp.topRows(num_eq);
        *Beq = Beq_tmp.head(num_eq);
        *Aineq = Aineq_tmp.topRows(num_ineq);
        *Bineq = Bineq_tmp.head(num_ineq);
    }
}  // namespace benchmark
