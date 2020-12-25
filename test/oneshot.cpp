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

            std::cout << qp_container.problem_.id_ << "  ";
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


    try
    {
        benchmark::Timer timer;
        std::size_t result_index = 0;

        for (benchmark::Problem<qp::QP> &qp_container : qps)
        {
            std::shared_ptr<qpOASES::QProblem> qp = std::make_shared<qpOASES::QProblem>(
                    qp_container.problem_.getNumberOfVariables(),
                    qp_container.problem_.getNumberOfConstraints(),
                    qpOASES::HST_POSDEF);

            qpOASES::Options options;
            options.setToMPC();
            qp->setOptions(options);


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

            qpOASES::returnValue qpoases_return_value;

            Eigen::VectorXd solution;
            solution.resize(qp_container.problem_.getNumberOfVariables());


            std::cout << qp_container.problem_.id_ << "  ";
            int number_of_iterations = 10000;
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
                    /*max_time_ptr=*/NULL);
            qp->getPrimalSolution(solution.data());
            results.qpoases_.durations_(result_index) = timer.stop();

            results.qpoases_.errors_(result_index) = (solution - qp_container.problem_.solution_.vector_).norm();
            if (results.qpoases_.errors_(result_index) < 5e-9)
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


    if (fail)
    {
        return (EXIT_FAILURE);
    }
    else
    {
        std::cout << "qpmad " << results.qpmad_.durations_.sum() << std::endl;
        std::cout << "qpOASES " << results.qpoases_.durations_.sum() << std::endl;
        std::cout << "eiquadprog " << results.eiquadprog_.durations_.sum() << std::endl;

        // ariles2::apply<ariles2::octave::Writer>(std::cout, results);
        ariles2::apply<ariles2::octave::Writer>("oneshot.m", results);

        std::ofstream octave_script;
        octave_script.open("oneshot.m", std::ofstream::out | std::ofstream::app);
/*
        octave_script <<
R"string(
figure('position',[0,0,1200,400])
hold on
plot(qpmad.durations, 'r')
plot(qpoases.durations, 'k')
plot(eiquadprog.durations, 'b')
legend(['qpmad, mean=', num2str(mean(qpmad.durations))], ...
       ['qpoases, mean=', num2str(mean(qpoases.durations))], ...
       ['eiquadprog, mean=', num2str(mean(eiquadprog.durations))], ...
        'Location', 'northwest')
axis tight
print ('oneshot.png', '-dpng', '-color', '-tight', '-F:12', '-S1200,400');
figure('position',[0,0,1200,400])
hold on
semilogy(qpmad.durations, 'r')
semilogy(qpoases.durations, 'k')
semilogy(eiquadprog.durations, 'b')
legend(['qpmad, mean=', num2str(mean(qpmad.durations))], ...
       ['qpoases, mean=', num2str(mean(qpoases.durations))], ...
       ['eiquadprog, mean=', num2str(mean(eiquadprog.durations))], ...
        'Location', 'northwest')
axis tight
print ('oneshot_log.png', '-dpng', '-color', '-tight', '-F:12', '-S1200,400');
)string" << std::endl;
*/
        octave_script <<
R"string(
figure('position',[0,0,1200,400])
hold on
bar([qpoases.durations, eiquadprog.durations, qpmad.durations], 'basevalue', 1e-7)
legend(['qpoases, mean=', num2str(mean(qpoases.durations))], ...
       ['eiquadprog, mean=', num2str(mean(eiquadprog.durations))], ...
       ['qpmad, mean=', num2str(mean(qpmad.durations))], ...
        'Location', 'northwest')
axis tight
print ('oneshot.png', '-dpng', '-color', '-tight', '-F:12', '-S1200,400');
figure('position',[0,0,1200,400])
hold on
bar([qpoases.durations, eiquadprog.durations, qpmad.durations], 'basevalue', 1e-7)
set(gca,'yscale','log')
legend(['qpoases'], ...
       ['eiquadprog'], ...
       ['qpmad'], ...
        'Location', 'north')
axis tight
print ('oneshot_log.png', '-dpng', '-color', '-tight', '-F:12', '-S1200,400');
)string" << std::endl;

        octave_script.close();

        return (EXIT_SUCCESS);
    }
}
