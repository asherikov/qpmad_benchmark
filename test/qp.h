/**
    @file
    @author  Alexander Sherikov

    @copyright 2020 Alexander Sherikov, Licensed under the Apache License, Version 2.0.
    (see @ref LICENSE or http://www.apache.org/licenses/LICENSE-2.0)

    @brief
*/


#include <ariles2/visitors/rapidjson.h>
#include <ariles2/visitors/octave.h>

#include <ariles2/adapters/basic.h>
#include <ariles2/adapters/eigen.h>
#include <ariles2/adapters/std_vector.h>

#include <ariles2/ariles.h>
#include <ariles2/extra.h>


namespace qp
{
    class Objective : public ariles2::NonFlatMatricesRelaxedSloppyBase
    {
#define ARILES2_ENTRIES(v)                                                                                             \
    ARILES2_TYPED_ENTRY_(v, hessian, Eigen::MatrixXd)                                                                  \
    ARILES2_TYPED_ENTRY_(v, vector, Eigen::VectorXd)                                                                   \
    ARILES2_TYPED_ENTRY_(v, positive_definite, bool)
#include ARILES2_INITIALIZE

    public:
        virtual ~Objective() = default;

        void arilesVisit(const ariles2::Defaults & /*visitor*/, const ariles2::Defaults::Parameters & /*param*/)
        {
            hessian_.resize(0, 0);
            vector_.resize(0);
            positive_definite_ = false;
        }
    };


    class Bounds : public ariles2::NonFlatMatricesRelaxedSloppyBase
    {
#define ARILES2_ENTRIES(v)                                                                                             \
    ARILES2_TYPED_ENTRY_(v, lower, Eigen::VectorXd)                                                                    \
    ARILES2_TYPED_ENTRY_(v, upper, Eigen::VectorXd)
#include ARILES2_INITIALIZE

    public:
        virtual ~Bounds() = default;

        void arilesVisit(const ariles2::Defaults & /*visitor*/, const ariles2::Defaults::Parameters & /*param*/)
        {
            lower_.resize(0);
            upper_.resize(0);
        }
    };


    class Solution : public ariles2::NonFlatMatricesRelaxedSloppyBase
    {
#define ARILES2_ENTRIES(v)                                                                                             \
    ARILES2_TYPED_ENTRY_(v, vector, Eigen::VectorXd)                                                                   \
    ARILES2_TYPED_ENTRY_(v, value, double)
#include ARILES2_INITIALIZE

    public:
        virtual ~Solution() = default;

        void arilesVisit(const ariles2::Defaults & /*visitor*/, const ariles2::Defaults::Parameters & /*param*/)
        {
            vector_.resize(0);
            value_ = std::numeric_limits<double>::signaling_NaN();
        }
    };


    class Constraints : public ariles2::NonFlatMatricesRelaxedSloppyBase
    {
#define ARILES2_ENTRIES(v)                                                                                             \
    ARILES2_TYPED_ENTRY_(v, lower, Eigen::VectorXd)                                                                    \
    ARILES2_TYPED_ENTRY_(v, upper, Eigen::VectorXd)                                                                    \
    ARILES2_TYPED_ENTRY_(v, matrix, Eigen::MatrixXd)
#include ARILES2_INITIALIZE

    public:
        virtual ~Constraints() = default;

        void arilesVisit(const ariles2::Defaults & /*visitor*/, const ariles2::Defaults::Parameters & /*param*/)
        {
            lower_.resize(0);
            upper_.resize(0);
            matrix_.resize(0, 0);
        }
    };


    class QP : public ariles2::NonFlatMatricesRelaxedSloppyBase
    {
#define ARILES2_DEFAULT_ID "qp"
#define ARILES2_ENTRIES(v)                                                                                             \
    ARILES2_TYPED_ENTRY_(v, objective, Objective)                                                                      \
    ARILES2_TYPED_ENTRY_(v, bounds, Bounds)                                                                            \
    ARILES2_TYPED_ENTRY_(v, constraints, Constraints)                                                                  \
    ARILES2_TYPED_ENTRY_(v, solution, Solution)                                                                        \
    ARILES2_TYPED_ENTRY_(v, id, std::string)
#include ARILES2_INITIALIZE

    public:
        virtual ~QP() = default;

        std::size_t getNumberOfVariables() const
        {
            return (objective_.hessian_.rows());
        }

        std::size_t getNumberOfConstraints() const
        {
            return (constraints_.matrix_.rows());
        }

        std::size_t hasBounds() const
        {
            return (bounds_.lower_.rows() > 0);
        }
    };


    class IterativeQP : public ariles2::NonFlatMatricesRelaxedSloppyBase
    {
#define ARILES2_DEFAULT_ID "qp_iterative"
#define ARILES2_ENTRIES(v)                                                                                             \
    ARILES2_TYPED_ENTRY_(v, common, QP)                                                                                \
    ARILES2_TYPED_ENTRY_(v, instances, std::vector<QP>)                                                                \
    ARILES2_TYPED_ENTRY_(v, id, std::string)
#include ARILES2_INITIALIZE

    public:
        virtual ~IterativeQP() = default;

        std::size_t getNumberOfVariables() const
        {
            return (common_.getNumberOfVariables());
        }

        std::size_t getNumberOfConstraints() const
        {
            if (instances_.size() > 0)
            {
                return (common_.getNumberOfConstraints() + instances_[0].getNumberOfConstraints());
            }
            else
            {
                return (common_.getNumberOfConstraints());
            }
        }

        std::size_t hasBounds() const
        {
            if (instances_.size() > 0)
            {
                return (common_.hasBounds() > 0 or instances_[0].hasBounds());
            }
            else
            {
                return (common_.hasBounds() > 0);
            }
        }
    };
}  // namespace qp
