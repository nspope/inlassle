#ifndef DLIB_MODIFIED_NEWTON_H
#define DLIB_MODIFIED_NEWTON_H

#include <dlib/optimization.h>

// modified Newton search strategy using generalized inverse of Hessian

// ----------------------------------------------------------------------------------------

namespace dlib
{
    template <typename inv_hessian_funct>
    class modified_newton_search_strategy_obj
    {
    public:
        explicit modified_newton_search_strategy_obj(
            const inv_hessian_funct& inv_hess
        ) : inv_hessian(inv_hess) {}

        double get_wolfe_rho (
        ) const { return 0.01; }

        double get_wolfe_sigma (
        ) const { return 0.9; }

        unsigned long get_max_line_search_iterations (
        ) const { return 100; }

        template <typename T>
        const matrix<double,0,1> get_next_direction (
            const T& x,
            const double ,
            const T& funct_derivative
        )
        {
            return -inv_hessian(x)*funct_derivative; 
        }

    private:
        inv_hessian_funct inv_hessian;
    };

    template <typename inv_hessian_funct>
    modified_newton_search_strategy_obj<inv_hessian_funct> modified_newton_search_strategy (
        inv_hessian_funct inv_hessian
    ) { return modified_newton_search_strategy_obj<inv_hessian_funct>(inv_hessian); }
} // namespace dlib

// ----------------------------------------------------------------------------------------

#endif // DLIB_MODIFIED_NEWTON_H
