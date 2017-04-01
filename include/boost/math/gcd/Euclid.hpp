//  (C) Copyright Jeremy William Murphy 2017.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_GCD_EUCLID_HPP
#define BOOST_MATH_GCD_EUCLID_HPP

#include <algorithm>

namespace boost { namespace math {
    /* The Euclidean algorithm
     * 
     * First described in Euclid's Elements (c. 300 BC), it is the earliest
     * known non-trivial algorithm. (Perhaps second only to Egyptian multiplication.)
     * 
     * The Euclidean algorithm is based on the principle that the greatest common
     * divisor of two numbers does not change if the smaller number is subtracted 
     * from the larger number.
     * The original formulation of the algorithm was repeated subtraction but this
     * was simplified computationally to taking the modulo with the advent of
     * positional notation.
     * 
     * This implementation is taken from FM2GP (www.fm2gp.com):
     *  Alexander A. Stepanov and Daniel E. Rose. 2014. From Mathematics to Generic Programming (1st ed.). Addison-Wesley Professional.
     */
    template <typename EuclideanDomain>
    inline EuclideanDomain Euclid_gcd(EuclideanDomain a, EuclideanDomain b)
    {
        using std::swap;
        while (b != EuclideanDomain(0))
        {
            a %= b;
            swap(a, b);
        }
        return a;
    }
    
    template <typename EuclideanDomain, typename ModuloAssignment>
    inline EuclideanDomain Euclid_gcd(EuclideanDomain a, EuclideanDomain b, ModuloAssignment m)
    {
        using std::swap;
        while (b != EuclideanDomain(0))
        {
            m(a, b);
            swap(a, b);
        }
        return a;
    }
    
}} // boost::math

#endif
