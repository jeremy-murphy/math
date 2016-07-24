//  (C) Copyright John Maddock and Jeremy William Murphy 2016.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_POLYNOMIAL_GCD_HPP
#define BOOST_MATH_TOOLS_POLYNOMIAL_GCD_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/polynomial.hpp>
#include <boost/math/common_factor_rt.hpp>

namespace boost { namespace math {

//
// Special handling for polynomials:
// Note that gcd_traits_polynomial is templated on the coefficient type T,
// not on polynomial<T>.
//
template <typename T, typename Enabled = void>
struct gcd_traits_polynomial 
{};


// gcd_traits for Z[x].
template <typename T>
struct gcd_traits_polynomial<T, typename enable_if_c< std::numeric_limits<T>::is_integer >::type > : public gcd_traits_defaults< boost::math::tools::polynomial<T> >
{
    typedef boost::math::tools::polynomial<T> polynomial_type;

    // Just do a no-op for abs() as we'll normalize() later.
    static polynomial_type
    abs(const polynomial_type &val)
    {
        return leading_coefficient(val) < T(0) ? -val : val;
    }
    
    inline static 
    void normalize(polynomial_type &x)
    {
        using boost::lambda::_1;
        
        // This does assume that the coefficients are totally ordered.
        if (x)
        {
            if (leading_coefficient(x) < T(0))
                x.negate();
            // Find the first non-zero, because we can't do gcd(0, 0).
            T const d = gcd_range(find_if(x.data().begin(), x.data().end(), _1 != T(0)), x.data().end()).first;
            x /= d;
        }
    }
    
    inline static int make_odd(boost::math::tools::polynomial<T> &x)
    {
        unsigned r = 0;
        while (even(x))
        {
            x >>= 1;
            r++;
        }
        return r;
    }
    
    inline static bool
    less(boost::math::tools::polynomial<T> const &a, boost::math::tools::polynomial<T> const &b)
    {
        return a.size() < b.size();
    }
    
    inline static void
    subtract(boost::math::tools::polynomial<T> &a, boost::math::tools::polynomial<T> const &b)
    {
        // Antoine Joux's implementation: huge coefficients.
        T const tmp = constant_coefficient(a);
        a *= constant_coefficient(b);
        a -= tmp * b;
        normalize(a);
    }
};


// TODO: Enable for multiprecision floating point types.
template <typename T>
struct gcd_traits_polynomial<T, typename disable_if_c< std::numeric_limits<T>::is_integer >::type > : public gcd_traits_defaults< boost::math::tools::polynomial<T> >
{
    typedef boost::math::tools::polynomial<T> polynomial_type;
    typedef gcd_traits_defaults< boost::math::tools::polynomial<T> > parent_class;
    using typename parent_class::method_type;
    static const method_type method = method_type::method_euclid;
    
    static polynomial_type
    abs(const polynomial_type &val)
    {
        return leading_coefficient(val) < T(0) ? -val : val;
    }
    
    
    inline static int make_odd(polynomial_type &x)
    {
        unsigned r = 0;
        while (even(x))
        {
            x >>= 1;
            r++;
        }
        return r;
    }
    
    inline static bool
    less(polynomial_type const &a, polynomial_type const &b)
    {
        return a.size() < b.size();
    }
    
    inline static void
    subtract(polynomial_type &a, polynomial_type const &b)
    {
#if 1
        // Stepanov's implementation; suffers from floating point inaccuracy.
        using std::floor;
        
        T const r = constant_coefficient(a) / constant_coefficient(b);
        a -= r * b;
        // normalize coefficients so that leading coefficient is whole
        if (a && floor(a.data().back()) != a.data().back())
            a /= a.data().back();
#else
        using boost::lambda::_1;
        
        // Antoine Joux's implementation: huge coefficients.
        T const tmp = constant_coefficient(a);
        a *= constant_coefficient(b);
        a -= tmp * b;
        if (std::numeric_limits<T>::has_infinity)
        {
        T const &inf(std::numeric_limits<T>::infinity());
        if (std::find_if(a.data().begin(), a.data().end(), _1 == inf || _1 == -inf) != a.data().end())
            throw std::domain_error("Floating point overflow.");
    }
#endif
    }
};


template <typename T>
struct gcd_traits< boost::math::tools::polynomial<T> > : public gcd_traits_polynomial<T>
{
    
};


namespace tools{
    
/* From Knuth, 4.6.1:
* 
* We may write any nonzero polynomial u(x) from R[x] where R is a UFD as
*
*      u(x) = cont(u) · pp(u(x))
*
* where cont(u), the content of u, is an element of S, and pp(u(x)), the primitive
* part of u(x), is a primitive polynomial over S. 
* When u(x) = 0, it is convenient to define cont(u) = pp(u(x)) = O.
*/

template <class T>
T content(polynomial<T> const &x)
{
    return x ? gcd_range(x.data().begin(), x.data().end()).first : T(0);
}

// Knuth, 4.6.1
template <class T>
polynomial<T> primitive_part(polynomial<T> const &x, T const &cont)
{
    return x ? x / cont : polynomial<T>();
}


template <class T>
polynomial<T> primitive_part(polynomial<T> const &x)
{
    return primitive_part(x, content(x));
}


namespace detail
{
    /* Reduce u and v to their primitive parts and return the gcd of their 
    * contents. Used in a couple of gcd algorithms.
    */
    template <class T>
    T reduce_to_primitive(polynomial<T> &u, polynomial<T> &v)
    {
        T const u_cont = content(u), v_cont = content(v);
        u /= u_cont;
        v /= v_cont;
        return gcd(u_cont, v_cont);
    }
}


/**
* Knuth, The Art of Computer Programming: Volume 2, Third edition, 1998
* Algorithm 4.6.1C: Greatest common divisor over a unique factorization domain.
* 
* The subresultant algorithm by George E. Collins [JACM 14 (1967), 128-142], 
* later improved by W. S. Brown and J. F. Traub [JACM 18 (1971), 505-514].
* 
* Although step C3 keeps the coefficients to a "reasonable" size, they are
* still potentially several binary orders of magnitude larger than the inputs.
* Thus, this algorithm should only be used where T is a multi-precision type.
* 
* @tparam  T   Polynomial coefficient type.
* @return      Greatest common divisor of polynomials u and v.
*/
template <class T>
typename enable_if_c< std::numeric_limits<T>::is_integer, polynomial<T> >::type
subresultant_gcd(polynomial<T> u, polynomial<T> v)
{
    BOOST_ASSERT(u || v);
    
    if (!u)
        return v;
    if (!v)
        return u;
    
    typedef typename polynomial<T>::size_type N;
    
    T const d = detail::reduce_to_primitive(u, v);
    T g = 1, h = 1;
    polynomial<T> r;
    while (true)
    {
        // Pseudo-division.
        r = u % v;
        if (!r)
            return d * primitive_part(v); // Attach the content.
        if (r.degree() == 0)
            return d * polynomial<T>(T(1)); // The content is the result.
        N const delta = u.degree() - v.degree();
        BOOST_ASSERT(delta < u.degree());
        // Adjust remainder.
        u = v;
        BOOST_ASSERT(g && h);
        auto const ghd = g * detail::integer_power(h, delta);
        v = r / ghd;
        BOOST_ASSERT(v);
        g = leading_coefficient(u);
        T const tmp = detail::integer_power(g, delta);
        BOOST_ASSERT(tmp);
        if (delta <= N(1))
            h = tmp * detail::integer_power(h, N(1) - delta);
        else
            h = tmp / detail::integer_power(h, delta - N(1));
    }
}
    
}}} // namespace boost::math::tools

#endif
