//  Copyright Jeremy Murphy 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning (disable : 4224)
#endif

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/std_tuple.hpp>
#include <boost/math/common_factor_rt.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/tools/polynomial_gcd.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>
#include <boost/random.hpp>
#include <boost/array.hpp>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>
#include <functional>
#include "fibonacci.hpp"
#include "../../test/table_type.hpp"
#include "table_helper.hpp"
#include "performance.hpp"


using namespace std;
using namespace boost::math::tools;
using boost::math::detail::Stein_gcd;
boost::multiprecision::cpp_int total_sum(0);

template <typename Func, class Table>
pair<double, typename Table::value_type::first_type>
exec_timed_test_foo(Func f, const Table& data, double min_elapsed = 0.5)
{
    double t = 0;
    unsigned repeats = 1;
    typename Table::value_type::first_type sum{0};
    stopwatch<boost::chrono::high_resolution_clock> w;
    do
    {
       for(unsigned count = 0; count < repeats; ++count)
       {
           for(auto i = begin(data); i != end(data); i++)
           {
               // cerr << "sum: " << sum << "\n";
                sum += f(i->first, i->second);
           }
       }

        t = boost::chrono::duration_cast<boost::chrono::duration<double>>(w.elapsed()).count();
        if(t < min_elapsed)
            repeats *= 2;
    }
    while(t < min_elapsed);
    return {t / repeats, sum};
}


template <typename T>
struct test_function_template
{
    vector<pair<T, T> > const *data;
    const char* data_name;
    
    test_function_template() {}
    
    test_function_template(vector<pair<T, T> > const &data, const char* name) : data(&data), data_name(name) {}
    
    template <typename Function>
    void operator()(pair<Function, string> const &f) const
    {
        auto result = exec_timed_test_foo(f.first, *data);
        
        auto table_name = string("gcd method comparison with ") + compiler_name() + string(" on ") + platform_name();

        report_execution_time(result.first, 
                            table_name,
                            string(data_name), 
                            string(f.second) + "\n" + boost_name());
    }
};

boost::random::mt19937 rng;
boost::random::uniform_real_distribution<> d_0_1(0, 1);
boost::random::uniform_int_distribution<> d_0_6(0, 6);
boost::random::uniform_int_distribution<> d_1_5(1, 5);
boost::random::uniform_int_distribution<> d_1_10(1, 10);
boost::random::uniform_int_distribution<> d_1_20(1, 20);

template <typename F, typename T = typename result_of<F()>::type>
polynomial<T> random_polynomial(size_t degree, double p0, F random_coefficient)
{
    polynomial<T> x;

    auto const t = [&](T &z){ z = d_0_1(rng) > p0 ? random_coefficient() : T(0); };
    do
    {
        x.data().resize(degree);
        for_each(begin(x.data()), end(x.data()), t);
        x.normalize();
    }
    while (!x);
    return x;
}


template <class T>
T get_prime_products()
{
    auto const d = d_1_5;
   int n_primes = d_0_6(rng);
   switch(n_primes)
   {
   case 0:
   case 2:
   case 4:
      // Generate a power of 2:
      return static_cast<T>(1u) << d(rng);
   case 1:
   case 3:
   case 5:
   case 6:
      // prime number:
      return boost::math::prime(d_1_20(rng) + 3);
   }
}

template <class T>
T get_uniform_random()
{
   static boost::random::uniform_int_distribution<T> minmax(numeric_limits<T>::min(), 1000);
   return minmax(rng);
}

template <class T>
inline bool even(T const& val)
{
   return !(val & 1u);
}

template <class Backend, boost::multiprecision::expression_template_option ExpressionTemplates>
inline bool even(boost::multiprecision::number<Backend, ExpressionTemplates> const& val)
{
   return !bit_test(val, 0);
}


template <class T>
void test_type(const char* name)
{
   using namespace boost::math::detail;
   typedef T int_type;
   
   using pf_test = function<polynomial<int_type>(polynomial<int_type>, polynomial<int_type>)>;
   auto test_functions = make_tuple(
      make_pair(pf_test(subresultant_gcd<int_type>), "subresultant gcd (Knuth)"s),
      make_pair(pf_test(Stein_gcd< polynomial<int_type> >), "Stein gcd (Stepanov)"s)
    );
   
   vector<pair<polynomial<int_type>, polynomial<int_type>> > data;
   data.resize(100);
   string row_name;
   test_function_template<polynomial<int_type>> tft;
   /*
   for (auto i = begin(data); i != end(data); ++i)
      *i = make_pair(random_polynomial(3, 0.5, get_prime_products<T>), 
                     random_polynomial(3, 0.5, get_prime_products<T>));
   string row_name("gcd<");
   row_name += name;
   row_name += "> (random prime number products)";
   
   tft = test_function_template<polynomial<int_type>>(data, row_name.c_str());
   boost::fusion::for_each(test_functions, tft);
    */

   for (auto i = begin(data); i != end(data); ++i)
       *i = make_pair(random_polynomial(5, 0.2, get_uniform_random<T>), 
                      random_polynomial(5, 0.2, get_uniform_random<T>));

   row_name.erase();
   row_name += "gcd<";
   row_name += name;
   row_name += "> (uniform random numbers)";
   tft = test_function_template<polynomial<int_type>>(data, row_name.c_str());
   boost::fusion::for_each(test_functions, tft);
   
   /*
   // Fibonacci number tests:
   row_name.erase();
   row_name += "gcd<";
   row_name += name;
   row_name += "> (adjacent Fibonacci numbers)";
   for_each(begin(test_functions), end(test_functions), test_function_template<int_type>(fibonacci_numbers_permution_1<T>(), row_name.c_str()));

   row_name.erase();
   row_name += "gcd<";
   row_name += name;
   row_name += "> (permutations of Fibonacci numbers)";
   for_each(begin(test_functions), end(test_functions), test_function_template<int_type>(fibonacci_numbers_permution_2<T>(), row_name.c_str()));

   row_name.erase();
   row_name += "gcd<";
   row_name += name;
   row_name += "> (Trivial cases)";
   for_each(begin(test_functions), end(test_functions), test_function_template<int_type>(trivial_gcd_test_cases<T>(), row_name.c_str()));
   */
}

/*******************************************************************************************************************/

template <class T>
T generate_random(unsigned bits_wanted)
{
   static boost::random::mt19937 gen;
   typedef boost::random::mt19937::result_type random_type;

   T max_val;
   unsigned digits;
   if(numeric_limits<T>::is_bounded && (bits_wanted == (unsigned)numeric_limits<T>::digits))
   {
      max_val = (numeric_limits<T>::max)();
      digits = numeric_limits<T>::digits;
   }
   else
   {
      max_val = T(1) << bits_wanted;
      digits = bits_wanted;
   }

   unsigned bits_per_r_val = numeric_limits<random_type>::digits - 1;
   while((random_type(1) << bits_per_r_val) > (gen.max)()) --bits_per_r_val;

   unsigned terms_needed = digits / bits_per_r_val + 1;

   T val = 0;
   for(unsigned i = 0; i < terms_needed; ++i)
   {
      val *= (gen.max)();
      val += gen();
   }
   val %= max_val;
   return val;
}


#if 0
void test_n_bits(unsigned n, string data_name, const vector<pair<boost::multiprecision::cpp_int, boost::multiprecision::cpp_int> >* p_data = 0)
{
   using namespace boost::math::detail;
   typedef boost::multiprecision::cpp_int int_type;
   vector<pair<int_type, int_type> > data, data2;

   for(unsigned i = 0; i < 1000; ++i)
   {
      data.push_back(make_pair(generate_random<int_type>(n), generate_random<int_type>(n)));
   }

   typedef pair< function<int_type(int_type, int_type)>, string> f_test;
   array<f_test, 2> test_functions{ { /*{ Stein_gcd<int_type>, "Stein_gcd" } ,{ Euclid_gcd<int_type>, "Euclid_gcd" },{ binary_textbook<int_type>, "Stein_gcd_textbook" },{ euclid_textbook<int_type>, "gcd_euclid_textbook" },{ mixed_binary_gcd<int_type>, "mixed_binary_gcd" },{ gcd_stein<int_type>, "gcd_stein" },*/{ big_gcd, "boost::multiprecision::gcd" },{ big_gcd_new, "big_gcd_new" } } };
   for_each(begin(test_functions), end(test_functions), test_function_template<int_type>(p_data ? *p_data : data, data_name.c_str(), true));
}
#endif

int main()
{
    // test_type<unsigned short>("unsigned short");
    // test_type<unsigned>("unsigned");
    // test_type<unsigned long>("unsigned long");
    // test_type<unsigned long long>("unsigned long long");
    test_type<boost::multiprecision::cpp_int>("boost::multiprecision::cpp_int");
    /*
    test_type<boost::multiprecision::uint256_t>("boost::multiprecision::uint256_t");
    test_type<boost::multiprecision::uint512_t>("boost::multiprecision::uint512_t");
    test_type<boost::multiprecision::uint1024_t>("boost::multiprecision::uint1024_t");
    */
    // test_type< boost::math::tools::polynomial<unsigned> >("polynomial<unsigned>");
    
    /*
    test_n_bits(16, "   16 bit random values");
    test_n_bits(32, "   32 bit random values");
    test_n_bits(64, "   64 bit random values");
    test_n_bits(125, "  125 bit random values");
    test_n_bits(250, "  250 bit random values");
    test_n_bits(500, "  500 bit random values");
    test_n_bits(1000, " 1000 bit random values");
    test_n_bits(5000, " 5000 bit random values");
    test_n_bits(10000, "10000 bit random values");
    //test_n_bits(100000);
    //test_n_bits(1000000);

    test_n_bits(0, "consecutive first 1000 fibonacci numbers", &fibonacci_numbers_cpp_int_permution_1());
    test_n_bits(0, "permutations of first 1000 fibonacci numbers", &fibonacci_numbers_cpp_int_permution_2());
    */
}
