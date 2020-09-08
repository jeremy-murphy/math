// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP

#include <boost/math/special_functions/prime.hpp>
#include <boost/math/special_functions/interval_sieve.hpp>
#include <boost/assert.hpp>
#include <vector>
#include <iterator>
#include <cmath>
#include <thread>
#include <memory>
#include <future>
#include <numeric>
#include <algorithm>
#include <execution>


namespace boost::math {

template <typename Integer>
Integer upper_bound_prime_count(Integer x)
{
    using std::floor;
    using std::log;
    constexpr
    auto c = 30 * log(113) / 113; // Magic numbers from Wikipedia.
    return floor(c * x / log(x));
}

namespace detail
{
// https://mathworld.wolfram.com/SieveofEratosthenes.html
// https://www.cs.utexas.edu/users/misra/scannedPdf.dir/linearSieve.pdf
template<class Integer, class OutputIterator>
OutputIterator linear_sieve(Integer upper_bound, OutputIterator resultant_primes)
{
    auto const first = resultant_primes;
    auto const least_divisors_size = upper_bound + 1;
    std::unique_ptr<Integer[]> least_divisors{new Integer[least_divisors_size]{0}};

    for (Integer i{2}; i < upper_bound; ++i)
    {
        if (least_divisors[i] == 0)
        {
            least_divisors[i] = i;
            *resultant_primes++ = i;
        }

        for (Integer j = 0; j < resultant_primes - first
                            && i * resultant_primes[j] <= upper_bound
                            && resultant_primes[j] <= least_divisors[i]
                            && j < least_divisors_size;                     ++j)
        {
            least_divisors[i * resultant_primes[j]] = resultant_primes[j];
        }
    }

    return resultant_primes;
}

// This wrapper function could possibly drop the _container suffix with the
// judicious use of SFINAE.
template<class Integer, class Container>
void linear_sieve_container(Integer upper_bound, Container &resultant_primes)
{
    resultant_primes.resize(upper_bound_prime_count(upper_bound));
    auto const first = std::begin(resultant_primes);
    auto const last = linear_sieve(upper_bound, first);
    resultant_primes.resize(last - first);
}

// 4096 is where benchmarked performance of linear_sieve begins to diverge
template<class Integer>
const Integer linear_sieve_limit = Integer(4096); // Constexpr does not work with boost::multiprecision types

template<class Integer, class PrimeContainer, class Container>
void mask_sieve(Integer lower_bound, Integer upper_bound, const PrimeContainer& primes, Container &resultant_primes)
{
    Integer limit {static_cast<Integer>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};

    std::size_t primes_size {};
    auto it{primes.begin()};
    while(it != primes.end() && *it < limit)
    {
        ++primes_size;
        ++it;
    }

    const std::size_t n {static_cast<std::size_t>(upper_bound - lower_bound + 1)};
    std::unique_ptr<bool[]> is_prime {new bool[n]};
    memset(is_prime.get(), true, sizeof(*is_prime.get()) * (n));
     
    // Enable use of thread pool, not SIMD compatible
    std::for_each(std::execution::par, primes.begin(), it, [&is_prime, lower_bound, upper_bound](auto prime){
        for(Integer j {std::max(static_cast<std::size_t>(prime * prime), static_cast<std::size_t>((lower_bound + prime - 1) / prime * prime))}; 
            j < upper_bound; j += prime)
        {
            is_prime[static_cast<std::size_t>(j - lower_bound)] = false;
        }
    });

    if(lower_bound == 1)
    {
        is_prime[0] = false;
    }

    for(Integer i{lower_bound}; i <= upper_bound; ++i)
    {
        if(is_prime[static_cast<std::size_t>(i - lower_bound)])
        {
            resultant_primes.emplace_back(i);
        }
    }
}

template<class Integer, class Container>
void mask_sieve(Integer lower_bound, Integer upper_bound, Container &resultant_primes)
{
    auto limit{std::floor(std::sqrt(static_cast<double>(upper_bound))) + 1};
    std::vector<Integer> primes;

    boost::math::detail::linear_sieve_container(limit, primes);

    boost::math::detail::mask_sieve(lower_bound, upper_bound, primes, resultant_primes);
}

template<class Integer, class Container>
constexpr void prime_table(std::size_t min_index, Integer upper_bound, Container &resultant_primes)
{
    std::size_t current_index {min_index};
    Integer current_prime {2};

    while(current_prime < upper_bound)
    {
        resultant_primes.emplace_back(current_prime);
        ++current_index;
        current_prime = prime(current_index);
    }
}

template<class Integer, class Container>
constexpr void prime_table(Integer upper_bound, Container &resultant_primes)
{
    prime_table(0, upper_bound, resultant_primes);
}

template<class Integer, class PrimesContainer, class Container>
void segmented_sieve(Integer lower_bound, Integer upper_bound, const PrimesContainer &primes, Container &resultant_primes)
{
    const Integer L1_SIZE {32768};
    const Integer interval {L1_SIZE * 8};
    Integer current_lower_bound{lower_bound};
    Integer current_upper_bound{current_lower_bound + interval};

    if(current_upper_bound > upper_bound)
    {
        current_upper_bound = upper_bound;
    }

    std::size_t ranges {static_cast<std::size_t>((upper_bound - lower_bound) / interval)};

    std::vector<std::vector<Integer>> prime_vectors(ranges + 1);
    std::vector<std::future<void>> future_manager(ranges);

    auto primes_in_range {static_cast<std::size_t>(static_cast<double>(current_upper_bound) / std::log(static_cast<double>(current_upper_bound)) -
                          static_cast<double>(current_lower_bound) / std::log(static_cast<double>(current_lower_bound)))};

    for(std::size_t i {}; i < ranges; ++i)
    {
        prime_vectors[i].reserve(primes_in_range);

        future_manager.emplace_back(std::async(std::launch::async, [current_lower_bound, current_upper_bound, &primes, &prime_vectors, i]{
            boost::math::detail::IntervalSieve sieve(current_lower_bound, current_upper_bound, primes, prime_vectors[i]);
        }));

        current_lower_bound = current_upper_bound;
        current_upper_bound += interval;
    }

    prime_vectors[ranges].reserve(primes_in_range);
    future_manager.emplace_back(std::async(std::launch::async, [current_lower_bound, upper_bound, &primes, &prime_vectors]{
        boost::math::detail::IntervalSieve sieve(current_lower_bound, upper_bound, primes, prime_vectors.back());
    }));

    for(auto &&future : future_manager)
    {
        if(future.valid())
        {
            future.get();
        }
    }

    for(auto &v : prime_vectors)
    {
        resultant_primes.insert(resultant_primes.end(), v.begin(), v.end());
    }
}

template<class Integer, class Container>
void segmented_sieve(Integer lower_bound, Integer upper_bound, Container &resultant_primes)
{
    using boost::math::detail::linear_sieve_limit;
    
    Integer limit{static_cast<Integer>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};
    std::vector<Integer> primes {};
    primes.reserve(static_cast<double>(limit) / std::log(static_cast<double>(limit)));

    // Prepare for max value so you do not have to calculate this again
    if(limit < linear_sieve_limit<Integer>)
    {
        boost::math::detail::linear_sieve_container(static_cast<Integer>(limit), primes);
    }

    else
    {
        boost::math::detail::linear_sieve_container(linear_sieve_limit<Integer>, primes);
        boost::math::detail::segmented_sieve(linear_sieve_limit<Integer>, limit, primes, primes);
    }

    boost::math::detail::segmented_sieve(lower_bound, upper_bound, primes, resultant_primes);
}

template<class Integer, class Container>
void sequential_segmented_sieve(Integer lower_bound, Integer upper_bound, Container &resultant_primes)
{
    const Integer L1_SIZE {32768};
    const Integer interval {L1_SIZE * 8};
    Integer current_lower_bound{lower_bound};
    Integer current_upper_bound{current_lower_bound + interval};

    if(current_upper_bound > upper_bound)
    {
        current_upper_bound = upper_bound;
    }

    std::size_t ranges {static_cast<std::size_t>((upper_bound - lower_bound) / interval)};

    boost::math::detail::IntervalSieve sieve(current_lower_bound, current_upper_bound, resultant_primes, resultant_primes);
    if(ranges == 0)
    {
        return;
    }

    for(std::size_t i {}; i < ranges; ++i)
    {
        current_lower_bound = current_upper_bound;
        current_upper_bound += interval;
        if(current_upper_bound > upper_bound)
        {
            current_upper_bound = upper_bound;
        }
        sieve.NewRange(current_lower_bound, current_upper_bound, resultant_primes);
    }
}
} // End namespace detail

template<class Integer>
constexpr void prime_reserve(Integer upper_bound, std::vector<Integer> &prime_container)
{
    prime_container.reserve(static_cast<double>(upper_bound) / std::log(static_cast<double>(upper_bound)));
}

template<class ExecutionPolicy, class Integer, class Container>
void prime_sieve(ExecutionPolicy&& policy, Integer upper_bound, Container &primes)
{
    using boost::math::detail::linear_sieve_limit;

    if(upper_bound == 2)
    {
        return;
    }

    if(upper_bound <= linear_sieve_limit<Integer>)
    {
        boost::math::detail::linear_sieve_container(static_cast<Integer>(upper_bound), primes);
    }

    else if(typeid(policy) == typeid(std::execution::seq))
    {
        boost::math::detail::linear_sieve_container(linear_sieve_limit<Integer>, primes);
        boost::math::detail::sequential_segmented_sieve(linear_sieve_limit<Integer>, upper_bound, primes);
    }

    else
    {
        std::vector<Integer> small_primes;

        std::thread t1([&small_primes] {
            boost::math::detail::linear_sieve_container(static_cast<Integer>(linear_sieve_limit<Integer> * 2), small_primes);
        });
        std::thread t2([upper_bound, &primes] {
            boost::math::detail::segmented_sieve(static_cast<Integer>(linear_sieve_limit<Integer> * 2), upper_bound, primes);
        });

        t1.join();
        t2.join();
        primes.insert(primes.begin(), small_primes.begin(), small_primes.end());
    }
}

template<class Integer, class Container>
void prime_sieve(Integer upper_bound, Container &primes)
{
    prime_sieve(std::execution::seq, upper_bound, primes);
}


template<class ExecutionPolicy, class Integer, class Container>
void prime_range(ExecutionPolicy&& policy, Integer lower_bound, Integer upper_bound, Container &primes)
{
    using boost::math::detail::linear_sieve_limit;
    Integer limit {static_cast<Integer>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};

    if(upper_bound == 2)
    {
        return;
    }

    if(upper_bound <= linear_sieve_limit<Integer>)
    {
        boost::math::detail::linear_sieve_container(static_cast<Integer>(upper_bound), primes);
    }

    else if(typeid(policy) == typeid(std::execution::seq))
    {
        if(limit <= linear_sieve_limit<Integer>)
        {   
            boost::math::detail::linear_sieve_container(limit, primes);
            
            if(lower_bound <= limit)
            {
                boost::math::detail::sequential_segmented_sieve(limit, upper_bound, primes);
            }
            else
            {
                boost::math::detail::sequential_segmented_sieve(lower_bound, upper_bound, primes);
            }
            
        }

        else
        {
            boost::math::detail::linear_sieve_container(linear_sieve_limit<Integer>, primes);
            boost::math::detail::sequential_segmented_sieve(linear_sieve_limit<Integer>, limit, primes);
            boost::math::detail::sequential_segmented_sieve(lower_bound, upper_bound, primes);
        }
    }

    else
    {
        std::vector<Integer> small_primes {};

        if(limit <= static_cast<Integer>(linear_sieve_limit<Integer> * 2))
        {
            small_primes.reserve(1028);

            std::thread t1([limit, &small_primes] {
                boost::math::detail::linear_sieve_container(limit, small_primes);
            });
            
            std::thread t2([lower_bound, limit, upper_bound, &primes] {
                if(lower_bound <= limit)
                {
                    boost::math::detail::segmented_sieve(limit, upper_bound, primes);
                }
                else
                {
                    boost::math::detail::segmented_sieve(lower_bound, upper_bound, primes);
                }
                
            });

            t1.join();
            t2.join();

            primes.insert(primes.begin(), small_primes.begin(), small_primes.end());
        }

        else
        {
            boost::math::prime_reserve(limit, small_primes);

            std::thread t1([&small_primes] {
                boost::math::detail::linear_sieve_container(static_cast<Integer>(linear_sieve_limit<Integer> * 2), small_primes);
            });

            std::thread t2([limit, &primes] {
                boost::math::detail::segmented_sieve(static_cast<Integer>(linear_sieve_limit<Integer> * 2), limit, primes);
            });

            t1.join();
            t2.join();

            primes.insert(primes.begin(), small_primes.begin(), small_primes.end());

            boost::math::detail::segmented_sieve(lower_bound, upper_bound, primes);
        }
    }

    auto it{primes.begin()};
    while(*it < lower_bound && it != primes.end())
    {
        ++it;
    }

    primes.erase(primes.begin(), it);
}

template<class Integer, class Container>
inline void prime_range(Integer lower_bound, Integer upper_bound, Container &primes)
{
    prime_range(std::execution::seq, lower_bound, upper_bound, primes);
}
}

#endif //BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP
