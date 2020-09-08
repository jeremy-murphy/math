// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/prime_sieve.hpp>
#include <boost/math/special_functions/interval_sieve.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <benchmark/benchmark.h>
#include <primesieve.hpp>
#include <vector>

// Individual Algos
template<class Integer>
inline auto linear_sieve_helper(Integer upper_bound, std::vector<Integer> primes) -> std::vector<Integer>
{
    boost::math::detail::linear_sieve(upper_bound, primes);
    return primes;
}

template<class Integer>
void linear_sieve(benchmark::State& state)
{
    Integer upper = static_cast<Integer>(state.range(0));
    for(auto _ : state)
    {
        std::vector<Integer> primes;
        benchmark::DoNotOptimize(linear_sieve_helper(upper, primes));
    }
    state.SetComplexityN(state.range(0));
}

template<class Integer>
inline auto mask_sieve_helper(Integer lower_bound, Integer upper_bound, std::vector<Integer> primes) -> std::vector<Integer>
{
    boost::math::detail::mask_sieve(lower_bound, upper_bound, primes);
    return primes;
}

template<class Integer>
void mask_sieve(benchmark::State& state)
{
    Integer lower = static_cast<Integer>(2);
    Integer upper = static_cast<Integer>(state.range(0));
    for(auto _ : state)
    {
        std::vector<Integer> primes;
        benchmark::DoNotOptimize(mask_sieve_helper(lower, upper, primes));
    }
    state.SetComplexityN(state.range(0));
}

template<class Integer>
inline auto interval_sieve_helper(Integer lower_bound, Integer upper_bound, std::vector<Integer> primes) -> std::vector<Integer>
{
    std::vector<Integer> pre_sieved_primes {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    boost::math::detail::IntervalSieve sieve(lower_bound, upper_bound, pre_sieved_primes, primes);
    return primes;
}

template<class Integer>
void interval_sieve(benchmark::State& state)
{
    Integer lower = static_cast<Integer>(2);
    Integer upper = static_cast<Integer>(state.range(0));
    for(auto _ : state)
    {
        std::vector<Integer> primes;
        benchmark::DoNotOptimize(interval_sieve_helper(lower, upper, primes));
    }
    state.SetComplexityN(state.range(0));
}

template <class Integer>
void prime_sieve(benchmark::State& state)
{
    Integer upper = static_cast<Integer>(state.range(0));
    std::vector<Integer> primes;

    for(auto _ : state)
    {
        primes.clear();
        boost::math::prime_sieve(std::execution::par, upper, primes);
    }
    state.SetComplexityN(state.range(0));
}

template <class Integer>
inline auto kimwalish_primes_helper(Integer upper, std::vector<Integer> primes) -> std::vector<Integer>
{
    primesieve::generate_primes(upper, &primes);
    return primes;
}

template <class Integer>
void kimwalish_primes(benchmark::State& state)
{
    Integer upper = static_cast<Integer>(state.range(0));
    for (auto _ : state)
    {
        std::vector<Integer> primes;
        benchmark::DoNotOptimize(kimwalish_primes_helper(upper, primes));
    }
    state.SetComplexityN(state.range(0));
}

// Invidiual Implementations
// Linear
//BENCHMARK_TEMPLATE(linear_sieve, int32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 16)->Complexity(benchmark::oN);
//BENCHMARK_TEMPLATE(linear_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 16)->Complexity(benchmark::oN);
//BENCHMARK_TEMPLATE(linear_sieve, uint32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 16)->Complexity(benchmark::oN);

// Segmented
//BENCHMARK_TEMPLATE(mask_sieve, int32_t)->RangeMultiplier(2)->Range(1 << 2, 2 << 22)->Complexity(benchmark::oNLogN);
//BENCHMARK_TEMPLATE(mask_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 14, 2 << 26)->Complexity(benchmark::oNLogN);
//BENCHMARK_TEMPLATE(interval_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 14, 2 << 26)->Complexity();
//BENCHMARK_TEMPLATE(mask_sieve, uint32_t)->RangeMultiplier(2)->Range(1 << 2, 2 << 22)->Complexity(benchmark::oNLogN);

// Complete Implemenations
//BENCHMARK_TEMPLATE(prime_sieve, int32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
//BENCHMARK_TEMPLATE(prime_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
//BENCHMARK_TEMPLATE(kimwalish_primes, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime(); // Benchmark
BENCHMARK_TEMPLATE(prime_sieve, uint32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
BENCHMARK_TEMPLATE(prime_sieve, boost::multiprecision::cpp_int)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
BENCHMARK_TEMPLATE(prime_sieve, boost::multiprecision::mpz_int)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();

BENCHMARK_MAIN();
