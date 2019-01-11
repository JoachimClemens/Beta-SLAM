/*
 * Software License Agreement (BSD License)
 *
 *  Beta-SLAM - Simultaneous localization and grid mapping with beta distributions
 *  Copyright (c) 2013-2019, Joachim Clemens, Thomas Reineking, Tobias Kluth
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 *  * Neither the name of BSLAM nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef BSLAM_FACTORIAL_H_
#define BSLAM_FACTORIAL_H_

namespace bslam {

// Factorial calculations according to http://stackoverflow.com/questions/5721796/how-do-you-implement-the-factorial-function-in-c

// Recursive template implementation
template<int N>
struct StaticFactorial {
	enum { value = N * StaticFactorial<N-1>::value };
};

template<>
struct StaticFactorial<0> {
	enum { value = 1 };
};


struct Factorial {
	static inline uint64_t value( uint32_t n ) {
		switch( n ) {
			case 0:		return StaticFactorial<0>::value;
			case 1:		return StaticFactorial<1>::value;
			case 2:		return StaticFactorial<2>::value;
			case 3:		return StaticFactorial<3>::value;
			case 4:		return StaticFactorial<4>::value;
			case 5:		return StaticFactorial<5>::value;
			case 6:		return StaticFactorial<6>::value;
			case 7:		return StaticFactorial<7>::value;
			case 8:		return StaticFactorial<8>::value;
			case 9:		return StaticFactorial<9>::value;
			case 10:	return StaticFactorial<10>::value;
			case 11:	return StaticFactorial<11>::value;
			case 12:	return StaticFactorial<12>::value;
			// Higher values are no integer constants anymore
			default:	return valueRecursive( n );
		}
	}

	// Recursive
	static constexpr uint64_t valueRecursive( uint32_t n ) {
		return n == 0 ? 1 : n * valueRecursive( n - 1 );
	}

	// Iterative
	static inline uint64_t valueIterative( uint32_t n ) {
		uint64_t res = 1;

		for( uint32_t i = 1; i <= n; ++i )
			res *= i;

		return res;
	}
};




} /* namespace bslam */

#endif /* BSLAM_FACTORIAL_H_ */
