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

#ifndef BSLAM_BETADISTRIBUTION_H_
#define BSLAM_BETADISTRIBUTION_H_

#include "bslam/utils/Convenience.h"

namespace bslam {

class BetaDistribution {
public:
			inline							BetaDistribution( float alpha = 1.0, float beta = 1.0 );

			inline				bool		operator==( const BetaDistribution &other ) const;

			inline float&					operator[]( size_t i ) 			{ return m_data[i]; };
			inline const float&				operator[]( size_t i ) const	{ return m_data[i]; };


	// Beta function
	static	inline 				double		B( double alpha, double beta );
	static 	inline 				double		Binv( double alpha, double beta );

	// n choose k
	static 	inline 				uint32_t	choose( uint32_t n, uint32_t k );

	// Gamma and digamma function
	static 	inline 				double		gamma( double x );
	static 	inline 				double		digamma( double x );

	// Beta distribution probability density function
			inline				double		pdf( double x ) const;
	static 	inline 				double		pdf( double x, double alpha, double beta );

	// Beta-binomial distribution probability mass function
	static 	inline 				double		pmf( uint32_t k, uint32_t n, double alpha, double beta );

	// Beta distribution CDF and complement of CDF
			inline				double		cdf( double x ) const;
	static 	inline 				double		cdf( double x, double alpha, double beta );
	static 	inline 				double		cdfComp( double x, double alpha, double beta );

	// Beta-binomial distribution CDF
	static 	inline 				double		cdf( uint32_t k, uint32_t n, double alpha, double beta );

	// Beta distribution expected value
			inline				double		mean() const;
	static 	inline constexpr	double		mean( double alpha, double beta );

	// Beta-binomial distribution expected value
	static 	inline constexpr 	double		mean( uint32_t n, double alpha, double beta );
	static	inline constexpr 	double		mean( uint32_t k, uint32_t n, double alpha, double beta );

	// Beta distribution mode
			inline				double		mode() const;
	static 	inline constexpr 	double		mode( double alpha, double beta );

	// Beta distribution variance, standard deviation and entropy
			inline				double		var() const;
	static 	inline constexpr 	double		var( double alpha, double beta );
			inline				double		stddev() const												{ return sqrt( var() );	}
	static 	inline constexpr 	double		stddev( double alpha, double beta )							{ return sqrt( var( alpha, beta ) ); }
			inline				double		entropy() const;
	static 	inline 				double		entropy( double alpha, double beta );

	// Beta-binomial distribution variance and standard deviation
	static 	inline constexpr 	double		var( uint32_t n, double alpha, double beta );
	static 	inline constexpr 	double		var( uint32_t k, uint32_t n, double alpha, double beta );
	static 	inline constexpr	double		stddev( uint32_t n, double alpha, double beta )				{ return sqrt( var( n, alpha, beta ) ); }
	static 	inline constexpr 	double		stddev( uint32_t k, uint32_t n, double alpha, double beta )	{ return sqrt( var( k, n, alpha, beta ) ); }

	// uncertainty measures
			inline				double		ignorance( double priorAlpha, double priorBeta ) const;
			inline				double		dissonance() const;
			inline				double		differentialEntropy() const									{ return entropy();	}
			inline				double		conflict( double epsilon = 0.25 ) const;


protected:
	union {
		float	m_data[2];
		struct {
			float	m_alpha,
					m_beta;
		};
	};

public:
	GETTER( alpha );
	SETTER( alpha );
	GETTER( beta );
	SETTER( beta );
	GETTER( data );
	SETTER( data );


	/*
	// Integer variants, because the Gamma function can be calculated for integers with the faculty
	// However, the for large alpha and beta values, the large numbers result in computation errors
	static inline double	B( int alpha, int beta );
	static inline double	Binv( int alpha, int beta );
	static inline double	Gamma( int x );
	static inline double	pdf( double x, int alpha, int beta );
	*/
};

} /* namespace bslam */


#include "BetaDistribution.hpp"

#endif /* BSLAM_BETADISTRIBUTION_H_ */
