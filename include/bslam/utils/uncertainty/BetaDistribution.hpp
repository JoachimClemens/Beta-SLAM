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

#include <cmath>
#include <assert.h>

#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

#include "bslam/utils/Factorial.h"

namespace bslam {

BetaDistribution::BetaDistribution( float alpha, float beta ) :
		m_alpha( alpha ),
		m_beta( beta )
{
	// Nothing else to do here
}


bool
BetaDistribution::operator==( const BetaDistribution &other ) const {
	return m_alpha == other.m_alpha && m_beta == other.m_beta;
}


double
BetaDistribution::B( double alpha, double beta ) {
	//return gamma( alpha ) * gamma( beta ) / gamma( alpha + beta );
	return boost::math::beta( alpha, beta );
}


double
BetaDistribution::Binv( double alpha, double beta ) {
	//return gamma( alpha + beta ) / gamma( alpha ) * gamma( beta );
	return 1.0 / B( alpha, beta );
}


uint32_t
BetaDistribution::choose( uint32_t n, uint32_t k ) {
	// iterative
	uint32_t res = 1;
	for( uint32_t i = 1; i <= k; i++ )
		res *= (n + 1 - i) / i;
	return res;

	// recursive
	//return k == 0 ? 1 : (n * choose( n - 1, k - 1 )) / k;
}


double
BetaDistribution::gamma( double x ) {
	//return std::tgamma( x );
	return boost::math::tgamma( x );
}

double
BetaDistribution::digamma( double x ) {
	/*
	// According to http://web.science.mq.edu.au/~mjohnson/code/digamma.c
	assert(x > 0);

	double 	result = 0,
			xx, xx2, xx4;

	for ( ; x < 7; ++x)
		result -= 1/x;

	x	-= 1.0/2.0;
	xx	= 1.0/x;
	xx2	= xx*xx;
	xx4	= xx2*xx2;

	result += log( x ) + (1. / 24.) * xx2 - (7.0 / 960.0) * xx4 + (31.0 / 8064.0) * xx4 * xx2 - (127.0 / 30720.0) * xx4 * xx4;

	return result;
	*/

	return boost::math::digamma( x );
}


double
BetaDistribution::pdf( double x ) const {
	return pdf( x, m_alpha, m_beta );
}


double
BetaDistribution::pdf( double x, double alpha, double beta ) {
	return Binv( alpha, beta ) * pow( x, alpha - 1 ) * pow( 1 - x, beta - 1 );
}


double
BetaDistribution::pmf( uint32_t k, uint32_t n, double alpha, double beta ) {
	return choose( n, k ) * B( k + alpha, n - k + beta ) * Binv( alpha, beta );
}


double
BetaDistribution::cdf( double x ) const {
	return cdf( x, m_alpha, m_beta );
}


double
BetaDistribution::cdf( double x, double alpha, double beta ) {
	return boost::math::ibeta( alpha, beta, x ); // regularized incomplete beta function
}


double
BetaDistribution::cdfComp( double x, double alpha, double beta ) {
	return boost::math::ibetac( alpha, beta, x ); // regularized incomplete beta function
}

double
BetaDistribution::cdf( uint32_t k, uint32_t n, double alpha, double beta ) {
	throw std::runtime_error( "Not implemented yet" );
	//return 1 - B( beta + n - k - 1, alpha + k + 1 ) * hyp3F2( 1, alpha + k + 1, -n + k + 1; k + 2, -beta - n + k + 2; 1 ) / (B( alpha, beta ) * B( n - k, k + 2 ) * (n + 1));
	return 0.0;
}


double
BetaDistribution::mean() const {
	return mean( m_alpha, m_beta );
}


constexpr double
BetaDistribution::mean( double alpha, double beta ) {
	return alpha / (alpha + beta);
}


constexpr double
BetaDistribution::mean( uint32_t n, double alpha, double beta ) {
	return n * alpha / (alpha + beta);
}


constexpr double
BetaDistribution::mean( uint32_t k, uint32_t n, double alpha, double beta ) {
	return (alpha + k) / (alpha + beta + n);
}


double
BetaDistribution::mode() const {
	return mode( m_alpha, m_beta );
}


constexpr double
BetaDistribution::mode( double alpha, double beta ) {
	/*
	assert( alpha >= 1 );
	assert( alpha + beta > 2 );
	*/
	return (alpha - 1) / (alpha + beta - 2);
}


double
BetaDistribution::var() const {
	return var( m_alpha, m_beta );
}


constexpr double
BetaDistribution::var( double alpha, double beta ) {
	return alpha * beta / ((alpha + beta + 1) * (alpha + beta) * (alpha + beta));
}


double
BetaDistribution::entropy() const {
	return entropy( m_alpha, m_beta );
}


double
BetaDistribution::entropy( double alpha, double beta ) {
	return log( B( alpha, beta ) ) - (alpha - 1)*digamma( alpha ) - (beta - 1)*digamma( beta ) + (alpha + beta - 2)*digamma( alpha + beta );
}


constexpr double
BetaDistribution::var( uint32_t n, double alpha, double beta ) {
	return n * alpha * beta * (alpha + beta + n) / ((alpha + beta + 1) * (alpha + beta) * (alpha + beta));
}


constexpr double
BetaDistribution::var( uint32_t k, uint32_t n, double alpha, double beta ) {
	return var( alpha + k, beta + n - k );
}

/*
double
BetaDistribution::pdf( double x, int alpha, int beta ) {
	assert( alpha > 0 );
	assert( beta > 0 );
	return Binv( alpha, beta ) * pow( x, alpha - 1 ) * pow( 1 - x, beta - 1 );
}


double
BetaDistribution::B( int alpha, int beta ) {
	return gamma( alpha ) * gamma( beta ) / gamma( alpha + beta );
}


double
BetaDistribution::Binv( int alpha, int beta ) {
	return gamma( alpha + beta ) / gamma( alpha ) * gamma( beta );
}


double
BetaDistribution::gamma( int x ) {
	return Factorial::value( x - 1 );
}
*/

double
BetaDistribution::ignorance( double priorAlpha, double priorBeta ) const {
	return ( priorAlpha + priorBeta + 1 ) / ( m_alpha + m_beta + 1 );
}


double
BetaDistribution::dissonance() const {
	double	curMean = mean();

	// Shannon entropy
	return -curMean * log( curMean ) - (1.0 - curMean) * log( 1.0 - curMean );
}


double
BetaDistribution::conflict( double epsilon ) const {
	return cdf( 0.5 + epsilon ) - cdf( 0.5 - epsilon );
}

} /* namespace bslam */
