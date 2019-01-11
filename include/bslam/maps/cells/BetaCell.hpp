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

namespace bslam {


template<int N>	float	BetaCell<N>::sm_priorAlpha	= 1.0;
template<int N>	float 	BetaCell<N>::sm_priorBeta	= 1.0;


template<int N>
BetaCell<N>::BetaCell( int i ) :
	BayesCell<N>( i ),
	m_distribution( sm_priorAlpha, sm_priorBeta )
{
	// Nothing else to do here
}


template<int N>
BetaCell<N>::BetaCell( const std::string &str ) {
	fromStr( str );
}


template<int N>
BetaCell<N>::~BetaCell() {
	// Nothing to do here
}


template<int N>
void
BetaCell<N>::updateNoHit( const BetaDistribution &dist, bool behind  ) {
	if( !behind )
		BayesCell<N>::updateNoHit();
	updateDistribution( dist );
}


template<int N>
void
BetaCell<N>::updateHit( const BetaDistribution &dist, const PointNw &p ) {
	BayesCell<N>::updateHit( p );
	updateDistribution( dist );
}


template<int N>
void
BetaCell<N>::integrate( const BetaCell<N> &c ) {
	BayesCell<N>::integrate( c );
	updateDistribution( c.m_distribution );
}


template<int N>
void
BetaCell<N>::updateDistribution( const BetaDistribution &dist ) {
	m_distribution.alpha()	+= dist.alpha() - sm_priorAlpha;
	m_distribution.beta()	+= dist.beta() - sm_priorBeta;
}


/*
Use the BayesCell fullness(), otherwise scan-matching will fail!

template<int N>
double
BetaCell<N>::fullness() const {
	if( !BayesCell<N>::visits() )
		return -1;

	// this is added in order that mean() is not retrieved from an unhit cell in scan-matcher
	if( !BayesCell<N>::hits() )
		return 0;

	//return BetaDistribution::mean( alpha(), beta() );

	// this should be faster to calculate
	return (sm_priorAlpha + k()) / (sm_priorAlpha + sm_priorBeta + n());
}
 */

template<int N>
double
BetaCell<N>::entropy() const {
	return m_distribution.entropy();
}


template<int N>
float
BetaCell<N>::alpha() const {
	//return sm_priorAlpha + k();	// beta bernoulli model
	return m_distribution.alpha();
}


template<int N>
float
BetaCell<N>::beta() const {
	//return sm_priorBeta + n() - k();	// beta bernoulli model
	return m_distribution.beta();
}

/*
template<int N>
float
BetaCell<N>::n() const {
	return BayesCell<N>::visits();
}


template<int N>
float
BetaCell<N>::k() const {
	return BayesCell<N>::hits();
}
*/


template<int N>
bool
BetaCell<N>::operator==( const BetaCell<N> &other ) const noexcept {
	return BayesCell<N>::operator==( other ) && m_distribution == other.m_distribution;
}


template<int N>
bool
BetaCell<N>::operator!=( const BetaCell<N> &other ) const noexcept {
	return !(*this == other);
}


template<int N>
bool
BetaCell<N>::prune() const noexcept {
	return BayesCell<N>::prune();
}


template<int N>
bool
BetaCell<N>::pruneEqual( const BetaCell<N> &other ) const noexcept {
	return m_distribution == other.m_distribution && BayesCell<N>::pruneEqual( other );
}


template<int N>
bool
BetaCell<N>::isDefault() const noexcept {
	return BayesCell<N>::isDefault() && m_distribution.alpha() == sm_priorAlpha && m_distribution.beta() == sm_priorBeta;
}


template<int N>
void
BetaCell<N>::fromStr( const std::string &str ) {
	auto splitted = split( str, ' ' );

	if( splitted.size() != N + 2 + 2 ) {
		throw std::invalid_argument( to_string( "Wrong number of values, expected " ) + to_string( N + 2 + 2 ) + ", got " + to_string( splitted.size() ) );
	}

	uint32_t hits = std::stol( splitted[N] );
	if( hits ) {
		if( !this->m_acc )
			this->m_acc = std::unique_ptr<Acc>( new Acc() );

		this->m_acc->hits = hits;
		for( int i = 0; i < N; i++ )
			this->m_acc->sum[i] = std::stod( splitted[i] );
	} else {
		this->m_acc = nullptr;
	}

	this->m_visits					= std::stol( splitted[N+1] );
	this->m_distribution.alpha() 	= std::stod( splitted[N+2] );
	this->m_distribution.beta()		= std::stod( splitted[N+3] );
}


template<int N>
std::string
BetaCell<N>::toStr() const {
	std::ostringstream	sstr;

	sstr << std::setprecision( 30 );

	if( this->m_acc )
		sstr << this->m_acc->sum.transpose() << " " << this->m_acc->hits << " ";
	else
		for( int i = 0; i < N + 1; i++ )
			sstr << "0 ";

	sstr << this->m_visits;

	sstr << " " << alpha();
	sstr << " " << beta();

	return sstr.str();
}


template<int N>
void
BetaCell<N>::fromBinary( const uint8_t *data ) {
	BayesCell<N>::fromBinary( data );
	data += BayesCell<N>::binarySize();
	for( size_t i = 0; i < 2; i++ ) {
		m_distribution[i] = CellBase::castPtr<float>( data );
		data += sizeof( float );
	}
}


template<int N>
void
BetaCell<N>::toBinary( uint8_t *data ) const {
	BayesCell<N>::toBinary( data );
	data += BayesCell<N>::binarySize();
	for( size_t i = 0; i < 2; i++ ) {
		CellBase::castPtr<float>( data ) = m_distribution[i];
		data += sizeof( float );
	}
}


template<int N>
size_t
BetaCell<N>::binarySize() {
	//									distribution
	return BayesCell<N>::binarySize() +	sizeof( float ) * 2;
}


template<int N>
size_t
BetaCell<N>::bytes() const noexcept {
	return sizeof( BetaCell<N> ) + (this->m_acc ? sizeof( Acc ) : 0);
}


} /* namespace bslam */
