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


/*************************************
 * Base Class
 *************************************/

template<typename T, size_t N>
NormalDistributionEstimatorBase<T, N>::NormalDistributionEstimatorBase() :
	m_cov( CovMatrix::Zero() ),
	m_scale( 0 ),
	m_changed( false )
{
	// nothing else to do here
}


template<typename T, size_t N>
void
NormalDistributionEstimatorBase<T, N>::addSample( const T &s, double weight ) {
	m_samples.emplace_back( std::make_pair( weight, s ) );
	m_changed 	=	true;
	m_scale		+=	weight;
}


template<typename T, size_t N>
const T&
NormalDistributionEstimatorBase<T, N>::mean() {
	this->calcMeanAndCov();
	return m_mean;
}


template<typename T, size_t N>
const typename NormalDistributionEstimatorBase<T, N>::CovMatrix&
NormalDistributionEstimatorBase<T, N>::cov() {
	this->calcMeanAndCov();
	return m_cov;
}


template<typename T, size_t N>
double
NormalDistributionEstimatorBase<T, N>::scale() {
	return m_scale;
}


template<typename T, size_t N>
void
NormalDistributionEstimatorBase<T, N>::reserve( size_t n ) {
	m_samples.reserve( n );
}


/*************************************
 * For PoseSE2
 *************************************/

void
NormalDistributionEstimator<PoseSE2>::calcMeanAndCov() {
	if( !m_changed ) {
		return;
	}

	// according to https://en.wikipedia.org/wiki/Mean_of_circular_quantities and http://stackoverflow.com/questions/1686994/weighted-average-of-angles
	// is similar to quaternions in 2D

	Eigen::Vector4d	mean( Eigen::Vector4d::Zero() );
	m_cov.setZero();

	// calculate mean
	for( const auto &sample : m_samples ) {
		mean.head<2>()	+= (sample.first / m_scale) * sample.second.translation().template cast<double>();
		mean[2]			+= (sample.first / m_scale) * sin( sample.second.phi() );
		mean[3]			+= (sample.first / m_scale) * cos( sample.second.phi() );
	}
	assert( mean[2] != 0 || mean[3] != 0 );
	m_mean.translation()	= mean.head<2>().template cast<world_t>();
	m_mean.phi()			= atan2( mean[2], mean[3] );

	// calculate cov
	for( const auto &sample : m_samples ) {
		Eigen::Vector3d diff = (sample.second.boxminus( m_mean )).template cast<double>();
		m_cov += (sample.first / m_scale) * diff * diff.transpose();
	}

	m_changed = false;
}


/*************************************
 * For Eigen::Vector
 *************************************/

template<size_t N>
NormalDistributionEstimatorVec<N>::NormalDistributionEstimatorVec() :
	BaseType(),
	m_u( VectorDoubleType::Zero() ),
	m_K( CovMatrix::Zero() ),
	m_size( 0 )
{
	// Nothing else to do here
}


template<size_t N>
void
NormalDistributionEstimatorVec<N>::addSample( const VectorType &s, double weight ) {
	VectorDoubleType	x = s.template cast<double>();
	m_u 			+= x * weight;
	m_K				+= x * x.transpose() * weight;
	this->m_scale	+= weight;
	m_size++;
	this->m_changed	= true;
}


template<size_t N>
void
NormalDistributionEstimatorVec<N>::calcMeanAndCov() {
	if( !this->m_changed )
		return;

	VectorDoubleType	mean		= m_u / this->m_scale;
	this->				m_cov		= m_K / this->m_scale - mean * mean.transpose();
	this->				m_mean		= mean.template cast<Scalar>();
	this->				m_changed 	= false;
}



} /* namespace bslam */

