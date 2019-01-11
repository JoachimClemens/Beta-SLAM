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

#ifndef BSLAM_NORMALDISTRIBUTIONESTIMATOR_H_
#define BSLAM_NORMALDISTRIBUTIONESTIMATOR_H_

#include "bslam/utils/geometry/PoseSE2.h"


namespace bslam {


template<typename T, size_t N>
class NormalDistributionEstimatorBase {
public:
	using CovMatrix = Eigen::Matrix<double, N, N>;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	inline						NormalDistributionEstimatorBase();
	inline	virtual				~NormalDistributionEstimatorBase() {};

	inline	void				addSample( const T &s, double weight = 1 );

	inline	const T&			mean();
	inline	const CovMatrix&	cov();
	inline	double				scale();

	inline	void				reserve( size_t n );
	inline	size_t				size() const	{ return m_samples.size(); }

protected:
	virtual	inline void			calcMeanAndCov() {};

	aligned_vector< std::pair<double, T> >	m_samples;
	T			m_mean;
	CovMatrix	m_cov;
	double		m_scale;
	bool		m_changed;
};


template<typename T>
class NormalDistributionEstimator {
};


template<>
class NormalDistributionEstimator<PoseSE2> : public NormalDistributionEstimatorBase<PoseSE2, PoseSE2::DOF> {
public:
	using	BaseType	= NormalDistributionEstimatorBase<PoseSE2, PoseSE2::DOF>;
	using 	CovMatrix	= typename BaseType::CovMatrix;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
	inline	void	calcMeanAndCov();
};


template<size_t N>
class NormalDistributionEstimatorVec : public NormalDistributionEstimatorBase<Eigen::Matrix<world_t, N, 1>, N> {
public:
	using	Scalar				= world_t;
	using 	VectorType			= Eigen::Matrix<Scalar, N, 1>;
	using 	VectorDoubleType	= Eigen::Matrix<double, N, 1>;
	using	BaseType			= NormalDistributionEstimatorBase<VectorType, N>;
	using 	CovMatrix			= typename BaseType::CovMatrix;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	inline			NormalDistributionEstimatorVec();

	inline	void	addSample( const VectorType &s, double weight = 1 );

	inline	void	reserve( size_t n )	{ /* we do not have to reserve anything here */ }
	inline	size_t	size() const		{ return m_size;	}

private:
	inline	void	calcMeanAndCov();

	VectorDoubleType	m_u;
	CovMatrix			m_K;
	size_t				m_size;
};


template<>
class NormalDistributionEstimator<Point2w> : public NormalDistributionEstimatorVec<2> {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<>
class NormalDistributionEstimator<Point3w> : public NormalDistributionEstimatorVec<3> {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


} /* namespace bslam */


#include "NormalDistributionEstimator.hpp"

#endif /* BSLAM_NORMALDISTRIBUTIONESTIMATOR_H_ */
