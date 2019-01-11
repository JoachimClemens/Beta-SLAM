/*
 * Software License Agreement (BSD License)
 *
 *  Beta-SLAM - Simultaneous localization and grid mapping with beta distributions
 *  Copyright (c) 2013-2019, Joachim Clemens, Thomas Reineking, Tobias Kluth
 *  All rights reserved.
 *
 *  This file is partially based on the ROS CovarianceEllipsoid class
 *  Copyright (c) 2009, Daniel Stonier
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

#ifndef BSLAM_COVARIANCEELLIPSE_H_
#define BSLAM_COVARIANCEELLIPSE_H_

#include "bslam/utils/geometry/Point.h"
#include "bslam/utils/geometry/Quaternion.h"
#include "bslam/utils/Convenience.h"

namespace bslam {

class CovarianceEllipsoid {
public:
	using Matrix3w	= Eigen::Matrix<world_t, 3, 3>;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	inline				CovarianceEllipsoid();
	inline				CovarianceEllipsoid( const Matrix3w &cov );

	inline	void		compute( const Matrix3w &cov );

protected:
	Matrix3w	m_axes;
	Point3w		m_lengths;
	Quaternion	m_rotation;

public:
	GETTER( axes );
	GETTER( lengths );
	GETTER( rotation );
};

} /* namespace bslam */


#include "CovarianceEllipsoid.hpp"

#endif /* BSLAM_COVARIANCEELLIPSE_H_ */
