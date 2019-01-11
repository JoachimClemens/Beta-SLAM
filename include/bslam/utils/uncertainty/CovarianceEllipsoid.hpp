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

namespace bslam {

CovarianceEllipsoid::CovarianceEllipsoid() :
	m_axes( Matrix3w::Identity() ),
	m_lengths( Point3w::Ones() )
{

}


CovarianceEllipsoid::CovarianceEllipsoid( const Matrix3w &cov ) {
	compute( cov );
}


void
CovarianceEllipsoid::compute( const Matrix3w &cov ) {
	// Based on ROS' ecl::CovarianceEllipsoid3f

	Eigen::EigenSolver<Matrix3w> esolver( cov );

	m_lengths[0]	= sqrt( (double) esolver.pseudoEigenvalueMatrix()( 0, 0 ) );
	m_lengths[1]	= sqrt( (double) esolver.pseudoEigenvalueMatrix()( 1, 1 ) );
	m_lengths[2]	= sqrt( (double) esolver.pseudoEigenvalueMatrix()( 2, 2 ) );
	m_axes			= esolver.pseudoEigenvectors();

	// Note that sorting of eigenvalues may end up with left-hand coordinate system.
	// So here we correctly sort it so that it does end up being righ-handed and normalised.
	Point3w 	c0 = m_axes.block<3,1>(0,0),
				c1 = m_axes.block<3,1>(0,1),
				c2 = m_axes.block<3,1>(0,2);
	c0.normalize();
	c1.normalize();
	c2.normalize();
	Point3w		cc = c0.cross(c1);
	if( cc.dot( c2 ) < 0 ) {
		m_axes << c1, c0, c2;
		float e 		= m_lengths[0];
		m_lengths[0]	= m_lengths[1];
		m_lengths[1]	= e;
	} else {
		m_axes << c0, c1, c2;
	}

	m_rotation.fromMatrix( m_axes );
}


} /* namespace bslam */

