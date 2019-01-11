/*
 * Software License Agreement (BSD License)
 *
 *  Beta-SLAM - Simultaneous localization and grid mapping with beta distributions
 *  Copyright (c) 2013-2019, Joachim Clemens, Thomas Reineking, Tobias Kluth
 *  All rights reserved.
 *
 *  This file is partially based on the GMapping MotionModel class
 *  Copyright (c) 2004-2007, Giorgio Grisetti, Cyrill Stachniss, Wolfram Burgard
 *  Originally licensed under the Creative Commons (Attribution-NonCommercial-ShareAlike).
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

#include "bslam/utils/Config.h"
#include "bslam/utils/uncertainty/NormalDistribution.h"


namespace bslam {

MotionModel<PoseSE2>::MotionModel() {
	setParamsConfig();
}


MotionModel<PoseSE2>::~MotionModel() {
	// Nothing to do here
}


void
MotionModel<PoseSE2>::setParamsConfig() {
	m_transTransSigma	= Config::get( "TRANS_TRANS_SIGMA", 0.1 );
	m_transRotSigma		= Config::get( "TRANS_ROT_SIGMA", 0.1 );
	m_rotRotSigma		= Config::get( "ROT_ROT_SIGMA", 0.1 );
	m_rotTransSigma		= Config::get( "ROT_TRANS_SIGMA", 0.1 );
}


PoseSE2
MotionModel<PoseSE2>::predict( const PoseSE2 &pose, double linearMove, double angularMove ) const {
	PoseSE2 res( pose );

	res.x() 	+= linearMove * cos( res.phi() );
	res.y() 	+= linearMove * sin( res.phi() );
	res.phi() 	+= angularMove;
	res.normalizeRotation();

	return res;
}


PoseSE2
MotionModel<PoseSE2>::predict( const PoseSE2 &pose, const PoseSE2 &oldPose, const PoseSE2 &newPose ) const {
	return pose.oplus( newPose.ominus( oldPose ) );
}


PoseSE2
MotionModel<PoseSE2>::draw( const PoseSE2 &pose, double linearMove, double angularMove ) const {
	PoseSE2 res( pose );

	double	lin = linearMove  	+ fabs( linearMove )  * NormalDistribution::drawStddev( m_transTransSigma )
								+ fabs( angularMove ) * NormalDistribution::drawStddev( m_transRotSigma ),
			ang = angularMove	+ fabs( linearMove )  * NormalDistribution::drawStddev( m_rotTransSigma )
								+ fabs( angularMove ) * NormalDistribution::drawStddev( m_rotRotSigma );

	res.x() 	+= lin * cos( res.phi() + .5 * ang );
	res.y() 	+= lin * sin( res.phi() + .5 * ang );
	res.phi() 	+= ang;
	res.normalizeRotation();

	return res;
}


PoseSE2
MotionModel<PoseSE2>::draw( const PoseSE2 &pose, const PoseSE2 &oldPose, const PoseSE2 &newPose ) const {
	PoseSE2 		delta 	= newPose.ominus( oldPose );
	Eigen::Vector3d	stddev 	= calcStddev( delta );

	delta.x()	+= NormalDistribution::drawStddev( (double) stddev[0] );
	delta.y()	+= NormalDistribution::drawStddev( (double) stddev[1] );
	delta.phi()	+= NormalDistribution::drawStddev( (double) stddev[2] );
	delta.normalizeRotation();

	return pose.oplus( delta );
}



double
MotionModel<PoseSE2>::logLikelihood( const PoseSE2 &pose, double linearMove, double angularMove ) const {
	throw std::runtime_error( "Not implemented yet" );
}


double
MotionModel<PoseSE2>::logLikelihood( const PoseSE2 &pose, const PoseSE2 &oldPose, const PoseSE2 &newPose ) const {
	PoseSE2 		deltaOdom 	= newPose.ominus( oldPose ),
					deltaPose	= pose.ominus( oldPose );
	Eigen::Vector3d	stddev = calcStddev( deltaOdom );

	return NormalDistribution::logPdfStddev( deltaOdom.x(), deltaPose.x(), (double) stddev[0] )
		+ NormalDistribution::logPdfStddev( deltaOdom.y(), deltaPose.y(), (double) stddev[1] )
		+ NormalDistribution::logPdfStddev( ANGLE_DIFF( deltaOdom.phi(), deltaPose.phi() ), (double) stddev[2] );
}


typename MotionModel<PoseSE2>::CovarianceMatrix
MotionModel<PoseSE2>::covariance( const PoseSE2 &oldPose, const PoseSE2 &newPose ) const {
	return covariance( newPose.ominus( oldPose ) );
}


typename MotionModel<PoseSE2>::CovarianceMatrix
MotionModel<PoseSE2>::covariance( const PoseSE2 &delta ) const {
	Eigen::Vector3d		stddev = calcStddev( delta );
	CovarianceMatrix	res( CovarianceMatrix::Zero() );

	// TODO: cross covariance?
	res.diagonal() = stddev.cwiseProduct( stddev );

	return res;
}


Eigen::Vector3d
MotionModel<PoseSE2>::calcStddev( const PoseSE2 &delta ) const {
	double xySigma = 0.3 * m_transTransSigma;

	return Eigen::Vector3d( m_transTransSigma * fabs( (world_t) delta.x() )	// x
								+ m_transRotSigma * fabs( delta.phi() )
								+ xySigma * fabs( (world_t) delta.y() ),
							m_transTransSigma * fabs( (world_t) delta.y() )	// y
								+ m_transRotSigma * fabs( delta.phi() )
								+ xySigma * fabs( (world_t) delta.x() ),
							m_rotRotSigma * fabs( delta.phi() )				// phi
								+ m_rotTransSigma * delta.pos().norm() );
}


} /* namespace bslam */
