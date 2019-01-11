/*
 * Software License Agreement (BSD License)
 *
 *  Beta-SLAM - Simultaneous localization and grid mapping with beta distributions
 *  Copyright (c) 2013-2019, Joachim Clemens, Thomas Reineking, Tobias Kluth
 *  All rights reserved.
 *
 *  This file is partially based on the GMapping GridSlamProcessor class
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

#ifndef BSLAM_FASTSLAM_H_
#define BSLAM_FASTSLAM_H_

#include <vector>
#include <list>

#include "SLAMTypes.h"
#include "Particle.h"
#include "TrajectoryNode.h"

#include "bslam/utils/ConditionalMutex.h"
#include "bslam/utils/Stopwatch.h"
#include "bslam/utils/Scan.h"
#include "bslam/utils/aligned_vector.h"

#include "bslam/registration/MapScanMatcher.h"
#include "bslam/models/MotionModel.h"


namespace bslam {

template<typename Pose, typename Map, bool ThreadSafe = true>
class FastSLAM {
public:
	using PoseType			= Pose;
	using MapType			= Map;
	using ParticleType		= Particle<Pose, Map>;
	using TrajectoryNodePtr = typename TrajectoryNode<Pose>::Ptr;
	using Trajectory		= aligned_vector< Pose >;
	using TrajectoryVector	= aligned_vector< Trajectory >;
	using ParticleVector	= aligned_vector< ParticleType >;
	using ParticlePtrVector	= std::vector< ParticleType* >;
	using PointType			= Pointw<Map::Dimension>;
	using ScanType			= Scan< PointType >;
	using ScanConstPtr		= typename ScanType::ConstPtr;
	using EdgeList			= std::list< std::pair<size_t, size_t> >;	// only for compatible interface with GraphSLAM
	using Covariance		= Eigen::Matrix<double, Pose::DOF, Pose::DOF>;
	using Covariances		= aligned_vector<Covariance>;
	using ScanMatcherType	= MapScanMatcher< ScanType, Map >;
	using MotionModelType	= MotionModel< Pose >;
	using SamplesType		= typename ScanMatcherType::template Samples<PoseType>;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

						FastSLAM();

	void				setParamsConfig();
	void				setSensorPose( const PoseSE2 &sensorPose, size_t sensorId = 0 );

	void				processOdom( double linearMove, double angularMove );
	void				processOdom( const Pose &oldPose, const Pose &newPose );

	bool				processScan( const std::vector<ScanConstPtr> &zv, bool forceProcessing = false ); // forceProcessing = true processes the scan, even if scanRequired() is false. This is useful when multiple sensors are used and the process handling is done by the surrounding task.
	bool				scanRequired() const; //  returns whether the next scan would be processed

	void				lock() const		{ m_mutex.lock(); 	}
	void				unlock() const		{ m_mutex.unlock();	}
	bool				try_lock() const	{ return m_mutex.try_lock(); }

	size_t				bestParticleIdx( bool lock = true ) const;
	const ParticleType&	particle( size_t idx ) const						{ return m_particles[idx]; 			}
	const Pose&			pose( size_t idx ) const							{ return m_particles[idx].pose;		}
	Covariances			covariances( const std::vector<size_t> &timePoints, bool lock = true ) const;
	const Map&			map( size_t idx, bool lock = true ) const			{ return m_particles[idx].map;		}
	EdgeList			graph( bool lock = true ) const						{ return EdgeList(); 				}	// return empty list
	double				weight( size_t idx ) const							{ return m_particles[idx].weight;	}
	TrajectoryVector	trajectories( bool lock = true ) const				{ return trajectories( 0, lock ); 	}
	TrajectoryVector	trajectories( size_t maxDepth, bool lock = true ) const;
	Trajectory			trajectory( size_t idx, bool lock = true ) const	{ return trajectory( idx, 0, lock ); }
	Trajectory			trajectory( size_t idx, size_t maxDepth, bool lock = true ) const;
	double				Neff() const 										{ return m_Neff; 					}
	size_t				numParticles() const								{ return m_particles.size(); 		}
	uint64_t			numMeasurements() const								{ return m_numMeasurements; 		}

#	ifdef HAVE_PCL
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr
						mapPointCloud( size_t idx, bool lock = true ) const;
#	endif

	static constexpr SLAMType type() {	return SLAM_FAST;	};

protected:
	void				registerScan( const std::vector<ScanConstPtr> &zv, ParticleVector &particles ) const;
	void				registerScan( const std::vector<ScanConstPtr> &zv, ParticlePtrVector &particles ) const;
	void				resampleAndRegister( const std::vector<ScanConstPtr> &zv );
	std::vector<size_t>	resampleIndexes() const;
	Pose				sampleImprovedProposal( const std::vector<ScanConstPtr> &zv, const ParticleType &particle, const Pose &scanMatchingPose, const SamplesType &samples, double *logLikelihood, bool *logLikelihoodValid ) const;


	void				createTreeNodes();
	void				updateNormalizedWeightsAndNeff();

	enum ImporvedProposalMode {
		IMPROVED_PROPOSAL_NO,
		IMPROVED_PROPOSAL_DURING_SCAN_MATCHING,
		IMPROVED_PROPOSAL_AFTER_SCAN_MACHTING
	};

	ScanMatcherType			m_scanMatcher;
	MotionModelType			m_motionModel;

    ImporvedProposalMode	m_improvedProposalMode;

	ParticleVector		m_particles;
	double				m_Neff;

    uint64_t	m_numMeasurements;

    world_t		m_linearDist,
    			m_angularDist,
    			m_linearDistThreshold,
    			m_angularDistThreshold;
    double		m_minScore,
    			m_resampleThreshold,
    			m_likelihoodGain,
    			m_improvedProposalLinearRange,
    			m_improvedProposalAngularRange,
				m_improvedProposalStepScale;
    size_t		m_improvedProposalNumSamples;
    bool		m_improvedProposalUseForwardModel;

    mutable Stopwatch	m_stopWatch;
    mutable ConditionalMutex<ThreadSafe>	m_mutex;
};

} /* namespace bslam */

#include "FastSLAM.hpp"

#endif /* BSLAM_FASTSLAM_H_ */
