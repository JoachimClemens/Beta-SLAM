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

#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#	include <omp.h>
#endif

#include "bslam/utils/Log.h"
#include "bslam/utils/Config.h"
#include "bslam/utils/Convenience.h"
#include "bslam/maps/CellSelectionFunction.h"
#include "bslam/viz/Color.h"

namespace bslam {


template<typename Pose, typename Map, bool ThreadSafe>
FastSLAM<Pose, Map, ThreadSafe>::FastSLAM() :
	m_numMeasurements( 0 ),
	m_linearDist( 0 ),
	m_angularDist( 0 ),
	m_likelihoodGain( 1 )
{
	l_inf( "Initializing FastSLAM" );
	setParamsConfig();
}


template<typename Pose, typename Map, bool ThreadSafe>
void
FastSLAM<Pose, Map, ThreadSafe>::setParamsConfig() {
	m_linearDistThreshold				= Config::getDouble( "LINEAR_DIST_THRESHOLD", 				1.0 	);
	m_angularDistThreshold				= Config::getDouble( "ANGULAR_DIST_THRESHOLD", 				30.0 	);
	m_angularDistThreshold				= BSLAM_DEG2RAD( m_angularDistThreshold );
	m_minScore							= Config::getDouble( "MIN_SCORE",							0.0		);
	m_resampleThreshold					= Config::getDouble( "RESAMPLE_THRESHOLD",					0.5		);
	m_likelihoodGain					= Config::getDouble( "LIKELIHOOD_GAIN",						1.0		);
	m_improvedProposalLinearRange		= Config::getDouble( "IMPROVED_PROPOSAL_LINEAR_RANGE",		0.025 	);
	m_improvedProposalAngularRange		= Config::getDouble( "IMPROVED_PROPOSAL_ANGULAR_RANGE",		1.5 	);
	m_improvedProposalAngularRange		= BSLAM_DEG2RAD( m_improvedProposalAngularRange );
	m_improvedProposalNumSamples		= Config::getInt( "IMPROVED_PROPOSAL_NUM_SAMPLES",			3		);
	m_improvedProposalStepScale			= 2.0 / ((double) m_improvedProposalNumSamples - 1.0);
	m_improvedProposalUseForwardModel	= Config::getInt( "IMPROVED_PROPOSAL_USE_FORWARD_MODEL",	0	);

	std::string improvedProposal 		= Config::getString( "IMPROVED_PROPOSAL",	"NO"	);
	std::string improvedProposalLower 	= to_lower( improvedProposal );

	if( improvedProposalLower == "no" ) {
		m_improvedProposalMode = IMPROVED_PROPOSAL_NO;
	} else if( improvedProposalLower == "after_scan_matching" ) {
		m_improvedProposalMode = IMPROVED_PROPOSAL_AFTER_SCAN_MACHTING;
	} else if( improvedProposalLower == "during_scan_matching" ) {
		m_improvedProposalMode = IMPROVED_PROPOSAL_DURING_SCAN_MATCHING;
	} else {
		throw std::invalid_argument( "Invalid value for IMPROVED_PROPOSAL: `" + improvedProposal + "'" );
	}

	typename Map::WorldPoint	min, max;
	min[0] 	= Config::getDouble( "MAP_MIN_X", -50 );
	max[0] 	= Config::getDouble( "MAP_MAX_X", 50 );
	min[1] 	= Config::getDouble( "MAP_MIN_Y", -50 );
	max[1] 	= Config::getDouble( "MAP_MAX_Y", 50 );
	if( Map::Dimension == 3 ) {
		min[2] 	= Config::getDouble( "MAP_MIN_Z", -5 );
		max[2] 	= Config::getDouble( "MAP_MAX_Z", 5 );
	}
	world_t delta = Config::getDouble( "MAP_DELTA", 0.05 );

	// init particles
	int numParticles = Config::getInt( "NUM_PARTICLES", 30 );
	m_particles.reserve( numParticles );
	m_Neff = numParticles;

	Pose	initialPose;
	Map		map( Map::WorldPoint::Zero(), min, max, delta );
	auto 	rootNode = aligned_allocate_shared< TrajectoryNode<Pose> >( initialPose );

	assert( numParticles >= 0 );

	for( size_t i = 0; i < (size_t) numParticles; i++ ) {
		m_particles.push_back( ParticleType( initialPose, map, 1.0 / numParticles ) );
		m_particles[i].node = rootNode;
	}

	// print some infos
	l_inf( "Number of particles: " << numParticles );
	l_inf( "World size:          " << map.worldSize().transpose() << " [m]" );
	l_inf( "Map size:            " << map.mapSize().transpose() );
	l_inf( "World center:        " << map.center().transpose() << " [m]" );
	l_inf( "Map center:          " << map.world2map( map.center() ).transpose() );
	l_inf( "Patch size:          " << map.patchSize().transpose() << " -> " << map.patchSize().transpose().template cast<world_t>() * map.delta() << " [m]" );
	l_inf( "Number of patches:   " << map.numPatches().transpose() );
	l_inf( "Cell size:           " << map.delta() << " [m]" );
	l_inf( "Initial pose:        " << initialPose );
	l_inf( "Improved proposal:   " << improvedProposalLower << (m_improvedProposalMode != IMPROVED_PROPOSAL_NO ? " (with forward model: " + to_string( m_improvedProposalUseForwardModel ) + ")" : "") );
#ifdef _OPENMP
	l_inf( "Running parallel with OpenMP and a maximum number of " << omp_get_max_threads() << " threads." );
#endif
}


template<typename Pose, typename Map, bool ThreadSafe>
void
FastSLAM<Pose, Map, ThreadSafe>::setSensorPose( const PoseSE2 &sensorPose, size_t sensorId ) {
	m_scanMatcher.setSensorPose( sensorPose, sensorId );
}


/**********************************
 * Processing
 **********************************/

template<typename Pose, typename Map, bool ThreadSafe>
bool
FastSLAM<Pose, Map, ThreadSafe>::scanRequired() const {
	return !m_numMeasurements || m_linearDist > m_linearDistThreshold || m_angularDist > m_angularDistThreshold;
}


#ifndef BSLAM_PRECOMPILED_HEADERS

template<typename Pose, typename Map, bool ThreadSafe>
void
FastSLAM<Pose, Map, ThreadSafe>::processOdom( double linearMove, double angularMove ) {
	m_mutex.lock();

	l_dbg( UNDERLINE "Processing odometry: " << linearMove << "m  " << BSLAM_RAD2DEG( angularMove ) << "deg" NORMAL );

	m_stopWatch.reset();
	if( m_improvedProposalMode == IMPROVED_PROPOSAL_NO ) {
#		ifdef _OPENMP
#		pragma omp parallel for
#		endif
		for( size_t i = 0; i < m_particles.size(); i++ ) {
			m_particles[i].pose = m_motionModel.draw( m_particles[i].pose, linearMove, angularMove );
		}
	} else {
#		ifdef _OPENMP
#		pragma omp parallel for
#		endif
		for( size_t i = 0; i < m_particles.size(); i++ ) {
			m_particles[i].pose = m_motionModel.predict( m_particles[i].pose, linearMove, angularMove );
		}
	}

	m_linearDist 	+= fabs( linearMove );
	m_angularDist	+= fabs( ANGLE_NORMALIZE( angularMove ) );

	l_dbg( " - Time needed:             " << m_stopWatch.timePast() );

	m_mutex.unlock();
}


template<typename Pose, typename Map, bool ThreadSafe>
void
FastSLAM<Pose, Map, ThreadSafe>::processOdom( const Pose &oldPose, const Pose &newPose ) {
	m_mutex.lock();

	l_dbg( UNDERLINE "Processing odometry: " << oldPose << " -> " << newPose << NORMAL );

	m_stopWatch.reset();

	if( m_improvedProposalMode == IMPROVED_PROPOSAL_NO ) {
#		ifdef _OPENMP
#		pragma omp parallel for
#		endif
		for( size_t i = 0; i < m_particles.size(); i++ ) {
			m_particles[i].pose = m_motionModel.draw( m_particles[i].pose, oldPose, newPose );
		}
	} else {
#		ifdef _OPENMP
#		pragma omp parallel for
#		endif
		for( size_t i = 0; i < m_particles.size(); i++ ) {
			m_particles[i].pose = m_motionModel.predict( m_particles[i].pose, oldPose, newPose );
		}
	}

	Pose diff = newPose - oldPose;
	m_linearDist 	+= diff.posNorm();
	m_angularDist	+= diff.angNorm();

	l_dbg( " - Time needed:             " << m_stopWatch.timePast() );

	m_mutex.unlock();
}


template<typename Pose, typename Map, bool ThreadSafe>
bool
FastSLAM<Pose, Map, ThreadSafe>::processScan( const std::vector<ScanConstPtr> &zv, bool forceProcessing ) {
	bool processed = false;

	if( forceProcessing || scanRequired() ) {
		m_mutex.lock();
		l_inf( UNDERLINE "Processing scan " << m_numMeasurements << ": " << zv.size() << " subscan(s)" NORMAL );
		auto startTime = m_stopWatch.now();

		// Not the first measurement?
		if( m_numMeasurements ) {
			double 				sumScore 	= 0;
			uint32_t			sumCount	= 0;

			// Scan matching
			l_inf( " - Scan matching" );

			m_stopWatch.reset();
#			ifdef _OPENMP
#			pragma omp parallel for reduction(+:sumScore,sumCount)
#			endif
			for( size_t i = 0; i < m_particles.size(); i++ ) {
				double		score,
							logLikelihood;
				uint32_t	count;
				bool		scanMatchingSuccessful 	= false,
							logLikelihoodValid		= false;

				SamplesType	samples,				// buffer for poses in optimization process
							*samplesPtr = nullptr;	// pass nullptr to m_scanMatcher.optimize() by default in order not to retrieve the poses used in optimization

				if( m_improvedProposalMode == IMPROVED_PROPOSAL_DURING_SCAN_MATCHING ) {
					samplesPtr = &samples;	// this will cause m_scanMatcher.optimize() to return the poses used optimization as well
				}

				// scan matching
				Pose scanMatchingPose 		= m_scanMatcher.optimize( zv, m_particles[i].pose, m_particles[i].map, &score, samplesPtr, m_improvedProposalUseForwardModel );
				scanMatchingSuccessful		= (score > m_minScore);

				// scan matching successful?
				if( scanMatchingSuccessful ) {
					if( m_improvedProposalMode != IMPROVED_PROPOSAL_NO ) {
						// sample from improved proposal
						m_particles[i].pose = sampleImprovedProposal( zv, m_particles[i], scanMatchingPose, samples, &logLikelihood, &logLikelihoodValid );
					} else {
						// just set corrected pose if improved proposal should not be used
						m_particles[i].pose = scanMatchingPose;
					}
				} else {
					l_wrn( "Scan matching failed for particle " << i << " with score " << score << " (must be >" << m_minScore << ")! Using odometry only." );

					if( m_improvedProposalMode != IMPROVED_PROPOSAL_NO ) {
						// scan matching failed and we where therefore not able to sample from the improved proposal, so we have to sample from the motion model here
						m_particles[i].pose = m_motionModel.draw( m_particles[i].node->pose, m_particles[i].node->pose, m_particles[i].pose );
					}
				}

				if( !logLikelihoodValid ) {
					// evaluate forward model to calculate likelihood if not given by improved proposal already
					m_scanMatcher.scoreAndLikelihood( zv, m_particles[i].pose, m_particles[i].map, &logLikelihood, &count );
				}

				m_particles[i].logWeight 	+= logLikelihood;
				m_particles[i].logWeightSum	+= logLikelihood;
				sumScore 					+= score;
				sumCount					+= count;
			}

			updateNormalizedWeightsAndNeff();

			l_inf( "   Average score:           " << sumScore / m_particles.size() );
			l_inf( "   Average used measurem.:  " << (double) sumCount / m_particles.size() );
			l_inf( "   Neff:                    " << m_Neff );
			l_inf( "   Time needed:             " << m_stopWatch.timePast() );

			resampleAndRegister( zv );

		// First measurement
		} else {
			l_inf( " - Registering first scan" );
			registerScan( zv, m_particles );
		}

		m_linearDist	= 0;
		m_angularDist	= 0;
		m_numMeasurements++;

		processed = true;
		l_inf( " - Total time needed:       " << m_stopWatch.timePast( startTime ) );

		m_mutex.unlock();
	}

	return processed;
}


template<typename Pose, typename Map, bool ThreadSafe>
Pose
FastSLAM<Pose, Map, ThreadSafe>::sampleImprovedProposal( const std::vector<ScanConstPtr> &z, const ParticleType &particle, const Pose &scanMatchingPose, const SamplesType &samples, double *logLikelihood, bool *logLikelihoodValid ) const {
	NormalDistributionEstimator<PoseSE2>	distEstimator;
	double	curLogLikelihood,
			curScore;
	PoseSE2 resPose;

	if( samples.size() ) {
		// estimate improved proposal based on scan matching samples
		assert( m_improvedProposalMode == IMPROVED_PROPOSAL_DURING_SCAN_MATCHING );

		distEstimator.reserve( samples.size() );
		for( const auto &s : samples ) {
			const PoseSE2 &curPose = s.second;

			// likelihood is provided by scan matcher (depending on m_improvedProposalUseForwardModel, this contains the scoreLogLikelihood or logLikelihood based on forward model)
			curLogLikelihood = s.first;

			curLogLikelihood += m_motionModel.logLikelihood( curPose, particle.node->pose, particle.pose ); // particle's trajectory node holds the last pose and particle's pose is the one predicted using odometry
			distEstimator.addSample( curPose, exp( curLogLikelihood ) );
		}
	} else {
		// estimate improved proposal based on samples around scan matching result
		assert( m_improvedProposalMode == IMPROVED_PROPOSAL_AFTER_SCAN_MACHTING );

		PoseSE2	curPose;
		distEstimator.reserve( m_improvedProposalNumSamples*m_improvedProposalNumSamples*m_improvedProposalNumSamples );
		for( world_t phi = -m_improvedProposalAngularRange; phi <= m_improvedProposalAngularRange; phi += m_improvedProposalAngularRange*m_improvedProposalStepScale ) {
			for( world_t x = -m_improvedProposalLinearRange; x <= m_improvedProposalLinearRange; x += m_improvedProposalLinearRange*m_improvedProposalStepScale ) {
				for( world_t y = -m_improvedProposalLinearRange; y <= m_improvedProposalLinearRange; y += m_improvedProposalLinearRange*m_improvedProposalStepScale ) {
					/*
					diffPose.x() 	= x;
					diffPose.y() 	= y;
					diffPose.phi()	= phi;

					curPose = scanMatchingPose.oplus( diffPose );
					*/
					curPose.x()		= scanMatchingPose.x() 		+ x;
					curPose.y()		= scanMatchingPose.y() 		+ y;
					curPose.phi()	= scanMatchingPose.phi()	+ phi;
					curPose.normalizeRotation();

					if( m_improvedProposalUseForwardModel ) {
						// estimate based on forward model
						curScore = m_scanMatcher.scoreAndLikelihood( z, curPose, particle.map, &curLogLikelihood );
					} else {
						// estimate based on score likelihood
						curScore = m_scanMatcher.score( z, curPose, particle.map, &curLogLikelihood );
					}

					if( curScore <= m_minScore ) {
						continue;
					}

					curLogLikelihood += m_motionModel.logLikelihood( curPose, particle.node->pose, particle.pose ); // particle's trajectory node holds the last pose and particle's pose is the one predicted using odometry

					double likelihood = exp( curLogLikelihood );
					if( likelihood == 0 || std::isinf( likelihood ) ) {
						continue;
					}

					distEstimator.addSample( curPose, likelihood );
				}
			}
		}
	}

	double scale = distEstimator.scale();
	if( distEstimator.size() > 0 && scale != 0 && !std::isinf( scale ) ) {
		// Gaussian approximation of proposal distribution
		Eigen::Vector3d	mean  	= distEstimator.mean().toVector().template cast<double>();
		Eigen::Matrix3d	cov		= distEstimator.cov();

		// multivariate sampling (might fail, if parameters are estimated based on bayes forward model)
		resPose 	= NormalDistribution::draw<PoseSE2::DOF>( mean, cov ).template cast<world_t>();

		if( resPose.pos().hasNaN() ) {
			l_wrn( "Multivariate sampling from improved proposal failed! Using uncorrelated sampling as fallback." );

			// uncorrelated sampling
			resPose.x()		= NormalDistribution::draw( mean[0], cov( 0, 0 ) );
			resPose.y()		= NormalDistribution::draw( mean[1], cov( 1, 1 ) );
			resPose.phi()	= NormalDistribution::draw( mean[2], cov( 2, 2 ) );
		}

		if( m_improvedProposalUseForwardModel ) {
			*logLikelihood 		= log( scale ); // likelihood is given by normalization term
			*logLikelihoodValid	= true;
		} else {
			*logLikelihoodValid = false;
		}

			resPose.normalizeRotation();
	} else {
		// use scan matcher pose, if sampling failed
		l_wrn( "Likelihood from normal distribution estimation is 0 or inf or no valid samples! Using scan matching pose as fallback." );
		resPose = scanMatchingPose;
		*logLikelihoodValid = false;
	}

	return resPose;
}


template<typename Pose, typename Map, bool ThreadSafe>
void
FastSLAM<Pose, Map, ThreadSafe>::registerScan( const std::vector<ScanConstPtr> &zv, ParticleVector &particles ) const {
	m_stopWatch.reset();
#	ifdef _OPENMP
#	pragma omp parallel for
#	endif
	for( size_t i = 0; i < particles.size(); i++ ) {
		m_scanMatcher.registerScan( zv, particles[i].pose, &particles[i].map );
	}
	l_inf( "   Time needed:             " << m_stopWatch.timePast() );
}


template<typename Pose, typename Map, bool ThreadSafe>
void
FastSLAM<Pose, Map, ThreadSafe>::registerScan( const std::vector<ScanConstPtr> &zv, ParticlePtrVector &particles ) const {
	m_stopWatch.reset();
#	ifdef _OPENMP
#	pragma omp parallel for
#	endif
	for( size_t i = 0; i < particles.size(); i++ ) {
		m_scanMatcher.registerScan( zv, particles[i]->pose, &particles[i]->map );
	}
	l_inf( "   Time needed:             " << m_stopWatch.timePast() );
}


template<typename Pose, typename Map, bool ThreadSafe>
void
FastSLAM<Pose, Map, ThreadSafe>::resampleAndRegister( const std::vector<ScanConstPtr> &zv ) {
	if( m_Neff < m_resampleThreshold * m_particles.size() ) {
		// Resampling required
		l_inf( " - Resampling" );
		m_stopWatch.reset();
		auto indexes = resampleIndexes();

		// Collect sampled and deleted particles
		ParticlePtrVector	uniqueParticles; // TODO: Use a list here?

		// get samples particles and reset weights
		for( size_t i = 0, size = indexes.size(); i < size; i++ ) {
			// each sampled particle only once
			if( i == 0 || indexes[i] != indexes[i-1] ) {
				ParticleType *p = &m_particles[indexes[i]];

				p->logWeight	= 0;
				p->weight		= 1.0 / m_particles.size();

				uniqueParticles.push_back( p );
			}
		}
		l_inf( "   Unique particles:        " << uniqueParticles.size() );
		l_inf( "   Time needed:             " << m_stopWatch.timePast() );

		// Registering scan only for unique particles
		l_inf( " - Registering scan" );
		registerScan( zv, uniqueParticles );

		// Pruning
		if( Map::Dimension > 2 ) {
			l_inf( " - Pruning maps" );
			m_stopWatch.reset();
			uint32_t	sumBranches = 0;
#			ifdef _OPENMP
#			pragma omp parallel for reduction(+:sumBranches)
#			endif
			for( size_t i = 0; i < uniqueParticles.size(); i++ ) {
				sumBranches += uniqueParticles[i]->map.prune();
			}
			l_inf( "   Average pruned branches: " << ((double) sumBranches / uniqueParticles.size()) );
			l_inf( "   Time needed:             " << m_stopWatch.timePast() );
		}

		// Copy particles according to there actual sample count
		ParticleVector	newParticles;
		newParticles.reserve( indexes.size() );
		size_t j = 0;
		for( size_t i = 0; i < indexes.size(); i++ ) {
			if( i != 0 && indexes[i] != indexes[i-1] )
				j++;

			newParticles.push_back( *uniqueParticles[j] );
		}

		// delete particles
		m_particles.clear();

		// Swap vectors
		m_particles.swap( newParticles );

		// Update Neff
		m_Neff = m_particles.size();
	} else {
		// No resampling required
		l_inf( " - Registering scan" );
		registerScan( zv, m_particles );
	}

	// create the new generation of tree nodes
	createTreeNodes();
}


template<typename Pose, typename Map, bool ThreadSafe>
std::vector<size_t>
FastSLAM<Pose, Map, ThreadSafe>::resampleIndexes() const {
	// Uniform resampling
    // Similar to Probabilistic Robotics, page 110, table 4.4

	static thread_local std::random_device						rd;
	static thread_local	std::mt19937							generator( 24 );
	static thread_local std::uniform_real_distribution<double>	distribution( 0.0, 1.0 );

	size_t	n		= m_particles.size(),
			j		= 0;
	double 	cweight	= 0;

	std::vector<size_t>	indexes( n );

	// cumulative weights (for interval computation, is set to zero again later)
	for( const auto &p : m_particles ) {
		cweight += p.weight;
	}

	// interval
	double interval = cweight / n;

	// initial target weight
	double target = interval * distribution( generator );

	// draw samples
	cweight	= 0;
	for( size_t i = 0; i < n; i++) {
		cweight += m_particles[i].weight;
		while( cweight > target ) {
			indexes[j++] =  i;
			target 		 += interval;
		}
	}

	assert( j == m_particles.size() );
	return indexes;
}


template<typename Pose, typename Map, bool ThreadSafe>
void
FastSLAM<Pose, Map, ThreadSafe>::updateNormalizedWeightsAndNeff() {
	double	gain	= 1. / (m_likelihoodGain * m_particles.size() ),
			lmax	= -std::numeric_limits<double>::max(),
			wcum	= 0;

	for( const auto &p : m_particles ) {
		if( p.logWeight > lmax ) {
			lmax = p.logWeight;
		}
	}

	for( auto &p : m_particles ) {
		p.weight	=	exp( gain * (p.logWeight - lmax) );
		wcum		+=	p.weight;
	}

	m_Neff = 0;
	for( auto &p : m_particles ) {
		p.weight	/= wcum;
		m_Neff		+= SQR( p.weight );
	}

	m_Neff = 1.0 / m_Neff;
}

#endif /* BSLAM_PRECOMPILED_HEADERS */


/**********************************
 * Getter
 **********************************/

#ifdef HAVE_PCL

template<typename Pose, typename Map, bool ThreadSafe>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr
FastSLAM<Pose, Map, ThreadSafe>::mapPointCloud( size_t idx, bool lock ) const {
	if( lock ) {
		m_mutex.lock();
	}

	CellSelectionFunction<typename Map::CellType> f = &selectHitCells<typename Map::CellType>;
	auto cells = m_particles[idx].map.template cells<typename Map::CellType>( f );

	if( lock ) {
		m_mutex.unlock();
	}

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr res( new pcl::PointCloud<pcl::PointXYZRGB> );
	res->reserve( cells.size() );
	pcl::PointXYZRGB p;
	for( const auto &c : cells ) {
		// get point position (we collect only hit cells, hence mean should be valid)
		const auto pos = c.second.mean();

		p.x = pos[0];
		p.y = pos[1];
		if( Map::Dimension == 3 ) {
			p.z = pos[2];
		}

		// color dependent of number if hits
		Color color = Color::hsv2rgb( (c.second.hits()*5 + 180) % 360, 1.0, 1.0 );
		p.r = color.r;
		p.g = color.g;
		p.b = color.b;

		res->push_back( p );
	}

	return res;
}

#endif


template<typename Pose, typename Map, bool ThreadSafe>
size_t
FastSLAM<Pose, Map, ThreadSafe>::bestParticleIdx( bool lock ) const {
	if( lock ) {
		m_mutex.lock();
	}

	size_t	bestIdx		= 0;
	double	bestWeight	= m_particles[0].logWeightSum;

	for( size_t i = 1; i < m_particles.size(); i++ ) {
		if( bestWeight < m_particles[i].logWeightSum ) {
			bestIdx		= i;
			bestWeight	= m_particles[i].logWeightSum;
		}
	}

	if( lock ) {
		m_mutex.unlock();
	}

	return bestIdx;
}


template<typename Pose, typename Map, bool ThreadSafe>
typename FastSLAM<Pose, Map, ThreadSafe>::Covariances
FastSLAM<Pose, Map, ThreadSafe>::covariances( const std::vector<size_t> &timePoints, bool lock ) const {
	if( !timePoints.size() ) {
		return Covariances();
	}

	if( lock ) {
		m_mutex.lock();
	}

	Covariances 		res;
	TrajectoryVector	trajVec = trajectories( false );	// TODO: calculation can be optimized by collecting the poses on the fly and take tree structure into account
	//size_t				bestIdx = bestParticleIdx( false );
	std::vector<double>	weights;

	res.resize( timePoints.size() );
	weights.resize( m_particles.size() );

	for( size_t i = 0; i < m_particles.size(); i++ ) {
		weights[i] = m_particles[i].weight;
	}

	// everything is copied now, so we can unlock the mutex
	if( lock ) {
		m_mutex.unlock();
	}

#	ifdef _OPENMP
#	pragma omp parallel for
#	endif
	for( size_t t = 0; t < timePoints.size(); t++ ) {
		NormalDistributionEstimator<PoseSE2>	estimator;
		estimator.reserve( m_particles.size() );

		for( size_t i = 0; i < m_particles.size(); i++ ) {
			estimator.addSample( trajVec[i][timePoints[t]], weights[i] );
		}

		res[t] = estimator.cov();
	}

	return res;
}


template<typename Pose, typename Map, bool ThreadSafe>
typename FastSLAM<Pose, Map, ThreadSafe>::TrajectoryVector
FastSLAM<Pose, Map, ThreadSafe>::trajectories( size_t maxDepth, bool lock ) const {
	TrajectoryVector res( m_particles.size() );

	if( lock ) {
		m_mutex.lock();
	}

#	ifdef _OPENMP
#	pragma omp parallel for
#	endif
	for( size_t i = 0; i < m_particles.size(); i++ ) {
		res[i] = trajectory( i, maxDepth, false );
	}

	if( lock ) {
		m_mutex.unlock();
	}

	return res;
}


template<typename Pose, typename Map, bool ThreadSafe>
typename FastSLAM<Pose, Map, ThreadSafe>::Trajectory
FastSLAM<Pose, Map, ThreadSafe>::trajectory( size_t idx, size_t maxDepth, bool lock ) const {
	if( lock ) {
		m_mutex.lock();
	}

	Trajectory	res;
	res.reserve( maxDepth ? maxDepth : m_numMeasurements );

	for( TrajectoryNode<Pose> *node = m_particles[idx].node.get(); node && (!maxDepth || res.size() < maxDepth); node = node->parent.get() ) {
		res.push_back( node->pose );
	}

	std::reverse( res.begin(), res.end() );

	if( lock ) {
		m_mutex.unlock();
	}

	return res;
}


template<typename Pose, typename Map, bool ThreadSafe>
void
FastSLAM<Pose, Map, ThreadSafe>::createTreeNodes() {
	// Do not run this in parallel, because multiple nodes may have the same parent
	for( auto &p : m_particles ) {
		p.newNode();
	}
}


} /* namespace bslam */

