/*
 * Software License Agreement (BSD License)
 *
 *  Beta-SLAM - Simultaneous localization and grid mapping with beta distributions
 *  Copyright (c) 2013-2019, Joachim Clemens, Thomas Reineking, Tobias Kluth
 *  All rights reserved.
 *
 *  This file is partially based on the GMapping ScanMatcher class
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

#include <unordered_map>

#include "bslam/maps/cells/BayesCell.h"
#include "bslam/maps/cells/BetaCell.h"

#include "bslam/utils/Log.h"
#include "bslam/utils/Config.h"
#include "bslam/utils/geometry/PointIterator.h"
#include "bslam/utils/geometry/PointCrossIterator.h"
#include "bslam/utils/uncertainty/NormalDistributionEstimator.h"


namespace bslam {


template<typename Scan, typename Map>
MapScanMatcher<Scan, Map>::MapScanMatcher( const PoseSE2 &sensorPose )
{
	m_sensor2robot.push_back( sensorPose );
	setParamsConfig();
}


template<typename Scan, typename Map>
void
MapScanMatcher<Scan, Map>::setParamsConfig() {
	m_freeCellRatio				= Config::getDouble( 	"FREE_CELL_RATIO", 				sqrt( 2.0 ) );
	m_usableRange				= Config::getDouble( 	"USABLE_RANGE", 				15.0 	);
	m_maxRange					= Config::getDouble( 	"MAX_RANGE", 					80.0 	);
	m_fullnessThreshold			= Config::getDouble( 	"FULLNESS_THRESHOLD", 			0.1 	);
	m_scoreSigmaSqr				= Config::getDouble(	"SCORE_SIGMA", 					0.25 	);
	m_scoreSigmaSqr				*= m_scoreSigmaSqr;
	m_scoreLogPlScale 			= log( NormalDistribution::pdfScale( m_scoreSigmaSqr ) );
	m_relevantDist				= Config::getDouble( 	"RELEVANT_DIST",				0.0		);
	m_optAngularStep			= Config::getDouble( 	"OPT_ANGULAR_STEP", 			3.0 	);
	m_optAngularStep			= BSLAM_DEG2RAD( m_optAngularStep );
	m_optLinearStep				= Config::getDouble( 	"OPT_LINEAR_STEP", 				.05		);
	m_optIterations				= Config::getInt( 		"OPT_ITERATIONS", 				5		);
	m_scoreRangeY				= Config::getDouble(	"SCORE_RANGE_Y",				0.05	);
	m_scoreRangeZ				= Config::getDouble(	"SCORE_RANGE_Z",				0.05	);
	m_scoreSkipPoints			= Config::getInt(		"SCORE_SKIP_POINTS",			1		);
	m_likelihoodNumMeasurements	= Config::getInt( 		"LIKELIHOOD_NUM_MEASUREMENTS",	0 		);
	m_likelihoodLikeScore		= Config::getInt( 		"LIKELIHOOD_LIKE_SCORE",		0		);

	int		kernelSize	= Config::getInt( "KERNEL_SIZE", 1 );
	bool	kernelCross = false;

	if( Dimension == 2 ) {
		kernelCross = Config::getInt( "KERNEL_CROSS_2D", 0 );
		m_scoreSkipPoints = 0;
	} else if( Dimension == 3 ) {
		kernelCross = Config::getInt( "KERNEL_CROSS_3D", 0 );	// TODO: This set to true is aprox. 4x faster, but performance has to be evaluated
	}

	// pre-compute points in kernel
	if( kernelCross ) {
		for( PointCrossIterator<Dimension> iter( MapPoint::Ones() * -kernelSize, MapPoint::Ones() * (kernelSize+1) ); iter; ++iter ) {
			m_kernelPoints.push_back( *iter );
		}
	} else {
		for( PointIterator<Dimension> iter( MapPoint::Ones() * -kernelSize, MapPoint::Ones() * (kernelSize+1) ); iter; ++iter ) {
			m_kernelPoints.push_back( *iter );
		}
	}

	m_sensorModel.setParamsConfig();
}


template<typename Scan, typename Map>
void
MapScanMatcher<Scan, Map>::setSensorPose( const PoseSE2 &sensorPose, size_t sensorId ) {
	if( sensorId >= m_sensor2robot.size() ) {
		m_sensor2robot.resize( sensorId + 1 );
	}

	m_sensor2robot[sensorId] = sensorPose;
}


template<typename Scan, typename Map>
const PoseSE2&
MapScanMatcher<Scan, Map>::getSensorPose( size_t sensorId ) const {
	return m_sensor2robot[sensorId];
}


template<typename Scan, typename Map>
double
MapScanMatcher<Scan, Map>::scoreAndLikelihood( const std::vector<ScanConstPtr> &zv, const PoseBase &pose, const Map &map, double *logLikelihood, uint32_t *count ) const {
#if USE_MAP_CACHE
	return scoreInternal( zv, pose, MapCacheRead<Map>( &map ), logLikelihood, nullptr, count );
#else
	return scoreInternal( zv, pose, map, logLikelihood, nullptr, count );
#endif
}


template<typename Scan, typename Map>
double
MapScanMatcher<Scan, Map>::score( const std::vector<ScanConstPtr> &zv, const PoseBase &pose, const Map &map, double *scoreLogLikelihood, uint32_t *count ) const {
#if USE_MAP_CACHE
	return scoreInternal( zv, pose, MapCacheRead<Map>( &map ), nullptr, scoreLogLikelihood, count );
#else
	return scoreInternal( zv, pose, map, nullptr, scoreLogLikelihood, count );
#endif
}


template<typename Scan, typename Map>
template<typename MapT>
double
MapScanMatcher<Scan, Map>::scoreInternal( const std::vector<ScanConstPtr> &zv, const PoseBase &pose, const MapT &map, double *logLikelihood, double *scoreLogLikelihood, uint32_t *count ) const {
	double 	score 			= 0,
			freeDelta		= map.delta() * m_freeCellRatio,
			noHit			= m_sensorModel.nullLogLikelihood();
	bool	excludeNoHit	= m_sensorModel.excludeNoHit();

	if( logLikelihood ) {
		*logLikelihood = 0;
	}

	if( scoreLogLikelihood ) {
		*scoreLogLikelihood = 0;
	}

	if( count ) {
		*count = 0;
	}

	for( const auto &z : zv ) {
		if( z->sensorId() >= m_sensor2robot.size() ) {
			l_wrn( "Sensor with id " << z->sensorId() << " is unkown. Skipping this scan." );
			continue;
		}

		auto sensor2world	=	pose.toPoseSE<Dimension>();
		sensor2world 		*=	m_sensor2robot[z->sensorId()];

		WorldPoint	sensorOrigin	= sensor2world * WorldPoint::Zero();	// the origin of our sensor is 0 0 (0) in sensor coordinates
		MapPoint 	sensorOriginMap = map.world2map( sensorOrigin );

		// Calculate number of points to skip
		int skipPoints = m_scoreSkipPoints;
		if( logLikelihood && !m_likelihoodLikeScore && m_likelihoodNumMeasurements ) {
			skipPoints = z->size() / m_likelihoodNumMeasurements - 1;
			if( skipPoints < 0 ) {
				skipPoints = 0;
			}
		}

		size_t i = 0;
		for( const auto &z_i : *z ) {
			// Skip m_scoreSkipPoints points or, in likelihood computation, so many that
			// the total number is m_likelihoodNumMeasurements before considering one for the calculation
			i++;
			if( (i % (skipPoints+1)) != 0 ) {
				continue;
			}

			// For score calculation in 3D maps, reduce number of considered measurements
			if( Dimension == 3 && (!logLikelihood || m_likelihoodLikeScore) ) {
				// TODO: Pre-filter measurements based on y and z coordinate and use a fixed scoreNumMeasurements (similar to likelihoodNumMeasurements) instead of scoreSkipPoints

				// For 2D pose consider a horizontal line
				if( pose.dimension() == 2 ) {
					if( fabs( z_i[2] ) > m_scoreRangeZ ) {
						continue;
					}

				// For 3D pose consider a cross
				} else if( pose.dimension() == 3 ) {
					if( fabs( z_i[1] ) > m_scoreRangeY && fabs( z_i[2] ) > m_scoreRangeZ ) {
						continue;
					}
				}
			}

			world_t dist = z_i.norm();
			if( dist == 0 || dist > m_usableRange || dist > m_maxRange ) {
				continue;
			}

			// Hit cell
			WorldPoint	pointHitSensor( z_i );
			WorldPoint	pointHit 		= sensor2world * pointHitSensor;
			MapPoint	pointHitMap		= map.world2map( pointHit );

			// Cell before hit cell
			WorldPoint	pointFree		= sensor2world * ( pointHitSensor * ((dist - freeDelta) / dist) );
			MapPoint	pointFreeMap	= map.world2map( pointFree );

			bool 		found 		= false;
			WorldPoint	bestMu		= WorldPoint::Zero();
			MapPoint	bestMuMap	= pointHitMap;

			/* TODO: If this iteration is to expensive for a 3D map, then try to get the neighboring
			 * cells more efficiently (see https://github.com/OctoMap -> OcTreeLUT.h)
			 * or evaluate if it is enough to consider only pHit and pFree or maybe a line in this direction
			 */
			for( const auto &point : m_kernelPoints ) {
				MapPoint		hitIdx 		= pointHitMap  + point,
								freeIdx		= pointFreeMap + point;

				const CellType	*hitCell	= &map.cell( hitIdx ),
								*freeCell	= &map.cell( freeIdx );

				if( hitCell->fullness() > m_fullnessThreshold && freeCell->fullness() < m_fullnessThreshold ) {
					WorldPoint mu = pointHit - hitCell->mean();
					if( !found || mu.dot( mu ) < bestMu.dot( bestMu ) ) {
						bestMu		= mu;
						bestMuMap	= hitIdx;
						found		= true;
					}
				}
			}

			if( found ) {
				double logScore = NormalDistribution::logPlFromSqr( bestMu.dot( bestMu ), m_scoreSigmaSqr );

				score += exp( logScore );

				if( count ) {
					(*count)++;
				}

				if( scoreLogLikelihood ) {
					*scoreLogLikelihood += m_scoreLogPlScale + logScore;
				}
			}

			if( logLikelihood ) {
				if( !excludeNoHit || found ) {
					*logLikelihood += m_sensorModel.forwardModel( dist, bestMu, bestMuMap, sensorOrigin, sensorOriginMap, map );
				} else {
					*logLikelihood += noHit;
				}
			}
		}
	}

	return score;
}


template<typename Scan, typename Map>
double
MapScanMatcher<Scan, Map>::registerScan( const std::vector<ScanConstPtr> &zv, const PoseBase &pose, Map *map_ ) const {
#if USE_MAP_CACHE
	MapCacheWrite<Map>	cachedMap( map_ );
	MapCacheWrite<Map>	*map( &cachedMap );
#else
	Map 				*map( map_ );
#endif

	double entropy = 0.0;

	typename SensorModelType::ScanMap	scanMap;

	for( const auto &z : zv ) {
		if( z->sensorId() >= m_sensor2robot.size() ) {
			l_wrn( "Sensor with id " << z->sensorId() << " is unkown. Skipping this scan." );
			continue;
		}

		auto sensor2world	=	pose.toPoseSE<Dimension>();
		sensor2world 		*=	m_sensor2robot[z->sensorId()];

		WorldPoint	sensorOrigin	= sensor2world * WorldPoint::Zero();	// the origin of our sensor is 0 0 (0) in sensor coordinates
		MapPoint 	sensorOriginMap = map->world2map( sensorOrigin );

		for( const auto &z_i : *z ) {
			world_t 	dist 	= z_i.norm();
			bool		usable	= (dist <= m_usableRange);

			if( !usable && m_relevantDist > 0.0 ) {
				continue;
			}

			if( dist > m_maxRange ) {
				continue;
			}

			// Hit cell
			WorldPoint	pointHitSensor( z_i );
			WorldPoint	pointHit	= sensor2world * pointHitSensor;
			WorldPoint	pointEnd	= usable ? pointHit : sensor2world * ( pointHitSensor * (m_usableRange / dist) );
			MapPoint	pointEndMap	= map->world2map( pointEnd );

			//l_dbg( "sensorOriginMap " << sensorOriginMap.transpose() << "   pHit " << pHit.transpose() << "   pEnd " << pEnd.transpose() << "   pEndMap " << pEndMap.transpose() );

			MapPoint	pointOriginMap( sensorOriginMap );
			if( m_relevantDist > 0.0 && dist > m_relevantDist ) {
				pointOriginMap = pointEndMap - ( (pointHit - sensorOrigin) * ( m_relevantDist / (dist * map->delta()) ) ).template cast<map_t>();
			}

			Line<Dimension>		scanLine( pointOriginMap, pointEndMap ); 		// all cells from laser pose to measurement
			size_t				hitIdx	= scanLine.size() + (usable ?  -1 : 0);	// if not usable, hitIdx is outside of scanLine

			// cells behind measurement (only for belief maps)
			if( usable && m_sensorModel.cellsBehind() ) {
				MapPoint pointBehindMap	= pointEndMap + ( (pointHit - sensorOrigin) * ( m_sensorModel.cellsBehind() / dist ) ).template cast<map_t>();
				scanLine.extend( pointBehindMap ); // extend the current scan line to the point behind the measurement
			}

			entropy += m_sensorModel.inverseModel( pointHit, scanLine, hitIdx, &scanMap, map );
		}

		entropy += m_sensorModel.integrate( &scanMap, map );
	}

	return entropy;
}


template<typename Scan, typename Map>
template<typename PoseType>
PoseType
MapScanMatcher<Scan, Map>::optimize( const std::vector<ScanConstPtr> &zv, const PoseType &pose, const Map &map_, double *bestScoreOut, Samples<PoseType> *samplesOut, bool samplesWithForwardModelLikelihood ) const {
#if USE_MAP_CACHE
	const MapCacheRead<Map>	map( &map_ );
#else
	const Map				&map( map_ );
#endif

	world_t		angStep 			= m_optAngularStep,
				linStep				= m_optLinearStep;
	PoseType	currentPose 		= pose;
	double		currentScore		= scoreInternal( zv, currentPose, map ),
				bestScore			= -1;
	int			refinement			= 0;

	double		logLikelihood			= 0,			// internal buffer for logLikelihood or scoreLogLikelihood;
				*logLikelihoodPtr		= nullptr,		// pass a nullptr to scoreInternal() as default in order to skip logLikelihood calculation
				*scoreLogLikelihoodPtr	= nullptr;		// pass a nullptr to scoreInternal() as default in order to skip scoreLogLikelihood calculation

	if( samplesOut ) {
		if( samplesWithForwardModelLikelihood ) {
 			logLikelihoodPtr		= &logLikelihood;	// this will cause scoreInternal() to calculate the logLikelihood as well
		} else {
			scoreLogLikelihoodPtr	= &logLikelihood;	// this will cause scoreInternal() to calculate the scoreLogLikelihood as well
		}
		samplesOut->clear();
	}


	assert( currentScore > bestScore );

	// while current score is better then last best score or the maximum number of refines is not reached
	while( currentScore > bestScore || refinement < m_optIterations ) {
		if( bestScore >= currentScore ) {
			refinement++;
			angStep *= .5;
			linStep *= .5;
		}

		bestScore					= currentScore;
		PoseType	bestLocalPose 	= currentPose,
					localPose;

		Move move = Front;
		while( move != Done ) {
			localPose = currentPose;

			move = movePose( &localPose, move, linStep, angStep );

			double odoGain = 1.0;
			/*
			// m_angularOdometryReliability and m_linearOdometryReliability are set to 0.0 in GMapping by default
			if (m_angularOdometryReliability > 0.) {
				world_t dphi 	= 	pose.phi - localPose.phi;
				dphi			= 	ANGLE_NORMALIZE( dphi );
				dphi			*= 	dphi;
				odoGain 		*= 	exp( -m_angularOdometryReliability * dphi );
			}
			if (m_linearOdometryReliability > 0.) {
				Point2w	dpos	= 	pose.pos - localPose.pos;
				world_t drho 	= 	dpos.squaredNorm();
				odoGain 		*=	exp( -m_linearOdometryReliability * drho );
			}
			*/

			double localScore = odoGain * scoreInternal( zv, localPose, map, logLikelihoodPtr, scoreLogLikelihoodPtr );

			if( samplesOut ) {
				// logLikelihood contains the logLikelihood or scoreLogLikelihood depending on samplesWithForwardModelLikelihood
				samplesOut->push_back( std::make_pair( logLikelihood, localPose ) );
			}

			if( localScore > currentScore ) {
				currentScore	= localScore;
				bestLocalPose	= localPose;
			}
		}

		currentPose = bestLocalPose;
	}

	if( bestScoreOut ) {
		*bestScoreOut = bestScore;
	}

	currentPose.normalizeRotation();
	return currentPose;
}


template<typename Scan, typename Map>
typename MapScanMatcher<Scan, Map>::Move
MapScanMatcher<Scan, Map>::movePose( PoseSE2 *p, Move move, world_t linStep, world_t angStep ) const {
	switch( move ) {
		case Front:
			p->x() += linStep;
			move = Back;
			break;
		case Back:
			p->x() -= linStep;
			move = Left;
			break;
		case Left:
			p->y() -= linStep;
			move = Right;
			break;
		case Right:
			p->y() += linStep;
			move = YawPos;
			break;
		case YawPos:
			p->phi() += angStep;
			move = YawNeg;
			break;
		case YawNeg:
			p->phi() -= angStep;
			move = Done;
			break;
		default:
			break;
	}

	return move;
}


} /* namespace bslam */

