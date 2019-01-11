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

#ifndef BSLAM_MAP_SCANMATCHER_H_
#define BSLAM_MAP_SCANMATCHER_H_

#include <cstdint>
#include <functional>

#ifndef CALC_ENTROPY
#	define CALC_ENTROPY 0
#endif

#ifndef USE_MAP_CACHE
#	define USE_MAP_CACHE 0	// Cache seams to be slower because of dynamic allocations  TODO: Test this again with larger maps
#endif


#include "bslam/models/SensorModelBayes.h"
#include "bslam/models/SensorModelBeta.h"

#if USE_MAP_CACHE
#	include "maps/MapCache.h"
#endif

#include "bslam/utils/geometry/PoseBase.h"
#include "bslam/utils/aligned_vector.h"


namespace bslam {


template<typename Scan, typename Map>
class MapScanMatcher {
public:
	static constexpr int 	Dimension		= Map::Dimension;
	using 					WorldPoint		= typename Map::WorldPoint;
	using 					MapPoint		= typename Map::MapPoint;
	using			 		CellType		= typename Map::CellType;
	using 					ScanConstPtr	= typename Scan::ConstPtr;
	template<typename PoseType>
	using					Samples			= aligned_vector< std::pair<double, PoseType> >;

							MapScanMatcher( const PoseSE2 &sensorPose = PoseSE2() );

			void			setParamsConfig();
	inline	void			setSensorPose( const PoseSE2 &sensorPose, size_t sensorId = 0 );
	inline	const PoseSE2&	getSensorPose(  size_t sensorId = 0 ) const;

	// logLikelihood is calculated based on forward model, while scoreLogLikelihood is calculated based on score model
	inline 	double			scoreAndLikelihood( const std::vector<ScanConstPtr> &zv, const PoseBase &pose, const Map &map, double *logLikelihood = nullptr, uint32_t *count = nullptr ) const;
	inline 	double			score( const std::vector<ScanConstPtr> &zv, const PoseBase &pose, const Map &map, double *scoreLogLikelihood = nullptr, uint32_t *count = nullptr ) const;
			double			registerScan( const std::vector<ScanConstPtr> &zv, const PoseBase &pose, Map *map ) const;

	template<typename PoseType>
			PoseType		optimize( const std::vector<ScanConstPtr> &z, const PoseType &pose, const Map &map, double *bestScore = nullptr, Samples<PoseType> *samples = nullptr, bool samplesWithForwardModelLikelihood = false ) const;

protected:
	enum Move {
		Front, Back, Left, Right, Up, Down, RollPos, RollNeg, PitchPos, PitchNeg, YawPos, YawNeg, Done
	};

	template<typename MapT>
			double			scoreInternal( const std::vector<ScanConstPtr> &zv, const PoseBase &pose, const MapT &map, double *logLikelihood = nullptr, double *scoreLogLikelihood = nullptr, uint32_t *count = nullptr ) const;
			Move			movePose( PoseSE2 *p, Move move, world_t linStep, world_t angStep ) const;


	using SensorModelType 	= SensorModel<Map>;

	SensorModelType 			m_sensorModel;
	aligned_vector<PoseSE2>		m_sensor2robot;
	aligned_vector<MapPoint>	m_kernelPoints;

	double	m_freeCellRatio,
			m_fullnessThreshold,
			m_scoreSigmaSqr,
			m_scoreLogPlScale;
	world_t	m_optAngularStep,
			m_optLinearStep,
			m_usableRange,
			m_maxRange,
			m_relevantDist,
			m_scoreRangeY,
			m_scoreRangeZ;
	int		m_optIterations,
			m_scoreSkipPoints,
			m_likelihoodNumMeasurements;
	bool	m_likelihoodLikeScore;
};


} /* namespace bslam */


#include "MapScanMatcher.hpp"


#endif /* BSLAM_MAP_SCANMATCHER_H_ */
