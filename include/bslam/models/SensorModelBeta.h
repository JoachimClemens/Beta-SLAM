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

#ifndef BSLAM_SENSORMODELBETA_H_
#define BSLAM_SENSORMODELBETA_H_

#include <vector>

#include "SensorModel.h"

#include "bslam/utils/aligned_boost_unordered_map.h"
#include "bslam/utils/geometry/Line.h"
#include "bslam/utils/uncertainty/BetaDistribution.h"
#include "bslam/maps/cells/BetaCell.h"


namespace bslam {


template<typename Map, int N>
class SensorModel< Map, BetaCell<N> > {
public:
	using	ScanMap = aligned_boost_unordered_map< Pointm<N>, BetaCell<N> >;

									SensorModel( bool unused = false );

						void		setParamsConfig();

	/**
	 * @param pHit		end point of measurement in world coordinates
	 * @param scanLine	line from sensor origin to hit cell (up to usable range) in map coordinates
	 * @param hitIdx	index of the hit cell in scanLine, number >= scanLine.size() if hit cell is not part of scanLine (e.g. not in usable range)
	 * @param scanMap	cells updated in the current scan (only used in Belief model)
	 * @param map		the map
	 *
	 * @return	information gained in this step (if CALC_ENTROPY = 1)
	 */
	inline				double 		inverseModel( const typename Map::WorldPoint &pHit, const Line<N> &scanLine, size_t hitIdx, ScanMap *scanMap, Map *map ) const;

	/**
	 * @param dist				distance between measurement and sensor origin (measured range)
	 * @param bestMu			distance between measurement and best matching cell
	 * @param bestMuMap			index of the best matching cell in the map
	 * @param sensorOrigin		origin of sensor in world coordinates
	 * @param sensorOriginMap 	origin of sensor in map coordinates
	 * @param map				the map
	 *
	 * @return log-likelihood
	 */
	inline				double 		forwardModel( world_t dist, const typename Map::WorldPoint &bestMu, const typename Map::MapPoint &bestMuMap, const typename Map::WorldPoint &sensorOrigin, const typename Map::MapPoint &sensorOriginMap, const Map &map ) const;

	/**
	 * @param scanMap	cells updated in the current scan
	 * @param map		the map
	 *
	 * @return	information gained in this step (if CALC_ENTROPY = 1)
	 */
	inline				double		integrate( ScanMap *scanMap, Map *map ) const;

	inline				double		sigmaSqr() const;
	inline				double		nullLogLikelihood() const;
	inline	constexpr	bool		excludeNoHit() const;
	inline	constexpr	int			cellsBehind() const;

private:
			void						computeBetaLUT();
			void						computePdfLUT();

	inline	const BetaDistribution&		getBeta( world_t dist ) const;
	inline	double						getPdf( world_t dist ) const;


	//inline double 	betaLikelihood( world_t dist, int numCells, const BetaCell<N> **cells, const world_t *distances, int i = -1 ) const;

	double	m_betaSigmaForward,
			m_betaSigmaForwardSqr,
			m_betaSigmaForwardFactor,
			m_betaSigmaInverse,
			m_betaSigmaInverseSqr,
			m_betaSigmaInverseFactor,
			m_betaRandom,
			m_betaPriorP,
			m_mapDelta;
	bool	m_betaExcludeNoHit,
			m_scanNoConflict;
	int64_t	m_betaForwardCellRange;


	std::vector<BetaDistribution>	m_betaLUT;
	std::vector<double>				m_pdfLUT;

	double							m_betaDistScale,	// scale factor for distance discretization for beta LUT (controls its resolution)
									m_pdfDistScale;		// scale factor for distance discretization for PDF LUT (controls its resolution)

	BetaDistribution				m_betaZero;

	static constexpr size_t			MAX_SIZE = 5000;	// maximum number if entries in LUT (the actual size will be 2*MAX_SIZE, because of the two LUTs)
};


} /* namespace bslam */


#include "SensorModelBeta.hpp"


#endif /* BSLAM_SENSORMODELBETA_H_ */
