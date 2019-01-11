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

#include "bslam/utils/Log.h"
#include "bslam/utils/Config.h"
#include "bslam/utils/uncertainty/NormalDistribution.h"


namespace bslam {


/************************************
 * Bayes sensor model
 ************************************/

template<typename Map, int N>
SensorModel< Map, BayesCell<N> >::SensorModel( bool unused )
{
	setParamsConfig();
}


template<typename Map, int N>
void
SensorModel< Map, BayesCell<N> >::setParamsConfig() {
	m_bayesSigmaSqr 	=	Config::getDouble( "BAYES_SIGMA",  	0.25 	);
	m_bayesSigmaSqr 	*=	m_bayesSigmaSqr;
	m_scanNoConflict	=	Config::getInt( "SCAN_NO_CONFLICT",	0		);
}


template<typename Map, int N>
double
SensorModel< Map, BayesCell<N> >::inverseModel( const typename Map::WorldPoint &pHit, const Line<N> &scanLine, size_t hitIdx, ScanMap *scanMap, Map *map ) const {
	static_assert( std::is_same< typename Map::CellType, BayesCell<N> >(), "Map cell type must be BayesCell<N>."  );

	double entropy = 0.0;

	for( size_t i = 0; i < scanLine.size(); i++ ) {
		BayesCell<N> *cell;

		if( m_scanNoConflict ) {
			cell = &(*scanMap)[scanLine[i]];
		} else {
			cell = &map->cell( scanLine[i] );
		}

#		if CALC_ENTROPY
		if( !m_scanNoConflict ) {
			entropy -= cell->entropy();
		}
#		endif

		if( i == hitIdx ) {
			if( m_scanNoConflict && cell->fullness() != 1 ) {
				(*cell) = BayesCell<N>();
			}
			cell->updateHit( pHit );
		} else {
			if( !m_scanNoConflict || cell->fullness() != 1 ) {
				cell->updateNoHit();
			}
		}

#		if CALC_ENTROPY
		if( !m_scanNoConflict ) {
			entropy -= cell->entropy();
		}
#		endif
	}

	return entropy;
}


template<typename Map, int N>
double
SensorModel< Map, BayesCell<N> >::integrate( ScanMap *scanMap, Map *map ) const {
	double entropy = 0.0;

	if( m_scanNoConflict ) {
		// Integrate scan map
		for( const auto &cell : *scanMap ) {
			BayesCell<N> &mapCell = map->cell( cell.first );

#			if CALC_ENTROPY
			entropy -= mapCell.entropy();
#			endif

			mapCell.integrate( cell.second );

#			if CALC_ENTROPY
			entropy += mapCell.entropy();
#			endif
		}
	}

	scanMap->clear();

	return entropy;
}


template<typename Map, int N>
double
SensorModel< Map, BayesCell<N> >::forwardModel( world_t dist, const typename Map::WorldPoint &bestMu, const typename Map::MapPoint &bestMuMap, const typename Map::WorldPoint &sensorOrigin, const typename Map::MapPoint &sensorOriginMap, const Map &map ) const {
	static_assert( std::is_same< typename Map::CellType, BayesCell<N> >(), "Map cell type must be BayesCell<N>."  );
	return NormalDistribution::logPlFromSqr( bestMu.dot( bestMu ), m_bayesSigmaSqr );	// TODO: Save this in a LUT as well? (see Belief sensor model)
}


template<typename Map, int N>
double
SensorModel< Map, BayesCell<N> >::sigmaSqr() const {
	return m_bayesSigmaSqr;
}


template<typename Map, int N>
double
SensorModel< Map, BayesCell<N> >::nullLogLikelihood() const {
	return -0.5 / (2.0 * m_bayesSigmaSqr);
}


template<typename Map, int N>
constexpr bool
SensorModel< Map, BayesCell<N> >::excludeNoHit() const {
	return true;	// always exclude for bayes model
}


template<typename Map, int N>
constexpr int
SensorModel< Map, BayesCell<N> >::cellsBehind() const {
	return 0;
}


} /* namespace bslam */
