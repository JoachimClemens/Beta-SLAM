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
 * Beta sensor model
 ************************************/

template<typename Map, int N>
SensorModel< Map, BetaCell<N> >::SensorModel( bool unused )
{
	setParamsConfig();
}


template<typename Map, int N>
void
SensorModel< Map, BetaCell<N> >::setParamsConfig() {
	m_mapDelta 						= Config::getDouble(	"MAP_DELTA", 					0.05 	);
	m_betaRandom					= Config::getDouble(	"BETA_RANDOM",					0.2 	);
	m_betaPriorP					= Config::getDouble(	"BETA_PRIOR_P",					0.5		);
	m_betaExcludeNoHit				= Config::getInt(		"BETA_EXCLUDE_NO_HIT",			0 		);
	m_scanNoConflict				= Config::getInt( 		"SCAN_NO_CONFLICT",				0		);

	m_betaSigmaForward 				= Config::getDouble(	"BETA_SIGMA_FORWARD", 			0.05 );
	m_betaSigmaForwardSqr 			= m_betaSigmaForward*m_betaSigmaForward;
	m_betaSigmaForwardFactor		= Config::getDouble(	"BETA_SIGMA_FORWARD_FACTOR", 	-1.0 );
	if( m_betaSigmaForwardFactor < 0 ) {
		m_betaSigmaForwardFactor = 0.5 * m_mapDelta / m_betaSigmaForward;	// results in an integration range of 0.5*MAP_DELTA
	} else {
		l_wrn( "BETA_SIGMA_FORWARD_FACTOR is >=0, but should be <0 (feature disabled and use MAP_DELTA/2 instead) for a correct theoretical foundation." );
	}

	m_betaSigmaInverse 				= Config::getDouble(	"BETA_SIGMA_INVERSE", 			0.0125 );
	m_betaSigmaInverseSqr 			= m_betaSigmaInverse*m_betaSigmaInverse;
	m_betaSigmaInverseFactor		= Config::getDouble(	"BETA_SIGMA_INVERSE_FACTOR", 	-1.0 );
	if( m_betaSigmaInverseFactor < 0 ) {
		m_betaSigmaInverseFactor = 0.5 * m_mapDelta / m_betaSigmaInverse;	// results in an integration range of 0.5*MAP_DELTA
	}

	if( m_betaSigmaForward != m_betaSigmaInverse ) {
		l_wrn( "BETA_SIGMA_FORWARD and BETA_SIGMA_INVERSE are different, but should be the same for a correct theoretical foundation." );
	}

	double betaForwardSigmaRange	= Config::getDouble(	"BETA_FORWARD_SIGMA_RANGE", 	3.0 );
	m_betaForwardCellRange			= std::ceil( betaForwardSigmaRange * m_betaSigmaForward / m_mapDelta );

	BetaCell<N>::priorAlpha()		= Config::getDouble( 	"BETA_PRIOR_ALPHA",				1.0 );
	BetaCell<N>::priorBeta()		= Config::getDouble( 	"BETA_PRIOR_BETA",				1.0 );

	computeBetaLUT();
	computePdfLUT();
}


template<typename Map, int N>
double
SensorModel< Map, BetaCell<N> >::inverseModel( const typename Map::WorldPoint &pHit, const Line<N> &scanLine, size_t hitIdx, ScanMap *scanMap, Map *map ) const {
	static_assert( std::is_same< typename Map::CellType, BetaCell<N> >(), "Map cell type must be BetaCell<N>."  );

	double entropy = 0.0;

	for( size_t i = 0; i < scanLine.size(); i++ ) {
		BetaCell<N> *cell;

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

		double		dist 				= i == hitIdx ? 0 : ( map->map2world( scanLine[i] ) - pHit ).norm();	// distance between end point and cell or zero if hit cell
		const auto 	&betaDistribution	= getBeta( dist );

		if( i == hitIdx ) {
			if( m_scanNoConflict && cell->distribution().beta() != cell->priorBeta() ) {
				// Reset cell, if it is not strictly occupied
				(*cell) = BetaCell<N>();
			}

			cell->updateHit( betaDistribution, pHit );
		} else {
			if( m_scanNoConflict && cell->distribution().beta() != cell->priorBeta() && betaDistribution.alpha() != cell->priorAlpha() ) {
				// Reset cell, if it is not strictly occupied and should receive mass on occupied
				(*cell) = BetaCell<N>();
			}

			if( !m_scanNoConflict || cell->distribution().alpha() == cell->priorAlpha() || betaDistribution.beta() == cell->priorBeta() ) {
				// Update only if cell is not occupied or it should receive only mass on occupied
				cell->updateNoHit( betaDistribution, (i > hitIdx) ); // i > hitIdx -> cell behind measurement
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
SensorModel< Map, BetaCell<N> >::integrate( ScanMap *scanMap, Map *map ) const {
	double entropy = 0.0;

	if( m_scanNoConflict ) {
		// Integrate scan map
		for( const auto &cell : *scanMap ) {
			BetaCell<N> &mapCell = map->cell( cell.first );

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
SensorModel< Map, BetaCell<N> >::forwardModel( world_t dist, const typename Map::WorldPoint &bestMu, const typename Map::MapPoint &bestMuMap, const typename Map::WorldPoint &sensorOrigin, const typename Map::MapPoint &sensorOriginMap, const Map &map ) const {
	static_assert( std::is_same< typename Map::CellType, BetaCell<N> >(), "Map cell type must be BetaCell<N>."  );

	// Direction from sensor origin to measurement in world coordinates
	typename Map::WorldPoint 	direction = bestMuMap.template cast<world_t>() - sensorOriginMap.template cast<world_t>();
	world_t directionNorm = direction.norm();
	direction /= directionNorm;

	// vector of cells to consider along the scan in map coordinates
	typename Map::MapPoint		cellDelta = ( direction * m_betaForwardCellRange * 2 ).template cast<map_t>();

	typename Map::MapPoint 		originMap( sensorOriginMap );
	if( m_betaForwardCellRange * 2 * m_mapDelta < directionNorm ) {
		originMap = bestMuMap - cellDelta;
	}

	// all relevant cells from laser pose to beam end point (or, if scan matching was successful, to best hit cell)
	Line<N>	line( originMap, bestMuMap );
	size_t hitIdx = line.size() - 1;

	// cells behind measurement
	line.extend( bestMuMap + cellDelta );

	// recalculate start and end with correct number of cells because we added a margin in line calculation
	size_t 	start	= std::max( (int64_t) hitIdx - m_betaForwardCellRange, 		(int64_t) 0 			),
			M		= std::min( (int64_t) hitIdx + m_betaForwardCellRange + 1, 	(int64_t) line.size() 	);

	double 	pSum = 0.0,
			cellDist;
	for( size_t i = start; i < M; i++ ) {
		const BetaCell<N> &cell = map.cell( line[i] );

		if( i == hitIdx && cell.hits() ) {
			cellDist = (sensorOrigin - cell.mean()).norm();
		} else {
			cellDist = (sensorOrigin - map.map2world( line[i] )).norm();
		}

		pSum += getPdf( dist - cellDist ) * BetaDistribution::mean( cell.alpha(), cell.beta() );

		//DEBUG_PRINT_COUT( dist << " <-> " << cellDist << ": " << "pdf " << NormalDistribution::pdf( dist, cellDist, m_betaSigmaSqr ) << ", beta " <<  BetaDistribution::mean( cell.alpha(), cell.beta() ) );
	}

	return log( pSum );
}


template<typename Map, int N>
double
SensorModel< Map, BetaCell<N> >::sigmaSqr() const {
	return m_betaSigmaForwardSqr;
}


template<typename Map, int N>
double
SensorModel< Map, BetaCell<N> >::nullLogLikelihood() const {
	return -0.5 / (2.0 * m_betaSigmaForwardSqr);
}


template<typename Map, int N>
constexpr bool
SensorModel< Map, BetaCell<N> >::excludeNoHit() const {
	return m_betaExcludeNoHit;
}


template<typename Map, int N>
constexpr int
SensorModel< Map, BetaCell<N> >::cellsBehind() const {
	return 0;
}


template<typename Map, int N>
void
SensorModel< Map, BetaCell<N> >::computeBetaLUT() {
	m_betaDistScale = 1000.0;							// default scaling factor

	world_t	maxDist		= m_betaSigmaInverse * 6;				// Calculate the values up to 6 sigma. Before and after this point, the values are aprox. equal to zero
	size_t	maxDistIdx 	= maxDist * m_betaDistScale;

	if( maxDistIdx > MAX_SIZE ) {
		maxDistIdx 			= MAX_SIZE;
		m_betaDistScale 	= maxDistIdx / maxDist;
		l_wrn( "Distance scale for beta LUT was reduced to " << m_betaDistScale << " in order to limit LUT size to " << MAX_SIZE << " entries." );
	}

	l_dbg( "Beta LUT size " << maxDistIdx << ", distance scale " << m_betaDistScale << "." );

	m_betaLUT.resize( maxDistIdx );

	const double	scale = MAX( BetaCell<N>::priorAlpha(), BetaCell<N>::priorBeta() );
	for( size_t i = 0; i < maxDistIdx; i++ ) {
		double	dist	= i / m_betaDistScale,
				pHit	= NormalDistribution::cdf( dist + m_betaSigmaInverseFactor*m_betaSigmaInverse, m_betaSigmaInverseSqr ) - NormalDistribution::cdf( dist - m_betaSigmaInverseFactor*m_betaSigmaInverse, m_betaSigmaInverseSqr ),
				pi		= m_betaRandom * m_betaPriorP + (1 - m_betaRandom) * pHit,
				c		= scale / MIN( pi, 1 - pi );
		m_betaLUT[i]	= BetaDistribution( c * pi, c * (1 - pi) );
	}

	double	piZero	= m_betaRandom * m_betaPriorP,
			cZero	= scale / MIN( piZero, 1 - piZero );

	m_betaZero.alpha()	= cZero * piZero;
	m_betaZero.beta()	= cZero * (1 - piZero);
}


template<typename Map, int N>
const BetaDistribution&
SensorModel< Map, BetaCell<N> >::getBeta( world_t dist ) const {
	size_t	distIdx = fabs( dist )*m_betaDistScale + 0.5;	// 0.5 is to round the value correctly

	if( distIdx >= m_betaLUT.size() ) {
		return m_betaZero;
	} else {
		return m_betaLUT[distIdx];
	}
}


template<typename Map, int N>
void
SensorModel< Map, BetaCell<N> >::computePdfLUT() {
	m_pdfDistScale = 1000.0;							// default scaling factor

	world_t	maxDist		= m_betaSigmaForward * 6;				// Calculate the values up to 6 sigma. Before and after this point, the values are aprox. equal to zero
	size_t	maxDistIdx 	= maxDist * m_pdfDistScale;

	if( maxDistIdx > MAX_SIZE ) {
		maxDistIdx 		= MAX_SIZE;
		m_pdfDistScale 	= maxDistIdx / maxDist;
		l_wrn( "Distance scale for PDF LUT was reduced to " << m_pdfDistScale << " in order to limit LUT size to " << MAX_SIZE << " entries." );
	}

	l_dbg( "PDF LUT size " << maxDistIdx << ", distance scale " << m_pdfDistScale << "." );

	m_pdfLUT.resize( maxDistIdx );

	for( size_t i = 0; i < maxDistIdx; i++ ) {
		double dist = i / m_pdfDistScale;
		m_pdfLUT[i] =	NormalDistribution::cdf( dist + m_betaSigmaForwardFactor*m_betaSigmaForward, m_betaSigmaForwardSqr ) - NormalDistribution::cdf( dist - m_betaSigmaForwardFactor*m_betaSigmaForward, m_betaSigmaForwardSqr );
		m_pdfLUT[i] /=	m_mapDelta;	// scale with 1/Delta r_i
	}
}


template<typename Map, int N>
double
SensorModel< Map, BetaCell<N> >::getPdf( world_t dist ) const {
	size_t	distIdx = fabs( dist )*m_pdfDistScale + 0.5;	// 0.5 is to round the value correctly

	if( distIdx >= m_pdfLUT.size() ) {
		return 0.0;
	} else {
		return m_pdfLUT[distIdx];
	}
}



} /* namespace bslam */
