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

#define BETA_UNCERTAINTY 	0  	// 0=gray, 1=blue, 2=black (only for BETA_CELL_COLOR=0)

namespace bslam {

/*********************************
 * BetaCell
 *********************************/

template<int N>
Color
CellColor< BetaCell<N> >::getColor( const BetaCell<N> &cell, bool normalize, uint8_t *colorAlpha, ColorModeE colorMode, double alphaScale ) {
	const BetaDistribution &betaDist = cell.distribution();

#if BETA_CELL_COLOR == 1
	double	r = 1 - exp( -betaDist.alpha() / 10 ),
			g = 1 - exp( -betaDist.beta() / 10 );

	if( colorAlpha ) {
		switch( colorMode ) {
			case COLOR_MODE_ALPHA_THETA:
				*colorAlpha = MAX( 0.0, (alphaScale - (2 - r - g) ) ) * 255;  	// 1 - (2 - occ - free) -> Alpha channel
				break;

			case COLOR_MODE_ALPHA_EMPTY:
				*colorAlpha = MAX( 0.0, (alphaScale - (1.0 - g) ) ) * 255;		// 1 - free -> Alpha channel
				break;

			case COLOR_MODE_ALPHA_NOT_OCCUPIED:
				*colorAlpha = r * 255;		 									// occupied -> Alpha channel
				break;

			case COLOR_MODE_DEFAULT:
			default:
				*colorAlpha = 255;            									// Full Alpha
				break;
		}
	}

	return {	(uint8_t)( r * 255 	),
				(uint8_t)( g * 255 	),
				(uint8_t)( 0 		)	};

#else // BETA_CELL_COLOR == 0
	static const double maxVar	= BetaDistribution::var( BetaCell<N>::priorAlpha(), BetaCell<N>::priorBeta() );
	static const double minLog	= log2( 0.0001 );

	double	mean	= betaDist.mean(),
			var		= MIN( 1.0, betaDist.var() / maxVar );	// scale from 0 to 1

	var = 1.0 - log2( 0.0001 + 0.9999 * var ) / minLog;

	/*
	// mean => occ -> free => red -> yellow -> green
	// var => uncertainty => gray -> color
	double	h		= (1.0 - mean) * 120,		// from red (occ) over yellow to green (free)
			s		= (1.0 - var),				// together with v from gray to full color
			v		= 0.5 + (1.0 - var) * 0.5;	// together with s from gray to full color

	return Color::hsv2rgb( h, s, v );
	*/

	// /*
	// mean => occ -> free => red -> yellow -> green
	double 	r = mean < 0.5 ? 2 * mean 			: 1.0,
			g = mean > 0.5 ? 2 * (1.0 - mean)	: 1.0;
	// */

	/*
	// mean => occ -> free => red -> green
	double 	r = mean,
			g = 1.0 - mean;
	*/

	if( colorAlpha ) {
		switch( colorMode ) {
			case COLOR_MODE_ALPHA_THETA:
				*colorAlpha = (alphaScale - var) * 255;  						// 1 - var -> Alpha channel
				break;

			case COLOR_MODE_ALPHA_EMPTY:
				*colorAlpha = MAX( 0.0, (alphaScale - (1.0 - mean) ) ) * 255;	// 1 - free -> Alpha channel
				break;

			case COLOR_MODE_ALPHA_NOT_OCCUPIED:
				*colorAlpha = mean * 255;		 								// occupied -> Alpha channel
				break;

			case COLOR_MODE_DEFAULT:
			default:
				*colorAlpha = 255;            									// Full Alpha
				break;
		}
	}

#if BETA_UNCERTAINTY == 0
	// var => uncertainty => gray -> color
	return {	(uint8_t)( ( r * (1.0 - var) + 	var/2.0 ) * 255 ),
				(uint8_t)( ( g * (1.0 - var) + 	var/2.0 ) * 255 ),
				(uint8_t)( (					var/2.0 ) * 255 )	};
#elif BETA_UNCERTAINTY == 1
	// var => uncertainty => blue -> color
	return {	(uint8_t)( ( r * (1.0 - var) 	) * 255 ),
				(uint8_t)( ( g * (1.0 - var) 	) * 255 ),
				(uint8_t)( ( var 				) * 255 )	};
#elif BETA_UNCERTAINTY == 2
	// var => uncertainty => black -> color
	return {	(uint8_t)( ( r * (1.0 - var) 	) * 255 ),
				(uint8_t)( ( g * (1.0 - var) 	) * 255 ),
				(uint8_t)( ( 0.0 				) * 255 )	};
#else
	static_assert( false, "Invalid value for BETA_UNCERTAINTY!" );
#endif

#endif // BETA_CELL_COLOR
}


template<int N>
constexpr Color
CellColor< BetaCell<N> >::defaultColor() {
#if BETA_CELL_COLOR == 1 or BETA_UNCERTAINTY == 2
	return {   0,   0,   0 };
#elif BETA_UNCERTAINTY == 0
	return { 127, 127, 127 };
#elif BETA_UNCERTAINTY == 1
	return {   0,   0, 255 };
#endif
}


template<int N>
constexpr Color
CellColor< BetaCell<N> >::pathColor() {
#if BETA_UNCERTAINTY == 1
	return { 255, 255, 0 };
#else
	return { 0, 127, 255 };
#endif
}


/*********************************
 * BayesCell
 *********************************/

template<int N>
Color
CellColor< BayesCell<N> >::getColor( const BayesCell<N> &cell, bool normalize, uint8_t *alpha, ColorModeE colorMode, double alphaScale ) {
	double v = cell.fullness();

	if( v < 0 ) {
		if( alpha )
			*alpha = 0;
		return defaultColor();
	}

	if( alpha ) {
		switch( colorMode ) {
			case COLOR_MODE_ALPHA_EMPTY:
				*alpha = MAX( 0.0, (alphaScale - 1.0+v) ) * 255;
				break;

			case COLOR_MODE_ALPHA_NOT_OCCUPIED:
				*alpha = v * 255;
				break;

			case COLOR_MODE_ALPHA_THETA:
			case COLOR_MODE_DEFAULT:
			default:
				*alpha = 255;
				break;
		}
	}

	uint8_t grayVal = (uint8_t) (255 - (uint8_t) (255 * v));
	return { grayVal, grayVal, grayVal };
}


template<int N>
constexpr Color
CellColor< BayesCell<N> >::defaultColor() {
	// return { 200, 200, 255 };	// bluish color as used by GMapping
	return { 125, 125, 125 };	// gray as correct for uniform distribution
}


template<int N>
constexpr Color
CellColor< BayesCell<N> >::pathColor() {
	return { 255, 0, 0 };
}

} /* namespace bslam */

