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

#include "bslam/utils/aligned_unordered_set.h"

namespace bslam {

// Declaration of specialized methods

template<>
inline void
Circle<2>::horn( map_t radius, const Point2m &center );

/*
template<>
inline void
Circle<3>::horn( const Pointm<3> &center, map_t radius ); // this has to draw a sphere, of course
*/


template<int N>
Circle<N>::Circle( map_t radius, const Pointm<N> &center ) {
	horn( radius, center );
}


template<int N>
Circle<N>::~Circle() {
	// Nothing to do here
}


template<int N>
inline void
Circle<N>::horn(  map_t radius, const Pointm<N> &center ) {
	static_assert( N != 2 && N != 3, "Specialization needed for this N." );
}


template<>
inline void
Circle<2>::horn(  map_t radius, const Point2m &center ) {
	// circle drawing method of Horn (http://de.wikipedia.org/wiki/Rasterung_von_Kreisen#Methode_von_Horn)

	map_t 	d = -radius,
			x = radius,
			y = 0;

	// use a set here, in order to have unique points only
	aligned_unordered_set<Point2m> points;

	while( y <= x ) {
		// reflect the current point in each direction
		points.emplace( Point2m(  x,  y ) + center );
		points.emplace( Point2m(  y,  x ) + center );
		points.emplace( Point2m( -x, -y ) + center );
		points.emplace( Point2m( -y, -x ) + center );
		points.emplace( Point2m( -x,  y ) + center );
		points.emplace( Point2m( -y,  x ) + center );
		points.emplace( Point2m(  x, -y ) + center );
		points.emplace( Point2m(  y, -x ) + center );

		d = d + 2*y + 1;
		y = y + 1;
		if( d > 0 ) {
			d = d - 2*x + 2;
			x = x - 1;
		}
	}

	// copy points from set to vector
	assign( points.begin(), points.end() );
}

} /* namespace bslam */

