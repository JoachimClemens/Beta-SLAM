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

#ifndef BSLAM_CELLBASE_H_
#define BSLAM_CELLBASE_H_

#include "bslam/utils/Convenience.h"

namespace bslam {

class CellBase {
public:
	enum TypeE {
		_CELL_TYPE_MIN_  = -1,
		CELL_POINT_ACCUMULATOR,
		CELL_BAYES,
		CELL_BETA,
		_CELL_TYPE_MAX_
	};

	virtual	inline 	constexpr	TypeE	type() const	{ return TypeE::_CELL_TYPE_MIN_; }
	static	inline	constexpr	TypeE	staticType()	{ return TypeE::_CELL_TYPE_MIN_; }

	static	inline				TypeE	getMapType( std::istream &f );
			inline	constexpr	bool	isDefault() const noexcept 		{ return false;			}

protected:
	template<typename T>
	static  inline				T&			castPtr( uint8_t *data )		{ return *((T*) data);	}
	template<typename T>
	static  inline				const T&	castPtr( const uint8_t *data )	{ return *((T*) data);	}
};


CellBase::TypeE
CellBase::getMapType( std::istream &f ) {
	int		lineNum = 1;

	while( f.good() ) {
		std::string line;
		std::getline( f, line );

		if( line.empty() || line[0] == '#' )
			continue;

		auto splitted = split( line, ':' );

		if( splitted.size() != 2 )
			throw std::runtime_error( to_string( "Wrong number of `:' in line " ) + to_string( lineNum ) );

		if( splitted[0] == "cell" ) {
			f.clear();
			f.seekg( 0, f.beg );
			return (CellBase::TypeE) std::stol( splitted[1] );
		}

		lineNum++;
	}

	f.clear();
	f.seekg( 0, f.beg );

	return CellBase::_CELL_TYPE_MIN_;
}


inline std::string
to_string( const CellBase::TypeE &t ) {
	switch( t ) {
		case CellBase::CELL_POINT_ACCUMULATOR:
			return "point accumulator";

		case CellBase::CELL_BAYES:
			return "bayes";

		case CellBase::CELL_BETA:
			return "beta";

		default:
			return "unknown";
	}
}

} /* namespace bslam */

#endif /* BSLAM_CELLBASE_H_ */
