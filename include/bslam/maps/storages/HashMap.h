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

#ifndef BSLAM_HASHMAP_H_
#define BSLAM_HASHMAP_H_

#include "bslam/utils/geometry/Point.h"
#include "bslam/utils/aligned_unordered_map.h"
#include "bslam/utils/aligned_list.h"
#include "bslam/maps/CellSelectionFunction.h"


namespace bslam {

template<typename Cell, int N = Cell::Dimension>
class HashMap {
public:
	static constexpr int	Dimension 	= N;
	template<typename CellOut = Cell>
	using					CellList	= aligned_list< std::pair<Pointm<N>, const CellOut &> >;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

								HashMap( const Pointm<N> &size );
								HashMap( const HashMap &other );
								HashMap( HashMap &&other );
								~HashMap();

	inline	void				clear();

	inline	HashMap&			operator=( const HashMap &other );
	inline	HashMap&			operator=( HashMap &&other );

	inline 	const Cell&			cell( const Pointm<N>& p ) const;
	inline 	Cell& 				cell( const Pointm<N>& p );

	template<typename CellOut = Cell>
	inline	CellList<CellOut>	cells( const CellSelectionFunction<Cell> &selectionFunction = &selectAllCells<Cell>, const Pointm<N> &offset = Pointm<N>::Zero() ) const;

	inline	uint32_t			prune()									{ return 0;			}	// nothing to prune

	inline	size_t				bytes() const;

	inline 	AccessibilityState 	cellState( const Pointm<N>& p ) const;
	inline 	bool 				isInside( const Pointm<N>& p ) const;
	inline 	bool 				isAllocated( const Pointm<N> &p ) const;

	inline 	Pointm<N>			size() const							{ return m_size;	}

protected:
	Pointm<N>								m_size;
	aligned_unordered_map<Pointm<N>, Cell>	m_map;
    static const Cell						sm_emptyValue;
};

} /* namespace bslam */

#include "HashMap.hpp"

#endif /* BSLAM_HASHMAP_H_ */
