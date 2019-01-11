/*
 * Software License Agreement (BSD License)
 *
 *  Beta-SLAM - Simultaneous localization and grid mapping with beta distributions
 *  Copyright (c) 2013-2019, Joachim Clemens, Thomas Reineking, Tobias Kluth
 *  All rights reserved.
 *
 *  This file is partially based on the GMapping Array2D class
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

#ifndef BSLAM_ARRAY2D_H_
#define BSLAM_ARRAY2D_H_

#include <iostream>

#include <assert.h>

#include "bslam/utils/geometry/Point.h"
#include "bslam/utils/aligned_list.h"
#include "bslam/maps/AccessState.h"
#include "bslam/maps/CellSelectionFunction.h"


namespace bslam {

template<class Cell>
class Array2D {
public:
	static constexpr int	Dimension 			= 2;
	template<typename CellOut = Cell>
	using					CellList			= aligned_list< std::pair<Point2m, const CellOut &> >;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

								Array2D( const Point2m &size );
								Array2D( map_t xsize = 0, map_t ysize = 0 );
								Array2D( const Array2D<Cell> &other );
								Array2D( Array2D<Cell> &&other );
								~Array2D();

			void 				clear();
			void 				resize( map_t xmin, map_t ymin, map_t xmax, map_t ymax );

			Array2D&			operator=( const Array2D &other );
			Array2D&			operator=( Array2D &&other );
			void				swap( Array2D &other );

	inline 	const Cell& 		cell( map_t x, map_t y ) const;
	inline 	Cell& 				cell( map_t x, map_t y );
	inline 	const Cell&			cell( const Point2m& p ) const 		{	return cell( (map_t) p[0], (map_t) p[1] );	}
	inline 	Cell& 				cell( const Point2m& p ) 			{	return cell( (map_t) p[0], (map_t) p[1] );	}

	template<typename CellOut = Cell>
	inline 	CellList<CellOut>	cells( const CellSelectionFunction<Cell> &selectionFunction = &selectAllCells<Cell>, const Point2m &offset = Point2m::Zero() ) const;

	inline	constexpr uint32_t	prune() const 						{ 	return 0; } // Array2D can't be pruned

	inline	size_t				bytes() const						{	return sizeof(*this) + sizeof(Cell*) * m_size[0] + sizeof(Cell) * m_size[0] * m_size[1]; }

	inline 	bool 				isInside( map_t x, map_t y ) const;
	inline 	bool 				isInside( const Point2m &p ) const 	{	return isInside( (map_t) p[0], (map_t) p[1] );	}

	inline 	AccessibilityState 	cellState( map_t x, map_t y ) const {	return (AccessibilityState) (isInside( x, y ) ? (AS_INSIDE | AS_ALLOCATED) : AS_OUTSIDE);	}
	inline 	AccessibilityState	cellState( const Point2m &p ) const	{	return cellState( (map_t) p[0], (map_t) p[1] );	}

	inline 	int 				getXSize() const 					{	return m_size[0];		}
	inline 	int 				getYSize() const 					{	return m_size[1];		}
	inline 	Point2m				size() const						{ 	return m_size; 			}
	inline	Point2m				patchSize() const					{ 	return size(); 			}
	inline	Point2m				numPatches() const					{ 	return Point3m::Ones();	}

	//inline Cell** 			cells() 										{	return m_cells;	}

protected:
	inline	constexpr	size_t	index( map_t x, map_t y ) const;

	Cell 	*m_cells;
	Point2m	m_size;
};

} /* namespace bslam */

#include "Array2D.hpp"

#endif /* BSLAM_ARRAY_2D_H_ */

