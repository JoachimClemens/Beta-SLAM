/*
 * Software License Agreement (BSD License)
 *
 *  Beta-SLAM - Simultaneous localization and grid mapping with beta distributions
 *  Copyright (c) 2013-2019, Joachim Clemens, Thomas Reineking, Tobias Kluth
 *  All rights reserved.
 *
 *  This file is partially based on the GMapping HArray2D class
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


namespace bslam {

template<class Cell>
Array2D<Cell>::Array2D( const Point2m &size ) :
	m_size( size )
{
	if( m_size[0] > 0 && m_size[1] > 0 ) {
		m_cells = new Cell[m_size[0] * m_size[1]];
	} else {
		m_size[0] = m_size[1] = 0;
		m_cells = nullptr;
	}
}


template<class Cell>
Array2D<Cell>::Array2D( map_t xsize, map_t ysize ) :
	m_size( xsize, ysize )
{
	if( m_size[0] > 0 && m_size[1] > 0 ) {
		m_cells = new Cell[m_size[0] * m_size[1]];
	} else {
		m_size[0] = m_size[1] = 0;
		m_cells = nullptr;
	}
}


template<class Cell>
Array2D<Cell> &
Array2D<Cell>::operator=( const Array2D<Cell> &g ) {
	if( m_size[0] != g.m_size[0] || m_size[1] != g.m_size[1] ) {
		if( m_cells )
			delete[] m_cells;

		m_size[0] = g.m_size[0];
		m_size[1] = g.m_size[1];

		m_cells = new Cell[m_size[0] * m_size[1]];
	}

	for( map_t i = 0; i < m_size[0]*m_size[1]; i++ )
		m_cells[i] = g.m_cells[i];

	return *this;
}


template<class Cell>
Array2D<Cell> &
Array2D<Cell>::operator=( Array2D<Cell> &&g ) {
	if( this != &g ) {
		if( m_cells )
			delete[] m_cells;

		m_size	 	= g.m_size;
		m_cells		= g.m_cells;
		g.m_cells	= nullptr;
	}
	return *this;
}


template<class Cell>
void
Array2D<Cell>::swap( Array2D<Cell> &other ) {
	if( this != &other ) {
		SWAP( m_cells, other.m_cells );
		SWAP( m_size, other.m_size );
	}
}


template<class Cell>
Array2D<Cell>::Array2D( const Array2D<Cell> &g ) :
	m_size( g.m_size )
{
	m_cells = new Cell[m_size[0] * m_size[1]];
	for( map_t i = 0; i < m_size[0]*m_size[1]; i++ )
		m_cells[i] = g.m_cells[i];
}


template<class Cell>
Array2D<Cell>::Array2D( Array2D<Cell> &&g ) :
	m_size( g.m_size ),
	m_cells( g.m_cells )
{
	g.m_cells = nullptr;
}


template<class Cell>
Array2D<Cell>::~Array2D() {
	if( m_cells ) {
		delete[] m_cells;
		m_cells = nullptr;
	}
}


template<class Cell>
void
Array2D<Cell>::clear() {
	for( map_t i = 0; i < m_size[0]*m_size[1]; i++ )
		m_cells[i] = Cell( 0 );
}


template<class Cell>
void
Array2D<Cell>::resize( map_t xmin, map_t ymin, map_t xmax, map_t ymax ) {
	map_t xsize = xmax - xmin;
	map_t ysize = ymax - ymin;
	Cell *newcells = new Cell[xsize * ysize];

	map_t dx = xmin < 0 ? 0 : xmin;
	map_t dy = ymin < 0 ? 0 : ymin;
	map_t Dx = xmax < m_size[0] ? xmax : m_size[0];
	map_t Dy = ymax < m_size[1] ? ymax : m_size[1];

	for (map_t x = dx; x < Dx; x++)
		for (map_t y = dy; y < Dy; y++)
			newcells[index( x - xmin, y - ymin )] = m_cells[index( x, y )];

	delete[] m_cells;
	m_cells = newcells;
	m_size[0] = xsize;
	m_size[1] = ysize;
}


template<class Cell>
bool
Array2D<Cell>::isInside( map_t x, map_t y ) const {
	return x >= 0 && y >= 0 && x < m_size[0] && y < m_size[1];
}


template<class Cell>
const Cell&
Array2D<Cell>::cell( map_t x, map_t y ) const {
	assert( isInside( x, y ) );
	return m_cells[index( x, y )];
}


template<class Cell>
Cell&
Array2D<Cell>::cell( map_t x, map_t y ) {
	assert( isInside( x, y ) );
	return m_cells[index( x, y )];
}


template<class Cell>
template<typename CellOut>
typename Array2D<Cell>::template CellList<CellOut>
Array2D<Cell>::cells( const CellSelectionFunction<Cell> &selectionFunction, const Point2m& offset ) const {
	CellList<CellOut> res;

	for( map_t x = 0; x < m_size[0]; x++ ) {
		for( map_t y = 0; y < m_size[1]; y++ ) {
			const Cell &c = m_cells[index( x, y )];

			if( selectionFunction( c ) && !c.isDefault() ) {
				res.push_back( std::pair<Point2m, const CellOut &>( offset + Point2m( x, y ), c ) );
			}
		}
	}

	return res;
}


template<class Cell>
constexpr size_t
Array2D<Cell>::index( map_t x, map_t y ) const {
	return x*m_size[1] + y;
}


} /* namespace bslam  */
