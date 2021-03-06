/*
 * Software License Agreement (BSD License)
 *
 *  Beta-SLAM - Simultaneous localization and grid mapping with beta distributions
 *  Copyright (c) 2013-2019, Joachim Clemens, Thomas Reineking, Tobias Kluth
 *  All rights reserved.
 *
 *  This file is partially based on the GMapping PointAccumulator class
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

#include <iostream>
namespace bslam {


template<int N>
PointAccumulator<N>::PointAccumulator( int i ) :
	m_acc( nullptr )
{
	// the i is just for compatibility with the map storages
}


template<int N>
PointAccumulator<N>::PointAccumulator( const PointAccumulator &p ) :
	m_acc( p.m_acc ? std::unique_ptr<Acc>( new Acc( *p.m_acc ) ) : nullptr )
{
	// Nothing else to do here
}


template<int N>
PointAccumulator<N>::PointAccumulator( PointAccumulator &&p ) :
	m_acc( std::move( p.m_acc ) )
{
	// Nothing else to do here
}


template<int N>
PointAccumulator<N>::~PointAccumulator() {
	// Nothing to do here
}


template<int N>
PointAccumulator<N>&
PointAccumulator<N>::operator=( const PointAccumulator &p ) {
	if( this != &p ) {
		if( p.m_acc ) {
			if( !m_acc ) {
				m_acc	=  std::unique_ptr<Acc>( new Acc( *p.m_acc ) );
			} else {
				*m_acc	= *p.m_acc;
			}
		} else {
			m_acc = nullptr;
		}
	}
	return *this;
}


template<int N>
PointAccumulator<N>&
PointAccumulator<N>::operator=( PointAccumulator &&p ) {
	if( this != &p ) {
		m_acc = std::move( p.m_acc );
	}
	return *this;
}


template<int N>
void
PointAccumulator<N>::update( const PointNw &p ) {
	if( !m_acc ) {
		m_acc	=  std::unique_ptr<Acc>( new Acc( p.template cast<float>() ) );
	} else {
		*m_acc	+= p.template cast<float>();
	}
}


template<int N>
void
PointAccumulator<N>::integrate( const PointAccumulator<N> &pa ) {
	if( pa.m_acc ) {
		if( !m_acc ) {
			m_acc	=  std::unique_ptr<Acc>( new Acc( *pa.m_acc ) );
		} else {
			*m_acc	+= *pa.m_acc;
		}
	}
}


template<int N>
typename PointAccumulator<N>::PointNw
PointAccumulator<N>::mean() const {
	assert( m_acc );
	assert( m_acc->hits );
	return m_acc->sum.template cast<world_t>() / m_acc->hits;
}


template<int N>
uint32_t
PointAccumulator<N>::hits() const {
	return m_acc ? m_acc->hits : 0;
}


template<int N>
bool
PointAccumulator<N>::operator==( const PointAccumulator<N> &other ) const noexcept {
	if( (bool) m_acc != (bool) other.m_acc )
		return false;

	return !m_acc || (m_acc->hits == other.m_acc->hits && m_acc->sum == other.m_acc->sum);
}


template<int N>
bool
PointAccumulator<N>::operator!=( const PointAccumulator<N> &other ) const noexcept {
	return !(*this == other);
}


template<int N>
bool
PointAccumulator<N>::prune() const noexcept {
	return !m_acc;
}


template<int N>
bool
PointAccumulator<N>::pruneEqual( const PointAccumulator<N> &other ) const noexcept {
	return !other.m_acc;
}


template<int N>
bool
PointAccumulator<N>::isDefault() const noexcept {
	return !m_acc;
}


template<int N>
void
PointAccumulator<N>::fromBinary( const uint8_t *data ) {
	uint32_t hits = CellBase::castPtr<uint32_t>( data );

	if( hits ) {
		if( !m_acc ) {
			m_acc = std::unique_ptr<Acc>( new Acc() );
		}

		m_acc->hits = hits;

		data += sizeof( uint32_t );
		for( size_t i = 0; i < N; i++ ) {
			m_acc->sum[i] = CellBase::castPtr<float>( data );
			data += sizeof( float );
		}
	} else {
		m_acc = nullptr;
	}
}


template<int N>
void
PointAccumulator<N>::toBinary( uint8_t *data ) const {
	if( m_acc ) {
		CellBase::castPtr<uint32_t>( data ) = m_acc->hits;

		data += sizeof( uint32_t );
		for( size_t i = 0; i < N; i++ ) {
			CellBase::castPtr<float>( data ) = m_acc->sum[i];
			data += sizeof( float );
		}
	} else {
		CellBase::castPtr<uint32_t>( data ) = 0;
	}
}


template<int N>
size_t
PointAccumulator<N>::binarySize() {
	//		sum						hits
	return	sizeof( float ) * N +	sizeof( uint32_t );
}


template<int N>
size_t
PointAccumulator<N>::bytes() const noexcept {
	return sizeof( PointAccumulator<N> ) + (m_acc ? sizeof( Acc ) : 0);
}


} /* namespace bslam */
