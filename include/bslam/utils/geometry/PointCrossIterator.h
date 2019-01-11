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

#ifndef BSLAM_POINTCROSSITERATOR_H_
#define BSLAM_POINTCROSSITERATOR_H_

#include "Point.h"

namespace bslam {

template<int N>
class PointCrossIterator {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

								PointCrossIterator( const Pointm<N> &start, const Pointm<N> &end );
								PointCrossIterator( const Pointm<N> &start, const Pointm<N> &end, const Pointm<N> &inc );

	inline						operator bool() const 	{ return !m_endReached; }
	inline	void				operator++(int);
	inline	void				operator++();
	inline	const Pointm<N>&	operator*() const		{ return m_cur; }

private:
	Pointm<N>	m_start,
				m_cur,
				m_end,
				m_inc,
				m_mid;
	bool		m_endReached;
	int			m_idx,
				m_midIdx;
};

} /* namespace bslam */

#include "PointCrossIterator.hpp"

#endif /* BSLAM_POINTCROSSITERATOR_H_ */
