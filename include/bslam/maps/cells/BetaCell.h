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

#ifndef BSLAM_BETA_CELL_H_
#define BSLAM_BETA_CELL_H_

#include "BayesCell.h"
#include "bslam/utils/uncertainty/BetaDistribution.h"


namespace bslam {


template<int N>
class BetaCell : public BayesCell<N> {
public:
	using typename BayesCell<N>::PointNw;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	inline					BetaCell( int i = 0 );
	inline					BetaCell( const std::string &str );
	inline virtual			~BetaCell();

	inline void				updateNoHit( const BetaDistribution &dist, bool behind = false );
	inline void				updateHit( const BetaDistribution &dist, const PointNw &p = PointNw::Zero() );
	inline void				integrate( const BetaCell &c );

	//inline double			fullness() const;
	inline double			entropy() const;

	inline float			alpha() const;
	inline float			beta() const;
	//inline float			n() const;
	//inline float			k() const;

	inline bool				operator==( const BetaCell &other ) const noexcept;
	inline bool				operator!=( const BetaCell &other ) const noexcept;

	inline bool				prune() const noexcept;
	inline bool				pruneEqual( const BetaCell &other ) const noexcept;
	inline bool				isDefault() const noexcept;

	inline void				fromStr( const std::string &str );
	inline std::string		toStr() const;

	inline void				fromBinary( const uint8_t *data );
	inline void				toBinary( uint8_t *data ) const;
	static inline size_t	binarySize();

	virtual inline size_t	bytes() const noexcept;

	virtual	inline 	constexpr CellBase::TypeE	type() const	{ return CellBase::TypeE::CELL_BETA; }
	static	inline	constexpr CellBase::TypeE	staticType()	{ return CellBase::TypeE::CELL_BETA; }

    friend std::ostream& operator<<( std::ostream& ostr, const BetaCell* c ){
        return operator<<( ostr,(*c) );
    }

    friend std::ostream& operator<<( std::ostream& ostr, const BetaCell& c ){
        ostr << "TODO"; // TODO
        return ostr;
    }


protected:
    using Acc = typename PointAccumulator<N>::Acc;

    inline	void			updateDistribution( const BetaDistribution &dist );

    // Prior beta distribution
    static	 float	sm_priorAlpha,
    				sm_priorBeta;

    BetaDistribution	m_distribution;

public:
    S_SETTER( priorAlpha );
    S_SETTER( priorBeta );

    GETTER( distribution );
    SETTER( distribution );
};


} /* namespace bslam */

#include "BetaCell.hpp"

#endif /* BSLAM_BETA_CELL_H_ */
