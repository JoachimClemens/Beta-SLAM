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

#ifndef BSLAM_POSESE2_HPP_
#define BSLAM_POSESE2_HPP_

#include "PoseBase.hpp"


namespace bslam {

PoseSE2::PoseSE2() :
	m_translation( Point2w::Zero() )
{
	// Nothing else to do here
}


PoseSE2::PoseSE2( PoseSE2 &&other ) :
	m_translation( std::move( other.m_translation ) ),
	m_rotation( std::move( other.m_rotation ) )
{
	// Nothing else to do here
}


PoseSE2::PoseSE2( const PoseSE2 &other ) :
	m_translation( other.m_translation ),
	m_rotation( other.m_rotation )
{
	// Nothing else to do here
}


PoseSE2::PoseSE2( world_t x, world_t y, double phi_ ) :
	m_translation( x, y ),
	m_rotation( phi_ )
{
	// Nothing else to do here
}


PoseSE2::PoseSE2( const Point2w &pos_, double phi_ ) :
	m_translation( pos_ ),
	m_rotation( phi_ )
{
	// Nothing else to do here
}


PoseSE2::PoseSE2( const EigenIsometry2w &other ) {
	*this = other;
}


PoseSE2::~PoseSE2() {
	// Nothing to do here
}


PoseSE2&
PoseSE2::operator=( const PoseSE2 &other ) {
	if( this != &other ) {
		m_translation	= other.m_translation;
		m_rotation 		= other.m_rotation;
	}
	return *this;
}


PoseSE2&
PoseSE2::operator=( PoseSE2 &&other ) {
	if( this != &other ) {
		m_translation	= std::move( other.m_translation );
		m_rotation 		= std::move( other.m_rotation );
	}
	return *this;
}


PoseSE2
PoseSE2::operator+( const PoseSE2 &p ) const {
	PoseSE2 res = *this;
	return res += p;
}


PoseSE2&
PoseSE2::operator+=( const PoseSE2 &p ) {
	m_translation	+= p.m_translation;
	m_rotation		*= p.m_rotation;

	return *this;
}


PoseSE2
PoseSE2::operator-( const PoseSE2 &p ) const {
	PoseSE2 res = *this;
	return res -= p;
}


PoseSE2&
PoseSE2::operator-=( const PoseSE2 &p ) {
	m_translation	-= p.m_translation;
	m_rotation		*= p.m_rotation.inverse();

	return *this;
}


PoseSE2
PoseSE2::oplus( const PoseSE2 &p ) const {
	/* Compounding operator according to
	 * 	 Smith, Randall, Matthew Self, and Peter Cheeseman. Estimating uncertain
	 * 	 spatial relationships in robotics. Autonomous robot vehicles. Springer
	 * 	 New York, 1990. 167-193.
	 * 	 (Sects 3.2.1)
	 * 	and
	 * 	 F. Lu and E. Milios. Globally consistent range scan alignment for
	 * 	 environment mapping. Journal of Autonomous Robots, 4:333–349, 1997.
	 *   (Eq. 17-20)
	 * This is basically equivalent to the concatenation of the corresponding transformations.
	 */

	return *this * p;

	/*
	double	sinPhi	= sin( phi() ),
			cosPhi	= cos( phi() );

	return PoseSE2(  x() + cosPhi*p.x() - sinPhi*p.y(),
					y() + sinPhi*p.x() + cosPhi*p.y(),
					phi()    + p.phi() );
	*/
}


PoseSE2&
PoseSE2::oplusAssign( const PoseSE2 &p ) {
	/* Compounding operator according to
	 * 	 Smith, Randall, Matthew Self, and Peter Cheeseman. Estimating uncertain
	 * 	 spatial relationships in robotics. Autonomous robot vehicles. Springer
	 * 	 New York, 1990. 167-193.
	 * 	 (Sects 3.2.1)
	 * 	and
	 * 	 F. Lu and E. Milios. Globally consistent range scan alignment for
	 * 	 environment mapping. Journal of Autonomous Robots, 4:333–349, 1997.
	 *   (Eq. 17-20)
	 * This is basically equivalent to the concatenation of the corresponding transformations.
	 */

	return *this *= p;
}


PoseSE2
PoseSE2::ominus( const PoseSE2 &p ) const {
	/* Inverse oplusing operator according to
	 * 	 Smith, Randall, Matthew Self, and Peter Cheeseman. Estimating uncertain
	 * 	 spatial relationships in robotics. Autonomous robot vehicles. Springer
	 * 	 New York, 1990. 167-193.
	 * 	 (Sects 3.2.2)
	 * 	and
	 * 	 F. Lu and E. Milios. Globally consistent range scan alignment for
	 * 	 environment mapping. Journal of Autonomous Robots, 4:333–349, 1997.
	 * 	 (Eq. 21-24)
	 * This is basically equivalent to the concatenation of the corresponding transformations.
	 */

	return p.inverse() *= *this;

	/*
	PoseSE2 	delta	= *this - p;
	double	sinPhi	= sin( p.phi() ),
			cosPhi	= cos( p.phi() );

	return PoseSE2(  cosPhi*delta.x() + sinPhi*delta.y(),
			       -sinPhi*delta.x() + cosPhi*delta.y(),
			       ANGLE_NORMALIZE( delta.phi() ) );
	 */
}


PoseSE2
PoseSE2::oplusBoxplus( const DOFVector &v ) const {
	PoseSE2 res( *this );

	return res.oplusBoxplusassign( v );
}


PoseSE2&
PoseSE2::oplusBoxplusassign( const DOFVector &v ) {
	return oplusAssign( PoseSE2( v[0], v[1], v[2] ) );
}


PoseSE2::DOFVector
PoseSE2::ominusBoxminus( const PoseSE2 &p ) const {
	return ominus( p ).toDOFVector();
}


PoseSE2
PoseSE2::boxplus( const DOFVector &v ) const {
	PoseSE2 res( *this );

	res.boxplusassign( v );

	return res;
}


PoseSE2&
PoseSE2::boxplusassign( const DOFVector &v ) {
	translation()	+= v.head<2>();
	m_rotation.boxplusassign( (double) v[2] );

	return *this;
}


PoseSE2::DOFVector
PoseSE2::boxminus( const PoseSE2 &p ) const {
	DOFVector	res;

	res.head<2>()	= translation() - p.translation();
	res[2]			= m_rotation.boxminus( p.m_rotation );

	return res;
}


PoseSE2
PoseSE2::inverse() const {
    PoseSE2 res;

    res.m_rotation		= m_rotation.inverse();
    res.normalizeRotation();
    res.m_translation	= res.m_rotation * (m_translation * -1.);

    return res;
}


PoseSE2
PoseSE2::operator*( const PoseSE2 &p ) const {
	PoseSE2 res( *this );
	return res *= p;
}


PoseSE2&
PoseSE2::operator*=( const PoseSE2 &p ) {
    m_translation	+=	m_rotation * p.m_translation;
    m_rotation		*= 	p.m_rotation;

    return *this;
}


template<typename Derived>
inline typename std::enable_if<Derived::IsVectorAtCompileTime && Derived::RowsAtCompileTime==2, Point2w>::type
PoseSE2::operator*( const Eigen::MatrixBase<Derived> &p ) const {
	return m_translation + m_rotation * p;
}


template<typename Derived>
inline typename std::enable_if<Derived::IsVectorAtCompileTime && Derived::RowsAtCompileTime==3, Point3w>::type
PoseSE2::operator*( const Eigen::MatrixBase<Derived> &p ) const {
	Derived res( p );

	res.template head<2>() = m_translation + m_rotation * p.template head<2>();

	return res;
}


void
PoseSE2::normalizeRotation() {
	m_rotation.normalize();
}


PoseSE2
PoseSE2::toPoseSE2() const {
	return *this;
}


PoseSE2::operator PoseBase::EigenIsometry2w() const {
	EigenIsometry2w res( m_rotation.toMatrix() );
    res.translation() = m_translation;
    return res;
}


PoseSE2&
PoseSE2::operator=( const EigenIsometry2w &other ) {
	m_translation = other.translation();
	m_rotation.fromMatrix( other.rotation() );

	return *this;
}


PoseSE2::EigenVector
PoseSE2::toVector() const {
	return EigenVector( x(), y(), phi() );
}


PoseSE2&
PoseSE2::operator=( const EigenVector &v ) {
	x() 	= v[0];
	y() 	= v[1];
	phi() 	= v[2];
	return *this;
}


PoseSE2::DOFVector
PoseSE2::toDOFVector() const {
	return toVector();
}


PoseSE2&
PoseSE2::fromDOFVector( const DOFVector &v ) {
	return operator=( v );
}


PoseSE2::TransformationMatrix
PoseSE2::toMatrix() const {
	TransformationMatrix	res( TransformationMatrix::Identity() );

	res.topLeftCorner<2,2>()	= m_rotation.toMatrix();
	res.topRightCorner<2,1>()	= m_translation;

	return res;
}


PoseSE2&
PoseSE2::fromMatrix( const TransformationMatrix &r ) {
	m_rotation.fromMatrix( r.topLeftCorner<2,2>() );
	m_translation = r.topRightCorner<2,1>();

	return *this;
}


PoseSE2::EigenDOFMatrix
PoseSE2::Adjoint() const {
	// According to E. Eade, "Lie Groups for 2D and 3D Transformations"

	EigenDOFMatrix	res( EigenDOFMatrix::Identity() );

	res.block<2, 2>( 0, 0 ) = m_rotation.toMatrix();
	res( 0, 2 )				= y();
	res( 1, 2 )				= -x();

	return res;
}


Point3w
PoseSE2::pos3D() const {
	return Point3w( x(), y(), 0 );
}


} /* namespace bslam */

#endif /* BSLAM_POSESE2_HPP_ */
