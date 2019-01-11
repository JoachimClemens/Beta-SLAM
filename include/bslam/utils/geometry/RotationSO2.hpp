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

namespace bslam {

RotationSO2::RotationSO2() :
	m_angle( 0 ),
	m_matrix( EigenMatrix2::Identity() ),
	m_matrixAngle( 0 )
{
	// Nothing else to do here
}
					

RotationSO2::RotationSO2( const RotationSO2 &other ) :
	m_angle( other.m_angle ),
	m_matrix( EigenMatrix2::Identity() ),
	m_matrixAngle( 0 )
{
	// Nothing else to do here
}
					

RotationSO2::RotationSO2( Scalar angle ) :
	m_angle( ANGLE_NORMALIZE( angle ) ),
	m_matrix( EigenMatrix2::Identity() ),
	m_matrixAngle( 0 )
{
	// Nothing else to do here
}


RotationSO2::RotationSO2( const EigenRotation &other ) :
	m_angle( ANGLE_NORMALIZE( other.angle() ) ),
	m_matrix( EigenMatrix2::Identity() ),
	m_matrixAngle( 0 )
{
	// Nothing else to do here
}


RotationSO2&
RotationSO2::operator=( const RotationSO2 &other ) {
	if( this != &other ) {
		m_angle = other.m_angle;
	}
	return *this;
}


RotationSO2&
RotationSO2::operator=( const EigenRotation &other ) {
	m_angle = ANGLE_NORMALIZE( other.angle() );
	return *this;
}


RotationSO2&
RotationSO2::boxplusassign( Scalar delta ) {
	m_angle += delta;
	normalize();
	return *this;
}


RotationSO2
RotationSO2::boxplus( Scalar delta ) const {
	RotationSO2 res( *this );
	res.boxplusassign( delta );
	return res;
}


RotationSO2::Scalar
RotationSO2::boxminus( const RotationSO2 &other ) const {
	return ANGLE_NORMALIZE( m_angle - other.m_angle );
}


RotationSO2&
RotationSO2::operator*=( const RotationSO2 &other ) {
	m_angle += other.m_angle;
	normalize();
	return *this;
}


RotationSO2
RotationSO2::operator*( const RotationSO2 &other ) const {
	return RotationSO2( m_angle + other.m_angle );
}


RotationSO2::EigenVector2
RotationSO2::operator*( const EigenVector2 &v ) const {
	return getMatrix() * v;
}


RotationSO2::EigenMatrix2
RotationSO2::toMatrix() const {
	return getMatrix();
}


void
RotationSO2::fromMatrix( const EigenMatrix2 &m ) {
	m_angle = acos( (double) m( 0, 0 ) );
}


void
RotationSO2::normalize() {
	ANGLE_NORMALIZE_INPLACE( m_angle );
}


RotationSO2
RotationSO2::inverse() const {
	RotationSO2 res( *this );
	res.invert();
	return res;
}


void
RotationSO2::invert() {
	m_angle *= -1;
}

const RotationSO2::EigenMatrix2&
RotationSO2::getMatrix() const {
	if( m_angle != m_matrixAngle ) {
		Scalar	cosAngle = cos( m_angle ),
				sinAngle = sin( m_angle );

		m_matrix <<	cosAngle, -sinAngle,
					sinAngle,  cosAngle;
	}

	return m_matrix;
}

} /* namespace bslam */
