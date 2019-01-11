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

#ifndef BSLAM_PCLPOINT_H_
#define BSLAM_PCLPOINT_H_

#include <pcl/point_types.h>

#include "bslam/utils/geometry/Point.h"


namespace bslam {

struct PclPointXYZ;

PCL_EXPORTS std::ostream& operator << (std::ostream& os, const PclPointXYZ& p);


// Wrapper for pcl::PointXYZ, that offers conversions regarding bslam point types
struct EIGEN_ALIGN16 PclPointXYZ : public pcl::_PointXYZ
{
	inline PclPointXYZ( const pcl::_PointXYZ &p ) {
		x = p.x;
		y = p.y;
		z = p.z;
		data[3] = 1.0f;
	}

	inline PclPointXYZ() {
		x = y = z = 0.0f;
		data[3] = 1.0f;
	}

	inline PclPointXYZ( float _x, float _y, float _z ) {
		x = _x;
		y = _y;
		z = _z;
		data[3] = 1.0f;
	}

	inline float squaredNorm() const noexcept {
		return x*x + y*y + z*z;
	}

	inline float norm() const noexcept {
		return sqrt( squaredNorm() );
	}

	inline operator Point2w() const {
		return Point2w( x, y );
	}

	inline operator Point3w() const {
		return Point3w( x, y, z );
	}

	inline operator pcl::PointXYZ() const {
		return pcl::PointXYZ( x, y, z );
	}

	inline float& operator[]( size_t i ) noexcept {
		return data[i];
	}

	inline const float& operator[]( size_t i ) const noexcept {
		return data[i];
	}

	friend std::ostream& operator << (std::ostream& os, const PclPointXYZ& p);
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


/*
template<typename InPointT, typename OutPointT>
inline typename pcl::PointCloud<OutPointT>::Ptr
convertPointCloud( const typename pcl::PointCloud<InPointT>::ConstPtr &inCloud ) {
	typename pcl::PointCloud<OutPointT>::Ptr res( new pcl::PointCloud<OutPointT> );

	res->reserve( inCloud.size() );
	for( const auto &p : *inCloud )
		res->push_back( OutPointT( p ) );

	return res;
}
*/


} /* namespace bslam */

/*
#define PCL_POINT_TYPES \
	PCL_POINT_TYPES \
	(bslam::PclPointXYZ)

#define PCL_XYZ_POINT_TYPES \
		PCL_XYZ_POINT_TYPES \
		(bslam::PclPointXYZ)
*/

// TODO: Check, if this is correct
// Alternative method: http://pointclouds.org/documentation/tutorials/adding_custom_ptype.php#how-to-add-a-new-pointt-type
POINT_CLOUD_REGISTER_POINT_WRAPPER(bslam::PclPointXYZ, pcl::_PointXYZ)

#endif /* BSLAM_PCLPOINT_H_ */
