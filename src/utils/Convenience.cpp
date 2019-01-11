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

#include <vector>
#include <sstream>

#include "bslam/utils/Convenience.h"

namespace bslam {


std::vector<std::string>
split( const std::string &s, char delim ) {
  // http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
  std::vector<std::string>  elems;
  std::istringstream        ss( s );
  std::string               item;

  while( std::getline( ss, item, delim ) )
	  if( !item.empty() )
		  elems.push_back( item );

  return elems;
}


std::string &
trim_inplace( std::string &s, std::function<bool( const char & )> f ) {
	// according to http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring and http://www.codeproject.com/Articles/10880/A-trim-implementation-for-std-string

	// trim left
	s.erase( s.begin(), std::find_if( s.begin(), s.end(), f ) );

	// trim right
	s.erase( std::find_if( s.rbegin(), s.rend(), f ).base(), s.end() );

	return s;
}


std::string
trim( const std::string &s, std::function<bool( const char & )> f ) {
	std::string res( s );
	trim_inplace( res, f );
	return res;
}


} /* namespace bslam */
