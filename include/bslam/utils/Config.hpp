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

#include "bslam/utils/Log.h"

#define BSLAM_CONFIG_CHECK_KEY( __key ) if( !isValidKey( __key ) ) throw std::invalid_argument( "invalid characters in key" );


namespace bslam {


std::string
Config::getString( const std::string &key ) {
	return sm_stringVals.at( key );
}


std::string
Config::getString( const std::string &key, const char *defaultVal ) {
	std::string val;

	try {
		val = getString( key );
	} catch( std::out_of_range &oor ) {
		val = std::string( defaultVal );
		l_wrn( "Unable to retrieve std::string value for `" << key.c_str() << "' from configuration, using default `" << defaultVal << "'." );
	}

	return val;
}


std::string
Config::getString( const std::string &key, const std::string defaultVal ) {
	return getString( key, defaultVal.c_str() );
}


double
Config::getDouble( const std::string &key ) {
	return sm_doubleVals.at( key );
}


double
Config::getDouble( const std::string &key, const double defaultVal ) {
	double val;

	try {
		val = getDouble( key );
	} catch( std::out_of_range &oor ) {
		val = defaultVal;
		l_wrn( "Unable to retrieve double value for `" << key.c_str() << "' from configuration, using default `" << defaultVal << "'." );
	}

	return val;
}


int
Config::getInt( const std::string &key ) {
	return sm_intVals.at( key );
}


int
Config::getInt( const std::string &key, const int defaultVal ) {
	int val;

	try {
		val = getInt( key );
	} catch( std::out_of_range &oor ) {
		val = defaultVal;
		l_wrn( "Unable to retrieve int value for `" << key.c_str() << "' from configuration, using default `" << defaultVal << "'." );
	}

	return val;
}


std::string
Config::get( const std::string &key, const std::string defaultVal ) {
	return getString( key, defaultVal );
}


double
Config::get( const std::string &key, const double defaultVal ) {
	return getDouble( key, defaultVal );
}


int
Config::get( const std::string &key, const int defaultVal ) {
	return getInt( key, defaultVal );
}


void
Config::set( const std::string &key, const std::string val ) {
	BSLAM_CONFIG_CHECK_KEY( key );
	sm_stringVals[key] = val;
}


void
Config::set( const std::string &key, const double val ) {
	BSLAM_CONFIG_CHECK_KEY( key );
	sm_doubleVals[key] = val;
}


void
Config::set( const std::string &key, int val ) {
	BSLAM_CONFIG_CHECK_KEY( key );
	sm_intVals[key] = val;
}


bool
Config::hasStringValue( const std::string &key ) {
  return sm_stringVals.find( key ) != sm_stringVals.end();
}


bool
Config::hasDoubleValue( const std::string &key ) {
  return sm_doubleVals.find( key ) != sm_doubleVals.end();
}


bool
Config::hasIntValue( const std::string &key ) {
	return sm_intVals.find( key ) != sm_intVals.end();
}


bool
Config::isValidKey( const std::string &key ) {
	for( size_t i = 0; i < key.size(); i++ )
		if( !isAllowedKeyChar( key[i] ) )
			return false;

	return true;
}


} /* namespace bslam */

