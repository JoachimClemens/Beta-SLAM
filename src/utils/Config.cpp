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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#ifndef WIN32
// Time functions
#	include <time.h>

// Dir access
#	include <sys/types.h>
#	include <sys/stat.h>
#	include <unistd.h>
#	include <errno.h>
#	include <string.h>
#endif

#include "bslam/utils/Config.h"
#include "bslam/utils/Convenience.h"

#if 0
#include "ConfigSqlAccessor.h"
#endif

namespace bslam {

#define KEY_WIDTH 32


std::unordered_map<std::string, std::string>	Config::sm_stringVals;
std::unordered_map<std::string, double>			Config::sm_doubleVals;
std::unordered_map<std::string, int>			Config::sm_intVals;

enum ParserState {
	SEARCHING_KEY,
	IN_KEY,
	SEARCHING_SEPARATOR,
	IN_VALUE,
	SEARCHING_VALUE,
	IN_COMMENT,
	IN_ESCAPE,
	SEARCHING_END,
	IN_QUOTE,
	ERROR
};


void
Config::loadDefaults( bool clearBeforeLoad ) {
	if( clearBeforeLoad )
		clear();

	// general
	set( "MAP_TYPE", 			"beta" 			);	// "beta" for Beta-SLAM or "bayes" for classical SLAM
	set( "FILENAME", 			"SET_ME" 		);
	set( "GT_FILENAME",			"" 				);
	set( "OUTPUT_BASEDIR",		"../results" 	);
	set( "OUTPUT_DIR",			""				);	// if this value remains empty, it will be overwritten by createOutputDir
	set( "NO_GUI",				0				); 	// only for BSlamSimpleGui
	set( "SLEEP_AFTER_INIT", 	2000			);  // [ms] sleeps this much after initialization, usefull to see warnings on the console
	set( "SAVE_EACH_STEP",		0				);  // saves the best map and path each X steps, 0 disables this feature (only implemented for BSlamCarmenGui, yet)

	// Laser scanner
	set( "LASER_START_ANGLE", 			-90.0 		); // [deg] Angle of first laser beam
	set( "LASER_ANGULAR_RES",   		1.0 		); // [deg] Angle between to laser beams


	//------------------------
	// Beta sensor model
	//------------------------
	set( "BETA_SIGMA_FORWARD",			0.0375	);	// [m]		standard deviation for forward model
	set( "BETA_SIGMA_INVERSE",			0.0375	);	// [m]		standard deviation for inverse model
	set( "BETA_SIGMA_FORWARD_FACTOR",	-1.0	);	// [0-inf]	sigma range (factor * BETA_SIGMA_FORWARD) to include in integration in forward model (controls the specificity of the model, set <0 for 0.5*MAP_DELTA)
	set( "BETA_SIGMA_INVERSE_FACTOR",	2.0		);	// [0-inf]	sigma range (factor * BETA_SIGMA_INVERSE) to include in integration in inverse model (controls the probability mass on occupied, set <0 for 0.5*MAP_DELTA)
	set( "BETA_FORWARD_SIGMA_RANGE",	3.0		);	// [0-inf]	include all cells in the range BETA_FORWARD_SIGMA_RANGE*BETA_SIGMA_FORWARD when computing the forward model
	set( "BETA_RANDOM",					0.3		);	// [0-1]	probability of an random measurement
	set( "BETA_PRIOR_P",				0.5		);	// [0-1]	prior for m_i, i.e., P(m_i|z,q)
	set( "BETA_PRIOR_ALPHA",			0.5		);	// [0-inf]	prior for alpha
	set( "BETA_PRIOR_BETA",				0.5		);	// [0-inf]	prior for beta
	set( "BETA_EXCLUDE_NO_HIT",			0		);	// do not use forward model if scan matching failed for this measurement


	//------------------------
	// Bayes sensor model
	//------------------------
	set( "BAYES_SIGMA",  			0.25	);		// [m]		standard deviation for forward model


	//--------------
	// Motion model
	//--------------
	set( "TRANS_TRANS_SIGMA",	0.1 	);	// [m^2]
	set( "TRANS_ROT_SIGMA",		0.1 	);	// [m*rad]
	set( "ROT_ROT_SIGMA",		0.1 	);	// [rad^2]
	set( "ROT_TRANS_SIGMA",		0.1 	);	// [rad*m]


	//--------------
	// SLAM
	//--------------
	set( "MAP_MIN_X", 				-100.0 	);	// [m]
	set( "MAP_MAX_X", 				100.0 	);	// [m]
	set( "MAP_MIN_Y", 				-100.0 	);	// [m]
	set( "MAP_MAX_Y", 				100.0	);	// [m]
	set( "MAP_MIN_Z", 				-5.0	);	// [m]
	set( "MAP_MAX_Z",				5.0		);	// [m]
	set( "MAP_DELTA",				0.05	);	// [m] 		size of a cell


	//--------------
	// FastSLAM
	//--------------
	set( "LINEAR_DIST_THRESHOLD", 				0.5 					);	// [m]		process the next scan after the robot traveled this distance
	set( "ANGULAR_DIST_THRESHOLD", 				10.0 					);	// [deg]	process the next scan after the robot turned by this angle
	set( "NUM_PARTICLES", 						30 						);
	set( "MIN_SCORE",							0.0						);	// [0-inf]
	set( "RESAMPLE_THRESHOLD",					0.5						);	// [0-1]	percentage of effective particles
	set( "LIKELIHOOD_GAIN",						20.0					);	// (0-inf]	used to compute weights from log weights, higher values reduces frequency of resampling
	set( "IMPROVED_PROPOSAL",					"after_scan_matching"	);	// "no" -> improved proposal not used, "after_scan_matching" -> estimate around corrected pose from scan matching, or "during_scan_matching" -> estimate during scan matching
	set( "IMPROVED_PROPOSAL_LINEAR_RANGE",		0.025 					);	// [m]		estimate improved proposal in +- this linear range (only for IMPROVED_PROPOSAL="after_scan_matcher")
	set( "IMPROVED_PROPOSAL_ANGULAR_RANGE",		1.5 					);	// [deg]	estimate improved proposal in +- this angular range (only for IMPROVED_PROPOSAL="after_scan_matcher")
	set( "IMPROVED_PROPOSAL_NUM_SAMPLES",		3						);	// [1-inf]	number of samples per DoF for estimating improved proposal, values <3 are usually a bad idea (only for IMPROVED_PROPOSAL="after_scan_matcher")
	set( "IMPROVED_PROPOSAL_USE_FORWARD_MODEL",	1						);	// [0,1]	use real forward model instead of score likelihood to compute improved proposal (only reasonable for IMPROVED_PROPOSAL="after_scan_matcher")


	//------------------
	// Map Scan matcher
	//------------------
	set( "FREE_CELL_RATIO", 	sqrt( 2.0 )	);	// Distance before the hit cell where an empty cell is expected (will be multiplied with MAP_DELTA)
	set( "USABLE_RANGE", 		15.0 		);	// [m] 	 Greater measurements will be reduced to this value
	set( "MAX_RANGE",			80.0 		);	// [m] 	 Measurements greater than this value will be skipped
	set( "FULLNESS_THRESHOLD", 	0.1 		);	// [0-1] A cell is considered to be occupied if the fullness is above this value and considered to be empty otherwise
	set( "RELEVANT_DIST",		0.0			);	// [m]   Only cells in the given distance from the measurement are considered in forward and inverse model, set to 0 to disable the feature
	set( "SCAN_NO_CONFLICT",	0			);	// [0,1] If there are conflicting values for a cell in one scan (i.e., due to grid resolution), prefer the occupied value instead of using a combination (parameter is actually used in sensor models)

	// Score calculation
	set( "SCORE_SIGMA", 				0.25 		);
	set( "KERNEL_SIZE", 				1 			);	// Half size of the kernel for score calculation, kernel will be 2*KERNEL_SIZE + 1
	set( "KERNEL_CROSS_2D", 			0 			);	// Use a cross instead of a square as kernel for score calculation in 2D
	set( "KERNEL_CROSS_3D", 			0 			);	// Use a cross instead of a cube as kernel for score calculation in 3D
	set( "SCORE_RANGE_Y",				0.05 		);	// [m] Only measurements in this y range are considered for score calculation (only used for 3D maps and PoseSE3, not for likelihood if LIKELIHOOD_LIKE_SCORE=0)
	set( "SCORE_RANGE_Z",				0.5 		);	// [m] Only measurements in this z range are considered for score calculation (only used for 3D maps, not for likelihood if LIKELIHOOD_LIKE_SCORE=0)
	set( "SCORE_SKIP_POINTS",			10			);	// Skip so many points before considering one for score and likelihood calculation (only used for 3D maps)
	set( "LIKELIHOOD_NUM_MEASUREMENTS",	500			);	// Fixed number of measurements to use to compute the likelihood. Note that additional measurements can be skipped when no matching cell is found. If set to 0 or if LIKELIHOOD_LIKE_SCORE=1, SCORE_SKIP_POINTS and SCORE_RANGE_* is used.
	set( "LIKELIHOOD_LIKE_SCORE",		1			);	// Select the points for likelihood calculation like for score calculation, i.e. use SCORE_RANGE_* and SCORE_SKIP_POINTS while LIKELIHOOD_NUM_MEASUREMENTS is ignored

	// Pose optimization
	set( "OPT_ANGULAR_STEP", 	3.0 		);	// [deg]
	set( "OPT_LINEAR_STEP", 	.05 		);	// [m]
	set( "OPT_ITERATIONS",		5 			);
}


void
Config::loadParams( int argc, const char *argv[], bool clearBeforeLoad, bool warnUnknown ) {
	if( argc < 3 )
		return;

	if( clearBeforeLoad )
		clear();

	int cur = 1;
	while( (argc - cur) >= 2 ) {
		// params start with "--"
		if( argv[cur][0] != '-' || argv[cur][1] != '-' ) {
			cur++;
			continue;
		}

		bool error = false;

		std::string	key 		= &argv[cur][2],
					value		= argv[cur+1];
		bool		wasQuoted	= false;

		// convert to upper case and '-' to '_'
		for( size_t i = 0; i < key.size(); i++ ) {
			key[i] = std::toupper( key[i] );
			if( key[i] == '-' )
				key[i] = '_';

			if( !isAllowedKeyChar( key[i] ) ) {
				l_wrn( "Unallowed character in command line argument #" << cur << " `" << argv[cur] << "': " << key[i] );
				error = true;
				break;
			}
		}

		if( error ) {
			cur++;
			continue;
		}

		// remove whitespaces at front and back
		trim_inplace( value );

		if( value.size() >= 2 ) {
			// detect double or single quote and erase them
			if( ( value.front() == '"' && value.back() == '"' )
				|| ( value.front() == '\'' && value.back() == '\'' ) ) {
				value.erase( 0, 1 );				// remove first character
				value.erase( value.size()-1, 1 );	// remove last character
				wasQuoted = true;
			}
		}

		addValue( key, value, wasQuoted, warnUnknown );

		cur += 2;
	}
}


void
Config::loadStream( std::istream &stream, bool clearBeforeLoad, bool warnUnknown ) {
	std::string line,
				key,
				value,
				error;
	int lineNum 		= 0;
	bool wasQuoted 		= false;
	ParserState state 	= SEARCHING_KEY;

	if( clearBeforeLoad )
		clear();

	// parse lines
	while( std::getline( stream, line ) ) {
		lineNum++;
		for( std::string::iterator c = line.begin(); c != line.end() && state != ERROR && state != IN_COMMENT; c++ ) {
			switch( state ) {
				case SEARCHING_KEY:
					if( *c == '#' ) {
						state = IN_COMMENT;
					} else if( isAllowedKeyChar( *c ) ) {
						key += *c;
						state = IN_KEY;
					}	else if( !isspace( *c ) ) {
						error = std::string( "invalid char for a key: `" ).append( 1, *c ).append( "\'" );
						state = ERROR;
					}
					break;

				case IN_KEY:
					if( *c == ':' ) {
						state = SEARCHING_VALUE;
					} else if( isspace( *c ) ) {
						state = SEARCHING_SEPARATOR;
					} else if( isAllowedKeyChar( *c ) ) {
						key += *c;
					} else {
						error = std::string( "invalid char for a key: `" ).append( 1, *c ).append( "'" );
						state = ERROR;
					}
					break;

				case SEARCHING_SEPARATOR:
					if( *c == ':' ) {
						state = SEARCHING_VALUE;
					} else if( !isspace( *c ) ) {
						error = std::string( "separator `:' expected, not `" ).append( 1, *c ).append( "'" );
						state = ERROR;
					}
					break;

				case SEARCHING_VALUE:
					if( *c == '#' ) {
						error = std::string( "value expected" );
						state = ERROR;
					} else if( *c == '\"' ) {
						state = IN_QUOTE;
					} else if( !isspace( *c ) ) {
						value += *c;
						state = IN_VALUE;
					}
					break;

				case IN_VALUE:
					if( *c == '#' ) {
						state = IN_COMMENT;
					} else if( isspace( *c ) ) {
						state = SEARCHING_END;
					} else {
						value += *c;
					}
					break;

				case IN_QUOTE:
					wasQuoted = true;
					if( *c == '"' ) {
						state = SEARCHING_END;
					} else if( *c == '\\' ) {
						state = IN_ESCAPE;
					} else {
						value += *c;
					}
					break;

				case IN_ESCAPE:
					if( *c == '"' ) {
						value += *c;
						state = IN_QUOTE;
					/* // this don't have to be escaped anymore
					} else if( *c == '#' ) {
						value += *c;
					*/
					} else if( *c == '\\' ) {
						value += '\\';
						state = IN_ESCAPE;
					} else {
						value += '\\';
						value += *c;
						state = IN_QUOTE;
					}
					break;

				case SEARCHING_END:
					if( *c == '#' ) {
						state = IN_COMMENT;
					} else if( !isspace( *c ) ) {
						error = std::string( "unexpected char: `" ).append( 1, *c ).append( "\' (missing quotes?)" );
						state = ERROR;
					}
					break;

				case IN_COMMENT:
				case ERROR:
					break;
				default:
					l_err( "Invalid parser state." );
					break;
			}
		}

		if( state == ERROR ) {
			l_wrn( "Syntax error in line " << lineNum << ": " << error.c_str() );
		} else if( state == IN_QUOTE ) {
			l_wrn( "Syntax error in line " << lineNum << ": missing `\"'" );
		} else if( key.empty() ) {
			// ignore
		} else if( !wasQuoted && value.empty() ) {
			l_wrn( "Syntax error in line " << lineNum << ": key `" << key << "' with empty value" );
		} else {
			addValue( key, value, wasQuoted );
		}

		error.clear();
		key.clear();
		value.clear(); 
		state = SEARCHING_KEY;
		wasQuoted = false;
	}

	//std::cout << toString() << std::endl;
}


void
Config::saveStream( std::ostream &stream ) {
	// write std::strings
	for( auto iter = sm_stringVals.begin(); iter != sm_stringVals.end(); iter++ ) {
		// escape double-quotes
		std::string outStr = iter->second;
		for( size_t pos = 0; pos < outStr.length(); pos++ ) {
			if( outStr[pos] == '"' ) {
				outStr.replace( pos, 1, "\\\"" );
				pos++;
			}
			/* //this don't have to escaped anymore
			else if( outStr[pos] == '#' ) {
				outStr.replace( pos, 1, "\\#" );
				pos++;
			}
			*/
		}

		stream << std::setw( KEY_WIDTH ) << std::left << iter->first + ":" << "\"" << outStr << "\"" << std::endl;
	}

	// write doubles
	for( auto iter = sm_doubleVals.begin(); iter != sm_doubleVals.end(); iter++ )
		stream << std::setw( KEY_WIDTH ) << std::left << iter->first + ":" << std::fixed << std::setprecision(20) << iter->second << std::endl;

	// write ints
	for( auto iter = sm_intVals.begin(); iter != sm_intVals.end(); iter++ )
		stream << std::setw( KEY_WIDTH ) << std::left << iter->first + ":" << iter->second << std::endl;
}


void
Config::loadFile( const char *filename, bool clearBeforeLoad, bool warnUnknown ) {
	std::ifstream filestream( filename );

	if( !filestream.is_open() )
		throw std::ios_base::failure( "Failed to open config file for reading." );

	loadStream( filestream, clearBeforeLoad, warnUnknown );

	filestream.close();
}


void
Config::loadFile( const std::string &filename, bool clearBeforeLoad, bool warnUnknown ) {
	loadFile( filename.c_str(), clearBeforeLoad, warnUnknown );
}


void
Config::loadDefaultFile( bool clearBeforeLoad, bool warnUnknown ) {
	loadFile( std::string( std::getenv( "HOME" ) ) + std::string( "/.bslamrc" ), clearBeforeLoad, warnUnknown );
}

void
Config::saveFile( char const *filename ) {
	// Note: This will overwrite comments, ordering, and formatting in existing files

	std::ofstream filestream( filename );

	if( !filestream.is_open() )
		throw std::ios_base::failure( "Failed to open config file for writing." );
	
	saveStream( filestream );

	filestream.close();
}


void
Config::saveFile( const std::string &filename ) {
	saveFile( filename.c_str() );
}


#ifndef WIN32
void
Config::createOutputDir( const std::string &suffix ) {
	std::string outputDir = get( "OUTPUT_DIR", "" );

	if( outputDir.empty() ) {
		time_t	curtime;
		char	prefix[200];

		time( &curtime );
		strftime( prefix, 199, "%F_%H-%M-%S_", localtime( &curtime ) );

		std::string	basedir	= get( "OUTPUT_BASEDIR", "." );

		outputDir = basedir + "/" + prefix + (!suffix.empty() ? suffix : std::string( "bslam" ) );
	}

	struct stat st;
	if( stat( outputDir.c_str(), &st ) != 0 ) // exists already?
		if( mkdir( outputDir.c_str(), 0755 ) != 0 && errno != EEXIST ) // EEXIST for race condition
			throw std::ios::failure( "Failed to create output directory `" + outputDir + "': " + std::string( strerror( errno ) ) );

	set( "OUTPUT_DIR", outputDir );
}
#endif

#if 0
void
Config::loadFromDb( SubsystemID source, bool clearBeforeLoad, bool overwriteExistingValues ) {
	if( clearBeforeLoad )
		clear();

	size_t  endIdx;

	// Load elements
	std::map<string, std::string>  keyValues = ConfigSqlAccessor::loadElements( source );

	// Sort the values into the correct lists (int, double, std::string)
	for( auto kv : keyValues ) {
	  // try to convert to int
	  int iVal = std::stoi( kv.second, &endIdx );

	  // is really an int value (everything was converted)?
	  if( endIdx == kv.second.size() ) {
	    if( overwriteExistingValues || !hasIntValue( kv.first ) )
	      sm_intVals[kv.first] = iVal;
	  } else {
      // try to convert to double
      double dVal = std::stod( kv.second, &endIdx );

      // is really a double value (everything was converted)?
      if( endIdx == kv.second.size() ) {
        if( overwriteExistingValues || !hasDoubleValue( kv.first ) )
          sm_doubleVals[kv.first] = dVal;
      } else {
        // if it is no int and no double, it has to be a std::string
        if( overwriteExistingValues || !hasStringValue( kv.first ) )
          sm_stringVals[kv.first] = kv.second;
      }
	  }
	}
}


void
Config::saveToDb() {
	// TODO
}
#endif

void
Config::clear() {
	sm_stringVals.clear();
	sm_intVals.clear();
	sm_doubleVals.clear();
}


std::string
Config::toString() {
	std::ostringstream strs;

	strs << "Strings:" << std::endl;
	for( auto iter = sm_stringVals.begin(); iter != sm_stringVals.end(); iter++ )
		strs << std::setw( KEY_WIDTH ) << std::left << iter->first + " ->" << "`" << iter->second << "'" << std::endl;

	strs << std::endl << "Doubles:" << std::endl;
	for( auto iter = sm_doubleVals.begin(); iter != sm_doubleVals.end(); iter++ )
		strs << std::setw( KEY_WIDTH ) << std::left << iter->first + " ->" << std::fixed << std::setprecision(20) << iter->second << std::endl;

	strs << std::endl << "Integers:" << std::endl;
	for( auto iter = sm_intVals.begin(); iter != sm_intVals.end(); iter++ )
		strs << std::setw( KEY_WIDTH ) << std::left << iter->first + " ->" << iter->second << std::endl;

	return strs.str();
}


bool
Config::isAllowedKeyChar( char c ) {
	if( isalnum( c ) || c == '_' || c == '-' )
		return true;
	else
		return false;
}


Config::ValueType
Config::getValueType( const std::string &str ) {
	if( str.empty() ) {
		return VALUE_TYPE_STRING;
	}

	ValueType type = VALUE_TYPE_INT;

	for( size_t i = 0; i < str.length(); i++ ) {
		//std::cout << i << " " << str[i] << std::endl;
		if( type == VALUE_TYPE_INT && str[i] == '.' ) { // found dot in number
			type = VALUE_TYPE_DOUBLE;
		} else if( (type == VALUE_TYPE_DOUBLE && str[i] == '.')					// found a second dot
					|| (!isdigit( str[i] ) && str[i] != '.' && str[i] != '-' ) 	// not a digit and not a point
					|| (str[i] == '-' && i > 0 )								// minus not as the first char
						) {
			type = VALUE_TYPE_STRING;
			break;
		}
	}

	//std::cout << "Result: " << type << std::endl;

	return type;
}


void
Config::addValue( const std::string &key, const std::string &value, bool wasQuoted, bool warnUnknown ) {
	ValueType type = getValueType( value );
	if( wasQuoted || type == VALUE_TYPE_STRING ) {
		if( warnUnknown && sm_stringVals.find( key ) == sm_stringVals.end() )
			l_wrn( "Key `" << key << "' was previously unknown for value type string!" );
		sm_stringVals[key]	= value;
	} else if (type == VALUE_TYPE_INT ) {
		if( warnUnknown && sm_intVals.find( key ) == sm_intVals.end() )
			l_wrn( "Key `" << key << "' was previously unknown for value type int!" );
		sm_intVals[key]		= std::stoi( value );
	} else {
		if( warnUnknown && sm_doubleVals.find( key ) == sm_doubleVals.end() )
			l_wrn( "Key `" << key << "' was previously unknown for value type double!" );
		sm_doubleVals[key]	= std::stof( value );
	}
}


std::string
Config::toLower( const std::string &in ) {
	std::string out( in );

	std::transform( out.begin(), out.end(), out.begin(), ::tolower );

	return out;
}


std::string
Config::toUpper( const std::string &in ) {
	std::string out( in );

	std::transform( out.begin(), out.end(), out.begin(), ::toupper );

	return out;
}

} /* namespace bslam */

