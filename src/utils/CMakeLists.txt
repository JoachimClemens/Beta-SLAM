SET( BSLAM_UTILS_SRSC
    Config.cpp
    Convenience.cpp
    Version.cpp
)

ADD_CUSTOM_COMMAND( OUTPUT VersionInfo.h
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../../scripts/create_version_info.sh ${CMAKE_CURRENT_SOURCE_DIR} Version.cpp
)
ADD_CUSTOM_TARGET( VersionInfo DEPENDS VersionInfo.h )

ADD_LIBRARY( bslam-utils SHARED ${BSLAM_UTILS_SRSC} )
#ADD_LIBRARY( bslam-utils-static STATIC ${BSLAM_UTILS_SRSC} )

ADD_DEPENDENCIES( bslam-utils VersionInfo )
#ADD_DEPENDENCIES( bslam-utils-static VersionInfo )

target_link_libraries( bslam-utils  ${BSLAM_ROS_LIBS} )
#target_link_libraries( bslam-utils-static  ${BSLAM_ROS_LIBS} )