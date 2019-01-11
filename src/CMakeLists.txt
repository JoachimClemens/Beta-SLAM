
ADD_SUBDIRECTORY( maps )
ADD_SUBDIRECTORY( utils )
ADD_SUBDIRECTORY( viz )
ADD_SUBDIRECTORY( io )
ADD_SUBDIRECTORY( slam )

ADD_LIBRARY( bslam SHARED empty.cpp )
SET_TARGET_PROPERTIES( bslam PROPERTIES
    VERSION ${BSLAM_VERSION}
    SOVERSION ${BSLAM_VERSION}
)
#ADD_LIBRARY( bslam-static STATIC empty.cpp )
#SET_TARGET_PROPERTIES( bslam-static PROPERTIES OUTPUT_NAME "bslam" )

TARGET_LINK_LIBRARIES( bslam bslam-slam2d bslam-utils bslam-io )
#TARGET_LINK_LIBRARIES( bslam-static bslam-slam-static bslam-utils-static bslam-io-static )
