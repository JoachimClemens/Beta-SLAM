SET( BSLAM_CARMEN_SRSC
    CarmenReader.cpp
    CarmenGtReader.cpp
)

ADD_LIBRARY( bslam-carmen SHARED ${BSLAM_CARMEN_SRSC} )
#ADD_LIBRARY( bslam-carmen-static STATIC ${BSLAM_CARMEN_SRSC} )

TARGET_LINK_LIBRARIES( bslam-carmen bslam-slam2d )
#TARGET_LINK_LIBRARIES( bslam-carmen-static bslam-slam2d-static )

ADD_EXECUTABLE( BSlamCarmenGui BSlamCarmenGui.cpp )
TARGET_LINK_LIBRARIES( BSlamCarmenGui bslam-carmen bslam-viz bslam-utils )
