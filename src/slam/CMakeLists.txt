SET( BSLAM_SLAM2D_SRSC
	BayesFastSLAM2D.cpp
	BetaFastSLAM2D.cpp
)


ADD_LIBRARY( bslam-slam2d SHARED ${BSLAM_SLAM2D_SRSC} )
#ADD_LIBRARY( bslam-slam2d-static STATIC ${BSLAM_SLAM2D_SRSC} )

target_link_libraries( bslam-slam2d bslam-utils )
#target_link_libraries( bslam-slam2d-static bslam-utils-static )

