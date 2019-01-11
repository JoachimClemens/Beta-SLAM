set(BSLAM_IO_SRSC
	../empty.cpp
)

ADD_SUBDIRECTORY( carmen )

ADD_LIBRARY( bslam-io SHARED ${BSLAM_IO_SRSC} )
#ADD_LIBRARY( bslam-io-static STATIC ${BSLAM_IO_SRSC} )

target_link_libraries(bslam-io bslam-carmen)
#target_link_libraries(bslam-io-static bslam-carmen-static)
