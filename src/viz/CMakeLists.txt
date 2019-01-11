SET( BSLAM_VIZ_SRSC
    MapView2D.cpp
    NeffView.cpp
)

SET( BSLAM_VIZ_INCLUDE_DIR ../../include/bslam/viz )

SET( BSLAM_VIZ_MOC_HDRS
    ${BSLAM_VIZ_INCLUDE_DIR}/MapView2D.h
    ${BSLAM_VIZ_INCLUDE_DIR}/NeffView.h
)

QT5_WRAP_CPP( BSLAM_VIZ_MOC_SRCS ${BSLAM_VIZ_MOC_HDRS} )

ADD_LIBRARY( bslam-viz SHARED ${BSLAM_VIZ_SRSC} ${BSLAM_VIZ_MOC_SRCS} )
#ADD_LIBRARY( bslam-viz-static STATIC ${BSLAM_VIZ_SRSC} ${BSLAM_VIZ_MOC_SRCS} )

target_link_libraries( bslam-viz Qt5::Widgets Qt5::Gui )
#target_link_libraries( bslam-viz-static ${QT_LIBRARIES} )

