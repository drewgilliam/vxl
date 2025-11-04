# brl/bseg/CMakeLists.txt
# main package for segmentation-related code and applications
# Mutual dependencies: brip --> sdet --> (vpro,segv)

# float image operations
add_subdirectory(brip)

# segmentation feature detection library
add_subdirectory(sdet)

# segmentation view library
add_subdirectory(segv)

# segmentation application executables
add_subdirectory(sbin)

# background modeling
add_subdirectory(bbgm)

# Lidar modeling (segmentation and meshing)
add_subdirectory(bmdl)

# octree
add_subdirectory(boct)

# semi-global matching
add_subdirectory(bsgm)

# Brown affine patch library, SIFT keypoint & descriptor
add_subdirectory(bapl)

#event trigger library
add_subdirectory(betr)

# VSI conflation tools
add_subdirectory(vcon)

# appearance modeling library
add_subdirectory(baml)
