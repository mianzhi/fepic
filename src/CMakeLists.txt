SET(fepic_src
  fepic.f90
  modFEPIC.f90
  modParticle.f90
  modPICGrid.f90
  modUglyFEM.f90
  modMGS.f90
  modPush.f90
  modImpact.f90
  modSource.f90
)

ADD_EXECUTABLE(fepic ${fepic_src})

TARGET_LINK_LIBRARIES(fepic
  fosolvers
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${SUNDIALS_LIBRARIES}
  ${UMFPACK_LIBRARIES}
  ${MATHEVAL_LIBRARIES}
  ${METIS_LIBRARIES}
)

INSTALL(TARGETS fepic DESTINATION bin)
