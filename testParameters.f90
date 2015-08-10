MODULE testParameters
  ! ========================================================================
  ! This module will contain information about the current physical domain and test parameters.
  ! Used to simplify passing of this information throughout subroutines and functions
  ! ========================================================================
  IMPLICIT NONE
  INTEGER :: limitingMeth,magCount
  DOUBLE PRECISION :: mu1,mu2,zsMinWeight
  SAVE

END MODULE testParameters
