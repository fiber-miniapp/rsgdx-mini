!C
!C    * add sliponly & deficit etc.  on 20120501
!C                         LAST CHANGE: 2012/05/01 18:00:57.
      module m_outputval

      USE precisions
      USE constant

      IMPLICIT NONE

      real   (dp) :: vsav  (0:NBMAX-1,0:NSAVMAX-1)
      real   (dp) :: thsav (0:NBMAX-1,0:NSAVMAX-1)
      real   (dp) :: ssav  (0:NBMAX-1,0:NSAVMAX-1)

      real   (dp), pointer :: vsav_gather  (:,:)=> null()
      real   (dp), pointer ::thsav_gather  (:,:)=> null()
      real   (dp), pointer :: ssav_gather  (:,:)=> null()
!>mh on 20120501
      real   (dp) :: sliponly(0:NBMAX-1,0:NSAVMAX-1)
      real   (dp) :: deficit (0:NBMAX-1,0:NSAVMAX-1)

      real   (dp), pointer ::  deficit_gather  (:,:)=> null()
      real   (dp), pointer :: sliponly_gather  (:,:)=> null()
!<mh on 20120501

      End MODULE m_outputval
