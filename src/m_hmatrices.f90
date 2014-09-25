!C***
!C*** Module : Hmatrices
!C***                          LAST CHANGE:2012/04/20 18:44:39.
!C***
      Module Hmatrices
      USE precisions
      IMPLICIT NONE
      Type fullmatrix
        Integer(kp)          :: row_offset_f
        Integer(kp)          :: col_offset_f
        Integer(kp)          :: rows_f
        Integer(kp)          :: cols_f
        Real   (dp), Pointer :: el_f(:)
      End Type fullmatrix

      Type rkmatrix
        Integer(kp)          :: row_offset_r
        Integer(kp)          :: col_offset_r
        Integer(kp)          :: rows_r
        Integer(kp)          :: cols_r
        Integer(kp)          :: k_r
        Real   (dp), Pointer :: el_a(:)
        Real   (dp), Pointer :: el_b(:)
      End Type rkmatrix

      Type(fullmatrix), pointer :: fulm    (:) => null()
      Type(rkmatrix)  , pointer ::  rkm    (:) => null()
      Type(fullmatrix), pointer :: fulm_dum(:) => null()
      Type(rkmatrix)  , pointer ::  rkm_dum(:) => null()

      Integer krmax,kfmax
      End Module Hmatrices
