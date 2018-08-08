      program demo
      implicit none

      call conrec
      write(*,*) 'Finished executing CONREC'
c
      end
c
c======================================================================

      SUBROUTINE conrec


C______________________________________________________________
c     global parameters - to be deactivated in Mistress
      REAL*4 C(20,20)
      integer*4 NX,NY
      parameter (NX = 20,NY = 20)
C______________________________________________________________
      
      INTEGER*4 NC
      PARAMETER (NC = 1)
c
      REAL*4 ZLVL(NC),C_test(NX*NY),GRIDLVL(NX,NY),VERTEX(NX+1,NY+1)
      INTEGER*4 I,J
C______________________________________________________________

C     Specify contour level
      ZLVL(1) = 0.5

C-----------------------------------------------------------------------  
C to be deactivated in mistress    
      open (unit=15, file="finger.txt", status='old')
      read (15, 120) C_test
      do 111 i=1, NX
          do 112 j=1,NY
            write (16,120) C_test((i-1)*NY+j)
            C(j,i) = C_test((i-1)*ny+j)
112   continue
111   continue
120   format (F7.2)

C     Print coordinates and concentration for sanity check
      OPEN (UNIT = 12, FILE = "print_finger.txt", STATUS='UNKNOWN')
      do 210 i=1,NX
         do 110 j=1,NY
            WRITE(12,1020) I*1.0, J*1.0, C(i,j)
110      continue
210   continue
1020  FORMAT(1X, 3(F6.3, 2X))

C------------------------------------------------------------------------
      OPEN (UNIT = 12, FILE = "contour.txt", STATUS='UNKNOWN')
C Set the threshold
      DO 2110 I=1,NX
         DO 1110 J=1,NY
            IF (C(I,J).GE.ZLVL(1)) THEN
                GRIDLVL(I,J) = 1.0
            ELSE
                GRIDLVL(I,J) = 0.0
            ENDIF
            WRITE(*,*) GRIDLVL(I,J)
1110      CONTINUE
2110   CONTINUE

C March through the grid 
C Case 0
      DO 2111 I=1,NX-1
         DO 1111 J=1,NY-1
            IF (GRIDLVL(I,J).EQ.0 .AND. GRIDLVL(I+1,J).EQ.0 .AND. 
     &          GRIDLVL(I,J+1).EQ.0 .AND. GRIDLVL(I+1,J+1).EQ.0) THEN
                VERTEX(I+1,J+1) = 0
C Case 1
            ELSEIF (GRIDLVL(I,J).EQ.0 .AND. GRIDLVL(I+1,J).EQ.1 .AND. 
     &          GRIDLVL(I,J+1).EQ.0 .AND. GRIDLVL(I+1,J+1).EQ.0) THEN
                VERTEX(I+1,J+1) = 1
                WRITE(12,1010) ZLVL(1), REAL(I)+0.5, REAL(J),
     &                                  REAL(I)+1.0, REAL(J)+0.5
C Case 2
            ELSEIF (GRIDLVL(I,J).EQ.0 .AND. GRIDLVL(I+1,J).EQ.0 .AND. 
     &          GRIDLVL(I,J+1).EQ.0 .AND. GRIDLVL(I+1,J+1).EQ.1) THEN
                VERTEX(I+1,J+1) = 2
                WRITE(12,1010) ZLVL(1), REAL(I)+1.0, REAL(J)+0.5,
     &                                  REAL(I)+0.5, REAL(J)+1.0
C Case 3
            ELSEIF (GRIDLVL(I,J).EQ.0 .AND. GRIDLVL(I+1,J).EQ.1 .AND. 
     &          GRIDLVL(I,J+1).EQ.0 .AND. GRIDLVL(I+1,J+1).EQ.1) THEN
                VERTEX(I+1,J+1) = 3
                WRITE(12,1010) ZLVL(1), REAL(I)+0.5, REAL(J),
     &                                  REAL(I)+0.5, REAL(J)+1.0
C Case 4
            ELSEIF (GRIDLVL(I,J).EQ.0 .AND. GRIDLVL(I+1,J).EQ.0 .AND. 
     &          GRIDLVL(I,J+1).EQ.1 .AND. GRIDLVL(I+1,J+1).EQ.0) THEN
                VERTEX(I+1,J+1) = 4
                WRITE(12,1010) ZLVL(1), REAL(I), REAL(J)+0.5,
     &                                  REAL(I)+0.5, REAL(J)+1.0
C Case 5
            ELSEIF (GRIDLVL(I,J).EQ.0 .AND. GRIDLVL(I+1,J).EQ.1 .AND. 
     &          GRIDLVL(I,J+1).EQ.1 .AND. GRIDLVL(I+1,J+1).EQ.0) THEN
                VERTEX(I+1,J+1) = 5
                WRITE(12,1010) ZLVL(1), REAL(I), REAL(J)+0.5,
     &                                  REAL(I)+0.5, REAL(J)
                WRITE(12,1010) ZLVL(1), REAL(I)+1.0, REAL(J)+0.5,
     &                                  REAL(I)+0.5, REAL(J)+1.0
C Case 6
            ELSEIF (GRIDLVL(I,J).EQ.0 .AND. GRIDLVL(I+1,J).EQ.0 .AND. 
     &          GRIDLVL(I,J+1).EQ.1 .AND. GRIDLVL(I+1,J+1).EQ.1) THEN
                VERTEX(I+1,J+1) = 6
                WRITE(12,1010) ZLVL(1), REAL(I), REAL(J)+0.5,
     &                                  REAL(I)+1.0, REAL(J)+0.5
C Case 7
            ELSEIF (GRIDLVL(I,J).EQ.0 .AND. GRIDLVL(I+1,J).EQ.1 .AND. 
     &          GRIDLVL(I,J+1).EQ.1 .AND. GRIDLVL(I+1,J+1).EQ.1) THEN
                VERTEX(I+1,J+1) = 7
                WRITE(12,1010) ZLVL(1), REAL(I), REAL(J)+0.5,
     &                                  REAL(I)+0.5, REAL(J)
C Case 8
            ELSEIF (GRIDLVL(I,J).EQ.1 .AND. GRIDLVL(I+1,J).EQ.0 .AND. 
     &          GRIDLVL(I,J+1).EQ.0 .AND. GRIDLVL(I+1,J+1).EQ.0) THEN
                VERTEX(I+1,J+1) = 8
                WRITE(12,1010) ZLVL(1), REAL(I), REAL(J)+0.5,
     &                                  REAL(I)+0.5, REAL(J)
C Case 9
            ELSEIF (GRIDLVL(I,J).EQ.1 .AND. GRIDLVL(I+1,J).EQ.1 .AND. 
     &          GRIDLVL(I,J+1).EQ.0 .AND. GRIDLVL(I+1,J+1).EQ.0) THEN
                VERTEX(I+1,J+1) = 9
                WRITE(12,1010) ZLVL(1), REAL(I), REAL(J)+0.5,
     &                                  REAL(I)+1.0, REAL(J)+0.5
C Case 10
            ELSEIF (GRIDLVL(I,J).EQ.1 .AND. GRIDLVL(I+1,J).EQ.0 .AND. 
     &          GRIDLVL(I,J+1).EQ.0 .AND. GRIDLVL(I+1,J+1).EQ.1) THEN
                VERTEX(I+1,J+1) = 10
                WRITE(12,1010) ZLVL(1), REAL(I)+0.5, REAL(J),
     &                                  REAL(I)+1.0, REAL(J)+0.5
                WRITE(12,1010) ZLVL(1), REAL(I), REAL(J)+0.5,
     &                                  REAL(I)+0.5, REAL(J)+1.0
C Case 11
            ELSEIF (GRIDLVL(I,J).EQ.1 .AND. GRIDLVL(I+1,J).EQ.1 .AND. 
     &          GRIDLVL(I,J+1).EQ.0 .AND. GRIDLVL(I+1,J+1).EQ.1) THEN
                VERTEX(I+1,J+1) = 11
                WRITE(12,1010) ZLVL(1), REAL(I)+0.5, REAL(J)+1.0,
     &                                  REAL(I), REAL(J)+0.5
C Case 12
            ELSEIF (GRIDLVL(I,J).EQ.1 .AND. GRIDLVL(I+1,J).EQ.0 .AND. 
     &          GRIDLVL(I,J+1).EQ.1 .AND. GRIDLVL(I+1,J+1).EQ.0) THEN
                VERTEX(I+1,J+1) = 12
                WRITE(12,1010) ZLVL(1), REAL(I)+0.5, REAL(J),
     &                                  REAL(I)+0.5, REAL(J)+1.0
C Case 13
            ELSEIF (GRIDLVL(I,J).EQ.1 .AND. GRIDLVL(I+1,J).EQ.1 .AND. 
     &          GRIDLVL(I,J+1).EQ.1 .AND. GRIDLVL(I+1,J+1).EQ.0) THEN
                VERTEX(I+1,J+1) = 13
                WRITE(12,1010) ZLVL(1), REAL(I)+1.0, REAL(J)+0.5,
     &                                  REAL(I)+0.5, REAL(J)+1.0
C Case 14
            ELSEIF (GRIDLVL(I,J).EQ.1 .AND. GRIDLVL(I+1,J).EQ.0 .AND. 
     &          GRIDLVL(I,J+1).EQ.1 .AND. GRIDLVL(I+1,J+1).EQ.1) THEN
                VERTEX(I+1,J+1) = 14
                WRITE(12,1010) ZLVL(1), REAL(I)+0.5, REAL(J),
     &                                  REAL(I)+1.0, REAL(J)+0.5
C Case 15
            ELSEIF (GRIDLVL(I,J).EQ.1 .AND. GRIDLVL(I+1,J).EQ.1 .AND. 
     &          GRIDLVL(I,J+1).EQ.1 .AND. GRIDLVL(I+1,J+1).EQ.1) THEN
                VERTEX(I+1,J+1) = 15

1010        FORMAT(5(F6.2, 1X))

            ENDIF

1111   CONTINUE
2111   CONTINUE

      CLOSE (12)

      END

