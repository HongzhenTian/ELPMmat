c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c       >This simple program to get the primivtive vectors 
c       $\delta$ strain, in order to calculate the independent
c       elastic constants of solids.
c       usage: C!!!!! Please first prepare the undeformed 
c       POSCAR in OLDPOS
c       >defevector.x
c       >Type defector.x > create new POSCAR in file fort.3
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        program defvector
        real*8 privect, strvect, strten, strain, pos, alat
        real*8 delta1, delta2, delta3, delta4, delta5, delta6
        dimension privect(3,3), strvect(3,3), strten(3,3), strain(6)
        dimension pos(50,3)
        character*10 bravlat, title, direct
        integer i, j, k, ntype, natomi, nn
        dimension natomi(10)

c%%%%%%%Read the undeformed primitive vector and atomic postion %%%%%%%%
        open(7,file='OLDPOS')
c%%% In first line of OLDPOS, please add the number
c%%%  of the type of atoms after the title

        read(7,*) title, ntype
        read(7,*) alat
        do i=1,3
                read(7,*) (privect(i,j),j=1,3)
                write(*,*) (privect(i,j),j=1,3)
        enddo
        read(7,*) (natomi(i),i=1,ntype)


        nn=0
        do i=1,ntype
                nn=nn+natomi(i)
        enddo


        read(7,*) direct
        do i=1,nn
                read(7,*) (pos(i,j),j=1,3)
        enddo

c%%%%%% Read the amti of strain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        read(*,*) delta1, delta2, delta3, delta4, delta5, delta6


c%%%%%% Define the strain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        strain(1)=delta1
        strain(2)=delta2
        strain(3)=delta3
        strain(4)=delta4
        strain(5)=delta5
        strain(6)=delta6
        
        write(*,*) (strain(i), i=1,6) 
        
c%%%%%% Define the strain tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        strten(1,1)=strain(1)+1.0
        strten(1,2)=0.5*strain(6)
        strten(1,3)=0.5*strain(5)
        strten(2,1)=0.5*strain(6)
        strten(2,2)=strain(2)+1.0
        strten(2,3)=0.5*strain(4)
        strten(3,1)=0.5*strain(5)
        strten(3,2)=0.5*strain(4)
        strten(3,3)=strain(3)+1.0
c%%%% Transfrom the privmitive vector to the new vector under strain %%%
c        strvect(i,j)=privect(i,j)*(I+strten(i,j))

        do k=1,3
          do i=1,3
            strvect(i,k)=0.0
              do j=1,3
                strvect(i,k)=strvect(i,k)+privect(i,j)*strten(j,k)
              enddo
          enddo
        enddo
c%%%%%% Write the new vector under strain %%%%%%
        do i=1,3
                write(*,100)(strvect(i,j),j=1,3)
        enddo
100     format(3f20.15)
c%%%%%% Create the POSCAR for total energy calculation %%%%%%
        write(3,'(A10)') title
        write(3,'(f15.10)') alat
        do i=1,3
                write(3,100) (strvect(i,j),j=1,3)
                write(*,*)(strvect(i,j), j=1,3)
        enddo
        write(3,'(10I4)') (natomi(i), i=1,ntype)
        write(3,'(A6)') Direct
        do i=1,nn
                write(3,100) (pos(i,j), j=1,3)
        enddo
c%%%%%%
        end

