      program getforces
c utility for analyzing the distances and forces
c uses the VASP OUTCAR file
c marcel sluiter Feb 19, 2002
      implicit none
      real(8), parameter:: distmax=3.2d0
      integer nions,i,j,in,n_dist1
      integer, allocatable:: n_dist(:)
      real(8) frc_average,dst_average,dst2,dr(3),
     &   mindst,maxdst,dst,f2,nn(3,8),tr(3,3),dl(3,8),
     &   frc_max
      real(8), allocatable:: xyz(:,:),frc(:,:)
      character*100 line,filein
      data nn(1:3,1)/0d0,0d0,0d0/
      data nn(1:3,2)/1d0,0d0,0d0/
      data nn(1:3,3)/0d0,1d0,0d0/
      data nn(1:3,4)/0d0,0d0,1d0/
      data nn(1:3,5)/1d0,1d0,0d0/
      data nn(1:3,6)/0d0,1d0,1d0/
      data nn(1:3,7)/1d0,0d0,1d0/
      data nn(1:3,8)/1d0,1d0,1d0/

      write(*,'(a)')'enter name of OUTCAR file'
      read *,filein
      open(10,file=trim(filein),status='old',form='formatted')
      open(20,file='distancesforces.txt')

10    read(10,'(a)') line       !to number of ions
      if( line(1:17).ne.'   number of dos ' ) goto 10
      read(line(66:80),*) nions
      write(*,'(a,i6)')'# ions:',nions
      allocate( xyz(3,nions),frc(3,nions),n_dist(nions) )

      write(*,'(a,f10.3)')'#cutoff distance set to:',distmax
      write(20,'(a,f10.3)')'#cutoff distance set to:',distmax
      write(*,'(a)') 'distance(avrg,min,max),force(avrg,max)'
      write(20,'(a)') 'distance(avrg,min,max),force(avrg,max)'

100   continue

15    read(10,'(a)',end=999,err=999) line       !to translations
      if( line(1:29).ne.'      direct lattice vectors ' ) goto 15
      read(10,'(3x,3f13.9)')tr(1:3,1)
      read(10,'(3x,3f13.9)')tr(1:3,2)
      read(10,'(3x,3f13.9)')tr(1:3,3)
      dl = matmul(tr,nn)

20    read(10,'(a)',end=999,err=999) line       !to positions
      if( line(1:10).ne.' POSITION ' ) goto 20
      read(10,'(a)') line       !line of dashes
      do i=1,nions
        read(10,*)xyz(:,i),frc(:,i)
      enddo

c actual calculations
      frc_average=0d0
      frc_max=0d0
      dst_average=0d0
      n_dist1=0
      mindst=huge(1d0)
      maxdst=0d0
      do i=1,nions
        f2=sqrt( dot_product(frc(:,i),frc(:,i)) )
        frc_average=frc_average+f2
        frc_max=max(frc_max,f2)
        do j=1,nions
          if(i.ne.j) then
            do in=1,8
             dr=xyz(:,i)-xyz(:,j) + dl(:,in)
             dst2=dot_product(dr,dr)
             if(dst2.le.distmax**2) then
               dst=sqrt(dst2)
               dst_average=dst_average+dst
               n_dist1=n_dist1+1
               mindst=min(mindst,dst)
               maxdst=max(maxdst,dst)
             endif
            enddo
          endif
        enddo
      enddo
      frc_average=frc_average/real(nions)
      dst_average=dst_average/real(n_dist1)
      write(*,'(5f10.3)')
     &    dst_average,mindst,maxdst,frc_average,frc_max
      write(20,'(5f10.3)')
     &    dst_average,mindst,maxdst,frc_average,frc_max

c once again
      goto 100

999   continue
      end program getforces

