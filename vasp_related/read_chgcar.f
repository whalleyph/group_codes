      program read_chgcar
c utility for analyzing the chargedensity from VASP as written
c in the CHGCAR file
c this utility computes local charges and local moments using
c Voronoi polyhedra
c it uses the assumption that VASP properly centers the cluster in the cell
c it also outputs the chargedensity and spindensity for vaspviewer
c (using smaller grids if necessary)
c marcel sluiter Sept 27-30 2001, Feb 19 2002
      implicit none
      integer, parameter:: ngrid=30
      integer natoms,atoms(100),iatoms,nx,ny,nz,ix,iy,iz,
     &   jx,jy,jz,nearest_atom,ia,in,n_neighborcells,ib,ic,
     &   nargs
      real(8) tr(3,3),dr(3),dt(3),dd,smallest_distance,scale,
     &   ddd,nn(3,27),rndm,centerofcluster(3),error,aux(100)
      real(8), allocatable:: density(:,:,:),xyz(:,:),xyzs(:,:),
     &   localcharge(:),localmoment(:),points(:)
      character*200 line,nxnynz
      character*255 fchgcar
      logical orthogonal
      data nn(1:3,1)/0d0,0d0,0d0/
      data nn(1:3,2)/1d0,0d0,0d0/
      data nn(1:3,3)/-1d0,0d0,0d0/
      data nn(1:3,4)/0d0,1d0,0d0/
      data nn(1:3,5)/0d0,-1d0,0d0/
      data nn(1:3,6)/0d0,0d0,1d0/
      data nn(1:3,7)/0d0,0d0,-1d0/
      data nn(1:3,8)/1d0,1d0,0d0/
      data nn(1:3,9)/-1d0,1d0,0d0/
      data nn(1:3,10)/1d0,-1d0,0d0/
      data nn(1:3,11)/-1d0,-1d0,0d0/
      data nn(1:3,12)/0d0,1d0,1d0/
      data nn(1:3,13)/0d0,-1d0,1d0/
      data nn(1:3,14)/0d0,1d0,-1d0/
      data nn(1:3,15)/0d0,-1d0,-1d0/
      data nn(1:3,16)/1d0,0d0,1d0/
      data nn(1:3,17)/-1d0,0d0,1d0/
      data nn(1:3,18)/1d0,0d0,-1d0/
      data nn(1:3,19)/-1d0,0d0,-1d0/
      data nn(1:3,20)/1d0,1d0,1d0/
      data nn(1:3,21)/-1d0,1d0,1d0/
      data nn(1:3,22)/1d0,-1d0,1d0/
      data nn(1:3,23)/-1d0,-1d0,1d0/
      data nn(1:3,24)/1d0,1d0,-1d0/
      data nn(1:3,25)/-1d0,1d0,-1d0/
      data nn(1:3,26)/1d0,-1d0,-1d0/
      data nn(1:3,27)/-1d0,-1d0,-1d0/


      nargs = iargc()

      if (nargs == 1) then
          call getarg(1, fchgcar)
      else
        fchgcar = 'CHGCAR'
      end if

      open(10,file=fchgcar,status='old',form='formatted')
      open(20,file='totalcharge.txt')
      open(30,file='spindensity.txt')
      open(40,file='localproperties.txt')
c read header
      read(10,'(a)') line        !title
      write(20,'(a)') trim(line)
      write(30,'(a)') trim(line)
      write(*,'(a)') trim(line)
      read(10,'(a)') line        !scalefactor
      write(20,'(a)') trim(line)
      write(30,'(a)') trim(line)
      write(*,'(a)') trim(line)
      read(10,'(a)') line        !translation 1
      read(line,*) tr(1:3,1)
      write(20,'(a)') trim(line)
      write(30,'(a)') trim(line)
      write(*,'(3f12.4)') tr(1:3,1)
      read(10,'(a)') line        !translation 2
      read(line,*) tr(1:3,2)
      write(20,'(a)') trim(line)
      write(30,'(a)') trim(line)
      write(*,'(3f12.4)') tr(1:3,2)
      read(10,'(a)') line        !translation 3
      read(line,*) tr(1:3,3)
      write(20,'(a)') trim(line)
      write(30,'(a)') trim(line)
      write(*,'(3f12.4)') tr(1:3,3)
      read(10,'(a)') line        !number of atoms for each species
      write(20,'(a)') trim(line)
      write(30,'(a)') trim(line)
      atoms=0
      read(line,*,end=10,err=10) atoms
10    continue
      natoms=sum(atoms)
      write(*,'(i6)') natoms
      read(10,'(a)') line        !should say "direct"
      write(20,'(a)') trim(line)
      write(30,'(a)') trim(line)
      write(*,'(a)') trim(line)
      allocate(xyz(3,natoms),xyzs(3,natoms),points(natoms),
     &         localmoment(natoms),localcharge(natoms))
      do iatoms=1,natoms
        read(10,'(a)') line        !xyz of each atom
        read(line,*) xyz(:,iatoms)
        write(20,'(a)') trim(line)
        write(30,'(a)') trim(line)
        write(*,'(3f12.4)') xyz(:,iatoms)
      enddo

c done with reading header, now read nx,ny,nz of mesh
20    read(10,'(a)') nxnynz        !nx,ny,nz of each atom
      read(nxnynz,*,end=20,err=20) nx,ny,nz
c rescale grid if it is too large
      jx=max(1,nx/ngrid)
      jy=max(1,ny/ngrid)
      jz=max(1,nz/ngrid)
      write(20,'(/,3i10)') nx/jx,ny/jy,nz/jz
      allocate(density(nx,ny,nz))
      read(10,*) density        !total charge at each xyz
      write(20,'(6(1x,e12.6))')
     &  (((density(ix,iy,iz),ix=1,nx,jx),iy=1,ny,jy),iz=1,nz,jz)
      close(20)
c a scale factor needs to be defined
      scale=1d0/real(nx*ny*nz)
c are translations orthogonal
      orthogonal=
     &     (abs( dot_product(tr(:,1),tr(:,2)) ).le.1d-6).and.
     &     (abs( dot_product(tr(:,2),tr(:,3)) ).le.1d-6).and.
     &     (abs( dot_product(tr(:,3),tr(:,1)) ).le.1d-6)
      n_neighborcells=27
c      if(orthogonal) n_neighborcells=1
cc check that center of cluster is at 1/2 1/2 1/2
c      centerofcluster(1) = sum(xyz(1,:))/real(natoms)
c      centerofcluster(2) = sum(xyz(2,:))/real(natoms)
c      centerofcluster(3) = sum(xyz(3,:))/real(natoms)
c      write(*,'(a,3f9.3)')'center of cluster',centerofcluster
c      if( abs(centerofcluster(1)-0.5d0).gt.1d-3  .or.
c     &    abs(centerofcluster(2)-0.5d0).gt.1d-3  .or.
c     &    abs(centerofcluster(3)-0.5d0).gt.1d-3) then
c        if(orthogonal) then
c          write(*,'(a,/,a)')'center of cluster deviates too much!',
c     &                 'considering all neighboring cells as well'
c          n_neighborcells=27
c        endif
c      endif

c***** determine approximate local charge *****
c scale atomic coordinates of that of the grid
      xyzs(1,:)=xyz(1,:)*real(nx)
      xyzs(2,:)=xyz(2,:)*real(ny)
      xyzs(3,:)=xyz(3,:)*real(nz)
c estimate errors in integrals of local charges
      iatoms=0
      do ix=1,100
       if(atoms(ix).ge.2) then
        do iy=1,atoms(ix)
         iatoms=iatoms+1
         ia=modulo(nint(xyzs(1,iatoms)),nx)+1
         ib=modulo(nint(xyzs(2,iatoms)),ny)+1
         ic=modulo(nint(xyzs(3,iatoms)),nz)+1
         aux(iy)=density(ia,ib,ic)
        write(40,'(a,i3,g14.6)')
     &  'max in local charge for atom:',iatoms,aux(iy)
        enddo
        error=0d0
        do iy=1,atoms(ix)-1
        do iz=iy+1,atoms(ix)
          error=max( error,abs(aux(iy)-aux(iz)) )
        enddo
        enddo
        write(40,'(a,i3,g11.3)')
     &  'error in local charge for atomtype:',ix,error*scale
       else
        iatoms=iatoms+atoms(ix)
       endif
      enddo
c scale coordinates of unitcell corners
      nn(1,:)=nn(1,:)*real(nx)
      nn(2,:)=nn(2,:)*real(ny)
      nn(3,:)=nn(3,:)*real(nz)
c for each point, determine which atom is nearest, add the local
c chargedensity to the nearest atom
      localcharge=0d0
      points=0d0
      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
        smallest_distance=huge(1d0)
        nearest_atom=0
        do ia=1,natoms
          dd=huge(1d0)
          do in=1,n_neighborcells
            dr(1)=real(ix)-xyzs(1,ia)-nn(1,in)
            dr(2)=real(iy)-xyzs(2,ia)-nn(2,in)
            dr(3)=real(iz)-xyzs(3,ia)-nn(3,in)
            dt = matmul(dr,tr)
            ddd= dot_product(dt,dt)
            dd = min(dd,ddd)
          enddo
          if(dd.lt.smallest_distance) then
            nearest_atom=ia
            smallest_distance=dd
          else if(dd.eq.smallest_distance) then
            call random_number(rndm)
            if(rndm.le.5d-1) nearest_atom=ia
          endif
c          write(*,'(3i3,7f8.2,i3)')ix,iy,iz,dr,dt,dd,nearest_atom
        enddo
        if(nearest_atom.eq.0) stop 'nearest atom not found'
        if(nearest_atom.gt.natoms) stop 'nearest atom bad'
        localcharge(nearest_atom)=localcharge(nearest_atom)+
     &                             density(ix,iy,iz)
        points(nearest_atom)=points(nearest_atom)+1d0
      enddo
      enddo
      enddo
      deallocate(density)
      localcharge=localcharge*scale
      write(40,'(a)')'local charges for each atom:'
      write(40,'(a)')'#atom, xyz, local charge, volume fraction'
      do ia=1,natoms
       write(40,'(i4,3(1x,f11.6),1x,2f12.6)')
     &   ia,xyz(:,ia),localcharge(ia),points(ia)/real(nx*ny*nz)
      enddo
      write(40,'(a,f12.6)')'total charge:',sum(localcharge)

c done with reading total chargedensity, now try spin-density
30    read(10,'(a)',end=999) line       !other crud
      if( trim(line).ne.trim(nxnynz) ) goto 30
      write(30,'(/,3i10)') nx/jx,ny/jy,nz/jz
      allocate(density(nx,ny,nz))
      read(10,*) density        !spin density at each xyz
      write(30,'(6(1x,e12.6))')
     &  (((density(ix,iy,iz),ix=1,nx,jx),iy=1,ny,jy),iz=1,nz,jz)

c***** determine approximate local moments *****
c for each point, determine which atom is nearest, add the local
c spindensity to the nearest atom
      localmoment=0d0
      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
        smallest_distance=huge(1d0)
        nearest_atom=0
        do ia=1,natoms
          dd=huge(1d0)
          do in=1,n_neighborcells
            dr(1)=real(ix)-xyzs(1,ia)-nn(1,in)
            dr(2)=real(iy)-xyzs(2,ia)-nn(2,in)
            dr(3)=real(iz)-xyzs(3,ia)-nn(3,in)
            dt = matmul(dr,tr)
            ddd= dot_product(dt,dt)
            dd = min(dd,ddd)
          enddo
          if(dd.lt.smallest_distance) then
            nearest_atom=ia
            smallest_distance=dd
          else if(dd.eq.smallest_distance) then
            call random_number(rndm)
            if(rndm.le.5d-1) nearest_atom=ia
          endif
        enddo
        if(nearest_atom.eq.0) stop 'nearest atom not found'
        localmoment(nearest_atom)=localmoment(nearest_atom)+
     &                             density(ix,iy,iz)
      enddo
      enddo
      enddo
      deallocate(density)
      localmoment=localmoment*scale
      write(40,'(/,a)')'local moments for each atom:'
      write(40,'(a)')'#atom, xyz, local charge, local moment'
      do ia=1,natoms
        write(40,'(i4,3(1x,f11.6),2(1x,f12.6))')
     &    ia,xyz(:,ia),localcharge(ia),localmoment(ia)
      enddo
      write(40,'(a,f12.6)')'total charge:',sum(localcharge)
      write(40,'(a,f12.6)')'total moment:',sum(localmoment)
      write(40,'(a,100(1x,f6.3))')'MAGMOM=',localmoment

999   continue
      end program read_chgcar

