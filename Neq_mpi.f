        program Ne
C parameter
C coords        xyz array of all molecules
C trajec        xyz array of one molecule
C boxlength     lower and upper limits of the sim box
C rprobe          probe radius used to measure the volume of a confining tube
C lee           end-to-end distance of a chain
C lpp           contour length of a primitive path
C griddel       grid size (distance between adjcent grid points)
C vol           tube volume
C chain         chain index
C numpts        total number of grid points
C grid          number of grid points along each axis
C s's           characters to read lines in dump file
C start         starting point in time
C finish        finishing point in time

        include 'mpif.h'
        include 'var.default'
        include 'var.options'

        real coords(num_atoms,6)
        real trajec(num_atoms_per_mol,3)
        real pts(max_rows,3)
        real boxlength(3,2)
        real rprobe
        real lee
        real lpp
        real griddel
        real vol
        integer chain
        integer*8 numpts
        integer*8 grid(3)
        character(len=10) s1(2), s2, s3(4), s4, s5(5), s9(8)
        real start, finish

C for mpi
        parameter (send_data_tag=9005, return_data_tag=9006)
        integer ierr
        integer status(MPI_STATUS_SIZE)
        integer my_id
        integer an_id
        integer num_procs
        integer root_process
        integer*8 start_row
        integer*8 end_row
        integer*8 avg_rows_per_process
        integer*8 num_rows
        integer*8 total_num
        integer*8 partial_num

        call cpu_time(start)

        call MPI_Init ( ierr )

        open(1,file='dump.PP') ! read in a dump file
        read(1,*) s1
        read(1,*) s2
        read(1,*) s3
        read(1,*) s4
        read(1,*) s5
        read(1,*) boxlength(1,1), boxlength(1,2)
        read(1,*) boxlength(2,1), boxlength(2,2)
        read(1,*) boxlength(3,1), boxlength(3,2)
        read(1,*) s9
        do j = 1, num_atoms, 1
           read(1,*) x1, x2, x3, x4, x5, x6
           coords(j,1)=x1
           coords(j,2)=x2
           coords(j,3)=x3
           coords(j,4)=x4
           coords(j,5)=x5
           coords(j,6)=x6
        enddo
        close(1)

        call sortrow(coords,num_atoms) ! sort rows with the first column ascending

        call MPI_Comm_Rank(MPI_Comm_World,my_id,ierr)
        call MPI_Comm_Size(MPI_Comm_World,num_procs,ierr)

        root_process=0
        open(2,file='chain.txt',status='old') ! read in chain number
        read(2,*) chain
        close(2)

        do i = chain, chain, 1
           trajec=coords(1+num_atoms_per_mol*(i-1):num_atoms_per_mol*i,
     &     4:6)
           call pbc(trajec,num_atoms_per_mol,boxlength) ! coords beyond pbc
           call gridn(trajec,num_atoms_per_mol,griddel,grid,vol)
           numpts=grid(1)*grid(2)*grid(3)

           if (my_id .eq. root_process) then
              write (*,*) 'Chain #', i, grid
              write (*,*) '* total grid pts num=', numpts
              write (*,*) '* volume=', vol
           endif

           avg_rows_per_process=numpts/num_procs

           do an_id = 0, num_procs-1
              if (my_id .eq. an_id) then
                 start_row=(an_id*avg_rows_per_process)+1
                 end_row=start_row+avg_rows_per_process-1
                 if (an_id .eq. (num_procs-1)) end_row=numpts

                 num_rows=end_row-start_row+1

                 call gridgen(trajec,num_atoms_per_mol,griddel,grid,
     &           start_row,num_rows,pts(1:num_rows,:))
              endif
           enddo

           if (my_id .eq. root_process) then

              total_num=0
              partial_num=tube(trajec,num_atoms_per_mol,rprobe,num_rows,
     &        pts(1:num_rows,:))
              total_num=total_num+partial_num

              do an_id = 1, num_procs-1
                 call MPI_Recv(partial_num,1,MPI_Int,MPI_Any_Source,
     &           MPI_Any_Tag,MPI_Comm_World,status,ierr)
                 total_num=total_num+partial_num
              enddo

              write (*,*) '* total num=', total_num
              lee=distcal(trajec(num_atoms_per_mol,:)-trajec(1,:))
              lpp=float(total_num)/float(numpts)*vol/(c_pi*rprobe**2)
              write (*,*) '* lee=', lee
              write (*,*) '* lpp=', lpp
           else
              partial_num=tube(trajec,num_atoms_per_mol,rprobe,num_rows,
     &        pts(1:num_rows,:))

              call MPI_Send(partial_num,1,MPI_Int,root_process,
     &        return_data_tag,MPI_Comm_World,ierr)
           endif
        enddo

        call MPI_Finalize(ierr)
        
        call cpu_time(finish)
        if (my_id .eq. root_process) then
           write (*,*) 'Time  =', finish-start, 'sec'
           write (*,*) lee, lpp
        endif
        stop
        end
 
C -----------------------------------------------------------

        subroutine gridn(xyz,num_atoms_per_mol,griddel,grid,vol)
C xyz           xyz array of a molecule
C griddel       grid size
C ul, ll        upper and lower limits along each axis
C lref          reference position 
C ltemp         temporary position
C l             box length along each axis
C vol           tube volume
C grid          number of grid points along each axis

        real xyz(num_atoms_per_mol,3)
        real griddel
        real ul(3)
        real ll(3)
        real lref
        real ltemp
        real l(3)
        real vol
        integer*8 grid(3)

        ll(1)=minval(xyz(:,1))
        ll(2)=minval(xyz(:,2))
        ll(3)=minval(xyz(:,3))
        ul(1)=maxval(xyz(:,1))
        ul(2)=maxval(xyz(:,2))
        ul(3)=maxval(xyz(:,3))
         
        do i = 1, 3, 1
           l(i)=float(ceiling(ul(i))-floor(ll(i)))
        enddo

        do i = 1, 3, 1
           lref=0
           ltemp=0
           grid(i)=0
           do while (ltemp .le. l(i))
              grid(i)=grid(i)+1
              ltemp=lref+griddel*(grid(i)-1)
           enddo
        enddo

        vol=l(1)*l(2)*l(3)

        return
        end

C -----------------------------------------------------------
C generates a grid

        subroutine gridgen(xyz,num_atoms_per_mol,griddel,grid,start_row,
     &  num_rows,pts)
C xyz           xyz array of a molecule
C griddel       grid size
C num_rows      number of grid points per each processor
C pts           xyz array of all grid points
C ul, ll        upper and lower limits along each axis
C lref          reference position 
C ltemp         temporary position
C l             box length along each axis
C grid          number of grid points along each axis
C row
C start_row
C end_row
C idx           
C i, j, k_ind's s

        real xyz(num_atoms_per_mol,3)
        real griddel
        integer*8 num_rows
        real pts(num_rows,3)
        real ul(3)
        real ll(3)
        real lref
        real ltemp
        real l(3)
        integer*8 grid(3)
        integer*8 row
        integer*8 start_row
        integer*8 end_row
        integer*8 idx
        integer*8 i_ind
        integer*8 j_ind
        integer*8 k_ind

        ll(1)=minval(xyz(:,1))
        ll(2)=minval(xyz(:,2))
        ll(3)=minval(xyz(:,3))
        ul(1)=maxval(xyz(:,1))
        ul(2)=maxval(xyz(:,2))
        ul(3)=maxval(xyz(:,3))
         
        do i = 1, 3, 1
           l(i)=float(ceiling(ul(i))-floor(ll(i)))
        enddo

        do i = 1, 3, 1
           lref=0
           ltemp=0
           grid(i)=0
           do while (ltemp .le. l(i))
              grid(i)=grid(i)+1
              ltemp=lref+griddel*(grid(i)-1)
           enddo
        enddo

        idx=1
        end_row=start_row+num_rows-1
        do row = start_row, end_row, 1
           i_ind=row/(grid(1)*grid(2))
           j_ind=(row-i_ind*(grid(1)*grid(2)))/grid(1)
           k_ind=row-j_ind*grid(1)-i_ind*(grid(1)*grid(2))
           
           dx=floor(ll(1))+griddel*(k_ind-1)
           dy=floor(ll(2))+griddel*(j_ind-1)
           dz=floor(ll(3))+griddel*(i_ind-1)

           pts(idx,1)=dx
           pts(idx,2)=dy
           pts(idx,3)=dz

           idx=idx+1
        enddo

        return
        end

C -----------------------------------------------------------
C counts the number of grid points falling withing the confining tube

        real function tube(xyz,num_atoms_per_mol,rprobe,numpts,pts)
C xyz           xyz array of a molecule
C rprobe        probe radius to measure volume 
C numpts        total number of grid points
C pts           xyz array of all grid points
C dist          distance between a grid point and a atom/bead
C dx, dy, dz    increment in each axis to generate grid points
C num           number of grid points falling within 'rprobe'
C flag          flag for while loop
C idx           index fo while loop

        real xyz(num_atoms_per_mol,3) 
        real rprobe
        integer*8 numpts
        real pts(numpts,3)
        real dist
        real dx
        real dy
        real dz
        integer*8 num 
        integer flag
        integer*8 idx

        num=0

        do j = 1, numpts, 1
           idx=1
           flag=0
           do while ((flag .eq. 0) .and. (idx .le. num_atoms_per_mol))
              dx=pts(j,1)-xyz(idx,1)
              dy=pts(j,2)-xyz(idx,2)
              dz=pts(j,3)-xyz(idx,3)
              dist=sqrt(dx**2+dy**2+dz**2)
              if (dist .le. rprobe) then
                  num=num+1
                  flag=1
              endif
              idx=idx+1
           enddo 
        enddo

        tube=num

        return
        end

C -----------------------------------------------------------
C sorts rows of an array

        subroutine sortrow(A,num_rows)
        real A(num_rows,6), buf(6)
        integer irow, krow

        do irow = 1, num_rows, 1
           krow=minloc(A(irow:num_rows,1),dim=1)+irow-1
           buf(:)=A(irow,:)
           A(irow,:)=A(krow,:)
           A(krow,:)=buf(:)
        enddo

        return
        end

C -----------------------------------------------------------
C relocates an atom beyond boundaries

        subroutine pbc(xyz,num_atoms_per_mol,boxlength)
        real xyz(num_atoms_per_mol,3)
        real temp(num_atoms_per_mol,3)
        real boxlength(3,2)
        real reflen(3)
        real vec(3)
        real bl(3)

        temp(1,:)=xyz(1,:)

        do j = 1, 3, 1
           bl(j)=boxlength(j,2)-boxlength(j,1)
           reflen(j)=bl(j)-10.0
        enddo

        do i = 1, (num_atoms_per_mol-1), 1
           do j = 1, 3, 1
              vec(j)=xyz(i+1,j)-xyz(i,j)
              
              if (vec(j) .lt. reflen(j)*(-1)) then
                 vec(j)=vec(j)+bl(j)
              endif
              if (vec(j) .gt. reflen(j)) then
                 vec(j)=vec(j)-bl(j)
              endif
           enddo
           temp(i+1,:)=temp(i,:)+vec(:)
        enddo

        xyz=temp

        return
        end

C -----------------------------------------------------------
C calculates the length of a vector made by two points

        real function distcal(a)
        real a(3)

        distcal=sqrt(a(1)**2+a(2)**2+a(3)**2)

        return
        end

C -----------------------------------------------------------
