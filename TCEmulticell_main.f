      program celldynamic
      implicit none

      real*8 l,dt
      integer nd,npart
      parameter (nd=6)
      parameter (l=240.0)      
      parameter (npart=nd*nd*nd)
      parameter (dt=0.02)
      real*8 rad
      parameter (rad=5.0)
      integer nsimu
      parameter (nsimu=5000000)
      real*8 ws ! adhesion energy
      real*8 ws_IN,ws_IT,ws_NT,ws_TT,ws_NN,ws_II ! I for immune cell; N for normal cell; T for tumor cell
      parameter (ws_IN=0.0382*1)
      parameter (ws_IT=0.0674*1)
      parameter (ws_NT=0.0)
      parameter (ws_II=0.0)
      parameter (ws_NN=0.0)
      parameter (ws_TT=0.0)	  
      real*8 pai
      parameter (pai=3.1415926)
      real*8 pois_ratio
      parameter (pois_ratio=1/3)
      real*8 Young
      parameter (Young=1)
      real*8 cm_frict_const
      parameter (cm_frict_const=0.4)
      real*8 cc_frict_const
      parameter (cc_frict_const=0.5)
      real*8 diffu_const
      parameter (diffu_const=100)
      integer niter
      parameter (niter=100)
      real*8 iter_tol
      parameter (iter_tol=0.0001)

      integer i,j,k
      real*8 x(npart),y(npart),z(npart)
      real*8 vx(npart),vy(npart),vz(npart)
      real*8 f_atr_x(npart)
      real*8 f_atr_y(npart)
      real*8 f_atr_z(npart)
      real*8 f_rps_x(npart)
      real*8 f_rps_y(npart)
      real*8 f_rps_z(npart)
      real*8 f_rdm_x(npart)
      real*8 f_rdm_y(npart)
      real*8 f_rdm_z(npart)
      real*8 vx_old(npart),vy_old(npart),vz_old(npart)
      real*8 vx_new(npart),vy_new(npart),vz_new(npart)
      real*8 ax(npart,npart),ay(npart,npart),az(npart,npart)
      real*8 bx(npart),by(npart),bz(npart)
      integer particle
      real*8 dl
      real*8 dij
      integer itime
      integer molecule,neighbor
      real*8 repulsive_parameter
      real*8 dir_x,dir_y,dir_z
      real*8 amp_dir
      real*8 scaled_dir_x,scaled_dir_y,scaled_dir_z
      real*8 RMSD
      character*3 cell_type_index(npart)
      character*3 cell_res_index(npart)
      real*8 cell_type_seed
      integer IN_num,IT_num,NT_num,TT_num,II_num,NN_num
      
      real rand3
      double precision r3
      real rand4
      double precision r4
      real rand5
      double precision r5
	  
      r3=5.0
      r4=5.0
      r5=5.0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   v3: the cells are confined in the box and there are two kinds
c       of different cells with different adhesive intensity
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccc
c  set initial positions
cccccccccccccccccccccccccccccccccccccc
c      print*,'*'
      particle=0
      dl=l/nd
      do i=1,nd
         do j=1,nd
            do k=1,nd
               particle=particle+1
c               print*,i,j,k,particle
               x(particle)=(real(i)-0.5)*dl+(dl-2*rad)*(rand3(r3)-0.5)
               y(particle)=(real(j)-0.5)*dl+(dl-2*rad)*(rand4(r4)-0.5) 
               z(particle)=(real(k)-0.5)*dl+(dl-2*rad)*(rand3(r3)-0.5)
            enddo
         enddo
      enddo
      do j=1,npart
         x(j)=x(j)-l*anint(x(j)/l)
         y(j)=y(j)-l*anint(y(j)/l)
         z(j)=z(j)-l*anint(z(j)/l)
      enddo

ccccccccccccccccccccccccccccccccccccccccccc
c initiate the types of cells
ccccccccccccccccccccccccccccccccccccccccccc
c      print*,'**'
      do i=1,npart
         cell_type_seed=rand4(r4)
         if(cell_type_seed.gt.0.8)then
            cell_type_index(i)='ICL' ! ICL means immune cells
         elseif((cell_type_seed.le.0.8)
     &           .and.(cell_type_seed.gt.0.4))then
            cell_type_index(i)='NCL'  ! NCL means normal cells
         else
            cell_type_index(i)='TCL' ! ECL means tumor cells
         endif
      enddo

cccccccccccccccccccccccccccccccccccccccccccc
c  main simulation
cccccccccccccccccccccccccccccccccccccccccccc

      do itime=1,nsimu
c         print*,itime
c>  calculate the attractive and repulsive force of each cell
         do molecule=1,npart
            f_atr_x(molecule)=0
            f_atr_y(molecule)=0
            f_atr_z(molecule)=0
            f_rps_x(molecule)=0
            f_rps_y(molecule)=0
            f_rps_z(molecule)=0
         enddo
         do molecule=1,npart-1
            do neighbor=molecule+1,npart
               dij=sqrt((x(molecule)-x(neighbor))**2+
     &              (y(molecule)-y(neighbor))**2+
     &              (z(molecule)-z(neighbor))**2)
               if(dij.lt.2*rad)then
                  if((cell_type_index(molecule).eq.'ICL').AND.
     &                 (cell_type_index(neighbor).eq.'NCL')
     &                 )then
                     ws=ws_IN
                  elseif((cell_type_index(molecule).eq.'NCL').AND.
     &                 (cell_type_index(neighbor).eq.'ICL')
     &                 )then
                     ws=ws_IN
                  elseif((cell_type_index(molecule).eq.'ICL').AND.
     &                 (cell_type_index(neighbor).eq.'TCL')
     &                 )then
                     ws=ws_IT
                  elseif((cell_type_index(molecule).eq.'TCL').AND.
     &                 (cell_type_index(neighbor).eq.'ICL')
     &                 )then
                     ws=ws_IT    
                  elseif((cell_type_index(molecule).eq.'NCL').AND.
     &                 (cell_type_index(neighbor).eq.'TCL')
     &                 )then
                     ws=ws_NT
                  elseif((cell_type_index(molecule).eq.'TCL').AND.
     &                 (cell_type_index(neighbor).eq.'NCL')
     &                 )then
                     ws=ws_NT         					 
                  elseif((cell_type_index(molecule).eq.'TCL').AND.
     &                 (cell_type_index(neighbor).eq.'TCL')
     &                 )then
                     ws=ws_TT
                  elseif((cell_type_index(molecule).eq.'ICL').AND.
     &                 (cell_type_index(neighbor).eq.'ICL')
     &                 )then
                     ws=ws_II
                  elseif((cell_type_index(molecule).eq.'NCL').AND.
     &                 (cell_type_index(neighbor).eq.'NCL')
     &                 )then
                     ws=ws_NN					 
                  endif
                  f_atr_x(molecule)=f_atr_x(molecule)+
     &                 ws*pai*(x(neighbor)-x(molecule))
                  f_atr_y(molecule)=f_atr_y(molecule)+
     &                 ws*pai*(y(neighbor)-y(molecule))
                  f_atr_z(molecule)=f_atr_z(molecule)+
     &                 ws*pai*(z(neighbor)-z(molecule))
                  f_atr_x(neighbor)=f_atr_x(neighbor)-
     &                 ws*pai*(x(neighbor)-x(molecule))
                  f_atr_y(neighbor)=f_atr_y(neighbor)-
     &                 ws*pai*(y(neighbor)-y(molecule))
                  f_atr_z(neighbor)=f_atr_z(neighbor)-
     &                 ws*pai*(z(neighbor)-z(molecule))
                  repulsive_parameter=(3*sqrt(2*rad)/8)*
     &                 ((1-pois_ratio**2)/Young)*
     &                 ((sqrt(2*rad-dij))**3)
                  f_rps_x(molecule)=f_rps_x(molecule)+
     &                 (repulsive_parameter/dij)*
     &                 (x(molecule)-x(neighbor))
                  f_rps_y(molecule)=f_rps_y(molecule)+
     &                 (repulsive_parameter/dij)*
     &                 (y(molecule)-y(neighbor))
                  f_rps_z(molecule)=f_rps_z(molecule)+
     &                 (repulsive_parameter/dij)*
     &                 (z(molecule)-z(neighbor))
                  f_rps_x(neighbor)=f_rps_x(neighbor)-
     &                 (repulsive_parameter/dij)*
     &                 (x(molecule)-x(neighbor))
                  f_rps_y(neighbor)=f_rps_y(neighbor)-
     &                 (repulsive_parameter/dij)*
     &                 (y(molecule)-y(neighbor))
                  f_rps_z(neighbor)=f_rps_z(neighbor)-
     &                 (repulsive_parameter/dij)*
     &                 (z(molecule)-z(neighbor))
               endif
            enddo
         enddo
c>> calculate the stochastic force of each cel
         do molecule=1,npart
            dir_x=2*rand3(r3)-1
            dir_y=2*rand5(r5)-1
            dir_z=2*rand4(r4)-1
            amp_dir=sqrt(dir_x**2+dir_y**2+dir_z**2)
            scaled_dir_x=dir_x/amp_dir
            scaled_dir_y=dir_y/amp_dir
            scaled_dir_z=dir_z/amp_dir
            f_rdm_x(molecule)=2*sqrt(diffu_const)*
     &           cm_frict_const*scaled_dir_x
            f_rdm_y(molecule)=2*sqrt(diffu_const)*
     &           cm_frict_const*scaled_dir_y
            f_rdm_z(molecule)=2*sqrt(diffu_const)*
     &           cm_frict_const*scaled_dir_z
         enddo
c>>> calculate the velocity of each cell

         do molecule=1,npart
            vx_old(molecule)=0
            vy_old(molecule)=0
            vz_old(molecule)=0
            do neighbor=1,npart
               ax(molecule,neighbor)=0
               ay(molecule,neighbor)=0
               az(molecule,neighbor)=0
            enddo
            bx(molecule)=0
            by(molecule)=0
            bz(molecule)=0
         enddo

         do molecule=1,npart
            ax(molecule,molecule)=
     &           ax(molecule,molecule)+cm_frict_const
            ay(molecule,molecule)=
     &           ay(molecule,molecule)+cm_frict_const
            az(molecule,molecule)=
     &           az(molecule,molecule)+cm_frict_const
            do neighbor=1,npart
               if(molecule.ne.neighbor)then
                  dij=sqrt((x(molecule)-x(neighbor))**2+
     &                 (y(molecule)-y(neighbor))**2+
     &                 (z(molecule)-z(neighbor))**2)
                  if(dij.lt.2*rad)then
                     ax(molecule,molecule)=
     &                    ax(molecule,molecule)+
     &                    cc_frict_const
                     ax(molecule,neighbor)=
     &                    ax(molecule,neighbor)-
     &                    cc_frict_const
                     ay(molecule,molecule)=
     &                    ay(molecule,molecule)+
     &                    cc_frict_const
                     ay(molecule,neighbor)=
     &                    ay(molecule,neighbor)-
     &                    cc_frict_const
                     az(molecule,molecule)=
     &                    az(molecule,molecule)+
     &                    cc_frict_const
                     az(molecule,neighbor)=
     &                    az(molecule,neighbor)-
     &                    cc_frict_const
                  endif
               endif
            enddo
         enddo
            
         do molecule=1,npart
            bx(molecule)=(f_atr_x(molecule)+
     &           f_rps_x(molecule)+f_rdm_x(molecule))
            by(molecule)=(f_atr_y(molecule)+
     &           f_rps_y(molecule)+f_rdm_y(molecule))
            bz(molecule)=(f_atr_z(molecule)+
     &           f_rps_z(molecule)+f_rdm_z(molecule))
         enddo
         
         do i=1,niter

            do molecule=1,npart
               vx_new(molecule)=bx(molecule)
               vy_new(molecule)=by(molecule)
               vz_new(molecule)=bz(molecule)
               do neighbor=1,npart
                  if(molecule.ne.neighbor)then
                     vx_new(molecule)=vx_new(molecule)-
     &                    ax(molecule,neighbor)*vx_old(neighbor)
                     vy_new(molecule)=vy_new(molecule)-
     &                    ay(molecule,neighbor)*vy_old(neighbor)
                     vz_new(molecule)=vz_new(molecule)-
     &                    az(molecule,neighbor)*vz_old(neighbor)
                  endif
               enddo
               vx_new(molecule)=vx_new(molecule)/ax(molecule,molecule)
               vy_new(molecule)=vy_new(molecule)/ay(molecule,molecule)
               vz_new(molecule)=vz_new(molecule)/az(molecule,molecule)
            enddo
            
            RMSD=0
            do molecule=1,npart
               RMSD=RMSD+(vx_new(molecule)-vx_old(molecule))**2+
     &              (vy_new(molecule)-vy_old(molecule))**2+
     &              (vz_new(molecule)-vz_old(molecule))**2
            enddo
            RMSD=sqrt(RMSD/real(npart))
c            print*,i,RMSD
            if(RMSD.lt.iter_tol)then
               goto 100
            endif
            do molecule=1,npart
               vx_old(molecule)=vx_new(molecule)
               vy_old(molecule)=vy_new(molecule)
               vz_old(molecule)=vz_new(molecule)
            enddo

         enddo
         
 100     continue
         do molecule=1,npart
            vx(molecule)=vx_new(molecule)
            vy(molecule)=vy_new(molecule)
            vz(molecule)=vz_new(molecule)
c            print*,vx(molecule),vy(molecule),vz(molecule)
         enddo
c>>>>  update the position of each cell
         do j=1,npart
            x(j)=x(j)+vx(j)*dt
            y(j)=y(j)+vy(j)*dt
            z(j)=z(j)+vz(j)*dt    
         enddo
         do j=1,npart
            x(j)=x(j)-l*anint(x(j)/l)
            y(j)=y(j)-l*anint(y(j)/l)
            z(j)=z(j)-l*anint(z(j)/l)
         enddo
		 
c>>>>> output data

         if(mod(itime,100000).eq.0)then
            open (unit=10,file=
     &           'TCE_MultiCellDynamic_traj_output_23.pdb',
     &           status='unknown',access='append')
            write(10,2101) itime
            do j=1,npart
               if(cell_type_index(j).eq.'ICL')then
                  cell_res_index(j)='ALA'
               elseif(cell_type_index(j).eq.'NCL')then
                  cell_res_index(j)='VAL'
               elseif(cell_type_index(j).eq.'TCL')then
                  cell_res_index(j)='ILE'
               endif
               write(10,2100) 'ATOM  ',j,' CA ',cell_res_index(j), 
     &              j,x(j),y(j),z(j)
            enddo
            write(10,2102) 'END'
            close(10)
         endif
         
         if(mod(itime,1000).eq.0)then
            IN_num=0
            IT_num=0
            NT_num=0
            TT_num=0
            II_num=0
            NN_num=0			
            do molecule=1,npart-1
               do neighbor=molecule+1,npart
                  dij=sqrt((x(molecule)-x(neighbor))**2+
     &                 (y(molecule)-y(neighbor))**2+
     &                 (z(molecule)-z(neighbor))**2)
                  if(dij.lt.2*rad)then
                     if((cell_type_index(molecule).eq.'ICL').AND.
     &                    (cell_type_index(neighbor).eq.'NCL')
     &                    )then
                        IN_num=IN_num+1
                     elseif((cell_type_index(molecule).eq.'NCL').AND.
     &                       (cell_type_index(neighbor).eq.'ICL')
     &                       )then
                        IN_num=IN_num+1		 
                     elseif((cell_type_index(molecule).eq.'ICL').AND.
     &                       (cell_type_index(neighbor).eq.'TCL')
     &                       )then
                        IT_num=IT_num+1
                     elseif((cell_type_index(molecule).eq.'TCL').AND.
     &                       (cell_type_index(neighbor).eq.'ICL')
     &                       )then
                        IT_num=IT_num+1
                     elseif((cell_type_index(molecule).eq.'NCL').AND.
     &                       (cell_type_index(neighbor).eq.'TCL')
     &                       )then
                        NT_num=NT_num+1
                     elseif((cell_type_index(molecule).eq.'TCL').AND.
     &                       (cell_type_index(neighbor).eq.'NCL')
     &                       )then
                        NT_num=NT_num+1        					 
                     elseif((cell_type_index(molecule).eq.'TCL').AND.
     &                       (cell_type_index(neighbor).eq.'TCL')
     &                       )then
                        TT_num=TT_num+1
                     elseif((cell_type_index(molecule).eq.'ICL').AND.
     &                       (cell_type_index(neighbor).eq.'ICL')
     &                       )then
                        II_num=II_num+1
                     elseif((cell_type_index(molecule).eq.'NCL').AND.
     &                       (cell_type_index(neighbor).eq.'NCL')
     &                       )then
                        NN_num=NN_num+1					 
                     endif						
                     
                  endif 
               enddo
            enddo
            open (unit=10,file=
     &           'TCE_MultiCellDynamic_red_output_23.dat',
     &           status='unknown',access='append')
            write(10,2103) itime,IN_num,IT_num,NT_num,
     &           TT_num,II_num,NN_num
            close(10)
         endif 
         
 2102    format(A3)   
 2101    format(I10)   
 2103    format(I10,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5) 
 2100    format(A6,I5,1x,A4,1x,A3,2X,I4,4x,3F8.3)
cccccccccccccccccccccccccccccccccccccccccc

      enddo

ccccccccccccccccccccccc

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand3(r3)
      double precision s,u,v,r3
      s=65536.0
      u=2053.0
      v=13849.0
      m=r3/s
      r3=r3-m*s
      r3=u*r3+v
      m=r3/s
      r3=r3-m*s
      rand3=r3/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand4(r4)
      double precision s,u,v,r4
      s=65536.0
      u=2053.0
      v=13849.0
      m=r4/s
      r4=r4-m*s
      r4=u*r4+v
      m=r4/s
      r4=r4-m*s
      rand4=r4/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc                      
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand5(r5)
      double precision s,u,v,r5
      s=65536.0
      u=2053.0
      v=13849.0
      m=r5/s
      r5=r5-m*s
      r5=u*r5+v
      m=r5/s
      r5=r5-m*s
      rand5=r5/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc                      
