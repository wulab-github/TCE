      program RBKMC_OneSpecies3D2D
      implicit none

      integer simu_step
      parameter (simu_step=5000000)
      real*8 reduce_phase
      parameter (reduce_phase=0.5)
      real*8 time_step
      parameter (time_step=0.1)
      real*8 distance_step
      parameter (distance_step=0.5)
      real*8 cell_range_x
      parameter (cell_range_x=500.0)
      real*8 cell_range_y
      parameter (cell_range_y=500.0)
      real*8 cell_range_z,cell_range_z_initial,cell_range_z_final
      parameter (cell_range_z_initial=70.0,cell_range_z_final=35.0)
c>>   the z-coordinate of 2D membrane surface is 0, the extracellular space is from 0 to cell_range_z
c>>   Periodic boundary is applied in x and y axis, but not z, the z coordinates of 3D molecules are not allowed to be smaller than 0 and larger than cell_range_z 
      real*8 TCR_radius
      parameter (TCR_radius=5.0)
      real*8 TMR_radius
      parameter (TMR_radius=5.0)	  
      real*8 TCR_Z_length
      parameter (TCR_Z_length=10.0)
      real*8 TMR_Z_length
      parameter (TMR_Z_length=10.0)
      real*8 TCE_linker_length
      parameter (TCE_linker_length=5.0)
      real*8 TCR_D
      parameter (TCR_D=0.25)
      real*8 TCR_rot_D
      parameter (TCR_rot_D=10.0)
      real*8 TMR_D
      parameter (TMR_D=0.25)
      real*8 TMR_rot_D
      parameter (TMR_rot_D=10.0)	  
      real*8 TCE_radius
      parameter (TCE_radius=5.0)
      real*8 TCE_D
      parameter (TCE_D=0.75)
      real*8 TCE_rot_D
      parameter (TCE_rot_D=30.0)
      real*8 complex_D
      parameter (complex_D=0.1)
      real*8 Complex_rot_D
      parameter (Complex_rot_D=2.0)
      real*8 Ass_Rate
      parameter (Ass_Rate=10.0)
      real*8 Binding_Affinity
      parameter (Binding_Affinity=-21.0)
      real*8 Diss_Rate
      parameter (Diss_Rate=Ass_Rate*DEXP(Binding_Affinity))
      integer TCR_tot_num
      parameter (TCR_tot_num=100)
      integer TCR_res_num
      parameter (TCR_res_num=4)
      integer TMR_tot_num
      parameter (TMR_tot_num=100)
      integer TMR_res_num
      parameter (TMR_res_num=4)	  
      integer TCE_tot_num
      parameter (TCE_tot_num=500)
      integer TCE_res_num
      parameter (TCE_res_num=8)
      integer num_trajec
      parameter (num_trajec=5)
      real*8 pai
      parameter (pai=3.1415926) 
      integer complex_nb_ctg_num
      parameter (complex_nb_ctg_num=2) ! the maximal number of subunits in complex can be TCR and TCE here, so it is 2.
      integer complex2_nb_ctg_num
      parameter (complex2_nb_ctg_num=2) ! the maximal number of subunits in complex can be TCR and TCE here, so it is 2.	  
      real*8 bond_dist_cutoff
      parameter (bond_dist_cutoff=2.0)
      real*8 bond_thetapd,bond_thetapd_cutoff
      parameter (bond_thetapd=180.0,bond_thetapd_cutoff=30.0)
      real*8 bond_thetaot,bond_thetaot_cutoff
      parameter (bond_thetaot=90.0,bond_thetaot_cutoff=30.0)

      real*8 TCR_x(TCR_tot_num,TCR_res_num)
      real*8 TCR_y(TCR_tot_num,TCR_res_num)
      real*8 TCR_z(TCR_tot_num,TCR_res_num)
      real*8 TCR_x_0(TCR_tot_num,TCR_res_num)
      real*8 TCR_y_0(TCR_tot_num,TCR_res_num)
      real*8 TCR_z_0(TCR_tot_num,TCR_res_num)
      integer TCR_status(TCR_tot_num)
      real*8 TCR_x_new(TCR_tot_num,TCR_res_num)
      real*8 TCR_y_new(TCR_tot_num,TCR_res_num)
      real*8 TCR_z_new(TCR_tot_num,TCR_res_num)
      real*8 TCR_x_new0(TCR_tot_num,TCR_res_num)
      real*8 TCR_y_new0(TCR_tot_num,TCR_res_num)
      real*8 TCR_z_new0(TCR_tot_num,TCR_res_num)
      integer TCR_status_new(TCR_tot_num)

      real*8 TMR_x(TMR_tot_num,TMR_res_num)
      real*8 TMR_y(TMR_tot_num,TMR_res_num)
      real*8 TMR_z(TMR_tot_num,TMR_res_num)
      real*8 TMR_x_0(TMR_tot_num,TMR_res_num)
      real*8 TMR_y_0(TMR_tot_num,TMR_res_num)
      real*8 TMR_z_0(TMR_tot_num,TMR_res_num)
      integer TMR_status(TMR_tot_num)
      real*8 TMR_x_new(TMR_tot_num,TMR_res_num)
      real*8 TMR_y_new(TMR_tot_num,TMR_res_num)
      real*8 TMR_z_new(TMR_tot_num,TMR_res_num)
      real*8 TMR_x_new0(TMR_tot_num,TMR_res_num)
      real*8 TMR_y_new0(TMR_tot_num,TMR_res_num)
      real*8 TMR_z_new0(TMR_tot_num,TMR_res_num)
      integer TMR_status_new(TMR_tot_num)	  
	  
      real*8 TCE_x(TCE_tot_num,TCE_res_num)
      real*8 TCE_y(TCE_tot_num,TCE_res_num)
      real*8 TCE_z(TCE_tot_num,TCE_res_num)
      real*8 TCE_x_0(TCE_tot_num,TCE_res_num)
      real*8 TCE_y_0(TCE_tot_num,TCE_res_num)
      real*8 TCE_z_0(TCE_tot_num,TCE_res_num)
      integer TCE_1_status(TCE_tot_num)
      integer TCE_2_status(TCE_tot_num)	  
      real*8 TCE_x_new(TCE_tot_num,TCE_res_num)
      real*8 TCE_y_new(TCE_tot_num,TCE_res_num)
      real*8 TCE_z_new(TCE_tot_num,TCE_res_num)
      real*8 TCE_x_new0(TCE_tot_num,TCE_res_num)
      real*8 TCE_y_new0(TCE_tot_num,TCE_res_num)
      real*8 TCE_z_new0(TCE_tot_num,TCE_res_num)
      integer TCE_1_status_new(TCE_tot_num)
      integer TCE_2_status_new(TCE_tot_num)	  
	  
      integer complex_num,complex_num_new ! Number of complex formed in simulation
      integer complex_nb_num(TCE_tot_num)
      integer complex_nb_num_new(TCE_tot_num) ! How many subunits in each complex. Here the max number is 2.
      integer complex_nb_ctg(TCE_tot_num,complex_nb_ctg_num)
      integer complex_nb_ctg_new(TCE_tot_num,complex_nb_ctg_num) ! e.g. for complex_nb_ctg(i,j), what is the molecule category of jth subunit of ith complex. (here 1 is TCR and 2 is TCE)
      integer complex_nb_ctg_idx(TCE_tot_num,complex_nb_ctg_num)
      integer complex_nb_ctg_idx_new(TCE_tot_num,complex_nb_ctg_num) ! e.g. for complex_nb_ctg(i,j), what is the molecule index of jth subunit of ith complex, given its category=complex_nb_ctg(i,j).

      integer complex2_num,complex2_num_new ! Number of complex formed in simulation
      integer complex2_nb_num(TCE_tot_num)
      integer complex2_nb_num_new(TCE_tot_num) ! How many subunits in each complex. Here the max number is 2.
      integer complex2_nb_ctg(TCE_tot_num,complex2_nb_ctg_num)
      integer complex2_nb_ctg_new(TCE_tot_num,complex2_nb_ctg_num) ! e.g. for complex_nb_ctg(i,j), what is the molecule category of jth subunit of ith complex. (here 1 is TCR and 2 is TCE)
      integer complex2_nb_ctg_idx(TCE_tot_num,complex2_nb_ctg_num)
      integer complex2_nb_ctg_idx_new(TCE_tot_num,complex2_nb_ctg_num) ! e.g. for complex_nb_ctg(i,j), what is the molecule index of jth subunit of ith complex, given its category=complex_nb_ctg(i,j).

      integer i,j,k,n_t
      real*8 temp_i,temp_j,temp_k
      real*8 theta,phi,psai,phai
      real*8 t(3,3)
      real*8 dist
      real*8 cm1_a_x
      real*8 cm1_a_y
      real*8 cm1_a_z
      real*8 CM_x,CM_y,CM_z
      real*8 PB_x,PB_y,PB_z
      real*8 initial_simu_time,current_simu_time
      integer mc_time_step
      integer iteration_mole_step
      integer selecting_mole_index
      integer TCR_index,TCE_index,TMR_index
      real*8 Prob_Diff
      integer collision_flag
      real*8 point_x(3),point_y(3),point_z(3)
      real*8 theta_pd,theta_ot
      real*8 Prob_Ass,Prob_Diss
      integer selected_complex,selected_TCR,selected_TCE,selected_TMR
      real*8 rot_angle
      real*8 x_axis_bg,y_axis_bg,z_axis_bg
      real*8 x_axis_ed,y_axis_ed,z_axis_ed
      real*8 x_bf_rot,y_bf_rot,z_bf_rot
      real*8 x_af_rot,y_af_rot,z_af_rot
      integer upper_bound_flag,lower_bound_flag
      real*8 prob
      real*8 temp
      integer part1,part2,part3
      integer part4,part5,part6
      integer part7,part8,part9
      integer xlink_num
      real*8 xlink_ave

      real rand3
      double precision r3
      real rand4
      double precision r4
      real rand5
      double precision r5
	  
      r3=5.0      
      r4=5.0   
      r5=5.0   
	  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc   generate multiple trajectories
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do n_t=1,num_trajec


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c   construct the initialized position and conformation of molecules in 3D and 2D
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   construct the position and conformation of TCR on 2D

c>>>>>   random position

         do i=1,TCR_tot_num
 100        continue
            temp_i=rand3(r3)*cell_range_x-cell_range_x/2  
            temp_j=rand3(r3)*cell_range_y-cell_range_y/2             
            do j=1,i-1
               dist=sqrt((temp_i-TCR_x(j,1))**2+
     &              (temp_j-TCR_y(j,1))**2)
               if(dist.le.
     &              (TCR_radius+TCR_radius))then
                  goto 100
               endif
            enddo
            TCR_x(i,1)=temp_i
            TCR_y(i,1)=temp_j
            TCR_z(i,1)=0
            TCR_x_0(i,2)=temp_i+TCR_radius
            TCR_y_0(i,2)=temp_j
            TCR_z_0(i,2)=0
            TCR_x_0(i,3)=temp_i
            TCR_y_0(i,3)=temp_j+TCR_radius
            TCR_z_0(i,3)=0
            TCR_x(i,4)=temp_i
            TCR_y(i,4)=temp_j
            TCR_z(i,4)=TCR_Z_length
            TCR_status(i)=0

c>>>>>   random orientation along membrane normal

            rot_angle=2*rand3(r3)*180-180
            x_axis_bg=TCR_x(i,1)
            y_axis_bg=TCR_y(i,1)
            z_axis_bg=-100.0
            x_axis_ed=TCR_x(i,TCR_res_num)
            y_axis_ed=TCR_y(i,TCR_res_num)
            z_axis_ed=100.0
            do j=2,TCR_res_num-1
               x_bf_rot=TCR_x_0(i,j)
               y_bf_rot=TCR_y_0(i,j)
               z_bf_rot=TCR_z_0(i,j)
               x_af_rot=0
               y_af_rot=0
               z_af_rot=0
               call rot_along_axis(rot_angle,
     &              x_axis_bg,y_axis_bg,z_axis_bg,
     &              x_axis_ed,y_axis_ed,z_axis_ed,
     &              x_bf_rot,y_bf_rot,z_bf_rot,
     &              x_af_rot,y_af_rot,z_af_rot)
               TCR_x(i,j)=x_af_rot
               TCR_y(i,j)=y_af_rot
               TCR_z(i,j)=z_af_rot
            enddo

         enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   construct the position and conformation of TMR on 2D

c>>>>>   random position

         do i=1,TMR_tot_num
 200        continue
            temp_i=rand3(r3)*cell_range_x-cell_range_x/2  
            temp_j=rand3(r3)*cell_range_y-cell_range_y/2             
            do j=1,i-1
               dist=sqrt((temp_i-TMR_x(j,1))**2+
     &              (temp_j-TMR_y(j,1))**2)
               if(dist.le.
     &              (TMR_radius+TMR_radius))then
                  goto 200
               endif
            enddo
            TMR_x_0(i,1)=temp_i
            TMR_y_0(i,1)=temp_j
            TMR_z_0(i,1)=cell_range_z_initial
            TMR_x_0(i,2)=temp_i+TMR_radius
            TMR_y_0(i,2)=temp_j
            TMR_z_0(i,2)=cell_range_z_initial
            TMR_x_0(i,3)=temp_i
            TMR_y_0(i,3)=temp_j+TMR_radius
            TMR_z_0(i,3)=cell_range_z_initial
            TMR_x_0(i,4)=temp_i
            TMR_y_0(i,4)=temp_j
            TMR_z_0(i,4)=cell_range_z_initial-TMR_Z_length
            TMR_status(i)=0

c>>>>>   random orientation along membrane normal

            rot_angle=2*rand3(r3)*180-180
            x_axis_bg=TMR_x_0(i,1)
            y_axis_bg=TMR_y_0(i,1)
            z_axis_bg=-100.0
            x_axis_ed=TMR_x_0(i,TMR_res_num)
            y_axis_ed=TMR_y_0(i,TMR_res_num)
            z_axis_ed=100.0
            do j=2,TMR_res_num-1
               x_bf_rot=TMR_x_0(i,j)
               y_bf_rot=TMR_y_0(i,j)
               z_bf_rot=TMR_z_0(i,j)
               x_af_rot=0
               y_af_rot=0
               z_af_rot=0
               call rot_along_axis(rot_angle,
     &              x_axis_bg,y_axis_bg,z_axis_bg,
     &              x_axis_ed,y_axis_ed,z_axis_ed,
     &              x_bf_rot,y_bf_rot,z_bf_rot,
     &              x_af_rot,y_af_rot,z_af_rot)
               TMR_x(i,j)=x_af_rot
               TMR_y(i,j)=y_af_rot
               TMR_z(i,j)=z_af_rot
            enddo
            TMR_x(i,1)=TMR_x_0(i,1)
            TMR_y(i,1)=TMR_y_0(i,1)
            TMR_z(i,1)=TMR_z_0(i,1)
            TMR_x(i,TMR_res_num)=TMR_x_0(i,TMR_res_num)
            TMR_y(i,TMR_res_num)=TMR_y_0(i,TMR_res_num)
            TMR_z(i,TMR_res_num)=TMR_z_0(i,TMR_res_num)


         enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   construct the position and conformation of TCE

c>>>>>   random position

         do i=1,TCE_tot_num
 300        continue
            temp_i=rand3(r3)*cell_range_x-cell_range_x/2  
            temp_j=rand3(r3)*cell_range_y-cell_range_y/2  
            temp_k=rand3(r3)*cell_range_z_initial           
            do j=1,TCR_tot_num
               dist=sqrt((temp_i-TCR_x(j,1))**2+
     &              (temp_j-TCR_y(j,1))**2+
     &              (temp_k-TCR_z(j,1))**2)
               if(dist.le.
     &              (TCE_radius+TCR_radius))then
                  goto 300
               endif
            enddo
            do j=1,TMR_tot_num
               dist=sqrt((temp_i-TMR_x(j,1))**2+
     &              (temp_j-TMR_y(j,1))**2+
     &              (temp_k-TMR_z(j,1))**2)
               if(dist.le.
     &              (TCE_radius+TMR_radius))then
                  goto 300
               endif
            enddo
            do j=1,i-1
               dist=sqrt((temp_i-TCE_x(j,1))**2+
     &              (temp_j-TCE_y(j,1))**2+
     &              (temp_k-TCE_z(j,1))**2)
               if(dist.le.
     &              (TCE_radius+TCE_radius))then
                  goto 300
               endif
            enddo
            TCE_x_0(i,1)=temp_i
            TCE_y_0(i,1)=temp_j
            TCE_z_0(i,1)=temp_k+TCE_linker_length*0.5
            TCE_x_0(i,2)=temp_i+TCE_radius
            TCE_y_0(i,2)=temp_j
            TCE_z_0(i,2)=temp_k+TCE_linker_length*0.5
            TCE_x_0(i,3)=temp_i
            TCE_y_0(i,3)=temp_j+TCE_radius
            TCE_z_0(i,3)=temp_k+TCE_linker_length*0.5
            TCE_x_0(i,4)=temp_i
            TCE_y_0(i,4)=temp_j
            TCE_z_0(i,4)=temp_k+TCE_radius+TCE_linker_length*0.5

            TCE_x_0(i,5)=temp_i
            TCE_y_0(i,5)=temp_j
            TCE_z_0(i,5)=temp_k-TCE_linker_length*0.5
            TCE_x_0(i,6)=temp_i+TCE_radius
            TCE_y_0(i,6)=temp_j
            TCE_z_0(i,6)=temp_k-TCE_linker_length*0.5
            TCE_x_0(i,7)=temp_i
            TCE_y_0(i,7)=temp_j+TCE_radius
            TCE_z_0(i,7)=temp_k-TCE_linker_length*0.5
            TCE_x_0(i,8)=temp_i
            TCE_y_0(i,8)=temp_j
            TCE_z_0(i,8)=temp_k-TCE_radius-TCE_linker_length*0.5

            TCE_1_status(i)=0
            TCE_2_status(i)=0

c>>>>>   random orientation

            theta=(2*rand3(r3)-1)*pai
            phi=(2*rand3(r3)-1)*pai
            psai=(2*rand3(r3)-1)*pai
            
            t(1,1)=cos(psai)*cos(phi)-cos(theta)*sin(phi)*sin(psai)
            t(1,2)=-sin(psai)*cos(phi)-cos(theta)*sin(phi)*cos(psai)
            t(1,3)=sin(theta)*sin(phi)
            
            t(2,1)=cos(psai)*sin(phi)+cos(theta)*cos(phi)*sin(psai)
            t(2,2)=-sin(psai)*sin(phi)+cos(theta)*cos(phi)*cos(psai)
            t(2,3)=-sin(theta)*cos(phi)
            
            t(3,1)=sin(psai)*sin(theta)
            t(3,2)=cos(psai)*sin(theta)
            t(3,3)=cos(theta)

            CM_x=(TCE_x_0(i,1)+TCE_x_0(i,5))*0.5
            CM_y=(TCE_y_0(i,1)+TCE_y_0(i,5))*0.5
            CM_z=(TCE_z_0(i,1)+TCE_z_0(i,5))*0.5
            
            do k=1,TCE_res_num
               TCE_x(i,k)=t(1,1)*(TCE_x_0(i,k)-CM_x)+
     &              t(1,2)*(TCE_y_0(i,k)-CM_y)
     &              +t(1,3)*(TCE_z_0(i,k)-CM_z)
     &              +CM_x
               TCE_y(i,k)=t(2,1)*(TCE_x_0(i,k)-CM_x)+
     &              t(2,2)*(TCE_y_0(i,k)-CM_y)
     &              +t(2,3)*(TCE_z_0(i,k)-CM_z)
     &              +CM_y
               TCE_z(i,k)=t(3,1)*(TCE_x_0(i,k)-CM_x)+
     &              t(3,2)*(TCE_y_0(i,k)-CM_y)
     &              +t(3,3)*(TCE_z_0(i,k)-CM_z)
     &              +CM_z

            enddo

         enddo

         complex_num=0
         do i=1,TCE_tot_num
            complex_nb_num(i)=0
            do j=1,complex_nb_ctg_num
               complex_nb_ctg(i,j)=0
               complex_nb_ctg_idx(i,j)=0
            enddo
         enddo
         complex2_num=0
         do i=1,TCE_tot_num
            complex2_nb_num(i)=0
            do j=1,complex2_nb_ctg_num
               complex2_nb_ctg(i,j)=0
               complex2_nb_ctg_idx(i,j)=0
            enddo
         enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      begin  main  loop of Diffusion-Reaction simulation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         initial_simu_time=0.0

         do mc_time_step=1,simu_step

            current_simu_time=current_simu_time+time_step

            if(mc_time_step.le.int(real(simu_step)*reduce_phase))then
               cell_range_z=cell_range_z_initial-
     &              ((cell_range_z_initial-cell_range_z_final)
     &              /(real(simu_step)*reduce_phase))
     &              *real(mc_time_step)
            else
               cell_range_z=cell_range_z_final
            endif



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   perform the diffusion for molecules
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

            do i=1,TCR_tot_num
               do j=1,TCR_res_num
                  TCR_x_new(i,j)=TCR_x(i,j)
                  TCR_y_new(i,j)=TCR_y(i,j)
                  TCR_z_new(i,j)=TCR_z(i,j)
               enddo
               TCR_status_new(i)=TCR_status(i)
            enddo
            do i=1,TMR_tot_num
               do j=1,TMR_res_num
                  TMR_x_new(i,j)=TMR_x(i,j)
                  TMR_y_new(i,j)=TMR_y(i,j)
                  TMR_z_new(i,j)=TMR_z(i,j)-TMR_z(i,1)+cell_range_z
               enddo
               TMR_status_new(i)=TMR_status(i)
            enddo
            do i=1,TCE_tot_num
               do j=1,TCE_res_num
                  TCE_x_new(i,j)=TCE_x(i,j)
                  TCE_y_new(i,j)=TCE_y(i,j)
                  if((TCE_2_status(i).eq.1).and.
     &                 (TCE_1_status(i).eq.0))then	 
                     do k=1,complex2_num
                        if(complex2_nb_ctg_idx(k,2)
     &                       .eq.i)then
                           TMR_index=
     &                          complex2_nb_ctg_idx(k,1)
                        endif
                     enddo					 
                     TCE_z_new(i,j)=TCE_z(i,j)-
     &                    TMR_z(TMR_index,1)+cell_range_z
                  else
                     TCE_z_new(i,j)=TCE_z(i,j)
                  endif
               enddo
               TCE_1_status_new(i)=TCE_1_status(i)
               TCE_2_status_new(i)=TCE_2_status(i)
            enddo
            complex_num_new=complex_num
            do i=1,complex_num
               complex_nb_num_new(i)=complex_nb_num(i)
               do j=1,complex_nb_ctg_num
                  complex_nb_ctg_new(i,j)=complex_nb_ctg(i,j)
                  complex_nb_ctg_idx_new(i,j)=complex_nb_ctg_idx(i,j)
               enddo
            enddo
            complex2_num_new=complex2_num
            do i=1,complex2_num
               complex2_nb_num_new(i)=complex2_nb_num(i)
               do j=1,complex2_nb_ctg_num
                  complex2_nb_ctg_new(i,j)=complex2_nb_ctg(i,j)
                  complex2_nb_ctg_idx_new(i,j)=complex2_nb_ctg_idx(i,j)
               enddo
            enddo

            do iteration_mole_step=1,
     &           TCR_tot_num+TMR_tot_num+TCE_tot_num

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   randomly select one molecule
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               selecting_mole_index=int(rand3(r3)*(TCR_tot_num+
     &              TMR_tot_num+TCE_tot_num))+1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   if the selected molecule is TCR

               if(selecting_mole_index.le.TCR_tot_num)then


                  TCR_index=int(rand3(r3)*TCR_tot_num)+1
cc>>>>   make random diffusion for selected molecule
cc>>>>   if it is unbound
                  if(TCR_status_new(TCR_index).eq.0)then
                     Prob_Diff=(6*TCR_D*time_step)/(distance_step**2)
                     if(rand3(r3).lt.Prob_Diff)then
cc>>     if smaller than the diffusion probability, move alone
                        phai=rand3(r3)*2*pai
                        do i=1,TCR_res_num
                           TCR_x_new0(TCR_index,i)=
     &                          TCR_x(TCR_index,i)+
     &                          distance_step*cos(phai)
                           TCR_y_new0(TCR_index,i)=
     &                          TCR_y(TCR_index,i)+
     &                          distance_step*sin(phai)
                           TCR_z_new0(TCR_index,i)=
     &                          TCR_z(TCR_index,i)
                        enddo
                        PB_x=cell_range_x*anint
     &                          (TCR_x_new0(TCR_index,1)/cell_range_x)
                        PB_y=cell_range_y*anint
     &                          (TCR_y_new0(TCR_index,1)/cell_range_y)
                        do i=1,TCR_res_num
                           TCR_x_new0(TCR_index,i)=
     &                          TCR_x_new0(TCR_index,i)-PB_x     
                           TCR_y_new0(TCR_index,i)=
     &                          TCR_y_new0(TCR_index,i)-PB_y
                        enddo
                     else
                        do i=1,TCR_res_num
                           TCR_x_new0(TCR_index,i)=
     &                          TCR_x(TCR_index,i)
                           TCR_y_new0(TCR_index,i)=
     &                          TCR_y(TCR_index,i)
                           TCR_z_new0(TCR_index,i)=
     &                          TCR_z(TCR_index,i)
                        enddo
                     endif
ccc>>>>   Check Collisions
                     collision_flag=0
                     do i=1,TCR_tot_num
                        if(TCR_index.ne.i)then
                           dist=sqrt((TCR_x_new(i,1)-
     &                          TCR_x_new0(TCR_index,1))**2+
     &                          (TCR_y_new(i,1)-
     &                          TCR_y_new0(TCR_index,1))**2+
     &                          (TCR_z_new(i,1)-
     &                          TCR_z_new0(TCR_index,1))**2)
                           if(dist.lt.(TCR_radius+TCR_radius))then
                              collision_flag=1
                           endif
                        endif
                     enddo
                     do i=1,TCE_tot_num
                        dist=sqrt((TCE_x_new(i,1)-
     &                       TCR_x_new0(TCR_index,4))**2+
     &                       (TCE_y_new(i,1)-
     &                       TCR_y_new0(TCR_index,4))**2+
     &                       (TCE_z_new(i,1)-
     &                       TCR_z_new0(TCR_index,4))**2)
                        if(dist.lt.TCE_radius)then
                           collision_flag=1
                        endif
                        dist=sqrt((TCE_x_new(i,5)-
     &                       TCR_x_new0(TCR_index,4))**2+
     &                       (TCE_y_new(i,5)-
     &                       TCR_y_new0(TCR_index,4))**2+
     &                       (TCE_z_new(i,5)-
     &                       TCR_z_new0(TCR_index,4))**2)
                        if(dist.lt.TCE_radius)then
                           collision_flag=1
                        endif						
                     enddo
                     if(collision_flag.eq.1)then
                        do i=1,TCR_res_num
                           TCR_x_new0(TCR_index,i)=
     &                          TCR_x(TCR_index,i)
                           TCR_y_new0(TCR_index,i)=
     &                          TCR_y(TCR_index,i)
                           TCR_z_new0(TCR_index,i)=
     &                          TCR_z(TCR_index,i)
                        enddo
                     endif
c>>>>>>>>> Rotate the molecule

                     TCR_x_new(TCR_index,1)=TCR_x_new0(TCR_index,1)
                     TCR_y_new(TCR_index,1)=TCR_y_new0(TCR_index,1)
                     TCR_z_new(TCR_index,1)=TCR_z_new0(TCR_index,1)

                     TCR_x_new(TCR_index,TCR_res_num)=
     &                    TCR_x_new0(TCR_index,TCR_res_num)
                     TCR_y_new(TCR_index,TCR_res_num)=
     &                    TCR_y_new0(TCR_index,TCR_res_num)
                     TCR_z_new(TCR_index,TCR_res_num)=
     &                    TCR_z_new0(TCR_index,TCR_res_num)

                     rot_angle=2*rand3(r3)*TCR_rot_D-TCR_rot_D
                     x_axis_bg=TCR_x_new(TCR_index,1)
                     y_axis_bg=TCR_y_new(TCR_index,1)
                     z_axis_bg=-100.0
                     x_axis_ed=TCR_x_new(TCR_index,TCR_res_num)
                     y_axis_ed=TCR_y_new(TCR_index,TCR_res_num)
                     z_axis_ed=100.0
                     do j=2,TCR_res_num-1
                        x_bf_rot=TCR_x_new0(TCR_index,j)
                        y_bf_rot=TCR_y_new0(TCR_index,j)
                        z_bf_rot=TCR_z_new0(TCR_index,j)
                        x_af_rot=0
                        y_af_rot=0
                        z_af_rot=0
                        call rot_along_axis(rot_angle,
     &                       x_axis_bg,y_axis_bg,z_axis_bg,
     &                       x_axis_ed,y_axis_ed,z_axis_ed,
     &                       x_bf_rot,y_bf_rot,z_bf_rot,
     &                       x_af_rot,y_af_rot,z_af_rot)
                        TCR_x_new(TCR_index,j)=x_af_rot
                        TCR_y_new(TCR_index,j)=y_af_rot
                        TCR_z_new(TCR_index,j)=z_af_rot
                     enddo

cc>>>>   elseif it is in the complex
                  elseif(TCR_status_new(TCR_index).eq.1)then
cc>>     find out the bound substrate
                     do j=1,complex_num_new
                        if(complex_nb_ctg_idx_new(j,1)
     &                       .eq.TCR_index)then
                           TCE_index=
     &                          complex_nb_ctg_idx_new(j,2)
                        endif
                     enddo
                     if(TCE_2_status_new(TCE_index).eq.0)then
                        Prob_Diff=(6*complex_D*time_step)/
     &                       (distance_step**2)
                        if(rand3(r3).lt.Prob_Diff)then
cc>>     if smaller than the diffusion probability, move the complex together
                           phai=rand3(r3)*2*pai
                           do i=1,TCR_res_num
                              TCR_x_new0(TCR_index,i)=
     &                             TCR_x(TCR_index,i)+
     &                             distance_step*cos(phai)
                              TCR_y_new0(TCR_index,i)=
     &                             TCR_y(TCR_index,i)+
     &                             distance_step*sin(phai)
                              TCR_z_new0(TCR_index,i)=
     &                             TCR_z(TCR_index,i)
                           enddo
                           do i=1,TCE_res_num
                              TCE_x_new0(TCE_index,i)=
     &                             TCE_x(TCE_index,i)+
     &                             distance_step*cos(phai)
                              TCE_y_new0(TCE_index,i)=
     &                             TCE_y(TCE_index,i)+
     &                             distance_step*sin(phai)
                              TCE_z_new0(TCE_index,i)=
     &                             TCE_z(TCE_index,i)
                           enddo
                           PB_x=cell_range_x*anint
     &                          (((TCR_x_new0(TCR_index,1)+
     &                          TCE_x_new0(TCE_index,1))/2)/
     &                          cell_range_x)
                           PB_y=cell_range_y*anint
     &                          (((TCR_y_new0(TCR_index,1)+
     &                          TCE_y_new0(TCE_index,1))/2)/
     &                          cell_range_y)
                           do i=1,TCR_res_num
                              TCR_x_new0(TCR_index,i)=
     &                             TCR_x_new0(TCR_index,i)-PB_x     
                              TCR_y_new0(TCR_index,i)=
     &                             TCR_y_new0(TCR_index,i)-PB_y     
                           enddo
                           do i=1,TCE_res_num
                              TCE_x_new0(TCE_index,i)=
     &                             TCE_x_new0(TCE_index,i)-PB_x     
                              TCE_y_new0(TCE_index,i)=
     &                             TCE_y_new0(TCE_index,i)-PB_y     
                           enddo
                        else
                           do i=1,TCR_res_num
                              TCR_x_new0(TCR_index,i)=
     &                             TCR_x(TCR_index,i)
                              TCR_y_new0(TCR_index,i)=
     &                             TCR_y(TCR_index,i)
                              TCR_z_new0(TCR_index,i)=
     &                             TCR_z(TCR_index,i)
                           enddo
                           do i=1,TCE_res_num
                              TCE_x_new0(TCE_index,i)=
     &                             TCE_x(TCE_index,i)
                              TCE_y_new0(TCE_index,i)=
     &                             TCE_y(TCE_index,i)
                              TCE_z_new0(TCE_index,i)=
     &                             TCE_z(TCE_index,i)
                           enddo
                        endif
ccc>>>>   Check Collisions
                        collision_flag=0
                        do i=1,TCR_tot_num
                           if(TCR_index.ne.i)then
                              dist=sqrt((TCR_x_new(i,1)-
     &                             TCR_x_new0(TCR_index,1))**2+
     &                             (TCR_y_new(i,1)-
     &                             TCR_y_new0(TCR_index,1))**2+
     &                             (TCR_z_new(i,1)-
     &                             TCR_z_new0(TCR_index,1))**2)
                              if(dist.lt.(TCR_radius+TCR_radius))then
                                 collision_flag=1
                              endif
                           endif
                        enddo
                        do i=1,TCE_tot_num
                           dist=sqrt((TCE_x_new(i,1)-
     &                          TCR_x_new0(TCR_index,4))**2+
     &                          (TCE_y_new(i,1)-
     &                          TCR_y_new0(TCR_index,4))**2+
     &                          (TCE_z_new(i,1)-
     &                          TCR_z_new0(TCR_index,4))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif
                           dist=sqrt((TCE_x_new(i,5)-
     &                          TCR_x_new0(TCR_index,4))**2+
     &                          (TCE_y_new(i,5)-
     &                          TCR_y_new0(TCR_index,4))**2+
     &                          (TCE_z_new(i,5)-
     &                          TCR_z_new0(TCR_index,4))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif						   
                        enddo
                        do i=1,TCR_tot_num
                           dist=sqrt((TCR_x_new(i,4)-
     &                          TCE_x_new0(TCE_index,1))**2+
     &                          (TCR_y_new(i,4)-
     &                          TCE_y_new0(TCE_index,1))**2+
     &                          (TCR_z_new(i,4)-
     &                          TCE_z_new0(TCE_index,1))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif
                           dist=sqrt((TCR_x_new(i,4)-
     &                          TCE_x_new0(TCE_index,5))**2+
     &                          (TCR_y_new(i,4)-
     &                          TCE_y_new0(TCE_index,5))**2+
     &                          (TCR_z_new(i,4)-
     &                          TCE_z_new0(TCE_index,5))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif						   
                        enddo
                        do i=1,TCE_tot_num
                           if(TCE_index.ne.i)then
                              dist=sqrt((TCE_x_new(i,1)-
     &                             TCE_x_new0(TCE_index,1))**2+
     &                             (TCE_y_new(i,1)-
     &                             TCE_y_new0(TCE_index,1))**2+
     &                             (TCE_z_new(i,1)-
     &                             TCE_z_new0(TCE_index,1))**2)
                              if(dist.lt.(TCE_radius+TCE_radius)
     &                             )then
                                 collision_flag=1
                              endif
                              dist=sqrt((TCE_x_new(i,5)-
     &                             TCE_x_new0(TCE_index,5))**2+
     &                             (TCE_y_new(i,5)-
     &                             TCE_y_new0(TCE_index,5))**2+
     &                             (TCE_z_new(i,5)-
     &                             TCE_z_new0(TCE_index,5))**2)
                              if(dist.lt.(TCE_radius+TCE_radius)
     &                             )then
                                 collision_flag=1
                              endif	
                              dist=sqrt((TCE_x_new(i,1)-
     &                             TCE_x_new0(TCE_index,5))**2+
     &                             (TCE_y_new(i,1)-
     &                             TCE_y_new0(TCE_index,5))**2+
     &                             (TCE_z_new(i,1)-
     &                             TCE_z_new0(TCE_index,5))**2)
                              if(dist.lt.(TCE_radius+TCE_radius)
     &                             )then
                                 collision_flag=1
                              endif
                              dist=sqrt((TCE_x_new(i,5)-
     &                             TCE_x_new0(TCE_index,1))**2+
     &                             (TCE_y_new(i,5)-
     &                             TCE_y_new0(TCE_index,1))**2+
     &                             (TCE_z_new(i,5)-
     &                             TCE_z_new0(TCE_index,1))**2)
                              if(dist.lt.(TCE_radius+TCE_radius)
     &                             )then
                                 collision_flag=1
                              endif							  
                           endif
                        enddo
                        do i=1,TMR_tot_num
                           dist=sqrt((TMR_x_new(i,4)-
     &                          TCE_x_new0(TCE_index,1))**2+
     &                          (TMR_y_new(i,4)-
     &                          TCE_y_new0(TCE_index,1))**2+
     &                          (TMR_z_new(i,4)-
     &                          TCE_z_new0(TCE_index,1))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif
                           dist=sqrt((TMR_x_new(i,4)-
     &                          TCE_x_new0(TCE_index,5))**2+
     &                          (TMR_y_new(i,4)-
     &                          TCE_y_new0(TCE_index,5))**2+
     &                          (TMR_z_new(i,4)-
     &                          TCE_z_new0(TCE_index,5))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif						   
                        enddo
                        if(collision_flag.eq.1)then
                           do i=1,TCR_res_num
                              TCR_x_new0(TCR_index,i)=
     &                             TCR_x(TCR_index,i)
                              TCR_y_new0(TCR_index,i)=
     &                             TCR_y(TCR_index,i)
                              TCR_z_new0(TCR_index,i)=
     &                             TCR_z(TCR_index,i)
                           enddo
                           do i=1,TCE_res_num
                              TCE_x_new0(TCE_index,i)=
     &                             TCE_x(TCE_index,i)
                              TCE_y_new0(TCE_index,i)=
     &                             TCE_y(TCE_index,i)
                              TCE_z_new0(TCE_index,i)=
     &                             TCE_z(TCE_index,i)
                           enddo
                        endif
c>>>>>>>>> Rotate the molecule
                       
                        cm1_a_x=0
                        cm1_a_y=0
                        cm1_a_z=0
                        do j=1,TCR_res_num
                           cm1_a_x=cm1_a_x+TCR_x_new0(TCR_index,j)
                           cm1_a_y=cm1_a_y+TCR_y_new0(TCR_index,j)
                           cm1_a_z=cm1_a_z+TCR_z_new0(TCR_index,j)
                        enddo
c                        do j=1,TCE_res_num
c                           cm1_a_x=cm1_a_x+TCE_x_new0(TCE_index,j)
c                           cm1_a_y=cm1_a_y+TCE_y_new0(TCE_index,j)
c                           cm1_a_z=cm1_a_z+TCE_z_new0(TCE_index,j)
c                        enddo
                        cm1_a_x=cm1_a_x/(real(TCR_res_num))
                        cm1_a_y=cm1_a_y/(real(TCR_res_num))
                        cm1_a_z=cm1_a_z/(real(TCR_res_num))
                        
                        rot_angle=2*rand3(r3)*Complex_rot_D-
     &                       Complex_rot_D
                        x_axis_bg=cm1_a_x
                        y_axis_bg=cm1_a_y
                        z_axis_bg=-100.0
                        x_axis_ed=cm1_a_x
                        y_axis_ed=cm1_a_y
                        z_axis_ed=100.0
                        do j=1,TCR_res_num
                           x_bf_rot=TCR_x_new0(TCR_index,j)
                           y_bf_rot=TCR_y_new0(TCR_index,j)
                           z_bf_rot=TCR_z_new0(TCR_index,j)
                           x_af_rot=0
                           y_af_rot=0
                           z_af_rot=0
                           call rot_along_axis(rot_angle,
     &                          x_axis_bg,y_axis_bg,z_axis_bg,
     &                          x_axis_ed,y_axis_ed,z_axis_ed,
     &                          x_bf_rot,y_bf_rot,z_bf_rot,
     &                          x_af_rot,y_af_rot,z_af_rot)
                           TCR_x_new(TCR_index,j)=x_af_rot
                           TCR_y_new(TCR_index,j)=y_af_rot
                           TCR_z_new(TCR_index,j)=z_af_rot
                        enddo
                        do j=1,TCE_res_num
                           x_bf_rot=TCE_x_new0(TCE_index,j)
                           y_bf_rot=TCE_y_new0(TCE_index,j)
                           z_bf_rot=TCE_z_new0(TCE_index,j)
                           x_af_rot=0
                           y_af_rot=0
                           z_af_rot=0
                           call rot_along_axis(rot_angle,
     &                          x_axis_bg,y_axis_bg,z_axis_bg,
     &                          x_axis_ed,y_axis_ed,z_axis_ed,
     &                          x_bf_rot,y_bf_rot,z_bf_rot,
     &                          x_af_rot,y_af_rot,z_af_rot)
                           TCE_x_new(TCE_index,j)=x_af_rot
                           TCE_y_new(TCE_index,j)=y_af_rot
                           TCE_z_new(TCE_index,j)=z_af_rot
                        enddo
                        
                     endif

                  endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   else if the selected molecule is TMR


               elseif((selecting_mole_index.gt.TCR_tot_num).AND.
     &                 (selecting_mole_index.le.(TCR_tot_num+
     &                 TMR_tot_num)))then


                  TMR_index=int(rand3(r3)*TMR_tot_num)+1
cc>>>>   make random diffusion for selected molecule
cc>>>>   if it is unbound
                  if(TMR_status_new(TMR_index).eq.0)then
                     Prob_Diff=(6*TMR_D*time_step)/(distance_step**2)
                     if(rand3(r3).lt.Prob_Diff)then
cc>>     if smaller than the diffusion probability, move alone
                        phai=rand3(r3)*2*pai
                        do i=1,TMR_res_num
                           TMR_x_new0(TMR_index,i)=
     &                          TMR_x(TMR_index,i)+
     &                          distance_step*cos(phai)
                           TMR_y_new0(TMR_index,i)=
     &                          TMR_y(TMR_index,i)+
     &                          distance_step*sin(phai)
                           TMR_z_new0(TMR_index,i)=
     &                          TMR_z(TMR_index,i)
                        enddo
                        PB_x=cell_range_x*anint
     &                          (TMR_x_new0(TMR_index,1)/cell_range_x)
                        PB_y=cell_range_y*anint
     &                          (TMR_y_new0(TMR_index,1)/cell_range_y)
                        do i=1,TMR_res_num
                           TMR_x_new0(TMR_index,i)=
     &                          TMR_x_new0(TMR_index,i)-PB_x     
                           TMR_y_new0(TMR_index,i)=
     &                          TMR_y_new0(TMR_index,i)-PB_y
                        enddo
                     else
                        do i=1,TMR_res_num
                           TMR_x_new0(TMR_index,i)=
     &                          TMR_x(TMR_index,i)
                           TMR_y_new0(TMR_index,i)=
     &                          TMR_y(TMR_index,i)
                           TMR_z_new0(TMR_index,i)=
     &                          TMR_z(TMR_index,i)
                        enddo
                     endif
ccc>>>>   Check Collisions
                     collision_flag=0
                     do i=1,TMR_tot_num
                        if(TMR_index.ne.i)then
                           dist=sqrt((TMR_x_new(i,1)-
     &                          TMR_x_new0(TMR_index,1))**2+
     &                          (TMR_y_new(i,1)-
     &                          TMR_y_new0(TMR_index,1))**2+
     &                          (TMR_z_new(i,1)-
     &                          TMR_z_new0(TMR_index,1))**2)
                           if(dist.lt.(TMR_radius+TMR_radius))then
                              collision_flag=1
                           endif
                        endif
                     enddo
                     do i=1,TCE_tot_num
                        dist=sqrt((TCE_x_new(i,1)-
     &                       TMR_x_new0(TMR_index,4))**2+
     &                       (TCE_y_new(i,1)-
     &                       TMR_y_new0(TMR_index,4))**2+
     &                       (TCE_z_new(i,1)-
     &                       TMR_z_new0(TMR_index,4))**2)
                        if(dist.lt.TCE_radius)then
                           collision_flag=1
                        endif
                        dist=sqrt((TCE_x_new(i,5)-
     &                       TMR_x_new0(TMR_index,4))**2+
     &                       (TCE_y_new(i,5)-
     &                       TMR_y_new0(TMR_index,4))**2+
     &                       (TCE_z_new(i,5)-
     &                       TMR_z_new0(TMR_index,4))**2)
                        if(dist.lt.TCE_radius)then
                           collision_flag=1
                        endif						
                     enddo
                     if(collision_flag.eq.1)then
                        do i=1,TCR_res_num
                           TMR_x_new0(TMR_index,i)=
     &                          TMR_x(TMR_index,i)
                           TMR_y_new0(TMR_index,i)=
     &                          TMR_y(TMR_index,i)
                           TMR_z_new0(TMR_index,i)=
     &                          TMR_z(TMR_index,i)
                        enddo
                     endif
c>>>>>>>>> Rotate the molecule

                     TMR_x_new(TMR_index,1)=TMR_x_new0(TMR_index,1)
                     TMR_y_new(TMR_index,1)=TMR_y_new0(TMR_index,1)
                     TMR_z_new(TMR_index,1)=TMR_z_new0(TMR_index,1)

                     TMR_x_new(TMR_index,TMR_res_num)=
     &                    TMR_x_new0(TMR_index,TMR_res_num)
                     TMR_y_new(TMR_index,TMR_res_num)=
     &                    TMR_y_new0(TMR_index,TMR_res_num)
                     TMR_z_new(TMR_index,TMR_res_num)=
     &                    TMR_z_new0(TMR_index,TMR_res_num)

                     rot_angle=2*rand3(r3)*TMR_rot_D-TMR_rot_D
                     x_axis_bg=TMR_x_new(TMR_index,1)
                     y_axis_bg=TMR_y_new(TMR_index,1)
                     z_axis_bg=-100.0
                     x_axis_ed=TMR_x_new(TMR_index,TMR_res_num)
                     y_axis_ed=TMR_y_new(TMR_index,TMR_res_num)
                     z_axis_ed=100.0
                     do j=2,TMR_res_num-1
                        x_bf_rot=TMR_x_new0(TMR_index,j)
                        y_bf_rot=TMR_y_new0(TMR_index,j)
                        z_bf_rot=TMR_z_new0(TMR_index,j)
                        x_af_rot=0
                        y_af_rot=0
                        z_af_rot=0
                        call rot_along_axis(rot_angle,
     &                       x_axis_bg,y_axis_bg,z_axis_bg,
     &                       x_axis_ed,y_axis_ed,z_axis_ed,
     &                       x_bf_rot,y_bf_rot,z_bf_rot,
     &                       x_af_rot,y_af_rot,z_af_rot)
                        TMR_x_new(TMR_index,j)=x_af_rot
                        TMR_y_new(TMR_index,j)=y_af_rot
                        TMR_z_new(TMR_index,j)=z_af_rot
                     enddo

cc>>>>   elseif it is in the complex
                  elseif(TMR_status_new(TMR_index).eq.1)then
cc>>     find out the bound substrate
                     do j=1,complex2_num_new
                        if(complex2_nb_ctg_idx_new(j,1)
     &                       .eq.TMR_index)then
                           TCE_index=
     &                          complex2_nb_ctg_idx_new(j,2)
                        endif
                     enddo
                    
                     if(TCE_1_status_new(TCE_index).eq.0)then
                        Prob_Diff=(6*complex_D*time_step)/
     &                       (distance_step**2)
                        if(rand3(r3).lt.Prob_Diff)then
cc>>     if smaller than the diffusion probability, move the complex together
                           phai=rand3(r3)*2*pai
                           do i=1,TMR_res_num
                              TMR_x_new0(TMR_index,i)=
     &                             TMR_x(TMR_index,i)+
     &                             distance_step*cos(phai)
                              TMR_y_new0(TMR_index,i)=
     &                             TMR_y(TMR_index,i)+
     &                             distance_step*sin(phai)
                              TMR_z_new0(TMR_index,i)=
     &                             TMR_z(TMR_index,i)
                           enddo
                           do i=1,TCE_res_num
                              TCE_x_new0(TCE_index,i)=
     &                             TCE_x(TCE_index,i)+
     &                             distance_step*cos(phai)
                              TCE_y_new0(TCE_index,i)=
     &                             TCE_y(TCE_index,i)+
     &                             distance_step*sin(phai)
                              TCE_z_new0(TCE_index,i)=
     &                             TCE_z(TCE_index,i)
                           enddo
                           PB_x=cell_range_x*anint
     &                          (((TMR_x_new0(TMR_index,1)+
     &                          TCE_x_new0(TCE_index,1))/2)/
     &                          cell_range_x)
                           PB_y=cell_range_y*anint
     &                          (((TMR_y_new0(TMR_index,1)+
     &                          TCE_y_new0(TCE_index,1))/2)/
     &                          cell_range_y)
                           do i=1,TMR_res_num
                              TMR_x_new0(TMR_index,i)=
     &                             TMR_x_new0(TMR_index,i)-PB_x     
                              TMR_y_new0(TMR_index,i)=
     &                             TMR_y_new0(TMR_index,i)-PB_y     
                           enddo
                           do i=1,TCE_res_num
                              TCE_x_new0(TCE_index,i)=
     &                             TCE_x_new0(TCE_index,i)-PB_x     
                              TCE_y_new0(TCE_index,i)=
     &                             TCE_y_new0(TCE_index,i)-PB_y     
                           enddo
                        else
                           do i=1,TMR_res_num
                              TMR_x_new0(TMR_index,i)=
     &                             TMR_x(TMR_index,i)
                              TMR_y_new0(TMR_index,i)=
     &                             TMR_y(TMR_index,i)
                              TMR_z_new0(TMR_index,i)=
     &                             TMR_z(TMR_index,i)
                           enddo
                           do i=1,TCE_res_num
                              TCE_x_new0(TCE_index,i)=
     &                             TCE_x(TCE_index,i)
                              TCE_y_new0(TCE_index,i)=
     &                             TCE_y(TCE_index,i)
                              TCE_z_new0(TCE_index,i)=
     &                             TCE_z(TCE_index,i)
                           enddo
                        endif
ccc>>>>   Check Collisions
                        collision_flag=0
                        do i=1,TMR_tot_num
                           if(TMR_index.ne.i)then
                              dist=sqrt((TMR_x_new(i,1)-
     &                             TMR_x_new0(TMR_index,1))**2+
     &                             (TMR_y_new(i,1)-
     &                             TMR_y_new0(TMR_index,1))**2+
     &                             (TMR_z_new(i,1)-
     &                             TMR_z_new0(TMR_index,1))**2)
                              if(dist.lt.(TMR_radius+TMR_radius))then
                                 collision_flag=1
                              endif
                           endif
                        enddo
                        do i=1,TCE_tot_num
                           dist=sqrt((TCE_x_new(i,1)-
     &                          TMR_x_new0(TMR_index,4))**2+
     &                          (TCE_y_new(i,1)-
     &                          TMR_y_new0(TMR_index,4))**2+
     &                          (TCE_z_new(i,1)-
     &                          TMR_z_new0(TMR_index,4))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif
                           dist=sqrt((TCE_x_new(i,5)-
     &                          TMR_x_new0(TMR_index,4))**2+
     &                          (TCE_y_new(i,5)-
     &                          TMR_y_new0(TMR_index,4))**2+
     &                          (TCE_z_new(i,5)-
     &                          TMR_z_new0(TMR_index,4))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif						   
                        enddo
                        do i=1,TMR_tot_num
                           dist=sqrt((TMR_x_new(i,4)-
     &                          TCE_x_new0(TCE_index,1))**2+
     &                          (TMR_y_new(i,4)-
     &                          TCE_y_new0(TCE_index,1))**2+
     &                          (TMR_z_new(i,4)-
     &                          TCE_z_new0(TCE_index,1))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif
                           dist=sqrt((TMR_x_new(i,4)-
     &                          TCE_x_new0(TCE_index,5))**2+
     &                          (TMR_y_new(i,4)-
     &                          TCE_y_new0(TCE_index,5))**2+
     &                          (TMR_z_new(i,4)-
     &                          TCE_z_new0(TCE_index,5))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif						   
                        enddo
                        do i=1,TCE_tot_num
                           if(TCE_index.ne.i)then
                              dist=sqrt((TCE_x_new(i,1)-
     &                             TCE_x_new0(TCE_index,1))**2+
     &                             (TCE_y_new(i,1)-
     &                             TCE_y_new0(TCE_index,1))**2+
     &                             (TCE_z_new(i,1)-
     &                             TCE_z_new0(TCE_index,1))**2)
                              if(dist.lt.(TCE_radius+TCE_radius)
     &                             )then
                                 collision_flag=1
                              endif
                              dist=sqrt((TCE_x_new(i,5)-
     &                             TCE_x_new0(TCE_index,5))**2+
     &                             (TCE_y_new(i,5)-
     &                             TCE_y_new0(TCE_index,5))**2+
     &                             (TCE_z_new(i,5)-
     &                             TCE_z_new0(TCE_index,5))**2)
                              if(dist.lt.(TCE_radius+TCE_radius)
     &                             )then
                                 collision_flag=1
                              endif
                              dist=sqrt((TCE_x_new(i,1)-
     &                             TCE_x_new0(TCE_index,5))**2+
     &                             (TCE_y_new(i,1)-
     &                             TCE_y_new0(TCE_index,5))**2+
     &                             (TCE_z_new(i,1)-
     &                             TCE_z_new0(TCE_index,5))**2)
                              if(dist.lt.(TCE_radius+TCE_radius)
     &                             )then
                                 collision_flag=1
                              endif
                              dist=sqrt((TCE_x_new(i,5)-
     &                             TCE_x_new0(TCE_index,1))**2+
     &                             (TCE_y_new(i,5)-
     &                             TCE_y_new0(TCE_index,1))**2+
     &                             (TCE_z_new(i,5)-
     &                             TCE_z_new0(TCE_index,1))**2)
                              if(dist.lt.(TCE_radius+TCE_radius)
     &                             )then
                                 collision_flag=1
                              endif							  
                           endif
                        enddo
                        do i=1,TCR_tot_num
                           dist=sqrt((TCR_x_new(i,4)-
     &                          TCE_x_new0(TCE_index,1))**2+
     &                          (TCR_y_new(i,4)-
     &                          TCE_y_new0(TCE_index,1))**2+
     &                          (TCR_z_new(i,4)-
     &                          TCE_z_new0(TCE_index,1))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif
                           dist=sqrt((TCR_x_new(i,4)-
     &                          TCE_x_new0(TCE_index,5))**2+
     &                          (TCR_y_new(i,4)-
     &                          TCE_y_new0(TCE_index,5))**2+
     &                          (TCR_z_new(i,4)-
     &                          TCE_z_new0(TCE_index,5))**2)
                           if(dist.lt.TCE_radius)then
                              collision_flag=1
                           endif						   
                        enddo
                        if(collision_flag.eq.1)then
                           do i=1,TMR_res_num
                              TMR_x_new0(TMR_index,i)=
     &                             TMR_x(TMR_index,i)
                              TMR_y_new0(TMR_index,i)=
     &                             TMR_y(TMR_index,i)
                              TMR_z_new0(TMR_index,i)=
     &                             TMR_z(TMR_index,i)
                           enddo
                           do i=1,TCE_res_num
                              TCE_x_new0(TCE_index,i)=
     &                             TCE_x(TCE_index,i)
                              TCE_y_new0(TCE_index,i)=
     &                             TCE_y(TCE_index,i)
                              TCE_z_new0(TCE_index,i)=
     &                             TCE_z(TCE_index,i)
                           enddo
                        endif
c>>>>>>>>> Rotate the molecule
                       
                        cm1_a_x=0
                        cm1_a_y=0
                        cm1_a_z=0
                        do j=1,TMR_res_num
                           cm1_a_x=cm1_a_x+TMR_x_new0(TMR_index,j)
                           cm1_a_y=cm1_a_y+TMR_y_new0(TMR_index,j)
                           cm1_a_z=cm1_a_z+TMR_z_new0(TMR_index,j)
                        enddo
c                        do j=1,TCE_res_num
c                           cm1_a_x=cm1_a_x+TCE_x_new0(TCE_index,j)
c                           cm1_a_y=cm1_a_y+TCE_y_new0(TCE_index,j)
c                           cm1_a_z=cm1_a_z+TCE_z_new0(TCE_index,j)
c                        enddo
                        cm1_a_x=cm1_a_x/(real(TMR_res_num))
                        cm1_a_y=cm1_a_y/(real(TMR_res_num))
                        cm1_a_z=cm1_a_z/(real(TMR_res_num))
                        
                        rot_angle=2*rand3(r3)*Complex_rot_D-
     &                       Complex_rot_D
                        x_axis_bg=cm1_a_x
                        y_axis_bg=cm1_a_y
                        z_axis_bg=-100.0
                        x_axis_ed=cm1_a_x
                        y_axis_ed=cm1_a_y
                        z_axis_ed=100.0
                        do j=1,TMR_res_num
                           x_bf_rot=TMR_x_new0(TMR_index,j)
                           y_bf_rot=TMR_y_new0(TMR_index,j)
                           z_bf_rot=TMR_z_new0(TMR_index,j)
                           x_af_rot=0
                           y_af_rot=0
                           z_af_rot=0
                           call rot_along_axis(rot_angle,
     &                          x_axis_bg,y_axis_bg,z_axis_bg,
     &                          x_axis_ed,y_axis_ed,z_axis_ed,
     &                          x_bf_rot,y_bf_rot,z_bf_rot,
     &                          x_af_rot,y_af_rot,z_af_rot)
                           TMR_x_new(TMR_index,j)=x_af_rot
                           TMR_y_new(TMR_index,j)=y_af_rot
                           TMR_z_new(TMR_index,j)=z_af_rot
                        enddo
                        do j=1,TCE_res_num
                           x_bf_rot=TCE_x_new0(TCE_index,j)
                           y_bf_rot=TCE_y_new0(TCE_index,j)
                           z_bf_rot=TCE_z_new0(TCE_index,j)
                           x_af_rot=0
                           y_af_rot=0
                           z_af_rot=0
                           call rot_along_axis(rot_angle,
     &                          x_axis_bg,y_axis_bg,z_axis_bg,
     &                          x_axis_ed,y_axis_ed,z_axis_ed,
     &                          x_bf_rot,y_bf_rot,z_bf_rot,
     &                          x_af_rot,y_af_rot,z_af_rot)
                           TCE_x_new(TCE_index,j)=x_af_rot
                           TCE_y_new(TCE_index,j)=y_af_rot
                           TCE_z_new(TCE_index,j)=z_af_rot
                        enddo
                        
                     endif

                  endif





cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   else if the selected molecule is a TCE

               elseif((selecting_mole_index.gt.TCR_tot_num+TMR_tot_num)
     &                 .AND.
     &                 (selecting_mole_index.le.(TCR_tot_num+
     &                 TMR_tot_num+TCE_tot_num)))then
                  TCE_index=int(rand3(r3)*
     &                 TCE_tot_num)+1


cc>>>>   make random diffusion for selected molecule
cc>>>>   if it is unbound
                  if((TCE_1_status_new(TCE_index).eq.0).AND.
     &                 (TCE_2_status_new(TCE_index).eq.0))then
                     Prob_Diff=(6*TCE_D*time_step)/(distance_step**2)
                     if(rand3(r3).lt.Prob_Diff)then
cc>>     if smaller than the diffusion probability, move alone
                        theta=rand3(r3)*pai
                        phai=rand3(r3)*2*pai
                        do i=1,TCE_res_num
                           TCE_x_new0(TCE_index,i)=
     &                          TCE_x(TCE_index,i)+
     &                          distance_step*sin(theta)*cos(phai)
                           TCE_y_new0(TCE_index,i)=
     &                          TCE_y(TCE_index,i)+
     &                          distance_step*sin(theta)*sin(phai)
                           TCE_z_new0(TCE_index,i)=
     &                          TCE_z(TCE_index,i)+
     &                          distance_step*cos(theta)
                        enddo
                        lower_bound_flag=0
                        do i=1,TCE_res_num
                           if(TCE_z_new0(TCE_index,i).le.0)then
                              lower_bound_flag=1
                           endif
                        enddo
                        upper_bound_flag=0
                        do i=1,TCE_res_num
                           if(TCE_z_new0(TCE_index,i)
     &                          .ge.cell_range_z)then
                              upper_bound_flag=1
                           endif
                        enddo
                        if(lower_bound_flag.eq.1)then
                           do i=1,TCE_res_num
                              TCE_z_new0(TCE_index,i)=
     &                             TCE_z(TCE_index,i)+
     &                             abs(2*distance_step*cos(theta))
                           enddo
                        elseif(upper_bound_flag.eq.1)then
                           do i=1,TCE_res_num
                              TCE_z_new0(TCE_index,i)=
     &                             TCE_z(TCE_index,i)-
     &                             abs(2*distance_step*cos(theta))
                           enddo
                        endif
                        PB_x=cell_range_x*anint
     &                          (TCE_x_new0(TCE_index,1)/cell_range_x)
                        PB_y=cell_range_y*anint
     &                          (TCE_y_new0(TCE_index,1)/cell_range_y)
                        do i=1,TCE_res_num
                           TCE_x_new0(TCE_index,i)=
     &                          TCE_x_new0(TCE_index,i)-PB_x     
                           TCE_y_new0(TCE_index,i)=
     &                          TCE_y_new0(TCE_index,i)-PB_y
                        enddo
                     else
                        do i=1,TCE_res_num
                           TCE_x_new0(TCE_index,i)=
     &                          TCE_x(TCE_index,i)
                           TCE_y_new0(TCE_index,i)=
     &                          TCE_y(TCE_index,i)
                           TCE_z_new0(TCE_index,i)=
     &                          TCE_z(TCE_index,i)
                        enddo
                     endif
ccc>>>>   Check Collisions
                     collision_flag=0
                     do i=1,TCE_tot_num
                        if(TCE_index.ne.i)then
                           dist=sqrt((TCE_x_new(i,1)-
     &                          TCE_x_new0(TCE_index,1))**2+
     &                          (TCE_y_new(i,1)-
     &                          TCE_y_new0(TCE_index,1))**2+
     &                          (TCE_z_new(i,1)-
     &                          TCE_z_new0(TCE_index,1))**2)
                           if(dist.lt.(TCE_radius+TCE_radius))then
                              collision_flag=1
                           endif
                           dist=sqrt((TCE_x_new(i,5)-
     &                          TCE_x_new0(TCE_index,5))**2+
     &                          (TCE_y_new(i,5)-
     &                          TCE_y_new0(TCE_index,5))**2+
     &                          (TCE_z_new(i,5)-
     &                          TCE_z_new0(TCE_index,5))**2)
                           if(dist.lt.(TCE_radius+TCE_radius))then
                              collision_flag=1
                           endif
                           dist=sqrt((TCE_x_new(i,1)-
     &                          TCE_x_new0(TCE_index,5))**2+
     &                          (TCE_y_new(i,1)-
     &                          TCE_y_new0(TCE_index,5))**2+
     &                          (TCE_z_new(i,1)-
     &                          TCE_z_new0(TCE_index,5))**2)
                           if(dist.lt.(TCE_radius+TCE_radius))then
                              collision_flag=1
                           endif
                           dist=sqrt((TCE_x_new(i,5)-
     &                          TCE_x_new0(TCE_index,1))**2+
     &                          (TCE_y_new(i,5)-
     &                          TCE_y_new0(TCE_index,1))**2+
     &                          (TCE_z_new(i,5)-
     &                          TCE_z_new0(TCE_index,1))**2)
                           if(dist.lt.(TCE_radius+TCE_radius))then
                              collision_flag=1
                           endif						   
                        endif
                     enddo
                     do i=1,TCR_tot_num
                        dist=sqrt((TCR_x_new(i,4)-
     &                       TCE_x_new0(TCE_index,1))**2+
     &                       (TCR_y_new(i,4)-
     &                       TCE_y_new0(TCE_index,1))**2+
     &                       (TCR_z_new(i,4)-
     &                       TCE_z_new0(TCE_index,1))**2)
                        if(dist.lt.TCE_radius)then
                           collision_flag=1
                        endif
                        dist=sqrt((TCR_x_new(i,4)-
     &                       TCE_x_new0(TCE_index,5))**2+
     &                       (TCR_y_new(i,4)-
     &                       TCE_y_new0(TCE_index,5))**2+
     &                       (TCR_z_new(i,4)-
     &                       TCE_z_new0(TCE_index,5))**2)
                        if(dist.lt.TCE_radius)then
                           collision_flag=1
                        endif						
                     enddo
                     do i=1,TMR_tot_num
                        dist=sqrt((TMR_x_new(i,4)-
     &                       TCE_x_new0(TCE_index,1))**2+
     &                       (TMR_y_new(i,4)-
     &                       TCE_y_new0(TCE_index,1))**2+
     &                       (TMR_z_new(i,4)-
     &                       TCE_z_new0(TCE_index,1))**2)
                        if(dist.lt.TCE_radius)then
                           collision_flag=1
                        endif
                        dist=sqrt((TMR_x_new(i,4)-
     &                       TCE_x_new0(TCE_index,5))**2+
     &                       (TMR_y_new(i,4)-
     &                       TCE_y_new0(TCE_index,5))**2+
     &                       (TMR_z_new(i,4)-
     &                       TCE_z_new0(TCE_index,5))**2)
                        if(dist.lt.TCE_radius)then
                           collision_flag=1
                        endif						
                     enddo

                     if(collision_flag.eq.1)then
                        do i=1,TCE_res_num
                           TCE_x_new0(TCE_index,i)=
     &                          TCE_x(TCE_index,i)
                           TCE_y_new0(TCE_index,i)=
     &                          TCE_y(TCE_index,i)
                           TCE_z_new0(TCE_index,i)=
     &                          TCE_z(TCE_index,i)
                        enddo
                     endif
c>>>>>>>>> Rotate the molecule
                     theta=(2*rand3(r3)-1)*TCE_rot_D*pai/180.0
                     phi=(2*rand3(r3)-1)*TCE_rot_D*pai/180.0 
                     psai=(2*rand3(r3)-1)*TCE_rot_D*pai/180.0
            
                     t(1,1)=cos(psai)*cos(phi)-
     &                    cos(theta)*sin(phi)*sin(psai)
                     t(1,2)=-sin(psai)*cos(phi)-
     &                    cos(theta)*sin(phi)*cos(psai)
                     t(1,3)=sin(theta)*sin(phi)
                     
                     t(2,1)=cos(psai)*sin(phi)+
     &                    cos(theta)*cos(phi)*sin(psai)
                     t(2,2)=-sin(psai)*sin(phi)+
     &                    cos(theta)*cos(phi)*cos(psai)
                     t(2,3)=-sin(theta)*cos(phi)
                     
                     t(3,1)=sin(psai)*sin(theta)
                     t(3,2)=cos(psai)*sin(theta)
                     t(3,3)=cos(theta)

                     
            
                     CM_x=(TCE_x_new0(TCE_index,1)+
     &                    TCE_x_new0(TCE_index,5))*0.5
                     CM_y=(TCE_y_new0(TCE_index,1)+
     &                    TCE_y_new0(TCE_index,5))*0.5
                     CM_z=(TCE_z_new0(TCE_index,1)+
     &                    TCE_z_new0(TCE_index,5))*0.5

                     do k=1,TCE_res_num
                        TCE_x_new(TCE_index,k)=t(1,1)*
     &                       (TCE_x_new0(TCE_index,k)
     &                       -CM_x)+
     &                       t(1,2)*(TCE_y_new0(TCE_index,k)-
     &                       CM_y)
     &                       +t(1,3)*(TCE_z_new0(TCE_index,k)-
     &                       CM_z)
     &                       +CM_x
                        TCE_y_new(TCE_index,k)=t(2,1)*
     &                       (TCE_x_new0(TCE_index,k)
     &                       -CM_x)+
     &                       t(2,2)*(TCE_y_new0(TCE_index,k)-
     &                       CM_y)
     &                       +t(2,3)*(TCE_z_new0(TCE_index,k)-
     &                       CM_z)
     &                       +CM_y
                        TCE_z_new(TCE_index,k)=t(3,1)*
     &                       (TCE_x_new0(TCE_index,k)
     &                       -CM_x)+
     &                       t(3,2)*(TCE_y_new0(TCE_index,k)-
     &                       CM_y)
     &                       +t(3,3)*(TCE_z_new0(TCE_index,k)-
     &                       CM_z)
     &                       +CM_z
                     enddo



                  endif


               endif

            enddo



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   endo of diffusion for molecules
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   perform the reaction for molecules
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  1) TCR and TCE associate into a complex
            do i=1,TCR_tot_num
               if(TCR_status_new(i).eq.0)then
                  do j=1,TCE_tot_num
                     if(TCE_1_status_new(j).eq.0)then
c>>>>>>>>>>>>>   calculate distance of two reaction sites between TCR and TCE
                        dist=sqrt((TCE_x_new(j,4)-
     &                       TCR_x_new(i,4))**2+
     &                       (TCE_y_new(j,4)-
     &                       TCR_y_new(i,4))**2+
     &                       (TCE_z_new(j,4)-
     &                       TCR_z_new(i,4))**2)
                        if(dist.lt.bond_dist_cutoff)then
c                           print*,'Hello',dist
c>>>>>>>>>>>> calculate two angle
c>> perpendicular angle
                           theta_pd=0
                           point_x(1)=TCR_x_new(i,1)-TCR_x_new(i,4)
                           point_y(1)=TCR_y_new(i,1)-TCR_y_new(i,4)
                           point_z(1)=TCR_z_new(i,1)-TCR_z_new(i,4)
                           point_x(2)=0
                           point_y(2)=0
                           point_z(2)=0
                           point_x(3)=TCE_x_new(j,1)-TCE_x_new(j,4)
                           point_y(3)=TCE_y_new(j,1)-TCE_y_new(j,4)
                           point_z(3)=TCE_z_new(j,1)-TCE_z_new(j,4)
                           call gettheta(point_x,point_y,point_z,
     &                          theta_pd)
c>> orientational angle
                           theta_ot=0
                           point_x(1)=TCR_x_new(i,1)-TCR_x_new(i,2)
                           point_y(1)=TCR_y_new(i,1)-TCR_y_new(i,2)
                           point_z(1)=TCR_z_new(i,1)-TCR_z_new(i,2)
                           point_x(2)=0
                           point_y(2)=0
                           point_z(2)=0
                           point_x(3)=TCE_x_new(j,1)-TCE_x_new(j,2)
                           point_y(3)=TCE_y_new(j,1)-TCE_y_new(j,2)
                           point_z(3)=TCE_z_new(j,1)-TCE_z_new(j,2)
                           call gettheta(point_x,point_y,point_z,
     &                          theta_ot)
c                           print*,theta_pd,theta_ot
                           if((abs(theta_ot-bond_thetaot)
     &                          .lt.bond_thetaot_cutoff).AND.
     &                          (abs(theta_pd-bond_thetapd)
     &                          .lt.bond_thetapd_cutoff))then                         

                              Prob_Ass=Ass_Rate*time_step ! ready to change
                              if(rand3(r3).lt.Prob_Ass)then
                                 TCR_status_new(i)=1
                                 TCE_1_status_new(j)=1
                                 complex_num_new=
     &                                complex_num_new+1
                                 complex_nb_num_new(complex_num_new)=2
                                 complex_nb_ctg_new(complex_num_new,1)=1
                                 complex_nb_ctg_new(complex_num_new,2)=2
                                 complex_nb_ctg_idx_new
     &                                (complex_num_new,1)=i
                                 complex_nb_ctg_idx_new
     &                                (complex_num_new,2)=j
                              endif
                           endif
                        endif
                     endif
                  enddo
               endif
            enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  2)  complex dissociate into TCR and TCE
            do i=1,TCR_tot_num
               if(TCR_status_new(i).eq.1)then
                  do j=1,complex_num_new
                     if(complex_nb_ctg_idx_new(j,1).eq.i)then
                        selected_complex=j
                        selected_TCR=complex_nb_ctg_idx_new(j,1)
                        selected_TCE=complex_nb_ctg_idx_new(j,2)
                     endif
                  enddo
                  Prob_Diss=Diss_Rate*time_step
                  part1=int(rand3(r3)*10)
                  part2=int(rand5(r5)*10)
                  part3=int(rand5(r5)*10)
                  temp=rand3(r3)
                  part4=int(rand3(r3)*10)
                  part5=int(rand4(r4)*10)
                  do j=1,int(rand4(r4)*10)
                     temp=rand5(r5)
                     temp=rand5(r5)
                  enddo
                  part6=int(rand5(r5)*10)
                  part7=int(rand5(r5)*10)
                  part8=int(rand4(r4)*10)
                  part9=int(rand5(r5)*10)
                  
                  prob=real(part1)/10.0+real(part2)/100.0
     &                 +real(part3)/1000.0+real(part4)/10000.0
     &                 +real(part5)/100000.0
     &                 +real(part6)/1000000.0
     &                 +real(part7)/10000000.0
     &                 +real(part8)/100000000.0
     &                 +real(part9)/1000000000.0
                  
                  if((prob.lt.Prob_Diss).AND.(prob.gt.0.000))then
                     TCR_status_new(selected_TCR)=0
                     TCE_1_status_new(selected_TCE)=0
                     do k=selected_complex+1,complex_num_new
                        complex_nb_num_new(k-1)=complex_nb_num_new(k)
                        complex_nb_ctg_idx_new(k-1,1)=
     &                       complex_nb_ctg_idx_new(k,1)
                        complex_nb_ctg_idx_new(k-1,2)=
     &                       complex_nb_ctg_idx_new(k,2)
                        complex_nb_ctg_new(k-1,1)=
     &                       complex_nb_ctg_new(k,1)
                        complex_nb_ctg_new(k-1,2)=
     &                       complex_nb_ctg_new(k,2)
                     enddo
                     complex_num_new=
     &                    complex_num_new-1
                  endif
               endif
            enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  3) TMR and TCE associate into a complex
            do i=1,TMR_tot_num
               if(TMR_status_new(i).eq.0)then
                  do j=1,TCE_tot_num
                     if(TCE_2_status_new(j).eq.0)then
c>>>>>>>>>>>>>   calculate distance of two reaction sites between TCR and TCE
                        dist=sqrt((TCE_x_new(j,8)-
     &                       TMR_x_new(i,4))**2+
     &                       (TCE_y_new(j,8)-
     &                       TMR_y_new(i,4))**2+
     &                       (TCE_z_new(j,8)-
     &                       TMR_z_new(i,4))**2)
                        if(dist.lt.bond_dist_cutoff)then
c                           print*,'Hello',dist
c>>>>>>>>>>>> calculate two angle
c>> perpendicular angle
                           theta_pd=0
                           point_x(1)=TMR_x_new(i,1)-TMR_x_new(i,4)
                           point_y(1)=TMR_y_new(i,1)-TMR_y_new(i,4)
                           point_z(1)=TMR_z_new(i,1)-TMR_z_new(i,4)
                           point_x(2)=0
                           point_y(2)=0
                           point_z(2)=0
                           point_x(3)=TCE_x_new(j,5)-TCE_x_new(j,8)
                           point_y(3)=TCE_y_new(j,5)-TCE_y_new(j,8)
                           point_z(3)=TCE_z_new(j,5)-TCE_z_new(j,8)
                           call gettheta(point_x,point_y,point_z,
     &                          theta_pd)
c>> orientational angle
                           theta_ot=0
                           point_x(1)=TMR_x_new(i,1)-TMR_x_new(i,2)
                           point_y(1)=TMR_y_new(i,1)-TMR_y_new(i,2)
                           point_z(1)=TMR_z_new(i,1)-TMR_z_new(i,2)
                           point_x(2)=0
                           point_y(2)=0
                           point_z(2)=0
                           point_x(3)=TCE_x_new(j,5)-TCE_x_new(j,6)
                           point_y(3)=TCE_y_new(j,5)-TCE_y_new(j,6)
                           point_z(3)=TCE_z_new(j,5)-TCE_z_new(j,6)
                           call gettheta(point_x,point_y,point_z,
     &                          theta_ot)
c                           print*,theta_pd,theta_ot
                           if((abs(theta_ot-bond_thetaot)
     &                          .lt.bond_thetaot_cutoff).AND.
     &                          (abs(theta_pd-bond_thetapd)
     &                          .lt.bond_thetapd_cutoff))then                         
                              Prob_Ass=Ass_Rate*time_step ! ready to change
                              if(rand3(r3).lt.Prob_Ass)then
                                 TMR_status_new(i)=1
                                 TCE_2_status_new(j)=1
                                 complex2_num_new=
     &                                complex2_num_new+1
                                 complex2_nb_num_new(complex2_num_new)=2
                                 complex2_nb_ctg_new
     &                                (complex2_num_new,1)=1
                                 complex2_nb_ctg_new
     &                                (complex2_num_new,2)=2
                                 complex2_nb_ctg_idx_new
     &                                (complex2_num_new,1)=i
                                 complex2_nb_ctg_idx_new
     &                                (complex2_num_new,2)=j
                              endif
                           endif
                        endif
                     endif
                  enddo
               endif
            enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  4)  complex dissociate into TMR and TCE
            do i=1,TMR_tot_num
               if(TMR_status_new(i).eq.1)then
                  do j=1,complex2_num_new
                     if(complex2_nb_ctg_idx_new(j,1).eq.i)then
                        selected_complex=j
                        selected_TMR=complex2_nb_ctg_idx_new(j,1)
                        selected_TCE=complex2_nb_ctg_idx_new(j,2)
                     endif
                  enddo
                  Prob_Diss=Diss_Rate*time_step

                  part1=int(rand3(r3)*10)
                  part2=int(rand5(r5)*10)
                  part3=int(rand5(r5)*10)
                  temp=rand3(r3)
                  part4=int(rand3(r3)*10)
                  part5=int(rand4(r4)*10)
                  do j=1,int(rand4(r4)*10)
                     temp=rand5(r5)
                     temp=rand5(r5)
                  enddo
                  part6=int(rand5(r5)*10)
                  part7=int(rand5(r5)*10)
                  part8=int(rand4(r4)*10)
                  part9=int(rand5(r5)*10)
                  
                  prob=real(part1)/10.0+real(part2)/100.0
     &                 +real(part3)/1000.0+real(part4)/10000.0
     &                 +real(part5)/100000.0
     &                 +real(part6)/1000000.0
     &                 +real(part7)/10000000.0
     &                 +real(part8)/100000000.0
     &                 +real(part9)/1000000000.0
                  
                  if((prob.lt.Prob_Diss).AND.(prob.gt.0.000))then
                     TMR_status_new(selected_TMR)=0
                     TCE_2_status_new(selected_TCE)=0
                     do k=selected_complex+1,complex2_num_new
                        complex2_nb_num_new(k-1)=complex2_nb_num_new(k)
                        complex2_nb_ctg_idx_new(k-1,1)=
     &                       complex2_nb_ctg_idx_new(k,1)
                        complex2_nb_ctg_idx_new(k-1,2)=
     &                       complex2_nb_ctg_idx_new(k,2)
                        complex2_nb_ctg_new(k-1,1)=
     &                       complex2_nb_ctg_new(k,1)
                        complex2_nb_ctg_new(k-1,2)=
     &                       complex2_nb_ctg_new(k,2)
                     enddo
                     complex2_num_new=
     &                    complex2_num_new-1
                  endif
               endif
            enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc         simulation update
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   update the coordinates for all molecules

            do i=1,TCR_tot_num
               do j=1,TCR_res_num
                  TCR_x(i,j)=TCR_x_new(i,j)
                  TCR_y(i,j)=TCR_y_new(i,j)
                  TCR_z(i,j)=TCR_z_new(i,j)
               enddo
                  TCR_status(i)=TCR_status_new(i)
            enddo
            do i=1,TMR_tot_num
               do j=1,TMR_res_num
                  TMR_x(i,j)=TMR_x_new(i,j)
                  TMR_y(i,j)=TMR_y_new(i,j)
                  TMR_z(i,j)=TMR_z_new(i,j)
               enddo
                  TMR_status(i)=TMR_status_new(i)
            enddo			
            do i=1,TCE_tot_num
               do j=1,TCE_res_num
                  TCE_x(i,j)=TCE_x_new(i,j)
                  TCE_y(i,j)=TCE_y_new(i,j)
                  TCE_z(i,j)=TCE_z_new(i,j)
               enddo
                  TCE_1_status(i)=TCE_1_status_new(i)
                  TCE_2_status(i)=TCE_2_status_new(i)				  
            enddo
            complex_num=complex_num_new
            do i=1,complex_num
               complex_nb_num(i)=complex_nb_num_new(i)
               do j=1,complex_nb_ctg_num
                  complex_nb_ctg(i,j)=complex_nb_ctg_new(i,j)
                  complex_nb_ctg_idx(i,j)=complex_nb_ctg_idx_new(i,j)
               enddo
            enddo
            complex2_num=complex2_num_new
            do i=1,complex2_num
               complex2_nb_num(i)=complex2_nb_num_new(i)
               do j=1,complex2_nb_ctg_num
                  complex2_nb_ctg(i,j)=complex2_nb_ctg_new(i,j)
                  complex2_nb_ctg_idx(i,j)=complex2_nb_ctg_idx_new(i,j)
               enddo
            enddo

            xlink_num=0
            do i=1,TCE_tot_num
               if((TCE_1_status(i).eq.1).AND.(TCE_2_status(i).eq.1))then
                  xlink_num=xlink_num+1
               endif
            enddo			   

            goto 400
c            print*,mc_time_step,complex_num
            if(MOD(mc_time_step,1000).eq.0)then
              open(unit=10,file='TCEoutputRIMD_ene_08242022_test3.dat',
     &              status='unknown',access='append')
               write(10,2103) mc_time_step,complex_num,
     &              complex2_num,xlink_num
 2103          format(4I10)
               close(10)
            endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c  output data along the trajectory
ccccccccccccccccccccccccccccccccccccccccccccccccccc           
           
            if(MOD(mc_time_step,10000).eq.0)then
               open (unit=10,file=
     &              'TCEoutputRIMD_trj_08242022_test3.pdb',
     &              status='unknown',access='append')
               do i=1,TCR_tot_num
                  write(10,2100) 'ATOM  ',i,' CA ','VAL'
     &                 ,i,TCR_x(i,1),TCR_y(i,1),TCR_z(i,1)
                  write(10,2100) 'ATOM  ',i,' CA ','VAL'
     &                 ,i,TCR_x(i,2),TCR_y(i,2),TCR_z(i,2)
                  write(10,2100) 'ATOM  ',i,' CA ','VAL'
     &                 ,i,TCR_x(i,3),TCR_y(i,3),TCR_z(i,3)
                  write(10,2100) 'ATOM  ',i,' CA ','ALA'
     &                 ,i,TCR_x(i,4),TCR_y(i,4),TCR_z(i,4)
               enddo
               write(10,2102) 'TER'    
               do i=1,TCE_tot_num
                  write(10,2100) 'ATOM  ',i,' CA ','GLY'
     &                 ,i,TCE_x(i,1),TCE_y(i,1),TCE_z(i,1)
                  write(10,2100) 'ATOM  ',i,' CA ','GLY'
     &                 ,i,TCE_x(i,2),TCE_y(i,2),TCE_z(i,2)
                  write(10,2100) 'ATOM  ',i,' CA ','GLY'
     &                 ,i,TCE_x(i,3),TCE_y(i,3),TCE_z(i,3)
                  write(10,2100) 'ATOM  ',i,' CA ','LEU'
     &                 ,i,TCE_x(i,4),TCE_y(i,4),TCE_z(i,4)
                  write(10,2100) 'ATOM  ',i,' CA ','GLY'
     &                 ,i,TCE_x(i,5),TCE_y(i,5),TCE_z(i,5)
                  write(10,2100) 'ATOM  ',i,' CA ','GLY'
     &                 ,i,TCE_x(i,6),TCE_y(i,6),TCE_z(i,6)
                  write(10,2100) 'ATOM  ',i,' CA ','GLY'
     &                 ,i,TCE_x(i,7),TCE_y(i,7),TCE_z(i,7)
                  write(10,2100) 'ATOM  ',i,' CA ','ILE'
     &                 ,i,TCE_x(i,8),TCE_y(i,8),TCE_z(i,8)	 
               enddo
               write(10,2102) 'TER'  
               do i=1,TMR_tot_num
                  write(10,2100) 'ATOM  ',i,' CA ','ARG'
     &                 ,i,TMR_x(i,1),TMR_y(i,1),TMR_z(i,1)
                  write(10,2100) 'ATOM  ',i,' CA ','ARG'
     &                 ,i,TMR_x(i,2),TMR_y(i,2),TMR_z(i,2)
                  write(10,2100) 'ATOM  ',i,' CA ','ARG'
     &                 ,i,TMR_x(i,3),TMR_y(i,3),TMR_z(i,3)
                  write(10,2100) 'ATOM  ',i,' CA ','LYS'
     &                 ,i,TMR_x(i,4),TMR_y(i,4),TMR_z(i,4)
               enddo			   
               write(10,2102) 'END'    
               close(10)
            endif

 400        continue

 2100       format(A6,I5,1x,A4,1x,A3,1x,I5,4x,3F8.3)
 2102       format(A3)          

ccccccccccccccccccccccccccccccccccc
cc  end current simulation step
ccccccccccccccccccccccccccccccccccc

         enddo

         open (unit=10,file=
     &        'TCERIMD_output_Aff21TCR100TMR100TCE500.dat',
     &        status='unknown',access='append')
         write(10,2103) n_t,complex_num,complex2_num,xlink_num
         close(10)
21043    format(4I10)

         xlink_ave=xlink_ave+real(xlink_num)

cccccccccccccccccccccccccccc
cc end current trajectory
cccccccccccccccccccccccccccc

      enddo

      open (unit=10,file=
     &     'TCERIMD_output_Aff21TCR100TMR100TCE500.dat',
     &     status='unknown',access='append')
      write(10,2104) xlink_ave/real(num_trajec)
      close(10)
 2104 format(F12.5)

cccccccccccccccccccc

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
ccccccccccccccccccccccccccccccccccccccccccccccccc            
            
      subroutine gettheta(point_x,point_y,point_z,test_theta)

      implicit none
      real*8 point_x(3),point_y(3),point_z(3),test_theta
ccc
      real*8 lx(2),ly(2),lz(2),lr(2)
      real*8 conv,doth1,doth2,conv2
      integer in

ccccccccccccccccccc
c>>  bond lengt
ccccccccccccccccccc

      do in=1,2
         lx(in)=0
         ly(in)=0
         lz(in)=0
         lr(in)=0
      enddo

      lx(1)=point_x(2)-point_x(1)
      ly(1)=point_y(2)-point_y(1)
      lz(1)=point_z(2)-point_z(1)
      lr(1)=sqrt(lx(1)**2+ly(1)**2+lz(1)**2)
      
      lx(2)=point_x(3)-point_x(2)
      ly(2)=point_y(3)-point_y(2)
      lz(2)=point_z(3)-point_z(2)
      lr(2)=sqrt(lx(2)**2+ly(2)**2+lz(2)**2)
 
cccccccccccccccccccc
c>>  theta value
cccccccccccccccccccc
      
      test_theta=0

      conv=180/3.14159 
      do in=1,1
         doth1=-1*(lx(in+1)*lx(in)+ly(in+1)*ly(in)+lz(in+1)*lz(in))
         doth2=doth1/(lr(in+1)*lr(in))
         if(doth2.gt.1)then
            doth2=1
         endif
         if(doth2.lt.-1)then
            doth2=-1
         endif
         test_theta=acos(doth2)*conv
      enddo

cccccccccccccccc

      return
      end

cccccccccccccccccccccccccccccccccccc
                        
