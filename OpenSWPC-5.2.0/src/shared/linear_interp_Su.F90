module linear_interp
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!https://blog.csdn.net/webzhuce/article/details/86585489
!!!data from Ilan_tomo.xyzv
!!!Model(lon_sample_num, lat_sample_num, dep_sample_num) = Model(k,j,i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

    subroutine interp(Model,x,y,z)
            real,intent(out) :: Model(96,82,78), x(96), y(82), z(78)
            integer :: aa
            real :: lon, lat, depth, p_vel
            integer :: i, j, k
            !!build the 3-D Vmodel matrix
            do i = 1,78
                z(i) = -99999.0
            enddo
            
            do j = 1,82
                y(j) = -99999.0
            enddo

            do k = 1,96
                x(k) = -99999.0
            enddo

            open(1010,file = 'Ilan_tomo.xyzv',status = 'OLD',action = 'read', access = 'sequential')
            do i = 1,78
                do j = 1,82
                    do k = 1,96                        
                        read(1010,*)lon, lat, depth, p_vel
                        do while (x(k) == -99999.0)
                            x(k) = lon
                        end do

                        do while (y(j) == -99999.0)
                            y(j) = lat
                        end do

                        do while (z(i) == -99999.0)
                            z(i) = depth
                        end do                       
                        
                        Model(k,j,i) = p_vel
                        
        
                    end do
                end do
            end do
           !write(*,'(F5.3,1x,F5.3,1x,F5.3)')Model(1,1,1),Model(10,2,1),Model(96,82,78) 
            close(1010)
            !write(*,*)z(78),y(3),x(50)
    end subroutine


    subroutine find_location(Model, x, y, z, lon, lat, dep, V_interp)

            real, intent(in) :: Model(96,82,78), x(96), y(82), z(78)
            real, intent(in) :: lon, lat, dep
            real, intent(out) :: V_interp
            integer :: i, j, k, x_0_index, x_1_index, y_0_index, y_1_index, z_0_index, z_1_index
            real :: V_000, V_001, V_011, V_111, V_100, V_110, V_101, V_010
            real :: x_d, y_d, z_d


            !find x-dir location
            do i = 1, 96
                if (x(i) > lon) then
                        x_0_index = i-1
                        x_1_index = i
                        exit
                else if (x(i) == lon) then
                        x_0_index = i
                        x_1_index = i+1
                        exit
                else
                        cycle
                end if            
            end do
            !find y-dir location(y-dir is decreasing)    
            do j = 1, 82
                if (y(j) < lat) then
                        y_0_index = j
                        y_1_index = j-1
                        exit
                else if (y(i) == lat) then
                        y_0_index = j
                        y_1_index = j-1
                        exit
                else
                        cycle
                end if          
            end do
            !find z-dir location
            do k = 1, 78
                if (z(k) > dep) then
                        z_0_index = k-1
                        z_1_index = k
                        exit
                else if (z(k) == dep) then
                        z_0_index = k
                        z_1_index = k+1
                        exit
                else
                        cycle
                end if            
            end do
            

            V_000 = Model(x_0_index, y_0_index, z_0_index)
            V_100 = Model(x_1_index, y_0_index, z_0_index)
            V_110 = Model(x_1_index, y_1_index, z_0_index)
            V_010 = Model(x_0_index, y_1_index, z_0_index)
            V_011 = Model(x_0_index, y_1_index, z_1_index)
            V_001 = Model(x_0_index, y_0_index, z_1_index)
            V_101 = Model(x_1_index, y_0_index, z_1_index)
            V_111 = Model(x_1_index, y_1_index, z_1_index)

            x_d = (lon - x(x_0_index)) / (x(x_1_index) - x(x_0_index))
            y_d = (lat - y(y_0_index)) / (y(y_1_index) - y(y_0_index))
            z_d = (dep - z(z_0_index)) / (z(z_1_index) - z(z_0_index))
          
            V_interp = V_000 * (1-x_d) * (1-y_d) * (1-z_d) + &
                       V_100 * x_d *(1-y_d) * (1-z_d) + &
                       V_001 * (1-x_d) * (1-y_d) * z_d + &
                       V_010 * (1-x_d) * y_d * (1-z_d) + &
                       V_101 * x_d * (1-y_d) * z_d + &
                       V_011 * (1-x_d) * y_d * z_d + &
                       V_110 * x_d * y_d *(1-z_d) + &
                       V_111 * x_d * y_d * z_d


    end subroutine





end module

