module mod_no_tree
use data_structure
use mod_interp
integer, parameter   :: n_no_tree = 400        ! the number of boxes in each direction (could be made different in R and Z)
integer, allocatable :: no_tree(:)             ! contains the element numbers overlapping with the bounding box (i,j) (i+1,j+1)
integer              :: no_tree_start(n_no_tree,n_no_tree)   ! contains the index in no_tree for the entries from (i,j)
real*8               :: no_tree_bounding_box(2,2)
real*8               :: no_tree_tolerance = 1.d-6
integer              :: n_entries(n_no_tree,n_no_tree)
real*8, allocatable  :: no_tree_r_min(:), no_tree_r_max(:), no_tree_z_min(:), no_tree_z_max(:)

contains

subroutine no_tree_init(node_list, element_list)
implicit none
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list

integer             :: offset(n_no_tree, n_no_tree)
integer             :: i, j, k, i_prev, j_prev, i_start_prev
integer             :: index_i_start, index_i_end, index_j_start, index_j_end

write(*,'(A)')      '*************************************'
write(*,'(A,i3,A)') '* no_tree_init, n : ',n_no_tree,'    *'
write(*,'(A)')      '*************************************'

if (element_list%n_elements .le. 0) return

no_tree_bounding_box(1,:) =  1.d33
no_tree_bounding_box(2,:) = -1.d33

if (allocated(no_tree_r_min)) deallocate(no_tree_r_min)
if (allocated(no_tree_r_max)) deallocate(no_tree_r_max)
if (allocated(no_tree_z_min)) deallocate(no_tree_z_min)
if (allocated(no_tree_z_max)) deallocate(no_tree_z_max)

allocate(no_tree_r_min(element_list%n_elements), no_tree_r_max(element_list%n_elements), &
         no_tree_z_min(element_list%n_elements), no_tree_z_max(element_list%n_elements))

do i=1, element_list%n_elements

  call RZ_minmax(node_list, element_list, i, no_tree_r_min(i), no_tree_r_max(i), no_tree_z_min(i), no_tree_z_max(i))

  no_tree_bounding_box(1,1) = min(no_tree_bounding_box(1,1), no_tree_r_min(i))
  no_tree_bounding_box(1,2) = min(no_tree_bounding_box(1,2), no_tree_z_min(i))
  no_tree_bounding_box(2,1) = max(no_tree_bounding_box(2,1), no_tree_r_max(i))
  no_tree_bounding_box(2,2) = max(no_tree_bounding_box(2,2), no_tree_z_max(i))

end do

no_tree_bounding_box(1,:) = no_tree_bounding_box(1,:) - no_tree_tolerance
no_tree_bounding_box(2,:) = no_tree_bounding_box(2,:) + no_tree_tolerance

write(*,'(A,6f12.6)') ' no_tree bounding box : ',no_tree_bounding_box, (no_tree_bounding_box(2,:) - no_tree_bounding_box(1,:))/n_no_tree

n_entries = 0

do k=1, element_list%n_elements       ! count the number of elements overlapping box(i,j)

! what if the element is smaller than the box, are the indices still valid???  

  index_i_start = (no_tree_r_min(k) - no_tree_bounding_box(1,1)) / (no_tree_bounding_box(2,1) - no_tree_bounding_box(1,1)) * n_no_tree + 1
  index_i_end   = (no_tree_r_max(k) - no_tree_bounding_box(1,1)) / (no_tree_bounding_box(2,1) - no_tree_bounding_box(1,1)) * n_no_tree + 1
  index_j_start = (no_tree_z_min(k) - no_tree_bounding_box(1,2)) / (no_tree_bounding_box(2,2) - no_tree_bounding_box(1,2)) * n_no_tree + 1
  index_j_end   = (no_tree_z_max(k) - no_tree_bounding_box(1,2)) / (no_tree_bounding_box(2,2) - no_tree_bounding_box(1,2)) * n_no_tree + 1

  do i = index_i_start, index_i_end
    do j= index_j_start, index_j_end
      n_entries(i,j)  = n_entries(i,j) + 1
    enddo
  enddo

enddo

write(*,*) ' total number of entries : ',sum(n_entries)
if (allocated(no_tree)) deallocate(no_tree)
allocate(no_tree(sum(n_entries)))

no_tree_start(1,1) = 1
i_start_prev       = no_tree_start(1,1)

i_prev = 1
j_prev = 1

do k = 2, n_no_tree * n_no_tree
 
  i = (k-1) / (n_no_tree) + 1
  j = mod(k-1, n_no_tree) + 1

  no_tree_start(i,j) = i_start_prev + n_entries(i_prev,j_prev)

  i_start_prev = no_tree_start(i,j)

  i_prev = i
  j_prev = j


enddo

offset = 0

do k=1, element_list%n_elements       ! count the number of elements overlapping box(i,j)

  index_i_start = (no_tree_r_min(k) - no_tree_bounding_box(1,1)) / (no_tree_bounding_box(2,1) - no_tree_bounding_box(1,1)) * n_no_tree + 1
  index_i_end   = (no_tree_r_max(k) - no_tree_bounding_box(1,1)) / (no_tree_bounding_box(2,1) - no_tree_bounding_box(1,1)) * n_no_tree + 1
  index_j_start = (no_tree_z_min(k) - no_tree_bounding_box(1,2)) / (no_tree_bounding_box(2,2) - no_tree_bounding_box(1,2)) * n_no_tree + 1
  index_j_end   = (no_tree_z_max(k) - no_tree_bounding_box(1,2)) / (no_tree_bounding_box(2,2) - no_tree_bounding_box(1,2)) * n_no_tree + 1

  do i = index_i_start, index_i_end
    do j= index_j_start, index_j_end
      no_tree(no_tree_start(i,j) + offset(i,j)) = k
      offset(i,j) = offset(i,j) + 1
    enddo
  enddo

enddo

end

!try to find all the elements in the box continaing i_elm
subroutine nearby_elements_no_tree(node_list, element_list, i_elm, i_nearby)
implicit none
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list
integer                             :: i_elm
integer, allocatable, intent(out)   :: i_nearby(:)
real*8                              :: x1(2), x2(2), x3(2), x4(2)
integer                             :: i, j, k, i_err
integer, dimension(:), allocatable  :: i_elms_1, i_elms_2, i_elms_3, i_elms_4


call interp_RZ(node_list, element_list, i_elm, 0.5d0, 0.0d0, x1(1), x1(2))
call interp_RZ(node_list, element_list, i_elm, 1.0d0, 0.5d0, x2(1), x2(2))
call interp_RZ(node_list, element_list, i_elm, 0.5d0, 1.0d0, x3(1), x3(2))
call interp_RZ(node_list, element_list, i_elm, 0.0d0, 0.5d0, x4(1), x4(2))

call elements_containing_point_no_tree(x1(1), x1(2), i_elms_1, .false.) 
call elements_containing_point_no_tree(x2(1), x2(2), i_elms_2, .false.) 
call elements_containing_point_no_tree(x3(1), x3(2), i_elms_3, .false.) 
call elements_containing_point_no_tree(x4(1), x4(2), i_elms_4, .false.) 

allocate(i_nearby(size(i_elms_1)+size(i_elms_2)+size(i_elms_3)+size(i_elms_4)))

if (i_elm .eq. 1) then
  i_nearby = [[i_elm+1], i_elms_1, i_elms_2, i_elms_3, i_elms_4]
else
  i_nearby = [[i_elm-1, i_elm+1], i_elms_1, i_elms_2, i_elms_3, i_elms_4]
endif
! i_nearby = [i_elms_1, i_elms_2, i_elms_3, i_elms_4]
end

!< Find the node containing the point
subroutine elements_containing_point_no_tree(R, Z, i_elms, filter)
implicit none
real*8, intent(in)                              :: R, Z
integer, dimension(:), allocatable, intent(out) :: i_elms
logical, intent(in) , optional                  :: filter
  
integer, dimension(:), allocatable              :: i_elms_full
real*8, dimension(:), allocatable               :: r_min_local, r_max_local, z_min_local, z_max_local
integer                                         :: i_err, i, j
logical                                         :: do_filter = .false.
real*8                                          :: tol = 1.d-6

if (present(filter)) then
  do_filter = filter
endif    

i = (R - no_tree_bounding_box(1,1)) / (no_tree_bounding_box(2,1) - no_tree_bounding_box(1,1)) * n_no_tree + 1
j = (Z - no_tree_bounding_box(1,2)) / (no_tree_bounding_box(2,2) - no_tree_bounding_box(1,2)) * n_no_tree + 1

if ((i .le. 0) .or. (j .le. 0) .or. (i .gt. n_no_tree) .or. (j .gt. n_no_tree)) return

if (.not. do_filter) then

  allocate(i_elms(n_entries(i,j)))
  i_elms = no_tree(no_tree_start(i,j):no_tree_start(i,j)+n_entries(i,j)-1)

else

  allocate(i_elms_full(n_entries(i,j)), r_min_local(n_entries(i,j)), r_max_local(n_entries(i,j)), &
                                        z_min_local(n_entries(i,j)), z_max_local(n_entries(i,j)))

  i_elms_full = no_tree(no_tree_start(i,j):no_tree_start(i,j)+n_entries(i,j)-1)

  r_min_local = no_tree_r_min(i_elms_full) - tol
  r_max_local = no_tree_r_max(i_elms_full) + tol
  z_min_local = no_tree_z_min(i_elms_full) - tol
  z_max_local = no_tree_z_max(i_elms_full) + tol
      
  i_elms = pack(i_elms_full, ((R .ge. r_min_local(:)) &
                        .and. (R .le. r_max_local(:)) & 
                        .and. (Z .ge. z_min_local(:)) &
                        .and. (Z .le. z_max_local(:))))
endif

end

end module
