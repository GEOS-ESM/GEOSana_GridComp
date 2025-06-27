program read_binary
  implicit none
  integer :: nr, i, nbytes
  real(kind=4), allocatable :: CNN_PBL_z(:), CNN_lat(:), CNN_lon(:)
  real(kind=8), allocatable :: CNN_jdy(:)
  integer :: unit_number
  character(len=100) :: filename

  ! Specify the filename
  !filename = "CNN_results_file"  ! Replace with the actual file name
  filename = "CATS-CNN_PBL_N-M7.2-V3.00.2015-09-16T21-11-58T21-57-12UTC.dat"
  unit_number = 10               ! File connection unit number
  
  ! Open the file
  open(unit=unit_number, file=filename, status="old", form="unformatted", access="stream")
  
  ! Read the first 4 bytes (nr: data length)
  read(unit_number) nr
  print *, "Number of elements (nr):", nr
  
  ! Allocate memory
  allocate(CNN_PBL_z(nr), CNN_lat(nr), CNN_lon(nr))
  allocate(CNN_jdy(nr))
  
  ! Read CNN_PBL_z (float32)
  read(unit_number) CNN_PBL_z
  print *, "CNN_PBL_z read successfully."

  ! Read CNN_lat (float32)
  read(unit_number) CNN_lat
  print *, "CNN_lat read successfully."

  ! Read CNN_lon (float32)
  read(unit_number) CNN_lon
  print *, "CNN_lon read successfully."

  ! Read CNN_jdy (float64)
  read(unit_number) CNN_jdy
  print *, "CNN_jdy read successfully."

  ! Print sample data (for verification)
  print *, "First element of CNN_PBL_z:", CNN_PBL_z(1)
  print *, "First element of CNN_lat:", CNN_lat(1)
  print *, "First element of CNN_lon:", CNN_lon(1)
  print *, "First element of CNN_jdy:", CNN_jdy(1)

  ! Deallocate memory and close the file
  deallocate(CNN_PBL_z, CNN_lat, CNN_lon, CNN_jdy)
  close(unit_number)
  
end program read_binary

