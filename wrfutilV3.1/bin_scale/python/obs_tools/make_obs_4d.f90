program main 
implicit real(a-h,o-z)

character(len=300) :: infile,outfile
character(len=20) :: cread

real(4)::wk(8)

call getarg(1,infile)
call getarg(2,outfile)
call getarg(3,cread)
read(cread,*) dif
call getarg(4,cread)
read(cread,*) iradar


open(11,file=trim(infile),form="unformatted",convert="little_endian")
if (iradar == 1) then
  read(11,iostat=ios) radarlon
  if ( ios /= 0 ) then
    write(*,*) ios
   stop
 end if
  read(11) radarlat
  read(11) radarz
end if

open(21,file=trim(outfile),form="unformatted",status="old",position="append",iostat=ierror,convert="little_endian")
if (ierror==0)then
  wk(8)=dif
  ios=0
  do while (ios==0)
    read(11,iostat=ios) wk(1:7)
    if (ios==0) write(21) wk(1:8)
  end do
else
  open(21,file=trim(outfile),form="unformatted",status="new",iostat=ierror,convert="little_endian")
  if (iradar==1) then
    write(21) radarlon
    write(21) radarlat
    write(21) radarz
  end if
  wk(8)=dif
  ios=0
  do while (ios==0)
    read(11,iostat=ios) wk(1:7)
    if (ios==0) write(21) wk(1:8)
  end do

end if
close(21)

stop
end program



