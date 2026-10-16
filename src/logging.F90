module logging
    !! Module to handle logging.
    !!
    !! ## Functions
    !! - log_info: Record a message to stdout and/or fort.16.
    !! - log_error: Record a message to stderr and fort.16, and then terminate.

#ifdef CMPI
    use messenger_elem, only: message_fini
#endif

    implicit none

    integer, parameter :: screen_unit = 6
    !! Unit for standard Fortran screen output.
    integer, parameter :: file_unit = 16
    !! Unit for fort.16 log file.

    contains

    subroutine log_info(msg, myproc, to_screen, to_file)
        !! Record logging info to the console and/or file.

        implicit none

        character(len=*), intent(in) :: msg
        !! The message to be recorded.
        integer, optional, intent(in) :: myproc
        !! The MPI process number from which the log was called.
        logical, optional, intent(in) :: to_screen
        !! Whether to record the message to screen. Default is true.
        logical, optional, intent(in) :: to_file
        !! Whether to record the message to file. Default is true.

        ! For defaulting
        logical :: the_to_screen !! Whether to record the message to screen
        logical :: the_to_file !! Whether to record the message to file

        if (present(to_screen)) then
            the_to_screen = to_screen
        else
            the_to_screen = .true.
        end if

        if (present(to_file)) then
            the_to_file = to_file
        else
            the_to_file = .true.
        end if

        if (the_to_screen) then
            if (present(myproc)) write(screen_unit, "(a, i4)") "Process ", myproc
            write(screen_unit, "(a)") msg
        end if

        if (the_to_file) then
            if (present(myproc)) write(file_unit, "(a, i4)") "Process ", myproc
            write(file_unit, "(a)") msg
        end if
    end subroutine log_info

    subroutine log_error(msg, myproc)
        !! Throw an error and terminate the program. Record logging info to the
        !! console and file.

        implicit none

        character(len=*), intent(in) :: msg
        !! The message to be recorded.
        integer, intent(in) :: myproc
        !! The MPI process number from which the log was called.

        call log_info(msg, myproc, to_screen=.false.)
#ifdef CMPI
        call message_fini()
#endif
        error stop msg !! Ensure that message goes to stderr, not stdout
    end subroutine log_error
end module logging
